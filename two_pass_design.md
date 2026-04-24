# Two-pass CSR build with dynamic num_reg — design plan

## Problems solved

| # | Problem | Root cause |
|---|---|---|
| 1 | Fatal error at runtime | `update_reg_list` aborts when any coarse box accumulates more than `max_reg_list_num` entries |
| 2 | ~5.5 GB peak build memory | `reg_list_lon/lat(nrx, nry, 800)` 3-D scratch array allocated upfront regardless of actual need |
| 3 | Wrong coarse grid size for large or small grids | Hardcoded step tiers cover only three cases; anything outside the documented range is silently mis-sized |

## What is NOT changed

- `quad_lon_lat_locate_ii` — unchanged
- `quad_lon_lat_evaluate_*` — unchanged
- All MPI usage — unchanged
- The type of `grid_lon_list`, `grid_lat_list`, `grid_num`, `grid_start` — unchanged
- Public API of `quad_utils_mod` — unchanged

---

## Part 1 — Dynamic `num_reg` sizing

### Current code (in `init_quad_interp`, `GRID_QUAD_FULLY_IRREGULAR` case)

```fortran
if (num_lats * num_lons > 6 * 1000 * 1000) then
   interp_handle%ii%num_reg_x = 900
   interp_handle%ii%num_reg_y = 900
   interp_handle%ii%max_reg_list_num = 800
else if (num_lats * num_lons > 250 * 1000) then
   interp_handle%ii%num_reg_x = 180
   interp_handle%ii%num_reg_y = 180
   interp_handle%ii%max_reg_list_num = 800
else
   interp_handle%ii%num_reg_x = 90
   interp_handle%ii%num_reg_y = 90
   interp_handle%ii%max_reg_list_num = 800
endif
```

### Problems

- Only three cases. A 5000×4500 grid (22.5 M points) hits the first tier: nrx=nry=900.
  Each coarse box covers ~22500/810000 = 2.8% of the domain. That is already far
  too coarse — many boxes will exceed 800 candidates.
- `max_reg_list_num = 800` is a hard limit: any grid where more than 800 target quads
  overlap a single coarse box causes an immediate fatal error.
- The tiers are tied to POP ocean grids (gx1, tx0.5, tx0.1) and have no basis for
  general irregular grids.

### Replacement

Choose `num_reg` so that each coarse box contains approximately
`TARGET_CANDIDATES` target quads on average.  The coarse grid has
`nrx × nry` boxes; the target grid has approximately `nx × ny` quads
(ignoring boundary effects).  Setting `nrx = nry = num_reg`:

```
num_reg^2 ≈ nx * ny / TARGET_CANDIDATES
num_reg   ≈ sqrt(nx * ny / TARGET_CANDIDATES)
```

Clamp to a minimum of 10 (avoid degenerate tiny coarse grids) and a maximum of
`min(nx, ny)` (coarse grid no finer than the target grid in either dimension).

Add a named constant to the module (or namelist — see note at end):

```fortran
integer, parameter :: TARGET_CANDIDATES = 8
```

Replace the entire `if / else if / else` block with:

```fortran
case(GRID_QUAD_FULLY_IRREGULAR)
   allocate(interp_handle%ii%lats_2D(num_lons, num_lats), &
            interp_handle%ii%lons_2D(num_lons, num_lats))
   interp_handle%ii%lats_2D = MISSING_R8
   interp_handle%ii%lons_2D = MISSING_R8

   interp_handle%ii%num_reg_x = max(10, min(num_lons, &
         int(sqrt(real(num_lons * num_lats, r8) / real(TARGET_CANDIDATES, r8)))))
   interp_handle%ii%num_reg_y = max(10, min(num_lats, &
         int(sqrt(real(num_lons * num_lats, r8) / real(TARGET_CANDIDATES, r8)))))

   if (debug > 0) then
      write(string1,'(A,I0,A,I0)') 'num_reg_x=', interp_handle%ii%num_reg_x, &
                                    ' num_reg_y=', interp_handle%ii%num_reg_y
      call log_it(string1)
   endif

   allocate(interp_handle%ii%grid_start(interp_handle%ii%num_reg_x, &
                                        interp_handle%ii%num_reg_y))
   allocate(interp_handle%ii%grid_num(  interp_handle%ii%num_reg_x, &
                                        interp_handle%ii%num_reg_y))
   interp_handle%ii%grid_num = 0
```

`max_reg_list_num` is removed from the type (see Part 3).  It is no longer needed.

### Resulting `num_reg` values

| Grid | nx × ny | Current | New (TARGET=8) |
|---|---|---|---|
| gx1v6 | 320 × 384 | 90 | 124 |
| tx0.5v1 | 720 × 480 | 180 | 207 |
| tx0.1v2 | 3600 × 2400 | 900 | 1039 |
| 5000 × 4500 | 22.5 M | 900 | 1677 |
| 50000 × 45000 | 2.25 G | 900 | 16771 |

---

## Part 2 — Two-pass CSR build

### Current algorithm in `init_irreg_interp`

```
1. allocate reg_list_lon(nrx, nry, maxlist)   ! PEAK: nrx*nry*maxlist*4 bytes
   allocate reg_list_lat(nrx, nry, maxlist)   ! same
2. loop over all quads → update_reg_list      ! counts AND stores in one pass
3. u_total = sum(grid_num)
4. allocate grid_lon_list(u_total), grid_lat_list(u_total)
5. copy reg_list_lon/lat → flat lists
6. deallocate reg_list_lon, reg_list_lat
```

Peak memory = `2 × nrx × nry × maxlist × 4` bytes.  At nrx=nry=900, maxlist=800:
`2 × 900 × 900 × 800 × 4 = 5.18 GB`.

The fatal error in `update_reg_list` fires when any box's count exceeds `maxlist`
before storage is exhausted.  It is not possible to pre-determine `maxlist` safely
without either over-allocating or risking the abort.

### Replacement algorithm

**Pass 1 — count only.**  Iterate over all quads; for each coarse box the quad
overlaps, increment `grid_num(ix, iy)`.  No list storage needed.  Peak memory:
`nrx × nry × 4` bytes (~11 MB at nrx=nry=1677).

**Prefix sum.**  Compute `grid_start` from `grid_num`.  Derive `u_total = sum(grid_num)`.
Allocate the exact-sized flat lists: `grid_lon_list(u_total)`, `grid_lat_list(u_total)`.

**Pass 2 — fill.**  Iterate over all quads again with an integer cursor array
`fill_pos(nrx, nry)` initialised from `grid_start`.  For each coarse box the quad
overlaps, write the quad's indices into the flat lists at `fill_pos(ix,iy)` and
advance the cursor.

No 3-D temporary is ever allocated.  No hard cap exists.

### New subroutines replacing `update_reg_list`

#### `count_reg_overlaps`

Increments `grid_num` for every coarse box overlapped by quad `(i, j)`.
Does not store anything else.

```fortran
subroutine count_reg_overlaps(grid_num, reg_lon_ind, reg_lat_ind, nrx, nry)
integer, intent(inout) :: grid_num(:,:)
integer, intent(inout) :: reg_lon_ind(2)
integer, intent(inout) :: reg_lat_ind(2)
integer, intent(in)    :: nrx, nry

integer :: ind_x, index_x, ind_y, index_y

if (reg_lon_ind(2) < reg_lon_ind(1)) reg_lon_ind(2) = reg_lon_ind(2) + nrx

do ind_x = reg_lon_ind(1), reg_lon_ind(2)
   index_x = ind_x
   if (index_x > nrx) index_x = index_x - nrx
   do ind_y = reg_lat_ind(1), reg_lat_ind(2)
      index_y = ind_y
      if (index_y > nry) index_y = index_y - nry
      grid_num(index_x, index_y) = grid_num(index_x, index_y) + 1
   enddo
enddo

end subroutine count_reg_overlaps
```

#### `fill_reg_lists`

Uses a cursor array to write the quad's indices into the flat lists.

```fortran
subroutine fill_reg_lists(fill_pos, grid_lon_list, grid_lat_list, &
                          reg_lon_ind, reg_lat_ind, nrx, nry, &
                          grid_lon_index, grid_lat_index)
integer, intent(inout) :: fill_pos(:,:)           ! cursor: current write position
integer, intent(inout) :: grid_lon_list(:)
integer, intent(inout) :: grid_lat_list(:)
integer, intent(inout) :: reg_lon_ind(2)
integer, intent(inout) :: reg_lat_ind(2)
integer, intent(in)    :: nrx, nry
integer, intent(in)    :: grid_lon_index, grid_lat_index

integer :: ind_x, index_x, ind_y, index_y, pos

if (reg_lon_ind(2) < reg_lon_ind(1)) reg_lon_ind(2) = reg_lon_ind(2) + nrx

do ind_x = reg_lon_ind(1), reg_lon_ind(2)
   index_x = ind_x
   if (index_x > nrx) index_x = index_x - nrx
   do ind_y = reg_lat_ind(1), reg_lat_ind(2)
      index_y = ind_y
      if (index_y > nry) index_y = index_y - nry
      pos = fill_pos(index_x, index_y)
      grid_lon_list(pos) = grid_lon_index
      grid_lat_list(pos) = grid_lat_index
      fill_pos(index_x, index_y) = pos + 1
   enddo
enddo

end subroutine fill_reg_lists
```

### New `init_irreg_interp` structure

```fortran
subroutine init_irreg_interp(h)
type(quad_interp_handle), intent(inout) :: h

integer, allocatable :: fill_pos(:,:)
real(r8) :: u_c_lons(4), u_c_lats(4), pole_row_lon
integer  :: i, j, nrx, nry, nx, ny, u_total, u_index, istatus
integer  :: reg_lon_ind(2), reg_lat_ind(2)
logical  :: cyclic, pole
integer  :: xlim

! --- setup (unchanged from current) ---
cyclic = h%opt%spans_lon_zero
pole   = h%opt%pole_wrap
nx = h%nlon;  ny = h%nlat
nrx = h%ii%num_reg_x;  nry = h%ii%num_reg_y

! ... min/max lon/lat setup, pole_x/u_pole_y/t_pole_y (unchanged) ...

! --- PASS 1: count ---
h%ii%grid_num = 0

do i = 1, xlim
   do j = 1, ny - 1
      if (all_corners_valid(h%opt, i, j, nx)) then
         call get_quad_corners(h%ii%lons_2d, i, j, cyclic, pole, nx, ny, u_c_lons, istatus)
         call get_quad_corners(h%ii%lats_2d, i, j, cyclic, pole, nx, ny, u_c_lats, istatus)
         call reg_box_overlap(h, u_c_lons, u_c_lats, .false., reg_lon_ind, reg_lat_ind)
         call count_reg_overlaps(h%ii%grid_num, reg_lon_ind, reg_lat_ind, nrx, nry)
      endif
   enddo
enddo

! --- prefix sum to build grid_start ---
u_total = 0
do j = 1, nry
   do i = 1, nrx
      h%ii%grid_start(i, j) = u_total + 1
      u_total = u_total + h%ii%grid_num(i, j)
   enddo
enddo

write(string1,'(A,I0)') 'total coarse-index entries = ', u_total
call error_handler(E_MSG, 'init_irreg_interp', string1)

! --- allocate exact-sized flat lists ---
allocate(h%ii%grid_lon_list(u_total))
allocate(h%ii%grid_lat_list(u_total))

! --- PASS 2: fill ---
allocate(fill_pos(nrx, nry))
fill_pos = h%ii%grid_start   ! cursor starts at each box's start position

do i = 1, xlim
   do j = 1, ny - 1
      if (all_corners_valid(h%opt, i, j, nx)) then
         call get_quad_corners(h%ii%lons_2d, i, j, cyclic, pole, nx, ny, u_c_lons, istatus)
         call get_quad_corners(h%ii%lats_2d, i, j, cyclic, pole, nx, ny, u_c_lats, istatus)
         call reg_box_overlap(h, u_c_lons, u_c_lats, .false., reg_lon_ind, reg_lat_ind)
         call fill_reg_lists(fill_pos, h%ii%grid_lon_list, h%ii%grid_lat_list, &
                             reg_lon_ind, reg_lat_ind, nrx, nry, i, j)
      endif
   enddo
enddo

deallocate(fill_pos)

! --- consistency check (replaces current u_index check) ---
u_index = sum(h%ii%grid_num) + 1
if (u_index /= u_total + 1) then
   call error_handler(E_ERR, 'init_irreg_interp', &
      'Storage indices did not balance: contact DART developers', &
      source, revision, revdate)
endif

end subroutine init_irreg_interp
```

Note: `get_quad_corners` is called twice per quad per pass.  The two passes together
call it 4× per quad, the same total as the current single-pass approach (2× per quad).
If profiling shows this is a bottleneck, the loop body can be factored into a helper
that is called from both passes.

---

## Part 3 — Type changes

### Remove `max_reg_list_num`

It has no purpose after the two-pass build.  Remove the field from
`quad_irreg_grid_coords`:

**Before:**
```fortran
integer  :: num_reg_x = 180
integer  :: num_reg_y = 180
integer  :: max_reg_list_num = 800
```

**After:**
```fortran
integer  :: num_reg_x = 180
integer  :: num_reg_y = 180
```

---

## Part 4 — Memory comparison

For a 5000×4500 grid with TARGET_CANDIDATES=8 (nrx=nry=1677):

| Phase | Current peak | Proposed peak |
|---|---|---|
| Build-time scratch | `2 × 900 × 900 × 800 × 4` = **5.18 GB** | `1677 × 1677 × 4` = **11 MB** (count array) |
| Persistent `grid_num` + `grid_start` | `2 × 900 × 900 × 4` = 6.5 MB | `2 × 1677 × 1677 × 4` = 22 MB |
| Flat lists (exact) | allocated after build | same — allocated after build |

The build-time peak drops from **5.18 GB to ~11 MB**.

---

## Part 5 — Summary of all file changes

One file only: `DART/models/utilities/quad_utils_mod.f90`

| Location | Change |
|---|---|
| `quad_irreg_grid_coords` type | Remove `max_reg_list_num` field |
| `init_quad_interp` — `GRID_QUAD_FULLY_IRREGULAR` case | Replace 3-tier `if/else` with single dynamic `num_reg` formula; remove `max_reg_list_num` assignments |
| `init_irreg_interp` | Replace single-pass + 3-D temporary with two-pass build using `count_reg_overlaps` + prefix sum + `fill_reg_lists` |
| `update_reg_list` | Replace with two new private subroutines: `count_reg_overlaps` and `fill_reg_lists` |
| Everything else | **No change** |

No changes to `mpi_utilities_mod.f90`, `null_mpi_utilities_mod.f90`, or
`mpif08_utilities_mod.f90`.

---

## Note: `TARGET_CANDIDATES` as a namelist parameter

`TARGET_CANDIDATES = 8` means on average 8 target quads per coarse box.
A smaller value makes the coarse grid finer (fewer candidates per box, faster search,
more memory for `grid_num`/`grid_start`/flat lists).  A larger value does the
opposite.  If users need to tune this, it can be added to the existing
`quad_utils_nml` namelist without any further API change.
