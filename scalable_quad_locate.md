# Scalable quad locate — design plan

## Problems

| # | Problem | Impact |
|---|---|---|
| 1 | `max_reg_list_num = 800` hard cap | Fatal error when any coarse box has > 800 candidate quads (correctness) |
| 2 | `reg_list_lon/lat(nrx, nry, 800)` 3-D temporary | ~5.5 GB peak at nrx=nry=900 (memory) |
| 3 | Hardcoded `num_reg` tiers | Undefined behaviour for grids larger than 1/10° |
| 4 | `lons_2d(nx,ny)` + `lats_2d` replicated on every rank | 360 MB × ranks/node. With round-robin placement and many ranks per node this is fatal |

Problems 1–3 are solved by **two-pass CSR + dynamic num_reg** (no MPI required).  
Problem 4 is solved by **distributing the 2D grid arrays across all ranks with MPI RMA**.

---

## Design overview

```
At init:
  - Each rank allocates a row slab of lons_2d and lats_2d  (nx × local_ny rows, +1 halo)
  - All ranks cooperatively build the coarse index (distributed count, MPI_Allreduce, replicated result)
  - Each rank exposes its slab as an MPI_Win_create window
  - Passive epoch opened with MPI_Win_lock_all

At search (quad_lon_lat_locate_ii):
  - Coarse box lookup — local (coarse index is replicated, tiny)
  - Candidate quad list — local (flat lists replicated)
  - Corner fetch — MPI_Get from owner rank's window, batched per observation
  - in_quad test — local after fetch

At finalize:
  - MPI_Win_unlock_all, MPI_Win_free, deallocate slabs
```

---

## Part 1 — Data distribution

### Row slab assignment

Divide the ny rows across P ranks using standard block distribution:

```
j_start(r) = r * (ny / P) + min(r,   ny mod P) + 1
j_end(r)   = j_start(r) + (ny / P) + merge(1, 0, r < ny mod P) - 1
local_ny(r) = j_end(r) - j_start(r) + 1
```

Each rank r allocates:

```fortran
allocate(h%ii%lons_slab(nx, local_ny + 1))   ! +1 = halo row
allocate(h%ii%lats_slab(nx, local_ny + 1))
```

The halo row is row `j_end(r) + 1` (first row of rank r+1's slab).  It ensures all
four corners of every quad `(i, j)` with `j` in `[j_start(r), j_end(r)]` are present
in rank r's local array without a remote fetch:

```
quad (i, j)  corners:  (i,j), (i+1,j)   → row j     = slab row j_local
                        (i,j+1),(i+1,j+1) → row j+1   = halo row j_local+1
```

**Edge case — cyclic j:** quad `(i, ny-1)` has corners at row ny and row 1 (wrap).
Row 1 is on rank 0.  This quad is handled as a remote fetch regardless of which
rank owns row ny-1.  This affects only the single wraparound row and is documented
in `fetch_candidate_corners` below.

**Longitude wrap (cyclic i):** column `i+1 = 1` wraps within the same rank's slab
(all nx columns are local).  No special treatment needed.

### What is replicated vs distributed

| Data | Replicated or distributed | Rationale |
|---|---|---|
| `grid_num(nrx, nry)` | Replicated | Small: ~11 MB at nrx=nry=1677 |
| `grid_start(nrx, nry)` | Replicated | Same |
| `grid_lon_list(:)` | Replicated | ~360 MB at 5000×4500; distribute in a further step if needed |
| `grid_lat_list(:)` | Replicated | Same |
| `lons_2d(nx, ny)` | Distributed (row slabs) | 180 MB total → 180/P MB per rank |
| `lats_2d(nx, ny)` | Distributed (row slabs) | Same |

For the current target grid (5000×4500) with P=128 ranks/node and round-robin
placement:

| Array | Per rank (current) | Per rank (proposed) |
|---|---|---|
| lons + lats | 360 MB | 360 MB / 128 = 2.8 MB |
| coarse index | 22 MB | 22 MB (replicated) |
| flat lists | 360 MB | 360 MB (replicated) |
| **Total** | **742 MB** | **385 MB** |

For P → very large (50000×45000 grid), flat lists also need distribution.
That is a further step; the architecture here supports it because the flat lists
are only read during search via the `grid_lon_list` / `grid_lat_list` arrays, which
can equally be exposed as RMA windows.

---

## Part 2 — MPI_Win_create windows

After each rank has populated its slab and halo, expose it:

```fortran
! disp_unit = 8 bytes (one real(r8) element)
call MPI_Win_create(h%ii%lons_slab, &
                    int(nx*(local_ny+1), MPI_ADDRESS_KIND) * 8_MPI_ADDRESS_KIND, &
                    8, MPI_INFO_NULL, get_dart_mpi_comm(), h%ii%win_lons, ierr)
call MPI_Win_create(h%ii%lats_slab, &
                    int(nx*(local_ny+1), MPI_ADDRESS_KIND) * 8_MPI_ADDRESS_KIND, &
                    8, MPI_INFO_NULL, get_dart_mpi_comm(), h%ii%win_lats, ierr)
```

Two windows: one for lons, one for lats.  Both over the full DART communicator
so any rank can fetch from any other.

### Passive epoch — opened once, closed at finalize

After window creation and slab population:

```fortran
call MPI_Win_lock_all(MPI_MODE_NOCHECK, h%ii%win_lons, ierr)
call MPI_Win_lock_all(MPI_MODE_NOCHECK, h%ii%win_lats, ierr)
```

`MPI_MODE_NOCHECK` asserts that no conflicting locks will be acquired — valid
because the slabs are **read-only** after init.

The epoch is closed in `finalize_quad_interp`:

```fortran
call MPI_Win_unlock_all(h%ii%win_lons, ierr)
call MPI_Win_unlock_all(h%ii%win_lats, ierr)
call MPI_Win_free(h%ii%win_lons, ierr)
call MPI_Win_free(h%ii%win_lats, ierr)
```

---

## Part 3 — Displacement calculation

Element `lons_slab(i, j_local)` of rank r is stored at byte offset:

```
disp = ((i - 1) * (local_ny_r + 1) + (j_local - 1))   ! in units of disp_unit = 8
```

Given quad `(i, j)`:

```fortran
r          = owner_rank(j, h%ii%j_starts, h%ii%num_tasks)   ! binary search on j_starts
j_local    = j - h%ii%j_starts(r) + 1
local_ny_r = h%ii%local_nys(r)

! Bottom corners  (row j)
disp_bl = (i      - 1) * (local_ny_r + 1) + (j_local - 1)   ! (i,   j  )
disp_br = (i_next - 1) * (local_ny_r + 1) + (j_local - 1)   ! (i+1, j  )

! Top corners  (row j+1 = j_local+1, in rank r's halo)
disp_tl = (i      - 1) * (local_ny_r + 1) + j_local          ! (i,   j+1)
disp_tr = (i_next - 1) * (local_ny_r + 1) + j_local          ! (i+1, j+1)
```

`owner_rank` is O(log P) binary search on a P-element array, called O(candidates)
times per observation — negligible cost.

---

## Part 4 — Batched corner fetch (`fetch_candidate_corners`)

Rather than one `MPI_Get` per element, issue all gets for all candidate quads
before any `in_quad` test, then flush once.

New private subroutine in `quad_utils_mod`:

```fortran
subroutine fetch_candidate_corners(h, start_ind, num_inds, &
                                   lon_corners, lat_corners)
type(quad_interp_handle), intent(in)  :: h
integer,                  intent(in)  :: start_ind, num_inds
real(r8),                 intent(out) :: lon_corners(4, num_inds)
real(r8),                 intent(out) :: lat_corners(4, num_inds)

integer :: k, r, i, j, i_next, j_local, local_ny_r
integer(dart_mpi_addr_kind) :: disps(4)

do k = 1, num_inds
   i      = h%ii%grid_lon_list(start_ind + k - 1)
   j      = h%ii%grid_lat_list(start_ind + k - 1)
   i_next = i + 1
   if (i_next > h%nlon) i_next = 1   ! cyclic lon wrap

   r          = owner_rank(j, h%ii%j_starts, h%ii%num_tasks)
   j_local    = j - h%ii%j_starts(r) + 1
   local_ny_r = h%ii%local_nys(r)

   disps(1) = (i      - 1) * (local_ny_r + 1) + (j_local - 1)
   disps(2) = (i_next - 1) * (local_ny_r + 1) + (j_local - 1)
   disps(3) = (i      - 1) * (local_ny_r + 1) + j_local
   disps(4) = (i_next - 1) * (local_ny_r + 1) + j_local

   if (r == h%ii%my_rank) then
      ! Local — direct array copy, no MPI
      lon_corners(1,k) = h%ii%lons_slab(i,      j_local    )
      lon_corners(2,k) = h%ii%lons_slab(i_next, j_local    )
      lon_corners(3,k) = h%ii%lons_slab(i,      j_local + 1)
      lon_corners(4,k) = h%ii%lons_slab(i_next, j_local + 1)
      lat_corners(1,k) = h%ii%lats_slab(i,      j_local    )
      lat_corners(2,k) = h%ii%lats_slab(i_next, j_local    )
      lat_corners(3,k) = h%ii%lats_slab(i,      j_local + 1)
      lat_corners(4,k) = h%ii%lats_slab(i_next, j_local + 1)
   else
      ! Remote — non-blocking MPI_Gets; completed by rma_flush_all below
      call rma_get_real8(h%ii%win_lons, r, disps(1), lon_corners(1,k), 1)
      call rma_get_real8(h%ii%win_lons, r, disps(2), lon_corners(2,k), 1)
      call rma_get_real8(h%ii%win_lons, r, disps(3), lon_corners(3,k), 1)
      call rma_get_real8(h%ii%win_lons, r, disps(4), lon_corners(4,k), 1)
      call rma_get_real8(h%ii%win_lats, r, disps(1), lat_corners(1,k), 1)
      call rma_get_real8(h%ii%win_lats, r, disps(2), lat_corners(2,k), 1)
      call rma_get_real8(h%ii%win_lats, r, disps(3), lat_corners(3,k), 1)
      call rma_get_real8(h%ii%win_lats, r, disps(4), lat_corners(4,k), 1)
   endif
enddo

! Complete all outstanding gets before caller tests in_quad
call rma_flush_all(h%ii%win_lons)
call rma_flush_all(h%ii%win_lats)

end subroutine fetch_candidate_corners
```

The calling code in `quad_lon_lat_locate_ii` then tests each candidate locally:

```fortran
call fetch_candidate_corners(interp_handle, start_ind, num_inds, &
                              lon_corners, lat_corners)
found = .false.
do k = 1, num_inds
   if (in_quad(lon, lat, lon_corners(:,k), lat_corners(:,k))) then
      lon_bot = interp_handle%ii%grid_lon_list(start_ind + k - 1)
      lat_bot = interp_handle%ii%grid_lat_list(start_ind + k - 1)
      found = .true.
      exit
   endif
enddo
```

The final `get_quad_corners` calls for the winning quad reuse the already-fetched
`lon_corners(:, k)` and `lat_corners(:, k)` — no additional RMA needed.

---

## Part 5 — Distributed coarse index build

Each rank builds only the entries for its own quads (j in [j_start, j_end]).
A global `MPI_Allreduce(MPI_SUM)` assembles the complete counts.  No rank ever
holds the full 2D grid during init.

### Build steps (init_irreg_interp)

1. **Dynamic num_reg:**
   ```fortran
   integer, parameter :: TARGET_CANDIDATES = 8
   nrx = max(10, min(nx, int(sqrt(real(nx*ny, r8) / TARGET_CANDIDATES))))
   nry = max(10, min(ny, int(sqrt(real(nx*ny, r8) / TARGET_CANDIDATES))))
   ```

2. **Slab allocation and fill:** each rank allocates `lons_slab(nx, local_ny+1)`,
   `lats_slab(nx, local_ny+1)` and populates rows [j_start, j_end] from the
   coordinate arrays passed in by `set_quad_coords`.  The halo row is exchanged
   with the next rank via a point-to-point (`send_to` / `receive_from`) using
   existing DART wrappers.

3. **Pass 1 — count (distributed):** each rank iterates only its quads
   (j in [j_start, j_end]), calls `reg_box_overlap`, accumulates into local
   `my_grid_num(nrx, nry)`.  Then global reduction:
   ```fortran
   call MPI_Allreduce(MPI_IN_PLACE, h%ii%grid_num, nrx*nry, &
                      MPI_INTEGER, MPI_SUM, get_dart_mpi_comm(), ierr)
   ```
   Every rank now holds the complete global `grid_num`.

4. **Prefix sum:** every rank computes `grid_start` identically — no communication.
   ```fortran
   u_total = 0
   do j = 1, nry
     do i = 1, nrx
       h%ii%grid_start(i,j) = u_total + 1
       u_total = u_total + h%ii%grid_num(i,j)
     enddo
   enddo
   ```

5. **Pass 2 — fill flat lists (distributed):** each rank allocates
   `grid_lon_list(u_total)` and `grid_lat_list(u_total)` initialised to 0.  It
   fills only the entries for its own quads using a cursor array `fill_pos(nrx,nry)`
   initialised from `grid_start`.  Then global reduction:
   ```fortran
   call MPI_Allreduce(MPI_IN_PLACE, h%ii%grid_lon_list, u_total, &
                      MPI_INTEGER, MPI_SUM, get_dart_mpi_comm(), ierr)
   call MPI_Allreduce(MPI_IN_PLACE, h%ii%grid_lat_list, u_total, &
                      MPI_INTEGER, MPI_SUM, get_dart_mpi_comm(), ierr)
   deallocate(fill_pos)
   ```
   This works because each entry is written by exactly one rank (the quad's owner)
   and is 0 on all others, so the sum is the correct global list.

6. **Create windows and open epoch** (Part 2).

---

## Part 6 — New wrappers in `mpi_utilities_mod`

All three files (`mpi_utilities_mod.f90`, `mpif08_utilities_mod.f90`,
`null_mpi_utilities_mod.f90`) must have identical public signatures.

### Handle-type strategy

All new wrapper signatures use plain `integer` for win handles.  The `mpif08`
version bridges internally via `.MPI_VAL` (the single integer component that
`mpi_f08` derived types are required to have by the standard).

### `dart_mpi_addr_kind` constant

```fortran
! mpi_utilities_mod.f90 and mpif08_utilities_mod.f90
integer, parameter, public :: dart_mpi_addr_kind = MPI_ADDRESS_KIND

! null_mpi_utilities_mod.f90
integer, parameter, public :: dart_mpi_addr_kind = 8
```

### New public routines (add to all three files)

```
create_rma_win_real8, open_rma_epoch, close_rma_epoch,
rma_get_real8, rma_flush_all, free_rma_win,
dart_mpi_addr_kind
```

---

#### `create_rma_win_real8`

```fortran
subroutine create_rma_win_real8(base, n_elem, win)
  real(r8), intent(in), target :: base(*)
  integer(dart_mpi_addr_kind), intent(in)  :: n_elem
  integer,                     intent(out) :: win
```

**`mpi_utilities_mod.f90`** (`use mpi`):
```fortran
call MPI_Win_create(base, n_elem * 8_MPI_ADDRESS_KIND, 8, &
                    MPI_INFO_NULL, my_local_comm, win, ierr)
```

**`mpif08_utilities_mod.f90`** (`use mpi_f08`):
```fortran
type(MPI_Win) :: w
call MPI_Win_create(base, n_elem * 8_MPI_ADDRESS_KIND, 8, &
                    MPI_INFO_NULL, my_local_comm, w)
win = w%MPI_VAL
```

**`null_mpi_utilities_mod.f90`**: `win = 0`  (1 rank, all data local, window never queried)

---

#### `open_rma_epoch`

```fortran
subroutine open_rma_epoch(win)
  integer, intent(in) :: win
```

**`mpi_utilities_mod.f90`:** `call MPI_Win_lock_all(MPI_MODE_NOCHECK, win, ierr)`  
**`mpif08_utilities_mod.f90`:** `w%MPI_VAL = win; call MPI_Win_lock_all(MPI_MODE_NOCHECK, w)`  
**`null_mpi_utilities_mod.f90`:** no-op

---

#### `close_rma_epoch`

```fortran
subroutine close_rma_epoch(win)
  integer, intent(in) :: win
```

**`mpi_utilities_mod.f90`:** `call MPI_Win_unlock_all(win, ierr)`  
**`mpif08_utilities_mod.f90`:** `w%MPI_VAL = win; call MPI_Win_unlock_all(w)`  
**`null_mpi_utilities_mod.f90`:** no-op

---

#### `rma_get_real8`

Non-blocking get of `count` `real(r8)` elements from target rank's window.
Completes only after a subsequent `rma_flush_all`.

```fortran
subroutine rma_get_real8(win, target_rank, target_disp, buf, count)
  integer,                     intent(in)  :: win
  integer,                     intent(in)  :: target_rank
  integer(dart_mpi_addr_kind), intent(in)  :: target_disp   ! in units of disp_unit
  real(r8),                    intent(out) :: buf(count)
  integer,                     intent(in)  :: count
```

**`mpi_utilities_mod.f90`:**
```fortran
call MPI_Get(buf, count, MPI_DOUBLE_PRECISION, target_rank, target_disp, &
             count, MPI_DOUBLE_PRECISION, win, ierr)
```

**`mpif08_utilities_mod.f90`:**
```fortran
type(MPI_Win) :: w
w%MPI_VAL = win
call MPI_Get(buf, count, MPI_DOUBLE_PRECISION, target_rank, target_disp, &
             count, MPI_DOUBLE_PRECISION, w)
```

**`null_mpi_utilities_mod.f90`:** `E_ERR` stub — never reached because the local
branch in `fetch_candidate_corners` is always taken when there is only 1 rank.

---

#### `rma_flush_all`

Completes all outstanding `MPI_Get` calls on a window.

```fortran
subroutine rma_flush_all(win)
  integer, intent(in) :: win
```

**`mpi_utilities_mod.f90`:** `call MPI_Win_flush_all(win, ierr)`  
**`mpif08_utilities_mod.f90`:** `w%MPI_VAL = win; call MPI_Win_flush_all(w)`  
**`null_mpi_utilities_mod.f90`:** no-op

---

#### `free_rma_win`

```fortran
subroutine free_rma_win(win)
  integer, intent(inout) :: win
```

**`mpi_utilities_mod.f90`:** `call MPI_Win_free(win, ierr)`  
**`mpif08_utilities_mod.f90`:** `w%MPI_VAL = win; call MPI_Win_free(w); win = MPI_WIN_NULL%MPI_VAL`  
**`null_mpi_utilities_mod.f90`:** no-op

---

## Part 7 — Changes to `quad_irreg_grid_coords` type

```fortran
type quad_irreg_grid_coords
   private

   ! ---- scalar metadata ----
   integer  :: num_reg_x = 180
   integer  :: num_reg_y = 180
   ! max_reg_list_num REMOVED — eliminated by two-pass CSR
   real(r8) :: min_lon, max_lon, lon_width, min_lat, max_lat, lat_width

   ! ---- distribution metadata ----
   integer :: my_rank    = -1
   integer :: num_tasks  = -1
   integer :: my_j_start = -1
   integer :: my_j_end   = -1
   integer :: my_local_ny = -1
   integer, allocatable :: j_starts(:)    ! (0:P-1) first j-index on each rank (1-based)
   integer, allocatable :: local_nys(:)   ! (0:P-1) local_ny for each rank

   ! ---- local row slabs (nx × (local_ny+1), last row = halo) ----
   real(r8), allocatable :: lons_slab(:,:)
   real(r8), allocatable :: lats_slab(:,:)

   ! ---- RMA window handles (integer regardless of mpi/mpi_f08 binding) ----
   integer :: win_lons = -1
   integer :: win_lats = -1

   ! ---- replicated coarse index ----
   integer, allocatable :: grid_start(:,:)
   integer, allocatable :: grid_num  (:,:)
   integer, allocatable :: grid_lon_list(:)
   integer, allocatable :: grid_lat_list(:)

end type quad_irreg_grid_coords
```

`lons_2D` and `lats_2D` (full grid arrays) are **removed**.  All coordinate access
goes through `lons_slab` / `lats_slab` (local) or `rma_get_real8` (remote).

---

## Part 8 — Changes to `quad_lon_lat_locate_ii`

The `get_grid_quad` call (which currently takes full 2D arrays) is replaced by the
`fetch_candidate_corners` + local `in_quad` loop (Part 4).  The existing
`get_quad_corners` calls for the winning quad's corner indices are replaced by
reuse of the already-fetched arrays.

The `quad_lon_lat_evaluate_*` routines are **unchanged** — they take integer
indices and field values from the model, not from the internal grid arrays.

---

## Part 9 — Summary of all file changes

| File | Change |
|---|---|
| `mpi_utilities_mod.f90` | Add 6 new routines + `dart_mpi_addr_kind` (integer handles, `use mpi`) |
| `mpif08_utilities_mod.f90` | Same 6 routines + constant — bridge `type(MPI_Win)` via `.MPI_VAL` |
| `null_mpi_utilities_mod.f90` | Same 6 stubs (no-op or E_ERR) + `dart_mpi_addr_kind = 8` |
| `quad_utils_mod.f90` — type | Remove `max_reg_list_num`, `lons_2D`, `lats_2D`; add slab arrays, window handles, distribution metadata |
| `quad_utils_mod.f90` — `init_quad_interp` | Replace hardcoded tiers with dynamic formula |
| `quad_utils_mod.f90` — `init_irreg_interp` | Rewrite: slab alloc, halo exchange, distributed two-pass CSR, `create_rma_win_real8`, `open_rma_epoch` |
| `quad_utils_mod.f90` — `update_reg_list` | Replace with `count_reg_overlaps` (pass 1) + `fill_reg_lists` (pass 2) |
| `quad_utils_mod.f90` — `quad_lon_lat_locate_ii` | Replace `get_grid_quad` + `get_quad_corners` with `fetch_candidate_corners` + local `in_quad` loop |
| `quad_utils_mod.f90` — `finalize_quad_interp` | Add `close_rma_epoch`, `free_rma_win`, deallocate slabs and distribution arrays |
| `quad_utils_mod.f90` — `quad_lon_lat_evaluate_*` | **No change** |
| `quad_utils_mod.f90` — `use` statements | Add `dart_mpi_addr_kind`, `create_rma_win_real8`, `open_rma_epoch`, `close_rma_epoch`, `rma_get_real8`, `rma_flush_all`, `free_rma_win` from `mpi_utilities_mod` |
