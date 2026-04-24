# `quad_utils_mod` — Summary

## Overview

`quad_utils_mod` provides initialisation, search, and bilinear interpolation for **logically-rectangular** lon/lat grids. The grid topology is hidden behind a single opaque handle (`quad_interp_handle`) so callers are independent of the underlying grid type.

---

## Grid types

| Constant | Storage | Description |
|---|---|---|
| `GRID_QUAD_FULLY_REGULAR` | origin + delta per axis | Uniform spacing; enclosing quad found by direct index arithmetic |
| `GRID_QUAD_IRREG_SPACED_REGULAR` | 1-D array per axis | Orthogonal but uneven; binary search along each 1-D axis |
| `GRID_QUAD_FULLY_IRREGULAR` | 2-D array `lons(nx,ny)`, `lats(nx,ny)` | Fully distorted (e.g. POP dipole grid); two-phase search described below |

---

## Locating an observation in a `GRID_QUAD_FULLY_IRREGULAR` grid

The fully-irregular case cannot use index arithmetic or 1-D searches because neither axis is monotone. Instead it uses a **coarse-grid pre-filter** to limit the number of candidates that need a precise geometric test.

### Phase 1 — Setup (`init_irreg_interp`)

Executed once when `set_quad_coords` is called.

```
┌─────────────────────────────────────────────────────┐
│  Overlay a coarse regular grid (num_reg_x × num_reg_y) │
│  that covers the full domain.                         │
│                                                        │
│  Resolution is chosen by grid size:                    │
│    nx·ny > 6M  →  900 × 900  (≈ 1/10° grids)         │
│    nx·ny > 250k →  180 × 180  (≈ 1/2° grids)          │
│    otherwise   →   90 ×  90                           │
└─────────────────────────────────────────────────────┘
```

For **every** target quad `(i, j)`:

1. **Get the 4 corners** of the quad from `lons_2d` / `lats_2d` via `get_quad_corners`.  
   Corners are stored counter-clockwise: `(i,j)`, `(i+1,j)`, `(i+1,j+1)`, `(i,j+1)` with wrap/pole handling.

2. **Find which coarse boxes it overlaps** (`reg_box_overlap`).  
   - For lat: trivial — min/max of corner lats mapped linearly.  
   - For lon: handle wrap-around by detecting `max-min > 180°` and separately finding the smallest value > 180° and largest < 180°.  
   - Pole-containing quads claim all longitude boxes.

3. **Store the mapping** (temporary 3-D list, then compacted into a CSR-style flat array):

```fortran
grid_num  (rx, ry)       ! number of target quads overlapping coarse box (rx,ry)
grid_start(rx, ry)       ! index into the flat list where that coarse box's candidates begin
grid_lon_list(u_index)   ! flat list of candidate i indices
grid_lat_list(u_index)   ! flat list of candidate j indices
```

### Phase 2 — Search (`quad_lon_lat_locate_ii`)

Called at assimilation time for each observation.

```
given: lon, lat

1. get_reg_box_indices             → (rx, ry)
   └── get_reg_lon_box: rx = int( num_reg_x * (lon - min_lon) / lon_width ) + 1
   └── get_reg_lat_box: ry = int( num_reg_y * (lat - min_lat) / lat_width ) + 1

   Non-global: if (rx,ry) outside coarse grid → istatus≠0, return

2. Retrieve candidate list for coarse box (rx, ry)
   num_inds  = grid_num  (rx, ry)
   start_ind = grid_start(rx, ry)

   If num_inds == 0 → no target quad can contain this point → istatus=1

3. get_grid_quad — loop over the num_inds candidates
   For each candidate quad (i,j):
     a. get_quad_corners  → x_corners(4), y_corners(4)
     b. in_quad(lon, lat, x_corners, y_corners)

4. in_quad uses a vertical ray-casting test:
   Cast a vertical line x = lon upward/downward from the point.
   For each of the 4 sides of the quad call line_intercept:
     - if colinear and point not on side → cant_be_in_box (early exit)
     - if point lies on a side           → in_box (immediate true)
     - otherwise count intersections above and below the point
   Point is inside iff at least one intersection above AND one below.
   (Lon wrap-around inside a side is handled by shifting values by +360°
    when max-min > 180°.)

5. On success: fill four_lon_indices, four_lat_indices
     corner 1: (lon_bot,   lat_bot)
     corner 2: (lon_top,   lat_bot)
     corner 3: (lon_top,   lat_top)
     corner 4: (lon_bot,   lat_top)
   where lon_top, lat_top come from quad_index_neighbors (handles cyclic/pole).
```

### Complexity

| Step | Cost |
|---|---|
| Setup | $O(N_x N_y)$ — once per run |
| Coarse box lookup | $O(1)$ |
| Candidate loop (`get_grid_quad`) | $O(K)$ where $K$ is the number of target quads overlapping that coarse box (typically a small constant) |
| `in_quad` per candidate | $O(1)$ — 4 sides |

The coarse-grid pre-filter reduces what would be an $O(N_x N_y)$ exhaustive search to an effectively $O(1)$ lookup for each observation.

---

## Key data structures inside `quad_interp_handle`

```
h%ii%lons_2d(nx,ny)       target grid longitudes
h%ii%lats_2d(nx,ny)       target grid latitudes
h%ii%num_reg_x/y          coarse grid dimensions
h%ii%grid_num (rx,ry)     how many target quads touch each coarse box
h%ii%grid_start(rx,ry)    offset into flat candidate lists
h%ii%grid_lon_list(:)     flat list: target quad i-indices
h%ii%grid_lat_list(:)     flat list: target quad j-indices
```
