"""
create_grid_files.py
--------------------
Generate a NetCDF file containing a logically-rectangular but fully irregular
(distorted) lon/lat grid for use by test_quad_search.f90.

Usage examples
--------------
Global grid, 36 x 18 cells:
    python create_grid_files.py --global --nx 36 --ny 18

Regional grid, 20 x 15 cells, lon 20-110 E, lat 40 S - 40 N:
    python create_grid_files.py --nx 20 --ny 15 \
        --lon-start 20 --lon-end 110 --lat-start -40 --lat-end 40

Optional perturbation magnitude (fraction of grid spacing, default 0.3):
    python create_grid_files.py --nx 20 --ny 15 --perturb 0.2

Random seed for reproducibility (default: 42):
    python create_grid_files.py --nx 20 --ny 15 --seed 123

Output file (default: irregular_grid.nc):
    python create_grid_files.py --nx 20 --ny 15 -o mygrid.nc
"""

import argparse
import numpy as np
import netCDF4 as nc


def _periodic_perturbation(nx, ny, rng, n_modes):
    """
    Return a perturbation field of shape (nx, ny) that is exactly periodic
    in the i direction with period nx, using random Fourier modes.
    Values are normalised to [-1, 1].
    """
    i_idx = np.arange(nx)
    j_idx = np.arange(ny)
    II, JJ = np.meshgrid(i_idx, j_idx, indexing='ij')

    field = np.zeros((nx, ny))
    for k in range(1, n_modes + 1):
        amp    = rng.uniform(0.0, 1.0)
        phases = rng.uniform(0.0, 2.0 * np.pi, ny)   # independent phase per row
        # sin(2π k i / nx + phase(j)) is periodic in i with period nx
        field += amp * np.sin(2.0 * np.pi * k * II / nx + phases[JJ])

    max_abs = np.abs(field).max()
    if max_abs > 0:
        field /= max_abs
    return field


def make_irregular_grid(nx, ny, lon_start, lon_end, lat_start, lat_end,
                        perturb, seed, cyclic_lon=False):
    """
    Build a fully irregular (nx, ny) grid by adding random perturbations
    to each grid point.  Every quad gets a different shape, giving a genuine
    stress-test for the locate algorithm.  Quads remain non-self-intersecting
    for perturb < 0.4.

    For global (cyclic) grids the longitude perturbation uses random Fourier
    modes so it is exactly periodic with period nx, keeping the wrap-around
    quad well-formed.  The grid spacing is 360/nx (no repeated endpoint).

    Parameters
    ----------
    cyclic_lon : bool
        True for global grids; enables periodic lon perturbation and
        cyclic spacing.
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    lons, lats : ndarray, shape (nx, ny)  – in degrees
    """
    if cyclic_lon:
        # Periodic spacing: nx points cover exactly 360°, no repeated endpoint
        dlon = 360.0 / nx
    else:
        dlon = (lon_end - lon_start) / (nx - 1) if nx > 1 else 0.0
    dlat = (lat_end - lat_start) / (ny - 1) if ny > 1 else 0.0

    i_idx = np.arange(nx)
    j_idx = np.arange(ny)
    II, JJ = np.meshgrid(i_idx, j_idx, indexing='ij')   # shape (nx, ny)

    rng = np.random.default_rng(seed)

    if cyclic_lon:
        n_modes = max(2, nx // 6)
        lon_perturb = _periodic_perturbation(nx, ny, rng, n_modes)
    else:
        lon_perturb = rng.uniform(-1, 1, (nx, ny))

    lat_perturb = rng.uniform(-1, 1, (nx, ny))

    lons = lon_start + II * dlon + perturb * dlon * lon_perturb
    lats = lat_start + JJ * dlat + perturb * dlat * lat_perturb

    # Clamp latitudes to valid range
    lats = np.clip(lats, -90.0, 90.0)

    return lons.astype(np.float64), lats.astype(np.float64)


def write_netcdf(filename, lons, lats, global_grid,
                 lon_start, lon_end, lat_start, lat_end, perturb, seed):
    nx, ny = lons.shape

    with nc.Dataset(filename, 'w', format='NETCDF4') as ds:
        # Global attributes
        ds.title = 'Irregular lon/lat grid for test_quad_search'
        ds.grid_type = 'global' if global_grid else 'regional'
        ds.nx = nx
        ds.ny = ny
        ds.lon_start = lon_start
        ds.lon_end = lon_end
        ds.lat_start = lat_start
        ds.lat_end = lat_end
        ds.perturbation_fraction = perturb
        ds.random_seed = seed
        ds.cyclic_lon = 'true' if global_grid else 'false'

        # Dimensions
        ds.createDimension('nx', nx)
        ds.createDimension('ny', ny)

        # Variables
        v_lon = ds.createVariable('longitudes', 'f8', ('nx', 'ny'))
        v_lon.units = 'degrees_east'
        v_lon.long_name = 'longitude of grid cell centres'
        v_lon[:, :] = lons

        v_lat = ds.createVariable('latitudes', 'f8', ('nx', 'ny'))
        v_lat.units = 'degrees_north'
        v_lat.long_name = 'latitude of grid cell centres'
        v_lat[:, :] = lats

    print(f'Written {nx} x {ny} irregular grid to {filename}')
    print(f'  lon range: {lons.min():.2f} to {lons.max():.2f} degrees_east')
    print(f'  lat range: {lats.min():.2f} to {lats.max():.2f} degrees_north')


def parse_args():
    p = argparse.ArgumentParser(
        description='Create a NetCDF file with an irregular lon/lat grid.')

    p.add_argument('--nx', type=int, default=36,
                   help='Number of grid points in the longitude direction (default: 36)')
    p.add_argument('--ny', type=int, default=18,
                   help='Number of grid points in the latitude direction (default: 18)')

    p.add_argument('--global', dest='global_grid', action='store_true',
                   help='Create a global grid (lon 0-360, lat -90 to 90). '
                        'Overrides --lon-start/end and --lat-start/end.')

    p.add_argument('--lon-start', type=float, default=20.0,
                   help='Regional: start longitude in degrees (default: 20)')
    p.add_argument('--lon-end',   type=float, default=110.0,
                   help='Regional: end longitude in degrees (default: 110)')
    p.add_argument('--lat-start', type=float, default=-40.0,
                   help='Regional: start latitude in degrees (default: -40)')
    p.add_argument('--lat-end',   type=float, default=40.0,
                   help='Regional: end latitude in degrees (default: 40)')

    p.add_argument('--perturb', type=float, default=0.3,
                   help='Perturbation magnitude as fraction of grid spacing '
                        '(0 = regular grid, <0.4 recommended; default: 0.3)')

    p.add_argument('--seed', type=int, default=42,
                   help='Random seed for reproducibility (default: 42)')

    p.add_argument('-o', '--output', default='irregular_grid.nc',
                   help='Output NetCDF filename (default: irregular_grid.nc)')

    return p.parse_args()


def main():
    args = parse_args()

    if args.perturb < 0 or args.perturb >= 0.5:
        raise ValueError('--perturb must be in [0, 0.5)')

    if args.global_grid:
        lon_start, lon_end = 0.0, 360.0
        lat_start, lat_end = -90.0, 90.0
    else:
        lon_start = args.lon_start
        lon_end   = args.lon_end
        lat_start = args.lat_start
        lat_end   = args.lat_end

        if lon_end <= lon_start:
            raise ValueError('--lon-end must be greater than --lon-start')
        if lat_end <= lat_start:
            raise ValueError('--lat-end must be greater than --lat-start')
        if lat_start < -90 or lat_end > 90:
            raise ValueError('Latitude bounds must be within [-90, 90]')

    lons, lats = make_irregular_grid(
        args.nx, args.ny, lon_start, lon_end, lat_start, lat_end, args.perturb, args.seed,
        cyclic_lon=args.global_grid)

    write_netcdf(args.output, lons, lats, args.global_grid,
                 lon_start, lon_end, lat_start, lat_end, args.perturb, args.seed)


if __name__ == '__main__':
    main()
