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

Output file (default: irregular_grid.nc):
    python create_grid_files.py --nx 20 --ny 15 -o mygrid.nc
"""

import argparse
import numpy as np
import netCDF4 as nc


def make_irregular_grid(nx, ny, lon_start, lon_end, lat_start, lat_end, perturb):
    """
    Build a fully irregular (nx, ny) grid by adding sinusoidal perturbations
    to a regular base grid.  The perturbation uses opposing sin/cos waves in
    i and j so that every row and column is shifted differently, producing
    genuinely non-rectangular quads while keeping quads non-self-intersecting
    for perturb < 0.4.

    Returns
    -------
    lons, lats : ndarray, shape (nx, ny)  – in degrees
    """
    # Base (regular) spacing
    if nx > 1:
        dlon = (lon_end - lon_start) / (nx - 1)
    else:
        dlon = 0.0
    if ny > 1:
        dlat = (lat_end - lat_start) / (ny - 1)
    else:
        dlat = 0.0

    i_idx = np.arange(nx)
    j_idx = np.arange(ny)
    II, JJ = np.meshgrid(i_idx, j_idx, indexing='ij')   # shape (nx, ny)

    # Normalised phase angles
    phi_i = np.pi * II / max(nx - 1, 1)
    phi_j = np.pi * JJ / max(ny - 1, 1)

    lons = lon_start + II * dlon + perturb * dlon * np.sin(phi_j)
    lats = lat_start + JJ * dlat + perturb * dlat * np.sin(phi_i)

    return lons.astype(np.float64), lats.astype(np.float64)


def write_netcdf(filename, lons, lats, global_grid,
                 lon_start, lon_end, lat_start, lat_end, perturb):
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
        args.nx, args.ny, lon_start, lon_end, lat_start, lat_end, args.perturb)

    write_netcdf(args.output, lons, lats, args.global_grid,
                 lon_start, lon_end, lat_start, lat_end, args.perturb)


if __name__ == '__main__':
    main()
