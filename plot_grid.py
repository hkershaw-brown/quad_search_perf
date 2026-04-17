"""
plot_grid.py
------------
Plot an irregular lon/lat grid (as produced by create_grid_files.py) on a map.

Grid points (cell centres) are drawn as dots; quad edges are drawn as lines
so the distortion of each cell is clearly visible.

Usage
-----
    python plot_grid.py                          # uses irregular_grid.nc
    python plot_grid.py mygrid.nc
    python plot_grid.py mygrid.nc --no-edges     # dots only, faster for large grids
    python plot_grid.py mygrid.nc -o mygrid.png  # save instead of display
"""

import argparse
import sys
import numpy as np
import netCDF4 as nc
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    HAS_CARTOPY = True
except ImportError:
    HAS_CARTOPY = False


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def load_grid(filename):
    with nc.Dataset(filename) as ds:
        lons = ds.variables['longitudes'][:]   # (nx, ny)
        lats = ds.variables['latitudes'][:]    # (nx, ny)
        global_grid = getattr(ds, 'grid_type', 'regional') == 'global'
        cyclic_lon  = getattr(ds, 'cyclic_lon', 'false') == 'true'
    return np.asarray(lons), np.asarray(lats), global_grid, cyclic_lon


def _normalize_lon(lon):
    """Wrap longitude(s) to [-180, 180]."""
    return ((np.asarray(lon) + 180.0) % 360.0) - 180.0


def _plot_segment(ax, lon0, lat0, lon1, lat1, transform, color, alpha, lw):
    """
    Plot a single edge between two points.  Longitudes are first normalised
    to [-180, 180].  After normalisation most wrap-around pairs (e.g. -2° and
    2°) have |Δlon| ≤ 180 and are drawn directly.  Only segments that
    genuinely straddle the antimeridian (e.g. 179° to -179°, a real grid
    crossing at 180°E) are split into two pieces.
    """
    lon0_n = float(_normalize_lon(lon0))
    lon1_n = float(_normalize_lon(lon1))

    if abs(lon1_n - lon0_n) <= 180.0:
        # Short path – no antimeridian crossing
        ax.plot([lon0_n, lon1_n], [lat0, lat1],
                '-', color=color, alpha=alpha, linewidth=lw, transform=transform)
        return

    # Segment genuinely straddles ±180° – split at the antimeridian
    # Extend lon1 to the same side as lon0 for linear interpolation
    lon1_ext = lon1_n + 360.0 if lon0_n > 0.0 else lon1_n - 360.0
    frac = (180.0 - abs(lon0_n)) / abs(lon1_ext - lon0_n)
    lat_cross = lat0 + frac * (lat1 - lat0)
    boundary = 180.0 if lon0_n > 0.0 else -180.0
    ax.plot([lon0_n,   boundary], [lat0,      lat_cross],
            '-', color=color, alpha=alpha, linewidth=lw, transform=transform)
    ax.plot([-boundary, lon1_n],  [lat_cross, lat1],
            '-', color=color, alpha=alpha, linewidth=lw, transform=transform)


def draw_quad_edges(ax, lons, lats, transform, cyclic_lon=False,
                    color='steelblue', alpha=0.6, lw=0.5):
    """Draw the edges of every quad as line segments.

    For cyclic grids the zonal ring is closed by connecting column nx-1
    back to column 0.  Segments that cross the antimeridian are split to
    avoid spurious lines across the map.
    """
    nx, ny = lons.shape
    kw = dict(transform=transform, color=color, alpha=alpha, lw=lw)

    # Normalise all longitudes to [-180, 180] once
    lons_n = _normalize_lon(lons)

    # Horizontal edges (constant j, vary i)
    for j in range(ny):
        for i in range(nx - 1):
            _plot_segment(ax, lons_n[i, j], lats[i, j],
                          lons_n[i+1, j], lats[i+1, j], **kw)
        # Closing wrap-around segment for cyclic grids
        if cyclic_lon:
            _plot_segment(ax, lons_n[nx-1, j], lats[nx-1, j],
                          lons_n[0, j],    lats[0, j],    **kw)

    # Vertical edges (constant i, vary j)
    for i in range(nx):
        ax.plot(lons_n[i, :], lats[i, :],
                '-', color=color, alpha=alpha, linewidth=lw, transform=transform)


def make_figure_cartopy(lons, lats, global_grid, cyclic_lon, draw_edges, output):
    if global_grid:
        projection = ccrs.Robinson()
    else:
        lon_c = 0.5 * (lons.min() + lons.max())
        lat_c = 0.5 * (lats.min() + lats.max())
        projection = ccrs.PlateCarree(central_longitude=lon_c)

    fig, ax = plt.subplots(subplot_kw={'projection': projection},
                           figsize=(12, 7))

    if global_grid:
        ax.set_global()
    else:
        pad_lon = max(5.0, 0.1 * (lons.max() - lons.min()))
        pad_lat = max(5.0, 0.1 * (lats.max() - lats.min()))
        ax.set_extent([lons.min() - pad_lon, lons.max() + pad_lon,
                       lats.min() - pad_lat, lats.max() + pad_lat],
                      crs=ccrs.PlateCarree())

    ax.add_feature(cfeature.LAND,       facecolor='wheat',  zorder=0)
    ax.add_feature(cfeature.OCEAN,      facecolor='lightcyan', zorder=0)
    ax.add_feature(cfeature.COASTLINE,  linewidth=0.5, zorder=1)
    ax.add_feature(cfeature.BORDERS,    linewidth=0.3, linestyle=':', zorder=1)
    ax.gridlines(draw_labels=True, linewidth=0.3, color='grey', alpha=0.5)

    transform = ccrs.PlateCarree()

    if draw_edges:
        draw_quad_edges(ax, lons, lats, transform=transform, cyclic_lon=cyclic_lon)

    ax.scatter(_normalize_lon(lons).ravel(), lats.ravel(), s=6, c='red', zorder=3,
               transform=transform, label='grid points')

    nx, ny = lons.shape
    ax.set_title(f'Irregular grid  ({nx} × {ny})  –  {lons.size} cell centres',
                 fontsize=12)

    _save_or_show(fig, output)


def make_figure_plain(lons, lats, cyclic_lon, draw_edges, output):
    """Fallback when cartopy is not available."""
    fig, ax = plt.subplots(figsize=(12, 7))

    ax.set_facecolor('lightcyan')

    if draw_edges:
        draw_quad_edges(ax, lons, lats, transform=ax.transData, cyclic_lon=cyclic_lon)

    ax.scatter(lons.ravel(), lats.ravel(), s=6, c='red', zorder=3,
               label='grid points')

    nx, ny = lons.shape
    ax.set_xlabel('Longitude (°E)')
    ax.set_ylabel('Latitude (°N)')
    ax.set_title(f'Irregular grid  ({nx} × {ny})  –  {lons.size} cell centres\n'
                 '(install cartopy for map background)',
                 fontsize=11)
    ax.legend(handles=[mlines.Line2D([], [], marker='o', color='red',
                                     markersize=4, linestyle='none',
                                     label='grid points')],
              loc='lower right')

    _save_or_show(fig, output)


def _save_or_show(fig, output):
    if output:
        fig.savefig(output, dpi=150, bbox_inches='tight')
        print(f'Figure saved to {output}')
    else:
        plt.tight_layout()
        plt.show()


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description='Plot an irregular lon/lat grid on a map.')
    p.add_argument('input', nargs='?', default='irregular_grid.nc',
                   help='NetCDF grid file (default: irregular_grid.nc)')
    p.add_argument('--no-edges', dest='draw_edges', action='store_false',
                   help='Do not draw quad edges (useful for very large grids)')
    p.add_argument('-o', '--output', default=None,
                   help='Save figure to this file instead of displaying it '
                        '(e.g. grid.png, grid.pdf)')
    return p.parse_args()


def main():
    args = parse_args()

    print(f'Loading {args.input} ...')
    lons, lats, global_grid, cyclic_lon = load_grid(args.input)
    nx, ny = lons.shape
    print(f'  Grid size : {nx} x {ny}')
    print(f'  Lon range : {lons.min():.2f} to {lons.max():.2f}')
    print(f'  Cyclic lon: {cyclic_lon}')
    print(f'  Lat range : {lats.min():.2f} to {lats.max():.2f}')
    print(f'  Type      : {"global" if global_grid else "regional"}')

    if not HAS_CARTOPY:
        print('cartopy not found – drawing without map background')
        make_figure_plain(lons, lats, cyclic_lon, args.draw_edges, args.output)
    else:
        make_figure_cartopy(lons, lats, global_grid, cyclic_lon, args.draw_edges, args.output)


if __name__ == '__main__':
    main()
