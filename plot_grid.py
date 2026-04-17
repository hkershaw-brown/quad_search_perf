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
    return np.asarray(lons), np.asarray(lats), global_grid


def draw_quad_edges(ax, lons, lats, transform, color='steelblue', alpha=0.6, lw=0.5):
    """Draw the edges of every quad as line segments."""
    nx, ny = lons.shape

    # Horizontal edges (constant j, vary i)
    for j in range(ny):
        ax.plot(lons[:, j], lats[:, j],
                '-', color=color, alpha=alpha, linewidth=lw, transform=transform)

    # Vertical edges (constant i, vary j)
    for i in range(nx):
        ax.plot(lons[i, :], lats[i, :],
                '-', color=color, alpha=alpha, linewidth=lw, transform=transform)


def make_figure_cartopy(lons, lats, global_grid, draw_edges, output):
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
        draw_quad_edges(ax, lons, lats, transform=transform)

    ax.scatter(lons.ravel(), lats.ravel(), s=6, c='red', zorder=3,
               transform=transform, label='grid points')

    nx, ny = lons.shape
    ax.set_title(f'Irregular grid  ({nx} × {ny})  –  {lons.size} cell centres',
                 fontsize=12)

    _save_or_show(fig, output)


def make_figure_plain(lons, lats, draw_edges, output):
    """Fallback when cartopy is not available."""
    fig, ax = plt.subplots(figsize=(12, 7))

    ax.set_facecolor('lightcyan')

    if draw_edges:
        draw_quad_edges(ax, lons, lats, transform=ax.transData)

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
    lons, lats, global_grid = load_grid(args.input)
    nx, ny = lons.shape
    print(f'  Grid size : {nx} x {ny}')
    print(f'  Lon range : {lons.min():.2f} to {lons.max():.2f}')
    print(f'  Lat range : {lats.min():.2f} to {lats.max():.2f}')
    print(f'  Type      : {"global" if global_grid else "regional"}')

    if not HAS_CARTOPY:
        print('cartopy not found – drawing without map background')
        make_figure_plain(lons, lats, args.draw_edges, args.output)
    else:
        make_figure_cartopy(lons, lats, global_grid, args.draw_edges, args.output)


if __name__ == '__main__':
    main()
