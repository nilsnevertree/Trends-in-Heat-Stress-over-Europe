
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
from matplotlib.colors import BoundaryNorm
import matplotlib.ticker as mticker
# from matplotlib.ticker import MaxNLocator
import matplotlib.colors as mplcolors
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

def area_plot(data = None, ax = None, levels = np.arange(20,38,1), contour_levels = np.arange(20,40,2), add_colorbar=False) :
    #create colormap
    cmap = plt.get_cmap("RdBu_r", len(levels))



    plot = data.plot.pcolormesh(ax = ax, cmap = cmap, \
        x = "longitude", y = "latitude",\
        yincrease = True, transform = ccrs.PlateCarree(), \
        norm = mplcolors.BoundaryNorm(levels, ncolors=len(levels) -1, clip=False) \
        )
    
    contourplot = data.plot.contour(ax= ax, levels= contour_levels, colors= 'k', linestyles = '-', linewidths = 0.2, x="longitude", y="latitude", \
        yincrease=True, transform=ccrs.PlateCarree(),\
        add_colorbar= add_colorbar)

    # plot gridlines and coastlines
    ax.coastlines(resolution = 'auto', color = 'k', zorder = 3)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                  linewidth=1, color='gray', alpha=0.5, linestyle=':')
    # gl.xlabels_top = False
    # gl.ylabels_right = False
    gl.xlines = True
    # gl.xlocator = mticker.FixedLocator(range(-180,181,30))
    # gl.xformatter = LONGITUDE_FORMATTER
    # gl.yformatter = LATITUDE_FORMATTER
    return plot
    # gl.xlabel_style = {'size': 15, 'color': 'gray'}
    # gl.xlabel_style = {'color': 'red', 'weight': 'bold'}
    #ax.add_feature(cart.feature.OCEAN, zorder=2, facecolor='w', edgecolor='k')