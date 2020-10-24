
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
from matplotlib.colors import BoundaryNorm
import matplotlib.ticker as mticker
# from matplotlib.ticker import MaxNLocator
import matplotlib.colors as mplcolors
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import make_axes_locatable

def area_plot(data = None, ax = None, 
        levels = np.arange(20,38,1), contour_levels = np.arange(20,40,2),
        add_colorbar=False, cmap = "RdBu",
        landcolor = "lightgrey", oceancolor = "lightskyblue",
        zorder_ocean = 2, zorder_land = 2,
        gridline_kwargs = {draw_labels : False, linewidth : 1, color :'k', alpha : 0.5, linestyle :':', zorder : 8}) :

    '''
    this function drwas the provided data and the gridliines and land and ocean to a geoaxes object provided.
    x and y dimensions need to be dimensions of the xarray dataarray
    INPUT:
    1. General
        data            - xarray DataArray with longitude and latitude as dimensions
        ax              - geoaxes object in which the data shall be drawn
    2. Plot OPtions
        levels          - levels of the pcolormesh plot
        contour_levels  - additional contour line levels
        add_colorbar    - False will add colorbar based on the add_cbar function (True will allow xarray to plot the colorbar)
        cmap            - colormap for the pcolormesh plot (default:"RdBu")
    Gridlines and further
        NOT WORKING NOW -> gridline_kwargs - dictonary with kwargs for ax.gridlines()
        landcolor       - color of land (default:"lightgrey")
        oceancolor      - color of ocean (default:"lightskyblue")
        zorder_ocean    - layer on which the land shall be plotted (default: behind data)
        zorder_land     - layer on which the ocean shall be plotted (default: behind data)
    OUTPUT: 
    - mapple objects of the pcolormesh (pmesh) and the contour plot (ct) 
    '''
    #create colormap
    cmap = plt.get_cmap(cmap, len(levels))

    # plot pcolormesh of the data with according colormap
    pmesh = data.plot.pcolormesh(ax = ax, cmap = cmap, \
        x = "longitude", y = "latitude",\
        yincrease = True, transform = ccrs.PlateCarree(), \
        norm = mplcolors.BoundaryNorm(levels, ncolors=len(levels) -1, clip=False), \
        zorder=3, add_colorbar= False)

    ct = contourplot = data.plot.contour(ax= ax, levels= contour_levels, colors= 'k',\
        linestyles = '-', linewidths = 0.5, x="longitude", y="latitude", \
        yincrease=True, transform=ccrs.PlateCarree(),\
        add_colorbar= False,\
        zorder=4)

    # plot gridlines and coastlines
    ax.coastlines(resolution = 'auto', color = 'k', )
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,\
                linewidth=1, color='k', alpha=0.5, linestyle=':',\
                zorder=8)
    # gl.xlabels_top = False
    # gl.ylabels_right = False
    gl.xlines = True
    # gl.xlocator = mticker.FixedLocator(range(-30,50,10))
    # gl.xformatter = LONGITUDE_FORMATTER
    # gl.yformatter = LATITUDE_FORMATTER

    # gl.xlabel_style = {'size': 15, 'color': 'gray'}
    # gl.xlabel_style = {'color': 'red', 'weight': 'bold'}
    ax.add_feature(cart.feature.OCEAN, facecolor= oceancolor, zorder=zorder_ocean)
    ax.add_feature(cart.feature.LAND, facecolor= landcolor, zorder=zorder_land)

    return pmesh, ct

def plot_30y_area(data = None, mask_res = None, 
                cmap= "Reds", cmap_days = "Reds", save_name = False, 
                mask_dict = {}):

    '''
    this function creates a figure with area_plot() function fro all heatwave days
    a) mean temperature 
    b) mean swbgt 
    c) days which are part of a heatwave 

    INPUT:
    1. data         - xarray DataSet with longitude and latitude as dimensions
                    - !! need to include variables t2m and swbgt and dimension time, latitude, longitude
                    - !! should be daily data
    2. mask_res     - mask which is True wherever there is a heatwave occuring
    3. cmap         - cmap for temperature and swbgt
    4. cmap_days    - cmap for days of heatwave plot
    5. save_name    - if False no save else save png and svg
    6. mask_dict    - "threshold" : what was the threshold for this mask,
                    - "duration" : what was the duration
                    - "quantile" : what was the quantile provided?
    OUTPUT: 
    - figure with the plots
    '''

    # create year range
    year_start = data.time[0].dt.strftime("%Y").values
    year_end = year = data.time[-1].dt.strftime("%Y").values
    years = int((data.time[1] - data.time[0]).dt.strftime("%Y").values)
    # get arguments for the mask (by default : NONE GIVEN)
    duration = mask_dict.setdefault("duration" , np.nan)
    threshold = mask_dict.setdefault("threshold" , np.nan)
    quantile = mask_dict.setdefault("quantile" , np.nan)
    #create figure title :
    figure_title = year_start + " - " + year_end + "\nheatwave similar to DWD definition \n {:.0f} days duration,\n {:.0f}°C daily max. temperature threshold,\n {:.0f}th percentile of temperature".format(duration, threshold, quantile*100)
    
    # ========== PLOT ==============
    # create the figure and the gridspecs
    fig = plt.figure(figsize = (9,27))
    gs = gridspec.GridSpec(3, 1)

    # plot the mean temperature for the heatwaves (where the mask_res is True)
    plot_data = data.t2m.where(mask_res, other = np.nan).mean("time")
    ax = fig.add_subplot(gs[0,0],\
            projection=ccrs.NearsidePerspective(
                            central_latitude=50.72,
                            central_longitude=10.53,
                            satellite_height=10000000.0))

    ax.coastlines(resolution='110m')
    ax.gridlines()
    pm, ct = area_plot(data= plot_data, ax = ax, \
        levels = np.arange(28, 42, 0.5), \
        contour_levels= np.arange(25,40,2), cmap = cmap,
        landcolor = [0.4,0.4,0.4])
    
    add_cbar(fig = fig, ax=ax, pm = pm, label = "Temperature in °C")
    ax.set_title("a) mean temperature during heat wave")

    # plot the mean temperature for the heatwaves (where the mask_res is True)
    plot_data = data.swbgt.where(mask_res, other = np.nan).mean("time")
    
    ax1 = fig.add_subplot(gs[1,0],\
            projection=ccrs.NearsidePerspective(
                            central_latitude=50.72,
                            central_longitude=10.53,
                            satellite_height=10000000.0))

    ax1.coastlines(resolution='110m')
    ax1.gridlines()
    pm, ct = area_plot(data= plot_data, ax = ax1, \
        levels = np.arange(22, 35, 0.5), \
        contour_levels= np.arange(15,35,2), \
        cmap = cmap,
        landcolor = [0.4,0.4,0.4])

    cbar = add_cbar(fig = fig, ax=ax1, pm = pm, label = "sWBGT")
    ax1.set_title("b) mean sWBGT during heat wave")

    # plot the days which are part of a heatwave
    plot_data = np.log(mask_res.sum("time")) # sum up all True values give the value for days 
    plot_data = plot_data.where(plot_data != 0)

    ax2 = fig.add_subplot(gs[2,0],\
            projection=ccrs.NearsidePerspective(
                            central_latitude=50.72,
                            central_longitude=10.53,
                            satellite_height=10000000.0))

    ax2.coastlines(resolution='110m')
    ax2.gridlines()

    pm, ct = area_plot(data= plot_data, ax = ax2, \
        levels =  np.arange(1, 10, 1), \
        contour_levels= np.arange(1,10,2), \
        cmap = cmap_days)
    ax2.set_title("heatwave days for the period - log scale")
    cbar = add_cbar(fig = fig, ax = ax2, pm = pm, label = "log(days)")

    # create figure title
    fig.suptitle(figure_title, fontsize= 15)

    # ============ SAVE ===============
    # save data if wanted
    year = data.time[0].dt.strftime("%Y-%m-%d").values
    if save_name:
        fig.savefig(plot_dict["savefolder"] + save_name + '.png')
        fig.savefig(plot_dict["savefolder"] + save_name + '.svg')
    return fig

def plot_city_data(data = None, data_max = None,
        mask_res = None, city_name = "None", 
        timeslice= slice(timestart, timeend), groupby= "time.year",
        color_bar = "tomato" , color_scatter = "tomato" , color_hist = "tomato",
        mask_dict = {}):
    
    duration = mask_dict.setdefault("duration" , "NONE GIVEN")
    threshold = mask_dict.setdefault("threshold" , "NONE GIVEN")
    quantile = mask_dict.setdefault("quantile" , "NONE GIVEN")

    try: 
            data = data.mean("longitude").mean("latitude")
            data_max = data_max.mean("longitude").mean("latitude")
    except: 
            pass

    fig = plt.figure(figsize = (15,15))
    gs = gridspec.GridSpec(3, 2)

    fig.suptitle(city_name + ": " + timeslice.start + " - " + timeslice.stop + "\nheatwave similar to DWD definition \n {:.0f} days duration, {:.0f}°C daily max. temperature threshold, {:.0f}th percentile of temperature".format(duration, threshold, quantile*100),\
            fontsize= 15)
    # plot temperature
    ax1 = fig.add_subplot(gs[0,:])
    data.t2m.sel(time = timeslice)\
            .plot(ax = ax1, color = "k", linewidth = "0.3",linestyle = ":", alpha= 0.75,\
            label= "a) temperature - 6 hourly")
    data_max.t2m.sel(time = timeslice)\
            .plot(ax = ax1, color = "tab:blue", linewidth = "0.5",\
            label= "b) temperature - daily max")
    data_max.t2m.where(mask_res).sel(time = timeslice)\
            .plot(ax = ax1, marker = 'x', color = "tab:orange",\
            label= "c) heatwaves days")
    ax1.set_ylim([-5,35])
    ax1.set_title("")
    ax1.set_ylabel("temperature in °C")
    leg = ax1.legend(loc = "lower left")
    ax1.grid()

    # plot bar plot
    ax2 = fig.add_subplot(gs[1,:], sharex = ax1)
    y = mask_res.sel(time = timeslice).groupby(groupby).sum('time')

    if "year" in groupby:
        # x = y.year
        x = data_max.time[data_max.time.dt.is_year_start].sel(time = timeslice)
        w = (x[1]-x[0])/2
        ax2.bar(x = x + w , height = y, width = w, align = "center")
    elif "month" in groupby:
        x = data_max.time[data_max.time.dt.is_month_start].sel(time = timeslice)
        w = (x[1]-x[0])/2
        ax2.bar(x = x + w , height = y, width = w/2, align = "center")
    elif "week" in groupby:
        x = x.week
        w = (x[1]-x[0])/2
    else:
        x = y.time
    
    
    # city(mask_org).sel(time = timeslice).groupby('time.month').sum('time').plot(ax = ax3,linestyle= ":", label = "xclim used 6 hourly")
    ax2.set_ylabel("heat wave days per year")
    
    # ax2.legend()
    ax2.grid()
    # ax1.set_xticklabels([])

    ax3 = fig.add_subplot(gs[2:3,0])
    x = data_max.t2m.where(mask_res).sel(time = timeslice)
    y = data_max.swbgt.where(mask_res).sel(time = timeslice)
    c = data_max.d2m.where(mask_res).sel(time = timeslice)
    # c = data_max["time.month"].values
    sct = ax3.scatter(x,y)
    ax3.set_xlabel("temperature in °C")
    ax3.set_ylabel("sWBGT")
    ax3.set_xlim([24,38])
    ax3.set_ylim([18,30])
    ax3.grid()

    ax4 = fig.add_subplot(gs[2:3,1])
    bins = np.arange(19,28)
    ax4.hist(x = y, bins = bins ,density = True)
    ax4.set_xlabel("sWBGT")
    ax4.set_ylabel("density")
    ax4.set_ylim([0,0.35])
    ax4.grid()

    return fig

def add_cbar(figure, axes, mapple_object, label = False, percentage = "5%", pad = 0.1) :
    # this function adds a colobar to the axes object given to it by dividing the 
    divider = make_axes_locatable(axes)
    axes_cb = divider.new_horizontal(size= percentage, pad= pad, axes_class= plt.Axes)
    figure.add_axes(axes_cb)
    cbar = plt.colorbar(mapple_object, cax= axes_cb, orientation='vertical')
    if label :
        cbar.ax.get_yaxis().labelpad = 15
        cbar.ax.set_ylabel(label, rotation=270)
    return cbar
