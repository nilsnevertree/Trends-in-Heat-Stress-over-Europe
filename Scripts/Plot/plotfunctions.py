
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import numpy as np
# from matplotlib.colors as mplcolors
import matplotlib.ticker as mticker
# from matplotlib.ticker import MaxNLocator
import matplotlib.colors as mplcolors
import matplotlib.gridspec as gridspec
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from mpl_toolkits.axes_grid1 import make_axes_locatable
from cartopy import feature
import Calculation.trend_func as trends
import cmocean as cmo
from memory_profiler import profile

rcParams_area = {
                'figure.titlesize': 17,
                'axes.titlesize': 17,
                'axes.labelsize': 17,
                'xtick.labelsize': 13,
                'ytick.labelsize': 13,
                'hatch.linewidth' : 1,
                'hatch.color' : 'k',
                'figure.facecolor':'w',
                'hatch.linewidth' : 1.5,
                'hatch.color' : 'w',
                'figure.facecolor':'w',
                'savefig.bbox': 'tight',
                'savefig.dpi': 300,
                'savefig.facecolor': 'w',
                'savefig.format': 'png',
                'savefig.pad_inches': 0.01,
                'savefig.transparent': False,
                'font.family': ['serif'],
                'font.sans-serif': ['Computer Modern Serif'],
                'legend.fontsize': 13,
                'font.weight': 'normal',
                'grid.alpha': 0.75,
                'grid.color': '#b0b0b0',
                'grid.linestyle': ':',

                'legend.fancybox': False,
                'legend.frameon': False,
                }


def area_plot(data = None, ax = None, 
        levels = [0,1], contour_levels = [0,1],
        add_colorbar=False, cmap = "RdBu",
        landcolor = False, oceancolor = "lightgrey", borders = True,
        pcolormesh = True, contourf = False, pcolor = False, contour = False,
        zorder_ocean = 2, zorder_land = 2,
        gridline_kwargs = dict(),
        pcolormesh_kwargs = dict(edgecolors= 'face', snap = True),
        contour_kwargs = dict(colors= 'k', linestyles = '--', linewidths = 1.5), 
        colorbar_kwargs = dict(),
        set_extent = True,
        axes_title = False):
        
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
    colorbar_kwargs.setdefault('divergent' , False)
    colorbar_kwargs.setdefault('add_colorbar' , True)
    colorbar_kwargs.setdefault('extend' , 'both')
    colorbar_kwargs.setdefault('contour_levels' , None)
    colorbar_kwargs.setdefault('norm' , False)
    contour_kwargs.setdefault('contour_levels' , None)

    # check how to extend the colorbar 

    if (data.max() > max(levels)) & (data.min() < min(levels)):
        colorbar_kwargs['extend'] = 'both'
    elif data.max() > max(levels) :
        colorbar_kwargs['extend'] = 'max'
    elif data.min() < min(levels) :
        colorbar_kwargs['extend'] = 'min'
    else :
        colorbar_kwargs['extend'] = 'neither'
    
    cmap_obj = plt.get_cmap(cmap, )
    
    if colorbar_kwargs['norm']:
        norm = colorbar_kwargs['norm']
    elif colorbar_kwargs['divergent'] :
        # create new levels so its evenly spaced
        length_div = sum(levels > 0) * 2 + 1
        max_div = np.max(np.abs(levels))
        levels = np.linspace(-max_div, max_div, length_div)
        # create coresponding norm
        norm = mplcolors.BoundaryNorm(levels, cmap_obj.N, extend = colorbar_kwargs['extend'])
    else : 
        norm = mplcolors.BoundaryNorm(levels, cmap_obj.N, extend = colorbar_kwargs['extend'])
            
    # plot the filling object
    if pcolormesh :
        mapple_object = ax.pcolormesh(data.longitude, data.latitude, data, cmap = cmap_obj, 
                transform = ccrs.PlateCarree(), 
                norm = norm, 
                zorder=3, **pcolormesh_kwargs)   
    elif contourf :
        mapple_object = ax.contourf(data.longitude, data.latitude, data, cmap = cmap_obj,
                yincrease = True, transform = ccrs.PlateCarree(), 
                norm = norm, \
                zorder=3, )
    elif pcolor :
        mapple_object = ax.pcolor(data.longitude, data.latitude, data, cmap = cmap_obj, 
                transform = ccrs.PlateCarree(),
                norm = norm, shading = 'nearest',
                zorder=3, **pcolormesh_kwargs)  
    else : 
        mapple_object = None

    if contour :
        ct = data.plot.contour(ax= ax, levels= contour_levels, 
            x="longitude", y="latitude",
            yincrease=True, transform=ccrs.PlateCarree(),
            add_colorbar= False,
            zorder=4, **contour_kwargs)
    # plot the colorbar
    if colorbar_kwargs['add_colorbar']:
            cbar = add_cbar(fig = ax.get_figure(), ax = ax, mapple_object= mapple_object, **colorbar_kwargs)

    # add ocean, land and borders
    if oceancolor :
        ax.add_feature(feature.OCEAN, facecolor= oceancolor, zorder=zorder_ocean)
    if landcolor :
        ax.add_feature(feature.LAND, facecolor= landcolor, zorder=zorder_land)
    if borders :
        ax.add_feature(feature.BORDERS, alpha = 0.3, zorder = 15)

        
    # plot gridlines and coastlines
    

    gridline_kwargs.setdefault('linewidth' , 1)
    gridline_kwargs.setdefault('color' , 'grey')
    gridline_kwargs.setdefault('alpha' , 0.5)
    gridline_kwargs.setdefault('linestyle' ,':')
    gridline_kwargs.setdefault('zorder' , 8)
    gridline_kwargs.setdefault('top_labels' , False)
    gridline_kwargs.setdefault('right_labels' , False)
    gridline_kwargs.setdefault('bottom_labels' , True)
    gridline_kwargs.setdefault('left_labels' , True)

    top_labels = gridline_kwargs.pop('top_labels')
    right_labels = gridline_kwargs.pop('right_labels')
    bottom_labels = gridline_kwargs.pop('bottom_labels')
    left_labels = gridline_kwargs.pop('left_labels')

    gl = ax.gridlines(crs=ccrs.PlateCarree(), ** gridline_kwargs)
    gl.xlines = True
    gl.top_labels = top_labels 
    gl.right_labels = right_labels 
    gl.bottom_labels = bottom_labels 
    gl.left_labels = left_labels
    gl.xlocator = mticker.FixedLocator(range(-10,51,10))
    gl.ylocator = mticker.FixedLocator(range(30,71,10))
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = dict(color = 'k', rotation = 'horizontal', fontsize  = rcParams_area['xtick.labelsize'])
    gl.ylabel_style = dict(color = 'k', rotation = 'horizontal', fontsize  = rcParams_area['ytick.labelsize'])
    ax.add_feature(feature.COASTLINE, zorder = 10, color = 'k')

#     gl.xlocator = mticker.FixedLocator(range(-30,71,10))
#     gl.ylocator = mticker.FixedLocator(range(30,80,10))
    
    # set axes title
    # ax.text(0.05, 0.945, '   ', ha='center', va='center', backgroundcolor = [1,1,1,0.85], transform=ax.transAxes, size = rcParams_area['axes.labelsize'], weight = 'normal', zorder = 15)
    if axes_title :
        ax.text(0.1, 0.94, ' ' + axes_title, ha='center', va='center', transform=ax.transAxes, size = rcParams_area['axes.labelsize'], weight = 'bold', zorder = 15)
    
    
    if set_extent :
        ax.set_extent([-15,51,29,71])
        
    return mapple_object, ct if contour else mapple_object


def add_cbar(fig, ax, mapple_object, divide = 'horizontal', label = False, percentage = "3.5%", pad = 0.075, labelpad = 10, aspect=25, extend = 'neither', ticks = None, orientation='vertical', **rest) :
    # this function adds a colobar to the ax object given to it by dividing the 
    # create colorbar axes objject with divider of axes object
    divider = make_axes_locatable(ax)
    if divide == 'horizontal':
        ax_cb = divider.new_horizontal(size= percentage, pad= pad, axes_class= plt.Axes)
    elif divide == 'vertical':
        ax_cb = divider.new_vertical(size= percentage, pad= pad, axes_class= plt.Axes)
    else :
        ax_cb = ax
    fig.add_axes(ax_cb)
    
    cbar = plt.colorbar(mapple_object, cax= ax_cb, extend = extend, orientation = orientation, aspect=aspect, ticks = ticks)
    if label and orientation == 'vertical':
        cbar.ax.get_yaxis().labelpad = labelpad
        cbar.ax.set_ylabel(label, rotation=90)
    elif label and orientation == 'horizontal':
        cbar.ax.get_xaxis().labelpad = labelpad
        cbar.ax.set_xlabel(label, rotation=0)
    return cbar



# ===============
# colormap creation by Kerry Halupka
def hex_to_rgb(value):
    '''
    Converts hex to rgb colours
    value: string of 6 characters representing a hex colour.
    Returns: list length 3 of RGB values'''
    value = value.strip("#") # removes hash symbol if present
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))


def rgb_to_dec(value):
    '''
    Converts rgb to decimal colours (i.e. divides each value by 256)
    value: list (length 3) of RGB values
    Returns: list (length 3) of decimal values'''
    return [v/256 for v in value]

def get_continuous_cmap(hex_list, float_list=None):
    ''' creates and returns a color map that can be used in heat map figures.
        If float_list is not provided, colour map graduates linearly between each color in hex_list.
        If float_list is provided, each color in hex_list is mapped to the respective location in float_list. 
        
        Parameters
        ----------
        hex_list: list of hex code strings
        float_list: list of floats between 0 and 1, same length as hex_list. Must start with 0 and end with 1.
        
        Returns
        ----------
        colour map'''
    rgb_list = [rgb_to_dec(hex_to_rgb(i)) for i in hex_list]
    if float_list:
        pass
    else:
        float_list = list(np.linspace(0,1,len(rgb_list)))
        
    cdict = dict()
    for num, col in enumerate(['red', 'green', 'blue']):
        col_list = [[float_list[i], rgb_list[i][num], rgb_list[i][num]] for i in range(len(float_list))]
        cdict[col] = col_list
    cmp = mplcolors.LinearSegmentedColormap('my_cmp', segmentdata=cdict, N=256)
    return cmp


# COLORMAPS SEFLMADE
hex_list = ['#ffffff', '#f4f4f9', '#f9dcc4', '#fcbf49', '#f77f00', '#d62828', '#6a040f', '#3c1642']
thermal = get_continuous_cmap(hex_list, float_list=[0, 0.02, 0.1, 0.3, 0.5, 0.75, 0.95, 1])

hex_list = ['#48bfe3', '#f3ffbd', '#ffb627', '#ff9505', '#d00000', '#6a040f', '#370617', '#3c1642']
thermal_high_extend = get_continuous_cmap(hex_list, float_list=[0, 0.1, 0.3, 0.5, 0.65, 0.85, 0.95, 1])

hex_list = ['#caf0f8', '#eaeaea', '#f9dcc4', '#fcbf49', '#f77f00', '#d62828', '#87031d','#720026', '#003049']
thermal_low_extend = get_continuous_cmap(hex_list, float_list=[0, 0.1, 0.25,0.4, 0.5, 0.75, 0.85, 0.95, 1])

hex_list = ['#ffffff', '#f4f4f9', '#f9dcc4', '#fcbf49', '#f77f00', '#d62828', '#6a040f', '#003049']
hwmid = get_continuous_cmap(hex_list, float_list=[0, 0.1, 0.3, 0.5, 0.75, 0.85, 0.5, 1])


def get_AR6_coords(url = 'https://raw.githubusercontent.com/IPCC-WG1/Atlas/master/reference-regions/IPCC-WGI-reference-regions-v4_coordinates.csv') :
    
    points = ['p0','p1','p2','p3','p4','p5','p6','p7','p8','p9','p10','p11','p12','p13']
    import pandas as pd
    df = pd.read_csv(url, nrows = 32, header = None,
                     names = ['continent','land/ocean','long_name','short_name'] + points)

    regions = [dict(name = 'NEU', region = 16), 
                   dict(name = 'EEU', region = 18), 
                   dict(name = 'MED', region = 19),
                   dict(name = 'WCE', region = 17), ]

    for region in regions : 
        region['latitude'] = []
        region['longitude'] = []
        region_row = df.loc[df['short_name'] == region['name']]
    #     print(region_row)
        for p in points[0:-1] :
    #         print(p)
    #         print(region_row[p].values[0])
            if isinstance(region_row[p].values[0], str):
                test = region_row[p].str.split('|').values[0]
                longitude = float(test[0])
                latitude = float(test[1])
                region['latitude'].append(latitude)
                region['longitude'].append(longitude)
        region['latitude'].append(region['latitude'][0])
        region['longitude'].append(region['longitude'][0])

    regions[3]['name'] = 'CEU'
    return regions

def plot_all_regions(ax, regions = get_AR6_coords(), colors = ['tab:red','tab:green','tab:blue','tab:orange'], 
                     linestyle = ['--','--','--','--'], 
                     linewidth = [0.75,0.75,0.75,0.75],
                     underline = False) :
    if not isinstance(linewidth,list) :
        linewidth = [linewidth, linewidth, linewidth, linewidth]
    if not isinstance(linestyle,list) :
        linestyle = [linestyle, linestyle, linestyle, linestyle]
    if not isinstance(colors,list) :
        colors = [colors, colors, colors, colors]
        
        
    for idx in range(len(regions)) :
        regions[idx]['plotkwargs'] = dict(color = colors[idx] , 
                                          linestyle = linestyle[idx], 
                                          linewidth = linewidth[idx])
    lines = []
    for region in regions :
        if underline :
            ax.plot(region['longitude'],region['latitude'], linestyle = '-', color = 'w', linewidth = region['plotkwargs']['linewidth'] + underline,transform = ccrs.PlateCarree(), zorder = 20)
        line = ax.plot(region['longitude'],region['latitude'], transform = ccrs.PlateCarree(), **region['plotkwargs'], zorder = 20)
        lines.append(line)

# REST
# def plot_30y_area(data = None, mask_res = None, 
#                 cmap= "Reds", cmap_days = "Reds", save_name = False, 
#                 mask_dict = {}):

#     '''
#     this function creates a figure with area_plot() function fro all heatwave days
#     a) mean temperature 
#     b) mean swbgt 
#     c) days which are part of a heatwave 

#     INPUT:
#     1. data         - xarray DataSet with longitude and latitude as dimensions
#                     - !! need to include variables t2m and swbgt and dimension time, latitude, longitude
#                     - !! should be daily data
#     2. mask_res     - mask which is True wherever there is a heatwave occuring
#     3. cmap         - cmap for temperature and swbgt
#     4. cmap_days    - cmap for days of heatwave plot
#     5. save_name    - if False no save else save png and svg
#     6. mask_dict    - "threshold" : what was the threshold for this mask,
#                     - "duration" : what was the duration
#                     - "quantile" : what was the quantile provided?
#     OUTPUT: 
#     - figure with the plots
#     '''

#     # create year range
#     year_start = data.time[0].dt.strftime("%Y").values
#     year_end = year = data.time[-1].dt.strftime("%Y").values
#     years = 31
#     # get arguments for the mask (by default : NONE GIVEN)
#     duration = mask_dict.setdefault("duration" , np.nan)
#     threshold = mask_dict.setdefault("threshold" , np.nan)
#     quantile = mask_dict.setdefault("quantile" , np.nan)
#     #create figure title :
#     figure_title = year_start + " - " + year_end + "\nheatwave similar to DWD definition \n {:.0f} days duration,\n {:.0f}°C daily max. temperature threshold,\n {:.0f}th percentile of temperature".format(duration, threshold, quantile*100)
    
#     # ========== PLOT ==============
#     # create the figure and the gridspecs
#     fig = plt.figure(figsize = (9,27))
#     gs = gridspec.GridSpec(3, 1)

#     # plot the mean temperature for the heatwaves (where the mask_res is True)
#     plot_data = data.t2m.where(mask_res, other = np.nan).mean("time")
#     ax = fig.add_subplot(gs[0,0],\
#             projection=ccrs.NearsidePerspective(
#                             central_latitude=50.72,
#                             central_longitude=10.53,
#                             satellite_height=10000000.0))

#     ax.coastlines(resolution='110m')
#     ax.gridlines()
#     pm, ct = area_plot(data= plot_data, ax = ax, \
#         levels = np.arange(28, 42, 0.5), \
#         contour_levels= np.arange(25,40,2), cmap = cmap,
#         landcolor = [0.4,0.4,0.4])
    
#     add_cbar(fig = fig, ax=ax, pm = pm, label = "Temperature in °C")
#     ax.set_title("a) mean temperature during heat wave")

#     # plot the mean temperature for the heatwaves (where the mask_res is True)
#     plot_data = data.swbgt.where(mask_res, other = np.nan).mean("time")
    
#     ax1 = fig.add_subplot(gs[1,0],\
#             projection=ccrs.NearsidePerspective(
#                             central_latitude=50.72,
#                             central_longitude=10.53,
#                             satellite_height=10000000.0))

#     ax1.coastlines(resolution='110m')
#     ax1.gridlines()
#     pm, ct = area_plot(data= plot_data, ax = ax1, \
#         levels = np.arange(22, 35, 0.5), \
#         contour_levels= np.arange(15,35,2), \
#         cmap = cmap,
#         landcolor = [0.4,0.4,0.4])

#     cbar = add_cbar(fig = fig, ax=ax1, pm = pm, label = "sWBGT")
#     ax1.set_title("b) mean sWBGT during heat wave")

#     # plot the days which are part of a heatwave
#     plot_data = np.log(mask_res.sum("time")) # sum up all True values give the value for days 
#     plot_data = plot_data.where(plot_data != 0)

#     ax2 = fig.add_subplot(gs[2,0],\
#             projection=ccrs.NearsidePerspective(
#                             central_latitude=50.72,
#                             central_longitude=10.53,
#                             satellite_height=10000000.0))

#     ax2.coastlines(resolution='110m')
#     ax2.gridlines()

#     pm, ct = area_plot(data= plot_data, ax = ax2, \
#         levels =  np.arange(1, 10, 1), \
#         contour_levels= np.arange(1,10,2), \
#         cmap = cmap_days)
#     ax2.set_title("heatwave days for the period - log scale")
#     cbar = add_cbar(fig = fig, ax = ax2, pm = pm, label = "log(days)")

#     # create figure title
#     fig.suptitle(figure_title, fontsize= 15)

#     # ============ SAVE ===============
#     # save data if wanted
#     year = data.time[0].dt.strftime("%Y-%m-%d").values
#     if save_name:
#         fig.savefig(plot_dict["savefolder"] + save_name + '.png')
#         fig.savefig(plot_dict["savefolder"] + save_name + '.svg')
#     return fig

# def plot_city_data(data = None, data_max = None,
#         mask_res = None, city_name = "None", 
#         timeslice= slice(np.nan,np.nan), groupby= "time.year",
#         color_bar = "tomato" , color_scatter = "tomato" , color_hist = "tomato",
#         mask_dict = {}):
    
#     duration = mask_dict.setdefault("duration" , np.nan)
#     threshold = mask_dict.setdefault("threshold" , np.nan)
#     quantile = mask_dict.setdefault("quantile" , np.nan)

#     try: 
#             data = data.mean("longitude").mean("latitude")
#             data_max = data_max.mean("longitude").mean("latitude")
#     except: 
#             pass

#     fig = plt.figure(figsize = (15,15))
#     gs = gridspec.GridSpec(3, 2)

#     fig.suptitle(city_name + ": " + timeslice.start + " - " + timeslice.stop + "\nheatwave similar to DWD definition \n {:.0f} days duration, {:.0f}°C daily max. temperature threshold, {:.0f}th percentile of temperature".format(duration, threshold, quantile*100),\
#             fontsize= 15)
#     # plot temperature
#     ax1 = fig.add_subplot(gs[0,:])
#     data.t2m.sel(time = timeslice)\
#             .plot(ax = ax1, color = "k", linewidth = "0.3",linestyle = ":", alpha= 0.75,\
#             label= "a) temperature - 6 hourly")
#     data_max.t2m.sel(time = timeslice)\
#             .plot(ax = ax1, color = "tab:blue", linewidth = "0.5",\
#             label= "b) temperature - daily max")
#     data_max.t2m.where(mask_res).sel(time = timeslice)\
#             .plot(ax = ax1, marker = 'x', color = "tab:orange",\
#             label= "c) heatwaves days")
#     ax1.set_ylim([-5,35])
#     ax1.set_title("")
#     ax1.set_ylabel("temperature in °C")
#     leg = ax1.legend(loc = "lower left")
#     ax1.grid()

#     # plot bar plot
#     ax2 = fig.add_subplot(gs[1,:], sharex = ax1)
#     y = mask_res.sel(time = timeslice).groupby(groupby).sum('time')

#     if "year" in groupby:
#         # x = y.year
#         x = data_max.time[data_max.time.dt.is_year_start].sel(time = timeslice)
#         w = (x[1]-x[0])/2
#         ax2.bar(x = x + w , height = y, width = w, align = "center")
#     elif "month" in groupby:
#         x = data_max.time[data_max.time.dt.is_month_start].sel(time = timeslice)
#         w = (x[1]-x[0])/2
#         ax2.bar(x = x + w , height = y, width = w/2, align = "center")
#     elif "week" in groupby:
#         x = x.week
#         w = (x[1]-x[0])/2
#     else:
#         x = y.time
    
    
#     # city(mask_org).sel(time = timeslice).groupby('time.month').sum('time').plot(ax = ax3,linestyle= ":", label = "xclim used 6 hourly")
#     ax2.set_ylabel("heat wave days per year")
    
#     # ax2.legend()
#     ax2.grid()
#     # ax1.set_xticklabels([])

#     ax3 = fig.add_subplot(gs[2:3,0])
#     x = data_max.t2m.where(mask_res).sel(time = timeslice)
#     y = data_max.swbgt.where(mask_res).sel(time = timeslice)
#     c = data_max.d2m.where(mask_res).sel(time = timeslice)
#     # c = data_max["time.month"].values
#     sct = ax3.scatter(x,y)
#     ax3.set_xlabel("temperature in °C")
#     ax3.set_ylabel("sWBGT")
#     ax3.set_xlim([24,38])
#     ax3.set_ylim([18,30])
#     ax3.grid()

#     ax4 = fig.add_subplot(gs[2:3,1])
#     bins = np.arange(19,28)
#     ax4.hist(x = y, bins = bins ,density = True)
#     ax4.set_xlabel("sWBGT")
#     ax4.set_ylabel("density")
#     ax4.set_ylim([0,0.35])
#     ax4.grid()

#     return fig

# def 