#%% HEAT WAVES IDENTIFICTION
import xarray as xr
import matplotlib.pyplot as plt
import time
import numpy as np
import xarray as xr
import time
import numpy as np
import cartopy.crs as ccrs
import cartopy as cart
import xclim as xc
import os
from memory_profiler import profile

from matplotlib.colors import BoundaryNorm
import matplotlib.ticker as mticker
import matplotlib.colors as mplcolors
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cmocean as cmo

import matplotlib.gridspec as gridspec
# own functions
import Plot.plotfunctions as plotfunc
import Calculation.humiditycalculation as humidcalc
import Calculation.trend_func as trends
import Calculation.save_detrend_deseason as sdd
import Calculation.heat_wave_indetification as hwi
#%%
# plt.style.use('ggplot')
plt.rcParams['figure.figsize'] = (15,10)

# output_dtype = float
# timestart = "1980-01-01"
# timeend = "2010-12-31"
# location = "europe_all"
# # rolling_window = 4*365
# timeslice = slice("2003-06-15", "2003-09-15")
# load_dict = {
#         't2m_path' : r"I:\Bachlor_thesis\Data\t2m_era20c_1900-2010.nc",
#         'd2m_path' : r"I:\Bachlor_thesis\Data\td2m_era20c_1900-2010.nc",
#         'time_slice' : slice(timestart, timeend),
#         "latitude_slice" : slice(53.5, 51.5),
#         "longitude_slice" : slice(12.5, 14.5),
#         }

# save_dict = {
        "savefolder" : r"C:\Users\Nils Niebaum\Documents\Uni\Bachlor_thesis\Prog\Data\\",
        "location" : location,
        "savename" : location + "_{}_{}.nc".format(timestart, timeend), # old: "savename" : r"swbgt" + "_{}_{}_{}.nc".format(location, timestart, timeend),
        "rechunked" : False,
        "quantiles" : False, #{"variables" : ["t2m"] , "quantiles" : np.array([0.98])} ,
        "deseason" : False, #{"variables" : ["t2m"] , "groupby" : "week"}, #"week",
        "detrend" : False, #{"variables" : ["t2m_deseason"]}, #False,
        "daily" : False #  please check this later !!! in the extern function
        }
plot_dict = {
        "doplot" : True,
        "levels" : np.arange(18, 35, 1),
        "savefolder" : r"C:\Users\Nils Niebaum\Documents\Uni\Bachlor_thesis\Prog\Plots\\"
        }

#%%b load data
data_1900= xr.open_dataset(r"C:\Users\Nils Niebaum\Documents\Uni\Bachlor_thesis\Prog\Data\europe_all_1900-01-01_1930-12-31.nc")
data_2010 = xr.open_dataset(r"C:\Users\Nils Niebaum\Documents\Uni\Bachlor_thesis\Prog\Data\europe_all_1980-01-01_2010-12-31.nc")
st = time.time()
data_1900_max = trends.calc_daily_max(data_1900, variable= "t2m")
print(time.time() - st)
st = time.time()
data_2010_max = trends.calc_daily_max(data_2010, variable= "t2m")
print(time.time() - st)

# data = sdd.save_calc_swbgt(save_dict= save_dict, load_dict = load_dict)
# print('loaded')
# st = time.time()
# data_max = trends.calc_daily_max(data, variable= "t2m")
# print(time.time() - st)
#%%
max_dict = {
        "quantiles" : False, #{"variables" : ["t2m"] , "quantiles" : np.array([0.98])} ,
        "deseason" : False, #{"variables" : ["t2m"] , "groupby" : "dayofyear"},#, "rolling" : 14}, #"week",
        "detrend" : False, #{"variables" : ["t2m_deseason"]}, #False,
        "daily" : False #  please check this later !!! in the extern function
        }
try :
        data_max = data_max.drop_vars("t2m_deseason")
        data = data.drop_vars("t2m_deseason")
except : pass
sdd.add_missing_variables(data= data_1900_max, save_dict = max_dict, load_dict = load_dict)
sdd.add_missing_variables(data= data_2010_max, save_dict = max_dict, load_dict = load_dict)
# sdd.add_missing_variables(data= data, save_dict = max_dict, load_dict = load_dict)
data_1900_max.t2m.attrs["units"] = "degC"
# data.t2m_deseason.attrs["units"] = "degC"
data_2010_max.t2m.attrs["units"] = "degC"
# data_max.t2m_deseason.attrs["units"] = "degC"


#%% calculate the mask 
# NOTE: for the quantile mask, it has to be above the general 98th percentile!!

threshold = 28
duration = 3
quantile = 0.98

data_1900_max["t2m_mask_quantiles"] = data_1900_max["t2m"] >= data_1900["t2m"].quantile(quantile)
data_2010_max["t2m_mask_quantiles"] = data_2010_max["t2m"] >= data_2010["t2m"].quantile(quantile)

mask_lambda = lambda x : (x.t2m >= threshold) & x.t2m_mask_quantiles

st = time.time()
mask_1900 = hwi.calc_heatwave_index(mask_original= mask_lambda(data_1900_max), duration = duration)
print(time.time() -st)

st = time.time()
mask_2010 = hwi.calc_heatwave_index(mask_original= mask_lambda(data_2010_max), duration = duration)
print(time.time() -st)

#%%

def area_plot(data = None, ax = None, 
        levels = np.arange(20,38,1), contour_levels = np.arange(20,40,2),
        add_colorbar=False, cmap = "Reds",
        landcolor = "lightgrey", oceancolor = "lightskyblue",
        zorder_ocean = 2, zorder_land = 2) :
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
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,\
                linewidth=1, color='w', alpha=0.5, linestyle='-',\
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

def add_cbar(fig, ax, pm, label = False) :
    divider = make_axes_locatable(ax)
    ax_cb = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
    fig.add_axes(ax_cb)
    cbar = plt.colorbar(pm, cax=ax_cb, orientation='vertical')
    if label :
        cbar.ax.get_yaxis().labelpad = 15
        cbar.ax.set_ylabel(label, rotation=270)
    return cbar

def plot_30y_area(data = None, mask_res = None, years = 31, cmap= "Reds", save_name = False):
    fig = plt.figure(figsize = (9,27))
    gs = gridspec.GridSpec(3, 1)
    ax = fig.add_subplot(gs[0,0],\
            projection=ccrs.NearsidePerspective(
                            central_latitude=50.72,
                            central_longitude=10.53,
                            satellite_height=10000000.0))

    ax.coastlines(resolution='110m')
    ax.gridlines()
    plot_data = data.t2m.where(mask_res, other = np.nan).max("time")
    pm, ct = area_plot(data= plot_data, ax = ax, \
        levels = np.arange(28, 42, 0.5), \
        contour_levels= np.arange(25,40,2), cmap = cmap,
        landcolor = [0.4,0.4,0.4])
    
    add_cbar(fig = fig, ax=ax, pm = pm, label = "Temperature in °C")
    # cax = divider.new_vertical(size="5%", pad=0.7, pack_start=True)
    # fig.add_axes(cax)
    # fig.colorbar(im, cax=cax, orientation="horizontal")

    ax.set_title("a) maximum temperature during heat wave")
    #
    ax1 = fig.add_subplot(gs[1,0],\
            projection=ccrs.NearsidePerspective(
                            central_latitude=50.72,
                            central_longitude=10.53,
                            satellite_height=10000000.0))

    ax1.coastlines(resolution='110m')
    ax1.gridlines()
    plot_data = data.swbgt.where(mask_res, other = np.nan).max("time")
    pm, ct = area_plot(data= plot_data, ax = ax1, \
        levels = np.arange(22, 35, 0.5), \
        contour_levels= np.arange(15,35,2), \
        cmap = cmap,
        landcolor = [0.4,0.4,0.4])

    cbar = add_cbar(fig = fig, ax=ax1, pm = pm, label = "sWBGT")
    ax1.set_title("b) maximum sWBGT during heat wave")

    ax2 = fig.add_subplot(gs[2,0],\
            projection=ccrs.NearsidePerspective(
                            central_latitude=50.72,
                            central_longitude=10.53,
                            satellite_height=10000000.0))

    ax2.coastlines(resolution='110m')
    ax2.gridlines()
    plot_data = np.log(mask_res.sum("time"))
    plot_data = plot_data.where(plot_data != 0)
    pm, ct = area_plot(data= plot_data, ax = ax2, \
        levels =  np.arange(1, 10, 1), \
        contour_levels= np.arange(1,10,2), \
        cmap = "Reds")
    ax2.set_title("heatwave days for the period - log scale")
    cbar = add_cbar(fig = fig, ax = ax2, pm = pm, label = "log(days)")

    year_start = data.time[0].dt.strftime("%Y").values
    year_end = year = data.time[-1].dt.strftime("%Y").values

    fig.suptitle(year_start + " - " + year_end + "\nheatwave similar to DWD definition \n {:.0f} days duration,\n {:.0f}°C daily max. temperature threshold,\n {:.0f}th percentile of temperature".format(duration, threshold, quantile*100),\
            fontsize= 15)
    # ax3 = fig.add_subplot(gs[0,2])
    # x = data_max.t2m.where(mask_res).sel(time = timeslice)
    # y = data_max.swbgt.where(mask_res).sel(time = timeslice)
    # ax3.scatter(x,y, marker = '.')
    plt.show()
    year = data.time[0].dt.strftime("%Y-%m-%d").values
    if save_name:
        fig.savefig(plot_dict["savefolder"] + save_name + '.png')
        fig.savefig(plot_dict["savefolder"] + save_name + '.svg')
    return fig

#%%
plot_30y_area(data = data_2010_max, mask_res = mask_2010, cmap = cmo.cm.thermal, save_name="europe_1980-2010_dwd_max")
plot_30y_area(data = data_1900_max, mask_res = mask_1900, cmap = cmo.cm.thermal, save_name="europe_1900-1930_dwd_max")
#%% plot the difference

cmap = "cmo.balance"

fig = plt.figure(figsize = (9,27))
gs = gridspec.GridSpec(3, 1)

#plot t2m
ax = fig.add_subplot(gs[0,0],\
        projection=ccrs.NearsidePerspective(
                        central_latitude=50.72,
                        central_longitude=10.53,
                        satellite_height=10000000.0))

ax.coastlines(resolution='110m')
ax.gridlines()

plot_lambda = lambda x, mask_res: x.t2m.where(mask_res, other = np.nan).mean("time")
plot_data = plot_lambda(data_2010_max, mask_2010) - plot_lambda(data_1900_max, mask_1900) 
pm, cf = area_plot(data= plot_data, ax = ax,
        levels = np.arange(-2, 2, 0.1),
        contour_levels= np.arange(-2,2,0.5), cmap = cmap)

add_cbar(fig = fig, ax=ax, pm = pm, label = "Temperature in °C")
ax.set_title("b) mean daily maximum temperature during heat wave - total difference")

# plot sWBGT
ax1 = fig.add_subplot(gs[1,0],\
        projection=ccrs.NearsidePerspective(
                        central_latitude=50.72,
                        central_longitude=10.53,
                        satellite_height=10000000.0))

ax1.coastlines(resolution='110m')
ax1.gridlines()

plot_lambda = lambda x, mask_res: x.swbgt.where(mask_res, other = np.nan).mean("time")
plot_data = plot_lambda(data_2010_max, mask_2010) - plot_lambda(data_1900_max, mask_1900) 

pm, ct = area_plot(data= plot_data, ax = ax1, levels = np.arange(-2, 2, 0.1), contour_levels= np.arange(-2,2,0.5), cmap = cmap)
ax1.set_title("b) mean daily maximum sWBGT during heat wave - total difference")

add_cbar(fig = fig, ax=ax1, pm = pm)

ax2 = fig.add_subplot(gs[2,0],\
        projection=ccrs.NearsidePerspective(
                        central_latitude=50.72,
                        central_longitude=10.53,
                        satellite_height=10000000.0))

ax2.coastlines(resolution='110m')
ax2.gridlines()
plot_lambda = lambda x, mask_res: mask_res.sum("time")
plot_mask_2010 = plot_lambda(data_2010_max, mask_2010)
plot_mask_1900 = plot_lambda(data_1900_max, mask_1900)
plot_data = plot_mask_2010 - plot_mask_1900

plot_data = plot_data.where((plot_mask_1900 != 0) + (plot_mask_2010 != 0))/31
pm, ct = area_plot(data= plot_data, ax = ax2, \
        levels =  np.arange(-20,21,1), \
        contour_levels= np.arange(-30,30,5), cmap = cmap)
add_cbar(fig = fig, ax=ax2, pm = pm)
ax2.set_title("heatwave days - change per year")
fig.suptitle("Change from 1900-1930 until 1980-2010" "\nheatwave similar to DWD defenition: \n {:.0f} days duration,\n {:.0f}°C daily max. temperature threshold,\n {:.0f}th percentile of temperature".format(duration, threshold, quantile*100),
            fontsize= 15)
# ax3 = fig.add_subplot(gs[0,2])
# x = data_max.t2m.where(mask_res).sel(time = timeslice)
# y = data_max.swbgt.where(mask_res).sel(time = timeslice)
# ax3.scatter(x,y, marker = '.')
plt.show()

fig.savefig(plot_dict["savefolder"] + "europe_change" + '.png')
fig.savefig(plot_dict["savefolder"] + "europe_change" + '.svg')