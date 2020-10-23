#%% HEAT WAVES IDENTIFICTION
import xarray as xr
import matplotlib.pyplot as plt
import time
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import time
import numpy as np
import cartopy.crs as ccrs
import cartopy as cart
import xclim as xc
from dask.diagnostics import ProgressBar
import os
import matplotlib.gridspec as gridspec
import Plot.plotfunctions as plotfunc
import Calculation.humiditycalculation as humidcalc
import Calculation.trend_func as trends
import Calculation.save_detrend_deseason as sdd
import Calculation.heat_wave_indetification as hwi
#%%
# plt.style.use('ggplot')
plt.rcParams['figure.figsize'] = (15,10)

output_dtype = float
timestart = "1980-01-01"
timeend = "2010-12-31"
location = "europe_all"
# rolling_window = 4*365
timeslice = slice("2003-06-15", "2003-09-15")
load_dict = {
        't2m_path' : r"I:\Bachlor_thesis\Data\t2m_era20c_1900-2010.nc",
        'd2m_path' : r"I:\Bachlor_thesis\Data\td2m_era20c_1900-2010.nc",
        'time_slice' : slice(timestart, timeend),
        "latitude_slice" : slice(53.5, 51.5),
        "longitude_slice" : slice(12.5, 14.5),
        }

save_dict = {
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
        "savefolder" : r"C:\Users\Nils Niebaum\Documents\Uni\Bachlor_thesis\Prog\Plots\Timeseries\\",
        "savename" : location + "_swbgt_detrend_timeseries" + "_{}_{}.png".format(timestart, timeend)
        }

#%%b load data
data_1900= xr.open_dataset(r"C:\Users\Nils Niebaum\Documents\Uni\Bachlor_thesis\Prog\Data\europe_1900-01-01_1930-12-31.nc")
data_2010 = xr.open_dataset(r"C:\Users\Nils Niebaum\Documents\Uni\Bachlor_thesis\Prog\Data\europe_1980-01-01_2010-12-31.nc")
st = time.time()
data_1900_max = trends.calc_daily_max(data, variable= "t2m")
print(time.time() - st)
st = time.time()
data_2010_max = trends.calc_daily_max(data, variable= "t2m")
print(time.time() - st)

# data = sdd.save_calc_swbgt(save_dict= save_dict, load_dict = load_dict)
# print('loaded')
# st = time.time()
# data_max = trends.calc_daily_max(data, variable= "t2m")
# print(time.time() - st)
#%%
max_dict = {
        "quantiles" : {"variables" : ["t2m"] , "quantiles" : np.array([0.98])} ,
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
data.t2m.attrs["units"] = "degC"
# data.t2m_deseason.attrs["units"] = "degC"
data_max.t2m.attrs["units"] = "degC"
# data_max.t2m_deseason.attrs["units"] = "degC"
#%% 
# data.t2m.sel(time = slice("2003-06", "2003-09")).plot()
# data.t2m.where(data.t2m_mask_quantiles.isel(quantile = 0)).sel(time = slice("2003-06", "2003-09")).plot()

# data_max.t2m.sel(time = slice("2003-06", "2003-09")).plot()
# data_max.t2m.where(data_max.t2m_mask_quantiles.isel(quantile = 0)).sel(time = slice("2003-06", "2003-09")).plot()

#%%
from memory_profiler import profile
threshold = 25
duration = 3

mask_lambda = lambda x : (x.t2m >= threshold) # & x.t2m_mask_quantiles.isel(quantile= 0)

mask_org = mask_lambda(data_max)
st = time.time()
mask_res = hwi.calc_heatwave_index(mask_original= mask_org, duration = duration)
print(time.time() -st)

#%%

from matplotlib.colors import BoundaryNorm
import matplotlib.ticker as mticker
# from matplotlib.ticker import MaxNLocator
import matplotlib.colors as mplcolors
import cmocean as cmo
cmap = cmo.cm.thermal

def area_plot(data = None, ax = None, levels = np.arange(20,38,1), contour_levels = np.arange(20,40,2), add_colorbar=False, cmap = "Reds") :
    #create colormap
    cmap = plt.get_cmap(cmap, len(levels))



    plot = data.plot.pcolormesh(ax = ax, cmap = cmap, \
        x = "longitude", y = "latitude",\
        yincrease = True, transform = ccrs.PlateCarree(), \
        norm = mplcolors.BoundaryNorm(levels, ncolors=len(levels) -1, clip=False), \
        )
    
    contourplot = data.plot.contour(ax= ax, levels= contour_levels, colors= 'k', linestyles = '-', linewidths = 0.2, x="longitude", y="latitude", \
        yincrease=True, transform=ccrs.PlateCarree(),\
        add_colorbar= add_colorbar)

    # plot gridlines and coastlines
    ax.coastlines(resolution = 'auto', color = 'k', )
    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,
                  linewidth=1, color='gray', alpha=0.5, linestyle=':')
    # gl.xlabels_top = False
    # gl.ylabels_right = False
    gl.xlines = True
    # gl.xlocator = mticker.FixedLocator(range(-30,50,10))
    # gl.xformatter = LONGITUDE_FORMATTER
    # gl.yformatter = LATITUDE_FORMATTER

    # gl.xlabel_style = {'size': 15, 'color': 'gray'}
    # gl.xlabel_style = {'color': 'red', 'weight': 'bold'}
    ax.add_feature(cart.feature.OCEAN,zorder=4)

    return plot

fig = plt.figure(figsize = (27,9))
gs = gridspec.GridSpec(1, 3)
ax = fig.add_subplot(gs[0,0],\
        projection=ccrs.NearsidePerspective(
                        central_latitude=50.72,
                        central_longitude=10.53,
                        satellite_height=10000000.0))

ax.coastlines(resolution='110m')
ax.gridlines()
plot_data = data_max.t2m.where(mask_res).sel(time = timeslice).mean("time")
area_plot(data= plot_data, ax = ax, levels = np.arange(int(plot_data.min()), int(plot_data.max()), 0.2), contour_levels= np.arange(15,35,5), cmap = cmap)

ax.set_title("mean of the daily maximum temperature during heat wave")
#
ax1 = fig.add_subplot(gs[0,1],\
        projection=ccrs.NearsidePerspective(
                        central_latitude=50.72,
                        central_longitude=10.53,
                        satellite_height=10000000.0))

ax1.coastlines(resolution='110m')
ax1.gridlines()
plot_data = data_max.swbgt.where(mask_res).sel(time = timeslice).mean("time")
area_plot(data= plot_data, ax = ax1, levels = np.arange(int(plot_data.min()), int(plot_data.max()), 0.2), contour_levels= np.arange(15,35,5), cmap = cmap)
ax1.set_title("mean daily sWBGT maximum during heat wave")

ax2 = fig.add_subplot(gs[0,2],\
        projection=ccrs.NearsidePerspective(
                        central_latitude=50.72,
                        central_longitude=10.53,
                        satellite_height=10000000.0))

ax2.coastlines(resolution='110m')
ax2.gridlines()
plot_data = mask_res.sel(time = timeslice).sum("time")
plot_data = plot_data.where(plot_data != 0)
area_plot(data= plot_data, ax = ax2, levels =  np.arange(1, 50, 5), contour_levels= np.arange(0,30,5), cmap = "Reds")
ax2.set_title("days of heatwave")




fig.suptitle(timeslice.start + timeslice.stop + " heatwave similar to DWD defenition \t {:.0f} days duration, {:.0f}°C threshold".format(duration, threshold))
# ax3 = fig.add_subplot(gs[0,2])
# x = data_max.t2m.where(mask_res).sel(time = timeslice)
# y = data_max.swbgt.where(mask_res).sel(time = timeslice)
# ax3.scatter(x,y, marker = '.')
plt.show()
fig.savefig("europe_own2_2003.png")
fig.savefig("europe_own2_2003.svg")

#%% plot area
def plot_hw():
        fig = plt.figure(figsize = (18,27))
        gs = gridspec.GridSpec(3, 2)
        swbgt_th = [15,10,15,10,15,15]
        window = [3,3,2,2,3,2]
        th = [3,3,3,3,2,2]
        idx = 0
        fig.suptitle('Days of heat waves from{} until {}'.format(timestart, timeend), fontsize=16)
        for ii in range(3):
                for jj in range(2) :
                        axss = fig.add_subplot(gs[ii,jj], projection=ccrs.EuroPP())
                        heatwave = xc.indices.heat_wave_index(data.swbgt_deseason.where(data.swbgt >= swbgt_th[idx]), thresh= str(th[idx]) + " degC", window= window[idx], freq= "30 YS")
                        print(idx)
                        plotfunc.area_plot(data= heatwave.isel(time = 0), ax = axss, levels = np.arange(0,200,2), contour_levels= [10,50,100,250])
                        # cb.remove()
                        axss.set_title("sWBGT >={:.0f}, days in a row >= {:.0f} days, sWBGT seasonal anomaly >= {}".format(swbgt_th[idx], window[idx], th[idx]))
                        del heatwave
                        idx += 1

        plt.tight_layout()
        fig.savefig("europe{}_{}.png".format(timestart, timeend))
        fig.savefig("europe{}_{}.pdf".format(timestart, timeend))
        plt.show()

#%% use swbgt for heat wave definition 
#fig, ax = plt.subplots(nrows= 3, ncols= 1, sharex = True, figsize= (24,10))

def swbgt_hw() :
        # plot swbgt
        data.swbgt.mean("latitude").mean("longitude").plot(ax = ax[0], color = "tab:blue", linewidth = "0.5",\
                label= "a) swbgt")
        (data.swbgt - data.swbgt_deseason).mean("latitude").mean("longitude").plot(ax = ax[0], linestyle = '--', color = "tab:red",\
                label= "b)seasonal swbgt")

        ax[0].set_title("(I) swbgt and seasonal swbgt\n(calculated from the mean of week per year)")

        # plot swbgt deseasonal
        data.swbgt_deseason.mean("latitude").mean("longitude").plot(ax = ax[1], marker = '.', markersize = 1, linestyle = '', linewidth = "0.5", color = "tab:blue",\
                label = "a) deseasonaled swbgt")

        # linreg, x = trends.calculate_linreg(data= data, variable= "swbgt_deseason", mutable= False)

        # (data.swbgt_deseason - data.swbgt_deseason_detrend).mean("latitude").mean("longitude").plot(ax = ax[1], color = "tab:red", linestyle = '--',\
        #         label = "b) linear regression of a) with slope = {:.2E}".format(linreg.slope*1e9*365*24*60*60*30) + r"$\frac{°C}{30 years}$")

        ax[1].set_title("(II) deseasonalised swbgt")

        # deseanonal and detrended data
        data.swbgt_deseason_detrend.where(data.swbgt_deseason_detrend >= 3).mean("latitude").mean("longitude").plot(ax= ax[2], marker= ".", markersize = 2, linestyle = "", color = "tab:blue" ,\
                label = "a) swbgt deseasonalised & detrended >= 3°C")
        data.swbgt_deseason_detrend.where(data.swbgt_deseason_detrend >= 3).mean("latitude").mean("longitude").where(data.swbgt >= 15).plot(ax= ax[2], marker= "x", markersize = 2, linestyle = "-" , color = "tab:orange",linewidth = "0.5",\
                label = "b) swbgt deseasonalised & detrended >= 3°C and swbgt >= 15°C")

        # linreg, x = trends.calculate_linreg(data = data.swbgt_deseason_detrend.where(data.swbgt_deseason_detrend >= 3).where(data.swbgt >= 15).mean("latitude").mean("longitude"), variable = "swbgt_deseason_detrend", mutable= False)
        # print(data.swbgt_deseason_detrend.where(data.swbgt_deseason_detrend >= 3).where(data.swbgt >= 15))
        # ax[2].plot(data.time, linreg.slope * x + linreg.intercept, color = "tab:red", linestyle = '--',\
        #         label = "linear regression of b) with slope = {:.2E}".format(linreg.slope*1e9*365*24*60*60*30) + r"$\frac{°C}{30 years}$")
        ax[2].set_title("(III) filtered data from plot (II)")


        for axs in ax :
                axs.grid()
                axs.legend()
                axs.set
                axs.set_ylabel("°C")

        fig.savefig(plot_dict["savefolder"] + plot_dict["savename"]+ "new.png")
        plt.show()

#%% use max temp for heat wave definition 

def plot_city_data(data = None, mask_data = None, data_max = None, mask_data_max = None, city_name = "None"):
        
        try: 
                data = data.mean("longitude").mean("latitude")
        except: 
                pass
        
        fig, ax = plt.subplots(nrows= 5, ncols= 1, sharex = False, figsize= (15,15))

        # plot t2m 6h
        data.t2m.sel(time = timeslice)\
                .plot(ax = ax[0], color = "tab:blue", linewidth = "0.5",\
                label= "a) t2m - 6 hourly")
        (data.t2m - data.t2m_deseason).sel(time = timeslice)\
                .plot(ax = ax[0], color = "tab:blue", linewidth = "1",linestyle = "--",\
                label= "b) t2m - 6 hourly - seasonal mean")
        data.t2m.where(data.t2m_mask_quantiles.isel(quantile= 0)).sel(time = timeslice)\
                .plot(ax = ax[0], marker = 'x', color = "tab:red",\
                label= "c) t2m - 6 hourly - 98th percentile")

        ax[0].set_title("(I) 2m temperature (t2m) and its 98th percentile")


        # plot t2m 6h
        data_max.t2m.sel(time = timeslice)\
                .plot(ax = ax[1], color = "tab:blue", linewidth = "0.5",\
                label= "a) t2m - daily max")
        (data_max.t2m -data_max.t2m_deseason).sel(time = timeslice)\
                .plot(ax = ax[1], color = "tab:blue", linewidth = "1",linestyle = "--",\
                label= "b) t2m - daily max - seasonal mean")
        data_max.t2m.where(data_max.t2m_mask_quantiles.isel(quantile= 0)).sel(time = timeslice)\
                .plot(ax = ax[1], marker = 'x', color = "tab:red",\
                label= "c) t2m - daily max - 98th percentile")

        ax[1].set_title("(I) daily maximum 2m temperature (t2m) and its 98th percentile")


        # DWD heat waves
        data.t2m.sel(time = timeslice)\
                .plot(ax = ax[2], color = "tab:blue", linewidth = "0.5",\
                label= "a) t2m - 6 hourly")
        data.t2m.where(city(mask_org)).sel(time = timeslice)\
                .plot(ax = ax[2], marker = 'x', color = "tab:orange",\
                label= "b) t2m - 6 hourly - 98th percentile and above {}°C".format(threshold))

        data_max.t2m.sel(time = timeslice)\
                .plot(ax = ax[2], color = "tab:purple", linewidth = "0.5",\
                label= "a) t2m - daily max")
        data_max.t2m.where(city(mask_res)).sel(time = timeslice)\
                .plot(ax = ax[2], marker = 'x', color = "tab:green",\
                label= "b) t2m - daily max - 98th percentile and above {}°C".format(threshold))

        data_max

        ax[2].set_title("(II) 98th percentile of 2m temperature (t2m) and above ")
        ax[2].set_ylim([23,33])
        ax[2].set_ylabel("°C")



        # heatwave = xc.indices.heat_wave_index(data.t2m.where(data.t2m_mask_quantiles.isel(quantile= 0)), thresh= str(threshold) + " degC", window= 2, freq= "MS").sel(time = timeslice)
        # heatwave_max = xc.indices.heat_wave_index(data_max.t2m.where(data_max.t2m_mask_quantiles.isel(quantile= 0)), thresh= str(threshold) + " degC", window= 2, freq= "MS").sel(time = timeslice)

        # heatwave_max.plot(ax = ax[3], label = "xclim used daily max")
        # heatwave.plot(ax = ax[3],linestyle= ":", label = "xclim used 6 hourly")
        
        city(mask_res).sel(time = timeslice).groupby('time.month').sum('time').plot(ax = ax[3], label = "xclim used daily max")
        city(mask_org).sel(time = timeslice).groupby('time.month').sum('time').plot(ax = ax[3],linestyle= ":", label = "xclim used 6 hourly")
        
        ax[3].set_ylabel("heat wave days")

        x = data_max.t2m.where(city(mask_res)) #.sel(time = timeslice)
        y = data_max.swbgt.where(city(mask_res)) #.sel(time = timeslice)
        c = data_max["time.month"].values
        sct = ax[4].scatter(x,y, c = c, cmap = "hsv")
        ax[4].set_xlabel("t2m in °C")
        ax[4].set_ylabel("sWBGT unitless")
        plt.colorbar(sct)

        for axs in ax[0:4] :
                axs.grid()
                axs.legend()

        # data_max.t2m.where(Paris(mask_org)).sel(time = timeslice).plot(ax = ax[3], marker = "+")
        # data_max.t2m.where(Paris(mask_res)).sel(time = timeslice).plot(ax = ax[3], marker = "x")

        fig.savefig(plot_dict["savefolder"] + plot_dict["savename"]+ "hw_{}_{}_{}.png".format(city_name,timeslice.start, timeslice.stop))
        plt.show()


mask_lambda = lambda x : (x.t2m >= threshold) # & x.t2m_mask_quantiles.isel(quantile= 0)

mask_org = mask_lambda(data_max)
st = time.time()
mask_res = hwi.calc_heatwave_index(mask_original= mask_org, duration = duration)
print(time.time() -st)
mask_org = mask_lambda(data)

city = lambda ds : ds.sel(method="nearest", latitude = 48.86, longitude = -2.35) 
plot_city_data(data = city(data), data_max = city(data_max), city_name = "Paris") 

city = lambda ds : ds.sel(method="nearest", latitude = 52.5, longitude = 13.4) 
plot_city_data(data = city(data), data_max = city(data_max), city_name = "Berlin")

city = lambda ds : ds.sel(method="nearest", latitude = 55.76, longitude = 37.62) 
plot_city_data(data = city(data), data_max = city(data_max), city_name = "Moscow") 

# #%%
# # def calc_hw_idx(data, duration):
# duration -= 1

# mask_new = (data_max.t2m_mask_quantiles.isel(quantile= 0) & (data_max.t2m >= threshold))
# mask_list = []
# max_length = int(data_max.time.size/duration) * duration
# for i in range(duration):
#         selection = np.arange(i,max_length, duration)
#         mask_list.append(mask_org.isel(time = selection).values)
#         mask_res = None
# for mask in mask_list:
#         if mask_res is None:
#                 mask_res = mask
#         mask_res = mask & mask_res
# mask_temp = mask_new == -99
# for i in range(duration):
#         selection = np.arange(i,max_length, duration)
#         mask_temp[selection,:,:] = mask_res


# data_max.t2m.where(mask_org).sel(time = timeslice).plot(ax = ax[3], marker = "X")
# data_max.t2m.where(mask_temp).sel(time = timeslice).plot(ax = ax[3], marker = "o")
#%%

city = lambda ds : ds.sel(method="nearest", latitude = 48.86, longitude = -2.35)
city = lambda ds : ds.sel(method="nearest", latitude = 52.5, longitude = 13.4) 
plt.hist(data_max.t2m.where(mask_res).values.flatten())