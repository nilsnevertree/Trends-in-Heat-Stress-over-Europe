#%% PERCEBTILE OF SWBGT
import xarray as xr
import matplotlib.pyplot as plt
import time
import numpy as np
# from scipy import stats
# import matplotlib.gridspec as gridspec
import Scripts.Plot.plotfunctions as plotfunc
import Scripts.Calculation.humiditycalculation as humidcalc
import Scripts.Calculation.trend_func as trends

# plt.style.use('ggplot')
plt.rcParams['figure.figsize'] = (24,10)

output_dtype = float
timestart = "1900-01-01"
timeend = "1929-12-31"
city = "test"
rolling_window = 4*365
latitude_slice = slice(55.0, 47.0)
longitude_slice = slice(5.0, 14.5)
# quantiles = np.array([0.95,0.98,0.99])
plot_dict = {"doplot" : True,
                "levels" : np.arange(18, 35, 1),
                "savefolder" : r"C:\Users\Nils Niebaum\Documents\Uni\Bachlor_thesis\Prog\Plots\Timeseries",
                "savename" : r"\swbgt_detrend_timeseries" + "_{}_{}_{}.png".format(city, timestart, timeend)
        }
save_dict = {"dosave" : True,
                "savefolder" : r"C:\Users\Nils Niebaum\Documents\Uni\Bachlor_thesis\Prog\Data",
                "savename" : r"\swbgt" + "_{}_{}_{}.nc".format(city, timestart, timeend),
                "savename_st" : r"\swbgt" + "_{}_st_{}_{}.nc".format(city, timestart, timeend)  
        }
#%% load data 
folderpath = r"I:\Bachlor_thesis\Data"
timeslice = slice(timestart, timeend)



#%% =============== LOAD DATA ================== 
found_st = False
found = False
# try to open deseason and detrend dataset
try : # if file exists load the data
    data = xr.open_dataset(save_dict["savefolder"] + save_dict["savename_st"])
    print("=======\nload deseason amd detrend data done")
    found_st = True
    found = True
except :
    pass

# if not found calc and save data 
while not found:
    try : # if file exists load the data
            data = xr.open_dataset(save_dict["savefolder"] + save_dict["savename"])
            print("=======\nload data done")
            try :
                data.t2m[0]
                print("t2m found" )
                data.d2m[0]
                print("d2m found" )
                data.swbgt[0]
                print("swbgt found" )

                found = True
            except :
                print("file not complete\n => calc new")
                found = False
    except: 
        pass
    
    if not found :
        #%% ==================== CALCULATION OF SWBGT ================
            print("=======\ncalculate and save data")
            folderpath = r"I:\Bachlor_thesis\Data"

            st = time.time()
            t2m = xr.open_dataset(folderpath + r"\t2m_era20c_1900-2010.nc")#, chunks={'time': 12*4*365})
            print("load temp took {:.2f}s".format(time.time()-st))
            st = time.time()
            d2m = xr.open_dataset(folderpath + r"\td2m_era20c_1900-2010.nc")#, chunks={'time': 12*4*365})
            print("load dew took {:.2f}s".format(time.time()-st))
            st = time.time()
            all_data = xr.merge([t2m, d2m])
            print("-------\ninput data\n",all_data)

            org_data = all_data.sel(time = timeslice, latitude= latitude_slice, longitude= longitude_slice)
            print("-------\nsliced data\n", org_data)
            input("-------\nsliced data ok?\n")
            # #%% ========== CHUNK DATA ==============
            # rechunk data for perfomance
            rechunk_dim = {
                    "time" : int(np.shape(org_data.t2m)[0] / 8),
                    }
            # !!!!!!!!!!!!! USE  °C !!!!!!!!!!!!!!!!!!!!
            org_data = org_data.chunk(rechunk_dim) -273.15
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            org_data.t2m.attrs = {"full name" : "temperature", "unit" : "°C"}
            org_data.d2m.attrs = {"full name" : "dew point temperature", "unit" : "°C"}
            swbgt = humidcalc.swbgt_func(org_data.t2m, org_data.d2m)
            swbgt.name = "swbgt"
            swbgt.attrs = {"place" : city, "full name" : "simple wet bulb globe temperature", "unit" : "degC"}
            # print("swbgt calc done")
            
            st = time.time()
            data = xr.merge([swbgt, org_data])

            print("swbgt compute {:.2f}s".format(time.time()- st))
            print(data)
            st = time.time()
            try :
                data.to_netcdf(save_dict["savefolder"] + save_dict["savename"])
            except :
                data.to_netcdf("test.nc")

                print('a')
                data.to_netcdf(save_dict["savefolder"] + save_dict["savename"], mode= "w")
            print("swbgt saved {:.2f}s".format(time.time()-st))
            
            # create quantile data
            # data["mask_swbgt_quantiles"] = data.swbgt >= data.swbgt.quantile(quantiles)
            # print("quanile mask done")
            # # create time_size_start 
            # data["time_since_start"] = (("time"), np.array(data.time.values - data.time[0].values, dtype=float))
            # print("time_since_start ")
            # data["t2m"] = org_data.t2m
            # data["d2m"] = org_data.d2m
            # print("temp and dew point added")

            st = time.time()
            # save all data
            data.to_netcdf(save_dict["savefolder"] + save_dict["savename"], mode = 'a')
            print("other stuff saved {:.2f}".format(time.time()-st))


print(data)
#
#%% calcuate deseason and detrend 
# if not found_st: 

# trends.deseason(data= data, variable= "swbgt", groupby= "week")
data = data.drop_vars("swbgt_deseason_detrend")
trends.detrend(data= data, variable= "swbgt_deseason")

# try :
#         data.to_netcdf(save_dict["savefolder"] + save_dict["savename_st"])
# except :
#         data.to_netcdf(save_dict["savefolder"] + save_dict["savename_st"], mode= "a")

# rechunk_dim = {
#                 "time" : int(np.shape(data.t2m)[0] / 8),
#                 }
# data = data.chunk({"time" : 1000})

# check if data is different over space and time 
fig, ax = plt.subplots(nrows=1, ncols=1)
(data.swbgt_deseason_detrend - data.swbgt_deseason).isel(longitude = 3, latitude = 3).plot(ax = ax, marker='.', markersize = 2, linestyle = "")
(data.swbgt_deseason_detrend - data.swbgt_deseason).isel(longitude = 6, latitude = 5).plot(ax = ax, marker='.', markersize = 1,linestyle = "")
plt.show()

input("ok?")
fig, ax = plt.subplots(nrows=1, ncols=1)
(data.swbgt_deseason).mean("longitude").mean("latitude").plot(ax = ax, marker='.', markersize = 2, linestyle = "", label="deseaon")
(data.swbgt_deseason_detrend).mean("longitude").mean("latitude").plot(ax = ax, marker='.', markersize = 1,linestyle = "", label = "detrend")
ax.legend()
plt.show()
input("ok?")

# calculate daily mean
print("calculate daily mean with trends.calc_daily_mean")

data = trends.calc_daily_mean(data= data)

# data =  data.resample(time="1d").mean().compute()
# data = data.rolling(time = 4).mean(keep_attrs = True)

#%% plot
import xclim as xc
data.swbgt_deseason.attrs["units"] = "degC"
data.swbgt_deseason_detrend.attrs["units"] = "degC"

fig, ax = plt.subplots(nrows= 2, ncols= 2, figsize= (15,15))

swbgt_th = [15,10,15,10]
window = [3,3,2,2]
th = "2 degC"
idx = 0
for axs in ax:
        for axss in axs :
                heatwave = xc.indices.heat_wave_index(data.swbgt_deseason.where(data.swbgt >= swbgt_th[idx]), thresh= th, window= window[idx], freq= "30 YS")
                print(idx)
                heatwave.isel(time = 0).plot(ax = axss, cmap = "Reds")
                axss.set_title("swbgt >={:.0f}, window = {:.0f} days, threshold = {}".format(swbgt_th[idx], window[idx], th))
                del heatwave
                idx += 1

fig.savefig("test.png")
plt.show()
input('ok?')


data.swbgt_deseason.mean("latitude").mean("longitude").plot
# print(heatwave)
#%%
fig, ax = plt.subplots(nrows= 3, ncols= 1, sharex = True, figsize= (24,10))


# plot swbgt
data.swbgt.mean("latitude").mean("longitude").plot(ax = ax[0], color = "tab:blue", linewidth = "0.5",\
        label= "a) swbgt")
(data.swbgt - data.swbgt_deseason).mean("latitude").mean("longitude").plot(ax = ax[0], linestyle = '--', color = "tab:red",\
        label= "b)seasonal swbgt")

ax[0].set_title("(I) swbgt and seasonal swbgt\n(calculated from the mean of week per year)")

# plot swbgt deseasonal
data.swbgt_deseason.mean("latitude").mean("longitude").plot(ax = ax[1], marker = '.', markersize = 1, linestyle = '', linewidth = "0.5", color = "tab:blue",\
        label = "a) deseasonaled swbgt")

linreg, x = trends.calculate_linreg(data= data, variable= "swbgt_deseason", mutable= False)

(data.swbgt_deseason - data.swbgt_deseason_detrend).mean("latitude").mean("longitude").plot(ax = ax[1], color = "tab:red", linestyle = '--',\
        label = "b) linear regression of a) with slope = {:.2E}".format(linreg.slope*1e9*365*24*60*60*30) + r"$\frac{°C}{30 years}$")

ax[1].set_title("(II) deseasonalised swbgt")

# deseanonal and detrended data
data.swbgt_deseason_detrend.where(data.swbgt_deseason_detrend >= 3).mean("latitude").mean("longitude").plot(ax= ax[2], marker= ".", markersize = 2, linestyle = "", color = "tab:blue" ,\
        label = "a) swbgt deseasonalised & detrended >= 3°C")
data.swbgt_deseason_detrend.where(data.swbgt_deseason_detrend >= 3).mean("latitude").mean("longitude").where(data.swbgt >= 15).plot(ax= ax[2], marker= "x", markersize = 2, linestyle = "-" , color = "tab:orange",linewidth = "0.5",\
        label = "b) swbgt deseasonalised & detrended >= 3°C and swbgt >= 15°C")

linreg, x = trends.calculate_linreg(data = data.swbgt_deseason_detrend.where(data.swbgt_deseason_detrend >= 3).where(data.swbgt >= 15).mean("latitude").mean("longitude"), variable = "swbgt_deseason_detrend", mutable= False)
# print(data.swbgt_deseason_detrend.where(data.swbgt_deseason_detrend >= 3).where(data.swbgt >= 15))
ax[2].plot(data.time, linreg.slope * x + linreg.intercept, color = "tab:red", linestyle = '--',\
        label = "linear regression of b) with slope = {:.2E}".format(linreg.slope*1e9*365*24*60*60*30) + r"$\frac{°C}{30 years}$")
ax[2].set_title("(III) filtered data from plot (II)")


for axs in ax :
    axs.grid()
    axs.legend()
    axs.set
    axs.set_ylabel("°C")

# fig.savefig(plot_dict["savefolder"] + plot_dict["savename"]+ "self.png")
plt.show()
# data = data.resample(time="1d").mean()
# print(data)
# data.to_netcdf("test.nc")
# data.swbgt_deseason_detrend.where(data.swbgt >= 20).plot(ax = ax[2], marker= "x", markersize = 4, linestyle = "-" )
# data.swbgt_deseason_detrend.where(data.swbgt_deseason_detrend >= 3).plot(ax = ax[2], marker= "+", markersize = 4, linestyle = "None" )


# # import Scripts.Calculation.trends as trends

# # trends.calculate_linreg(data= data)

# # plt.show()
# # #%%
# # print(np.shape(data.time["time.month"].values))
# # fig = plt.figure()
# # plt.scatter(data.time[0:500:1], data.swbgt_deseason_detrend.isel(latitude=0, longitude = 0)[0:500:1], c = data.swbgt.values[0:500:1], cmap = "RdBu") #, cmap = "HSV")

# # plt.scatter(ds = data, x = "time", y = "swbgt", hue = data.time["time.month"])
# # plt.show()
# # add temp to data 
# print("load temp")
# data["t2m"] = xr.open_dataset(folderpath + r"\t2m_era20c_1900-2010.nc").sel(time = timeslice, latitude= latitude_slice, longitude= longitude_slice).t2m -273.15
# data["t2m_rolling"] = data.t2m.rolling(time= rolling_window).mean()
# # create rolling mean
# print("swbgt rolling")
# data["swbgt_rolling"] = data.swbgt.rolling(time= rolling_window).mean()
# # create quantile data
# print("swbgt quantile")
# data["swbgt_quantile"] = data.swbgt.where(data.swbgt >= data.swbgt.quantile(quantiles))
# # create time_size_start
# print("time_size_start")
# data["time_size_start"] = (("time"), np.array(data.time.values - data.time[0].values, dtype=float))

# save all data
# data.to_netcdf(save_dict["savefolder"] + save_dict["savename"], mode = 'a')
# print("other stuff saved {:.2f}".format(time.time()-st))
# data = xr.open_dataset(save_dict["savefolder"] + save_dict["savename"] + 'test.nc')
#%% calculate other stuff

# # create rolling mean  
# data["swbgt_rolling"] = data.swbgt.rolling(time= rolling_window).mean()
# # create quantile data
# data["swbgt_quantile"] = data.swbgt.where(data.swbgt >= data.swbgt.quantile(quantiles))
#%%




# #%%
# fig = plt.figure(figsize=(15,6),tight_layout=True)
# gs = gridspec.GridSpec(3,1)

# ax1 = fig.add_subplot(gs[0:2,:])
# data.swbgt.isel(latitude = 0, longitude = 0).plot.line(ax = ax1, \
#         color = 'tab:blue', linestyle= ':', linewidth = 0.5, label = "a) swbgt")
# data.swbgt_rolling.isel(latitude = 0, longitude = 0).plot.line(ax = ax1, \
#         color = 'k', linestyle= '-', label= "b) swbgt with rolling window {} days".format(data.swbgt_rolling.attrs["rolling_window_length"]/4))

# def calculate_linreg(data, variable = "swbgt") : 
#         x = np.array(data.time.values - data.time[0].values, dtype=float)
#         y = data[variable].isel(latitude=0, longitude= 0).values
#         idx = (np.isfinite(x)) & (np.isfinite(y))
#         return stats.linregress(x[idx], y[idx])

# slope, intercept, r_value, p_value, std_err = calculate_linreg(data, variable = "swbgt_rolling")
# lin_reg = lambda x : slope*x + intercept

# lin_reg(data.time_size_start).plot(ax = ax1, \
#         color = "tab:red", linestyle= ':', linewidth = 3, label = "linear regression of b) with:\nR = {:.2f}".format(r_value))

# data.swbgt_quantile.isel(quantile = 2,latitude=0, longitude=0).plot.line('+', color= "tab:orange", label="0.99 quantile values")
# ax2.set_ylabel("swbgt in °C")

# ax2 = fig.add_subplot(gs[2,:])

# data.t2m.isel(latitude = 0, longitude = 0).plot.line(ax = ax2, \
#         color = 'tab:blue', linestyle= ':', linewidth = 0.5, label = "a) temperature")
# data.t2m_rolling.isel(latitude = 0, longitude = 0).plot.line(ax = ax2, \
#         color = 'k', linestyle= '-', label= "b) temperature with rolling window {} days".format(data.swbgt_rolling.attrs["rolling_window_length"]/4))
# ax2.set_ylabel("temperature in °C")

# # plt.show()
# plt.ylabel("swbgt")
# plt.legend()
# plt.grid()
# plt.xlim([data.time[0].values,data.time[-1].values])
# plt.title("swbgt timeseries for {}".format(city))
# plt.savefig(plot_dict["savefolder"] + plot_dict["savename"])
# # %%
