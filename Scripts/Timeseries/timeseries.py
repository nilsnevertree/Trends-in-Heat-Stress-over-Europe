#%% PERCEBTILE OF SWBGT
import xarray as xr
import matplotlib.pyplot as plt
import time
import numpy as np
from scipy import stats
import plotfunctions as plotfunc

import matplotlib.gridspec as gridspec

output_dtype = float
timestart = "1900-01-01"
timeend = "2010-12-31"
city = "kiel"
rolling_window = 4*365
latitude_slice = slice(55.0, 53.0)
longitude_slice = slice(10.0, 11.0)
quantiles = np.array([0.95,0.98,0.99])
plot_dict = {"doplot" : True,
                "levels" : np.arange(18, 35, 1),
                "savefolder" : r"C:\Users\Nils Niebaum\Documents\Uni\Bachlor_thesis\Prog\Plots",
                "savename" : r"\swbgt_timeseries" + "{}_{}_{}.png".format(city, timestart, timeend)
        }
save_dict = {"dosave" : True,
                "savefolder" : r"C:\Users\Nils Niebaum\Documents\Uni\Bachlor_thesis\Prog\Data",
                "savename" : r"\swbgt_kiel" + "_{}_{}.nc".format(timestart, timeend)
        }
#%% load data 
folderpath = r"I:\Bachlor_thesis\Data"
timeslice = slice(timestart, timeend)
#%% functions

# input('press any key to keep on going')
def e_gufunc(t_c) :
        return 6.1094 * np.exp(17.625 * t_c / (t_c + 243.04))

def e_func(t_c) : 
        return xr.apply_ufunc(
                e_gufunc,
                t_c,
                input_core_dims=[[]],
                dask="parallelized",
                output_dtypes=[output_dtype],
        )

def e_rh_gufunc(t_c, t_dew) :
        return (rh_func(t_c, t_dew) / 100) * e_func(t_dew)

def e_rh_func(t_c,t_dew, dim = "time"):
        return xr.apply_ufunc(
                e_rh_gufunc,
                t_c,
                t_dew,
                input_core_dims=[[], []],
                dask="parallelized",
                output_dtypes=[output_dtype],
        )

def rh_gufunc2(t_c,t_dew) :
        return e_gufunc(t_c) / e_gufunc(t_dew)

def rh_gufunc(t_c,t_dew) :
        return e_func(t_c) / e_func(t_dew)

def rh_func(t_c,t_dew, dim= 'time') :
        return xr.apply_ufunc(
                rh_gufunc,
                t_c,
                t_dew,
                input_core_dims=[[], []],
                dask="parallelized",
                output_dtypes=[output_dtype],
        )

def swbgt_gufunc(t_c, t_dew) :
        return 0.56*t_c + 0.393 * e_rh_func(t_c, t_dew) / 100 + 3.94

def swbgt_func(t_c, t_dew, dim = "time") :
        return xr.apply_ufunc(
                swbgt_gufunc,
                t_c,
                t_dew,
                input_core_dims=[[], []],
                dask="parallelized",
                output_dtypes=[output_dtype],
        )

def t_ws_gufunc(t_c, t_dew) :
        rh = rh_func(t_c,t_dew)
        return t_c * np.arctan(0.151977 * np.sqrt(rh + 8.313659)) + np.arctan(t_c + rh) - np.arctan(rh - 1.676331) + 0.00391838 * rh **(3/2) * np.arctan(0.023101*rh) - 4.68035

def t_ws_func(t_c,t_dew, dim = "time"):
        return xr.apply_ufunc(
                t_ws_gufunc,
                t_c,
                t_dew,
                input_core_dims=[[], []],
                dask="parallelized",
                output_dtypes=[output_dtype],
        )

#%% =============== LOAD DATA ================== 

try : # if file exists load the data
        data = xr.open_dataset(save_dict["savefolder"] + save_dict["savename"])
        print("=======\nload data done")
except: #%% ==================== CALCULATION OF SWBGT ================
        print("=======\nfile not found")
        folderpath = r"I:\Bachlor_thesis\Data"

        st = time.time()
        temp = xr.open_dataset(folderpath + r"\t2m_era20c_1900-2010.nc")#, chunks={'time': 12*4*365})
        print("load temp took {:.2f}".format(time.time()-st))
        st = time.time()
        temp_dew = xr.open_dataset(folderpath + r"\td2m_era20c_1900-2010.nc")#, chunks={'time': 12*4*365})
        print("load dew took {:.2f}".format(time.time()-st))
        st = time.time()
        all_data = xr.merge([temp, temp_dew])
        print("merge took {:.2f}".format(time.time()-st))
        print(all_data)

        st = time.time()
        data = all_data.sel(time = timeslice, latitude= latitude_slice, longitude= longitude_slice)
        print("select took {:.2f}".format(time.time()-st))

        # #%% ========== CHUNK DATA ==============
        # rechunk data for perfomance
        rechunk_dim = {
                "time" : int(np.shape(data.t2m)[0] / 8),
                }
        rechunked = data.chunk(rechunk_dim) -273.15

        swbgt = swbgt_func(rechunked.t2m, rechunked.d2m)
        swbgt.name = "swbgt"
        swbgt.attrs = {"place" : city}
        st = time.time()
        xr.merge([swbgt]).to_netcdf(save_dict["savefolder"] + save_dict["savename"])
        print("swbgt saved {:.2f}".format(time.time()-st))
        data = xr.open_dataset(save_dict["savefolder"] + save_dict["savename"])
        # create rolling mean 
        data["swbgt_rolling"] = data.swbgt.rolling(time= rolling_window).mean()
        # create quantile data
        data["swbgt_quantile"] = data.swbgt.where(data.swbgt >= data.swbgt.quantile(quantiles))
        # create time_size_start 
        data["time_size_start"] = (("time"), np.array(data.time.values - data.time[0].values, dtype=float))
        # add temp to data 
        data["t2m"] = rechunked.t2m
        data["t2m_rolling"] = data.t2m.rolling(time= rolling_window).mean()
        # save all data
        data.to_netcdf(save_dict["savefolder"] + save_dict["savename"], mode = 'a')
        print("other stuff saved {:.2f}".format(time.time()-st))
        data = xr.open_dataset(save_dict["savefolder"] + save_dict["savename"])


# add temp to data 
print("load temp")
data["t2m"] = xr.open_dataset(folderpath + r"\t2m_era20c_1900-2010.nc").sel(time = timeslice, latitude= latitude_slice, longitude= longitude_slice).t2m -273.15
data["t2m_rolling"] = data.t2m.rolling(time= rolling_window).mean()
# create rolling mean
print("swbgt rolling")
data["swbgt_rolling"] = data.swbgt.rolling(time= rolling_window).mean()
# create quantile data
print("swbgt quantile")
data["swbgt_quantile"] = data.swbgt.where(data.swbgt >= data.swbgt.quantile(quantiles))
# create time_size_start
print("time_size_start")
data["time_size_start"] = (("time"), np.array(data.time.values - data.time[0].values, dtype=float))

# save all data
data.to_netcdf(save_dict["savefolder"] + save_dict["savename"], mode = 'a')
print("other stuff saved {:.2f}".format(time.time()-st))
data = xr.open_dataset(save_dict["savefolder"] + save_dict["savename"] + 'test.nc')
#%% calculate other stuff

# # create rolling mean  
# data["swbgt_rolling"] = data.swbgt.rolling(time= rolling_window).mean()
# # create quantile data
# data["swbgt_quantile"] = data.swbgt.where(data.swbgt >= data.swbgt.quantile(quantiles))
#%%




#%%
fig = plt.figure(figsize=(15,6),tight_layout=True)
gs = gridspec.GridSpec(3,1)

ax1 = fig.add_subplot(gs[0:2,:])
data.swbgt.isel(latitude = 0, longitude = 0).plot.line(ax = ax1, \
        color = 'tab:blue', linestyle= ':', linewidth = 0.5, label = "a) swbgt")
data.swbgt_rolling.isel(latitude = 0, longitude = 0).plot.line(ax = ax1, \
        color = 'k', linestyle= '-', label= "b) swbgt with rolling window {} days".format(data.swbgt_rolling.attrs["rolling_window_length"]/4))

def calculate_linreg(data, variable = "swbgt") : 
        x = np.array(data.time.values - data.time[0].values, dtype=float)
        y = data[variable].isel(latitude=0, longitude= 0).values
        idx = (np.isfinite(x)) & (np.isfinite(y))
        return stats.linregress(x[idx], y[idx])

slope, intercept, r_value, p_value, std_err = calculate_linreg(data, variable = "swbgt_rolling")
lin_reg = lambda x : slope*x + intercept

lin_reg(data.time_size_start).plot(ax = ax1, \
        color = "tab:red", linestyle= ':', linewidth = 3, label = "linear regression of b) with:\nR = {:.2f}".format(r_value))

data.swbgt_quantile.isel(quantile = 2,latitude=0, longitude=0).plot.line('+', color= "tab:orange", label="0.99 quantile values")
ax2.set_ylabel("swbgt in °C")

ax2 = fig.add_subplot(gs[2,:])

data.t2m.isel(latitude = 0, longitude = 0).plot.line(ax = ax2, \
        color = 'tab:blue', linestyle= ':', linewidth = 0.5, label = "a) temperature")
data.t2m_rolling.isel(latitude = 0, longitude = 0).plot.line(ax = ax2, \
        color = 'k', linestyle= '-', label= "b) temperature with rolling window {} days".format(data.swbgt_rolling.attrs["rolling_window_length"]/4))
ax2.set_ylabel("temperature in °C")

# plt.show()
plt.ylabel("swbgt")
plt.legend()
plt.grid()
plt.xlim([data.time[0].values,data.time[-1].values])
plt.title("swbgt timeseries for {}".format(city))
plt.savefig(plot_dict["savefolder"] + plot_dict["savename"])
# %%
