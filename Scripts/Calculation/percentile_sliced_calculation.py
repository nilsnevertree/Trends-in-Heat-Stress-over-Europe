#%% PERCEBTILE OF SWBGT
import xarray as xr
import time
import numpy as np
import os
import gc
from memory_profiler import profile

output_dtype = np.float16
timestart = "1900-01-01"
timeend = "1929-12-31"
quantiles = np.array([0.95,0.98,0.99])
plot_dict = {"doplot" : True,
                "levels" : np.arange(18, 35, 1),
                "savefolder" : r"C:\Users\Nils Niebaum\Documents\Uni\Bachlor_thesis\Prog\Plots\Percentiles",
        }
save_dict = {"dosave" : True,
                "savefolder" : r"C:\Users\Nils Niebaum\Documents\Uni\Bachlor_thesis\Prog\Data",
                "savename" : r"\swbgt_mean_" + str(timestart) + '_' + str(timeend) + '_percentiles.nc'
        }
n = 8
#%% load data 
folderpath = r"I:\Bachlor_thesis\Data"

st = time.time()
temp = xr.open_dataset(folderpath + r"\t2m_era20c_1900-2010.nc")#, chunks={'time': 12*4*365})
print("=======\nload temp took {:.2f}".format(time.time()-st))

st = time.time()
temp_dew = xr.open_dataset(folderpath + r"\td2m_era20c_1900-2010.nc")#, chunks={'time': 12*4*365})
print("=======\nload dew took {:.2f}".format(time.time()-st))

st = time.time()
all_data = xr.merge([temp, temp_dew])
print("=======\nmerge took {:.2f}".format(time.time()-st))

st = time.time()
data = all_data.sel(time= slice(timestart, timeend))
print("=======\nselect took {:.2f}".format(time.time()-st))


#%% ========== CHUNK DATA ==============
# rechunk data for perfomance
rechunk_dim = {
        "time" : int(np.shape(data.t2m)[0] / 8), 
        # "latitude" : 20,
        # "longitude" : 40,
        }
rechunked = data.chunk(rechunk_dim) - 273.15
print(rechunked.chunks)

# input('press any key to keep on going')
def e_gufunc(t_c) :
        return 6.1094 * np.exp(17.625 * t_c / (t_c + 243.04))

def e_func(t_c) : 
        return xr.apply_ufunc(
                e_gufunc,
                t_c,
                input_core_dims=[[]],
                dask="parallelized",
                output_dtypes=[np.float16],
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
                output_dtypes=[np.float16],
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
                output_dtypes=[np.float16],
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
                output_dtypes=[np.float16],
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
                output_dtypes=[np.float16],
        )

#%% ==================== CALCULATION OF SWBGT ================

@profile
def calc_percentile(idx, lat):

        rechunk_dim = {
                "time" : int(np.shape(data.t2m)[0] / 8), 
                # "latitude" : 20,
                # "longitude" : 40,
                }

        print("=======\n{} of {}".format(idx+1, n-1))
        temporary = rechunked.sel(latitude= slice(lat[idx+1], lat[idx])).chunk(rechunk_dim)
        st = time.time()
        swbgt = swbgt_func(temporary.t2m, temporary.d2m).compute()
        processing_time = time.time() - st
        print("calc swbgt: {:.1f}s".format(processing_time))

        del temporary
        gc.collect()
        
        # calc mean with chunked data
        st = time.time()
        rechunk_dim = {
                # "time" : int(np.shape(data.t2m)[0] / 8), 
                "latitude" : 45,
                "longitude" : 90,
                }
        swbgt_quant_mean = swbgt.where(swbgt >= swbgt.quantile(q = quantiles, keep_attrs = True, dim = "time")).chunk(rechunk_dim).mean(dim= "time")
        
        del swbgt
        gc.collect()
        
        swbgt_quant_mean.name = "swbgt_mean"
        swbgt_quant_mean.attrs = {"long_name" : "simple wet bulb glob temperature metrics temporal mean"}
        processing_time = time.time() - st
        print("calc quantile and combine: {:.1f}s".format(processing_time))
        swbgt_quant_mean.to_netcdf(save_dict["savefolder"] + save_dict["savename"] + str(idx) + '.nc')
        del swbgt_quant_mean
        gc.collect()
        print("=======TIME: {:.1f}========".format(time.time() - st_all))

st_all = time.time()
@profile
def main(n, start=0):
        lat = np.linspace(-92,100,n)
        for idx in np.arange(start,n -1,1):
                calc_percentile(idx, lat)

main(n)