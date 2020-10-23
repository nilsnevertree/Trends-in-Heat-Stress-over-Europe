#%% PERCEBTILE OF SWBGT
import xarray as xr
import matplotlib.pyplot as plt
import time
import numpy as np
import cartopy.crs as ccrs
import cartopy as cart
from dask.diagnostics import ProgressBar
import os
# import humiditycalculation as hc
import plotfunctions as plotfunc
output_dtype = np.float16
timestart = "2009-01-01"
timeend = "2010-12-31"
quantiles = np.array([0.95,0.98,0.99])
plot_dict = {"doplot" : True,
                "levels" : np.arange(18, 35, 1),
                "savefolder" : r"C:\Users\Nils Niebaum\Documents\Uni\Bachlor_thesis\Prog\Plots\Percentiles",
        }
save_dict = {"dosave" : True,
                "savefolder" : r"C:\Users\Nils Niebaum\Documents\Uni\Bachlor_thesis\Prog\Data",
        }
#%% load data 
folderpath = r"I:\Bachlor_thesis\Data"
timeslice = slice(timestart, timeend)
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
print(all_data)

st = time.time()
data = all_data.sel(time= slice(timestart, timeend))
print("=======\nselect took {:.2f}".format(time.time()-st))
print(data, '\n')

# #%% ========== CHUNK DATA ==============
# rechunk data for perfomance
rechunk_dim = {
        "time" : int(np.shape(data.t2m)[0] / 8), 
        # "latitude" : 20,
        # "longitude" : 40,
        }
rechunked = data.chunk(rechunk_dim) -273.15
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
#calculate swbgt
st = time.time()
swbgt = swbgt_func(rechunked.t2m, rechunked.d2m).compute()
processing_time = time.time()-st
print("=======\ncalc sWGBT:\n{:.1f}s".format(processing_time))
print(swbgt)
# input('press any key')
#%% ================== QUANTILES ==================
for quant in quantiles :
        print("-------\n=========\ncalc " + str(quant) + "quantile of swbgt" )
        # # input('ok?\n')
        percentile = int(quant*100)
        
        # set savenames
        plot_dict["savename"] = r"\swbgt_mean_" + str(timestart) + '_' + str(timeend) + '_' + str(percentile) + ".png"
        save_dict["savename"] = r"\swbgt_mean_" + str(timestart) + '_' + str(timeend) + '_' + str(percentile) + ".nc"
        plot_dict["title"] = 'sWBGT {}th metric mean - from {} until {}'.format(percentile, timestart, timeend)
        # # calc masks for percentiles
        # st = time.time()
        # swbgt_quant = swbgt.where(swbgt >= swbgt.quantile(quant, dim = "time"))
        # processing_time = time.time()-st
        # print("=======\ncalc sWGBT quant {}:\n{:.1f}s".format(percentile, processing_time))

        # calc mean with chunked data
        st = time.time()
        rechunk_dim = {
                # "time" : int(np.shape(data.t2m)[0] / 8), 
                "latitude" : 45,
                "longitude" : 90,
                }
        swbgt_quant_mean = swbgt.where(swbgt >= swbgt.quantile(quant, dim = "time")).chunk(rechunk_dim).mean(dim= "time")
        print("=======\ncalc sWGBT {} mean :\n{:.1f}s".format(percentile, processing_time))
        
        # =============== SAVE DATA =============
        if save_dict["dosave"]:
                st = time.time()
                swbgt_quant_mean.to_netcdf(save_dict["savefolder"] + save_dict["savename"])
                processing_time = time.time()-st
                print("=======\nsave sWGBT {} mean:\n{:.1f}s".format(percentile, processing_time))

        # =================== PLOTING =======================
        if plot_dict["doplot"]:
                fig = plt.figure(figsize = (12,6))
                ax1 = fig.add_subplot(111, projection=ccrs.PlateCarree())
                st = time.time()
                plotfunc.area_plot(swbgt_quant_mean, ax= ax1, levels = plot_dict["levels"])
                ax1.title.set_text(plot_dict["title"])
                processing_time = time.time()-st
                print("=======\nplot sWGBT {} mean:\n{:.1f}s".format(percentile, processing_time))

                # fig.suptitle('timeslice ' + timestart + ' till '  + timeend)
                #plt.savefig(r"C:\Users\Nils Niebaum\Documents\Uni\Bachlor_thesis\Prog\sWBGT_both.svg")
                plt.savefig(plot_dict["savefolder"] + plot_dict["savename"], dpi = 400)
