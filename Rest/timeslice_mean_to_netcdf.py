#%%
import xarray as xr
import matplotlib.pyplot as plt
import time
import numpy as np
import cartopy.crs as ccrs
import cartopy as cart
import os


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
data = all_data.sel(time= slice("2005-12-31", "2010-12-31"))
print("=======\nselect took {:.2f}".format(time.time()-st))
print(data, '\n')

rechunked = data.chunk({"time": 4*30*4}) -273.15
print(rechunked.chunks)

def e_gufunc(T_C) :
        return 6.1094 * np.exp(17.625 * T_C / (T_C + 243.04))

def e_func(T_C) : 
        return xr.apply_ufunc(
                e_gufunc,
                T_C,
                input_core_dims=[[]],
                dask="parallelized",
                output_dtypes=[float],
        )

def e_RH_gufunc(T_C, T_dew) :
        return (RH_func(T_C, T_dew) / 100) * e_func(T_dew)

def e_RH_func(T_C,T_dew, dim = "time"):
        return xr.apply_ufunc(
                e_RH_gufunc,
                T_C,
                T_dew,
                input_core_dims=[[], []],
                dask="parallelized",
                output_dtypes=[float],
        )

def RH_gufunc2(T_C,T_dew) :
        return e_gufunc(T_C) / e_gufunc(T_dew)

def RH_gufunc(T_C,T_dew) :
        return e_func(T_C) / e_func(T_dew)

def RH_func(T_C,T_dew, dim= 'time') :
        return xr.apply_ufunc(
                RH_gufunc,
                T_C,
                T_dew,
                input_core_dims=[[], []],
                dask="parallelized",
                output_dtypes=[float],
        )

def sWBGT_gufunc(T_C, T_dew) :
        return 0.56*T_C + 0.393 * e_RH_func(T_C, T_dew) / 100 + 3.94

def sWBGT_func(T_C, T_dew, dim = "time") :
        return xr.apply_ufunc(
                sWBGT_gufunc,
                T_C,
                T_dew,
                input_core_dims=[[], []],
                dask="parallelized",
                output_dtypes=[float],
        )

def T_wS_gufunc(T_C, T_dew) :
        RH = RH_func(T_C,T_dew)
        return T_C * np.arctan(0.151977 * np.sqrt(RH + 8.313659)) + np.arctan(T_C + RH) - np.arctan(RH - 1.676331) + 0.00391838 * RH **(3/2) * np.arctan(0.023101*RH) - 4.68035

def T_wS_func(T_C,T_dew, dim = "time"):
        return xr.apply_ufunc(
                T_wS_gufunc,
                T_C,
                T_dew,
                input_core_dims=[[], []],
                dask="parallelized",
                output_dtypes=[float],
        )

input('stop?')

st = time.time()
mean = rechunked.mean(dim= 'time')
processing_time = time.time()-st
print("=======\nrechunked mean took {:.10f}".format(processing_time))

input('stop?')
st = time.time()
sWBGT_mean = swbgt_func(mean.t2m, mean.d2m)
processing_time = time.time()-st
print("=======\ncalc mean sWGBT took {:.10f}".format(processing_time))

st = time.time()
sWBGT = swbgt_func(rechunked.t2m, rechunked.d2m)
# sWBGT_high = sWBGT.where(sWBGT > sWBGT.mean(dim = "time") + sWBGT.std(dim = "time")).mean(dim= "time")
processing_time = time.time()-st
print("=======\ncalc sWGBT high took {:.10f}".format(processing_time))


st = time.time()
print(np.shape(sWBGT_mean))
print(np.shape(sWBGT))
print("=======\n took {:.10f}".format(processing_time))



def area_plot(data = None, ax = None, levels = np.arange(16,30,0.2)) :
        ax.coastlines(resolution='auto', color='k', zorder=3)
        ax.gridlines(color='grey', linestyle=':', zorder=3)
        plot = data.plot.pcolormesh(ax = ax, cmap= "RdBu_r", x="longitude", y="latitude", \
                levels= levels, yincrease=True, transform=ccrs.PlateCarree())
        #ax.add_feature(cart.feature.OCEAN, zorder=2, facecolor='w', edgecolor='k')


fig = plt.figure(figsize = (12,21))

ax1 = fig.add_subplot(311, projection=ccrs.Robinson())

st = time.time()
area_plot(sWBGT_mean, ax= ax1)
ax1.title.set_text('sWBGT mean')
processing_time = time.time()-st
print("=======\nplot sWGBT mean took {:.10f}".format(processing_time))

ax2 = fig.add_subplot(312, projection=ccrs.Robinson())

st = time.time()
area_plot(sWBGT.where(sWBGT > sWBGT.mean(dim = "time") + sWBGT.std(dim = "time")).mean(dim= "time"), ax= ax2)
ax2.title.set_text('sWBGT - mean of values higher than mean + stddev')
processing_time = time.time()-st
print("=======\nplot sWGBT high took {:.10f}".format(processing_time))

ax3 = fig.add_subplot(313, projection=ccrs.Robinson())

st = time.time()
area_plot(sWBGT.max(dim= "time"), ax= ax3)
ax3.title.set_text('sWBGT max')
processing_time = time.time()-st
print("=======\nplot sWGBT max took {:.10f}".format(processing_time))

fig.suptitle('timeslice ' + "2005-12-31" +' / ' "2010-12-31")
#plt.savefig(r"C:\Users\Nils Niebaum\Documents\Uni\Bachlor_thesis\Prog\sWBGT_both.svg")
plt.savefig(r"C:\Users\Nils Niebaum\Documents\Uni\Bachlor_thesis\Prog\sWBGT_three_2005-2010_RdBu.png", dpi = 1000)

# plt.show()
# input("ok?")
# # st = time.time()
# # sWBGT_func(data.t2m, data.d2m).mean(dim= 'time').to_netcdf(r'Data\'first10_year_mean_unchunked.nc')
# # processing_time = time.time()-st
# # print("original took {:.10f}".format(processing_time))