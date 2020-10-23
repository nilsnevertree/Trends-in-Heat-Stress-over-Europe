#%%
import xarray as xr
import matplotlib.pyplot as plt
import time
import numpy as np
import cartopy.crs as ccrs
import cartopy as cart
folderpath = r"I:\Bachlor_thesis\Data"

st = time.time()
temp = xr.open_dataset(r'Data\temp_firstmonth.nc', chunks={'time': 10})
print("took {:.2f}".format(time.time()-st))
temp_dew = xr.open_dataset(r'Data\temp_dew_firstmonth.nc', chunks={'time': 10})
print("took {:.2f}".format(time.time()-st))

#%% FUNCTIONS
e_lambda = lambda T_C : 6.1094 * np.exp(17.625 * T_C / (T_C + 243.04)) 
RH_lambda = lambda T_C, T_dew : e_lambda(T_C) / e_lambda(T_dew)
e_RH_lambda = lambda T_C, T_dew : (RH_lambda(T_C, T_dew) / 100) * e_lambda(T_dew)
sWBGT_lambda = lambda T_C, T_dew : 0.56*T_C + 0.393 * e_RH_lambda(T_C, T_dew) / 100 + 3.94
eq_8_lambda = lambda T_C, RH   : T_C * np.arctan(0.151977*np.sqrt(RH + 8.313659)) + np.arctan(T_C + RH) - np.arctan(RH - 1.676331) + 0.00391838 * RH **(3/2) * np.arctan(0.023101*RH) - 4.68035

def area_plot(data = None, ax = None, levels = np.linspace(-3,3,10)) :
        ax.coastlines(resolution='auto', color='k', zorder=3)
        ax.gridlines(color='grey', linestyle=':', zorder=3)
        plot = data.plot.pcolormesh(ax = ax, levels= levels, x="longitude", y="latitude", \
                yincrease=True, transform=ccrs.PlateCarree())
        ax.add_feature(cart.feature.OCEAN, zorder=2, facecolor='w', edgecolor='k')

#%%
starttime = time.time()
sWBGT = sWBGT_lambda(temp.t2m - 273.15, temp_dew.d2m - 273.15)
print(time.time() - starttime)

levels = np.arange(16,28,0.5)

#%%
starttime = time.time()
fig = plt.figure(figsize= (12,8))
gs = fig.add_gridspec(2,3)
ax1 = fig.add_subplot(gs[0,0:2], projection=ccrs.Robinson())
data = sWBGT.where(sWBGT > sWBGT.mean(dim = "time") + sWBGT.std(dim = "time")).mean(dim= "time")
area_plot(data, ax = ax1, levels = levels)
data.mean(dim = 'longitude').plot(ax=ax1)
ax2 = fig.add_subplot(gs[1, 0:2], projection=ccrs.Robinson())
area_plot(sWBGT.mean(dim= "time"), ax=ax2, levels = levels)
sWBGT.mean(dim= "time").mean(dim = 'longitude').plot(ax=ax1)
print(time.time() - starttime)
plt.savefig('sWBGT_firstmonth.png')
# %%
