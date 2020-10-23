#%% PERCEBTILE OF SWBGT
import xarray as xr
import matplotlib.pyplot as plt
import time
import numpy as np
import cartopy.crs as ccrs
import cartopy as cart
from dask.diagnostics import ProgressBar
import os
import matplotlib.gridspec as gridspec
# import humiditycalculation as hc
import plotfunctions as plotfunc
output_dtype = np.float16
from percentile_conditions import timestart, timeend

quantiles = [0.95, 0.98, 0.99]
plot_dict = {"doplot" : True,
                "levels" : np.arange(18, 33, 1),
                "quantile" : quantiles,
                "savefolder" : r"C:\Users\Nils Niebaum\Documents\Uni\Bachlor_thesis\Prog\Plots\Percentiles",
        }
save_dict = {"dosave" : True,
                "savefolder" : r"C:\Users\Nils Niebaum\Documents\Uni\Bachlor_thesis\Prog\Data",
        }
#%% load data 
filepath = r"C:\Users\Nils Niebaum\Documents\Uni\Bachlor_thesis\Prog\Data" + r"\swbgt_mean_" + str(timestart) + '_' + str(timeend) + '_percentiles.nc'

#%% =============== LOAD DATA ================== 

try :
        data = xr.open_dataset(filepath)
        print("=======\nload data done")
except :
        import slice_combine
        slice_combine.main()
        print("=======\ncombine data done")
        data = xr.open_dataset(filepath)
        print("=======\nload data done")


data.swbgt_mean.attrs = {"unit" : "°C","long_name" : "swbgt in °C"}
print(data.swbgt_mean)
#%% =================== PLOTING =======================
fig = plt.figure(figsize = (12,18))
gs = gridspec.GridSpec(3, 1)
idx = 0
for quant in quantiles:
        plot_dict["title"] = 'sWBGT {:.0f}th metric mean - from {} until {}'.format(quant*100, timestart, timeend)

        ax = fig.add_subplot(gs[idx,0], projection=ccrs.Robinson())
        st = time.time()
        PC = plotfunc.area_plot(data.swbgt_mean.sel(quantile= quant), ax= ax, levels = plot_dict["levels"])
        ax.title.set_text(plot_dict["title"])
        processing_time = time.time()-st
        print("=======\nplot sWGBT {} mean:\n{:.1f}s".format(quant, processing_time))
        idx += 1

# fig.suptitle('timeslice ' + timestart + ' till '  + timeend)
#plt.savefig(r"C:\Users\Nils Niebaum\Documents\Uni\Bachlor_thesis\Prog\sWBGT_both.svg")
plot_dict["savename"] = r"\swbgt_mean_" + str(timestart) + '_' + str(timeend) + '_percentile.png'
plt.savefig(plot_dict["savefolder"] + plot_dict["savename"], dpi = 400)

# %%
