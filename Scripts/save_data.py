#%% HEAT WAVES IDENTIFICTION
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import time
import numpy as np
import Calculation.humiditycalculation as humidcalc
import Calculation.trend_func as trends
import Calculation.save_detrend_deseason as sdd

output_dtype = float
timestart = "1900-01-01"
timeend = "1930-12-31"
location = "europe_all"
# rolling_window = 4*365

load_dict = {
        't2m_path' : r"I:\Bachlor_thesis\Data\t2m_era20c_1900-2010.nc",
        'd2m_path' : r"I:\Bachlor_thesis\Data\td2m_era20c_1900-2010.nc",
        'time_slice' : slice(timestart, timeend),
        "latitude_slice" : slice(75.0, 31.0),
        "longitude_slice" : slice(-14.5, 44),
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

data = sdd.save_calc_swbgt(save_dict= save_dict, load_dict = load_dict)
print(data)

# data = xr.open_dataset(save_dict["savefolder"] + save_dict["savename"])
# print(data)
data.swbgt.plot()
plt.show()

# data_west = xr.open_dataset(r"C:\Users\Nils Niebaum\Documents\Uni\Bachlor_thesis\Prog\Data\europe_west_1900-01-01_1930-12-31.nc")
# data_east = xr.open_dataset(r"C:\Users\Nils Niebaum\Documents\Uni\Bachlor_thesis\Prog\Data\europe_1900-01-01_1930-12-31.nc")
# print("loaded")
# st = time.time()
# xr.merge([data_west,data_east]).to_netcdf(r"C:\Users\Nils Niebaum\Documents\Uni\Bachlor_thesis\Prog\Data\europe_all_1900-01-01_1930-12-31.nc")
# print("done {:.2f} s".format(time.time() - st))