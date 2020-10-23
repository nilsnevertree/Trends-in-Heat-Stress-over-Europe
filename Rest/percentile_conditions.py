#%%

import numpy as np
# import humiditycalculation as hc
output_dtype = np.float16
timestart = "1900-01-01"
timeend = "1929-12-31"
quantiles = np.array([0.95,0.98,0.99])
latitude_slice = slice(55.0, 53.0)
longitude_slice = slice(10.0, 11.0)

plot_dict = {"doplot" : True,
                "levels" : np.arange(18, 33, 1),
                "savefolder" : r"C:\Users\Nils Niebaum\Documents\Uni\Bachlor_thesis\Prog\Plots\Percentiles",
        }
save_dict = {"dosave" : True,
                "savefolder" : r"C:\Users\Nils Niebaum\Documents\Uni\Bachlor_thesis\Prog\Data",
                "savename" : r"\swbgt_mean_" + str(timestart) + '_' + str(timeend) + '_percentiles.nc'
        }

n = 8

filepath = r"C:\Users\Nils Niebaum\Documents\Uni\Bachlor_thesis\Prog\Data" + r"\swbgt_mean_" + str(timestart) + '_' + str(timeend) + '_percentiles.nc'
folderpath = r"I:\Bachlor_thesis\Data"