#%%

def main():
    import xarray as xr
    import time
    import numpy as np
    from percentile_conditions import n, timestart, timeend
    # timestart = "1900-01-01"
    # timeend = "1929-12-31"
    # n = 8

    save_dict = {"savefolder" : r"C:\Users\Nils Niebaum\Documents\Uni\Bachlor_thesis\Prog\Data",
                    "savename" : r"\swbgt_mean_" + str(timestart) + '_' + str(timeend) + '_percentiles.nc'
            }

    def combine_mean(n):
            combined_data = False
            for idx in range(n-1):
                    swbgt = xr.open_dataset(save_dict["savefolder"] + save_dict["savename"] + str(idx) + '.nc')
                    print(swbgt)
                    if combined_data :
                            combined_data = xr.merge([combined_data, swbgt])
                    else :
                            combined_data = xr.merge([swbgt])

            return combined_data
    combine_mean(n).to_netcdf(save_dict["savefolder"] + save_dict["savename"])

if __name__ == "__main__":
    main()