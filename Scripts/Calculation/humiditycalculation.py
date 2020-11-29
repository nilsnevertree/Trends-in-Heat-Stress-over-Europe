import xarray as xr
import numpy as np

def e_gufunc(t_c) :
        return 6.1094 * np.exp(17.625 * t_c / (t_c + 243.04))

def e_func(t_c) : 
        return xr.apply_ufunc(
                e_gufunc,
                t_c,
                input_core_dims=[[]],
                dask="parallelized",
                output_dtypes=[float],
        )

def e_rh_gufunc(t_c, t_dew) :
        return (rh_func(t_c, t_dew)) * e_func(t_c)

def e_rh_func(t_c,t_dew, dim = "time"):
        return xr.apply_ufunc(
                e_rh_gufunc,
                t_c,
                t_dew,
                input_core_dims=[[], []],
                dask="parallelized",
                output_dtypes=[float],
        )

def rh_gufunc2(t_c,t_dew) :
        return e_gufunc(t_dew) / e_gufunc(t_c)

def rh_gufunc(t_c,t_dew) :
        return e_func(t_dew) / e_func(t_c)

def rh_func(t_c,t_dew, dim= 'time') :
        return xr.apply_ufunc(
                rh_gufunc,
                t_c,
                t_dew,
                input_core_dims=[[], []],
                dask="parallelized",
                output_dtypes=[float],
        )

def swbgt_gufunc(t_c, t_dew) :
#         return 0.56*t_c + 0.393 * e_rh_func(t_c,t_dew) * e_func(t_c) + 3.94
        return 0.56*t_c + 0.393 * e_func(t_dew) + 3.94

def swbgt_func(t_c, t_dew, dim = "time") :
        return xr.apply_ufunc(
                swbgt_gufunc,
                t_c,
                t_dew,
                input_core_dims=[[], []],
                dask="parallelized",
                output_dtypes=[float],
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
                output_dtypes=[float],
        )
    
def hwmid_func(t_c, mask):
    '''
    This fuction calculates the HWMId as in Russo et al. 2014.
    INPUT :
    t_c : temperature Dataarray 
    mask : mask of heatwave days (you can use calc_heatwave_index)
    OUTPUT:
    hwmid values 
    '''
    
    # 1. calculate the 25th and 75th qunatile of maximum yearly temperature
    t_25p = t_c.groupby("time.year").max('time').quantile(0.25, dim = "year")
    t_75p = t_c.groupby("time.year").max('time').quantile(0.75, dim = "year")

    # mask the temperature values (has to be part of a heatwave and hsa to be above the 25th percentile)
    t_masked = t_c.where(mask)
    mask_allowed = t_masked > t_25p
    t_allowed = t_masked.where(mask_allowed)
    
    # return the HWMId
    return (t_allowed - t_25p) / (t_75p -t_25p)