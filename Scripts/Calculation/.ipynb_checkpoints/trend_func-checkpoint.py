import numpy as np
import scipy.stats as stats
import xarray as xr
import warnings

def calculate_linreg(data, variable = "swbgt", mutable = True) : 
    '''
    Calculates temporal linear regression of variable with scipy.stats.linregress\n
    INPUT: 
    - data : xarray DataSet containing the "variable"
    - variable : name of variable you want to calculate the linear regression of (default: "swbgt")
    - latitude_idx and longitude_idx define whith isel() which spacial position shall used.
    - mutable = True : detrended data will be added to DataSet
    OUTPUT:
    - linreg : output of scipy.stats.linregress
    - x : x values used to calc lin reg.
    NOTE:
    - only 1D data can be used!!!!
    '''
    # create x values from time dimesion with timedelta and chage to float
    x = np.array(data.time.values - data.time[0].values, dtype=float)
    try : # try to select the latitude and longitude
        try:
            y = data[variable].values #isel(latitude= latitude_idx, longitude= longitude_idx).values
        except :
            y = data[variable].values
    except :
        try:
            y = data.values # isel(latitude= latitude_idx, longitude= longitude_idx).values
        except :
            y = data.values
    # values need to be finite
    print(x)
    idx = (np.isfinite(x)) & (np.isfinite(y))
    linreg = stats.linregress(x[idx], y[idx])
    if mutable:
        data[variable + "_trend"] = linreg.slope * x + linreg.intercept
        print(variable + "_trend added to DataSet")
        return linreg, x
    else :
        return linreg, x

def detrend_1d(data, variable = "swbgt" , latitude_idx = 0, longitude_idx = 0, mutable = True) :
    '''
    !!!!!!!!!!!!!!!!! OLD - CAN ONLY COMPUTE 1D !!!!!!!!!!!!!!!!!!!!!!!!!!!
    Calculates detrended values of "varaible" using calculate_linreg
    INPUT: 
    - data : xarray DataSet containing the "variable"
    - variable : name of variable you want to detrend (default: "swbgt")
    - latitude_idx and longitude_idx define whith isel() which spacial position shall used.
    - mutable = True : detrended data will be added to DataSet
    OUTPUT:
        mutable = True : detrended data will be added to DataSet , calculate_linreg(mutable= False) will be returned
        mutable 0 False : detrended data will be returned
    '''
    linreg, x = calculate_linreg(data= data, variable= variable, mutable= False)

    trend = linreg.slope* x + linreg.intercept
    print(data[variable].isel(latitude= latitude_idx, longitude= longitude_idx) - trend)
    if mutable:
        data[variable + "_detrend"] = data[variable].isel(latitude= latitude_idx, longitude= longitude_idx) - trend
        print(variable + "_detrend added to DataSet")
        return linreg, x
    else :
        return data[variable].isel(latitude= latitude_idx, longitude= longitude_idx) - trend

def detrend(data, variable = "swbgt" , area_mean = True, mutable = True) :
    '''
    Calculates detrended values of "varaible" using calculate_linreg
    INPUT: 
    - data : xarray DataSet containing the "variable"
    - variable : name of variable you want to detrend (default: "swbgt")
    - latitude_idx and longitude_idx define whith isel() which spacial position shall used.
    - mutable = True : detrended data will be added to DataSet
    OUTPUT:
        mutable = True : detrended data will be added to DataSet , calculate_linreg(mutable= False) will be returned
        mutable 0 False : detrended data will be returned
    '''
    temporary = None
    if area_mean:
        linreg, x = calculate_linreg(data= data.mean("latitude").mean('longitude'), variable= variable, mutable= False)
        trend = linreg.slope* x + linreg.intercept
        for latitude in data.latitude :
            for longitude in data.longitude :
                if temporary is None:
                    temporary = data[variable].copy(deep = True) * np.nan
                temporary = temporary.combine_first(data[variable].sel(latitude= latitude, longitude= longitude) - trend)#, concat_dim=["time"])
        if mutable:
            data[variable + "_detrend"] = temporary
            print(variable + "_detrend added to DataSet")
            return None
        return temporary

    else :
        temporary = None
        for latitude in data.latitude :
            for longitude in data.longitude :
                print(latitude.values, longitude.values)
                linreg, x = calculate_linreg(data= data.sel(latitude= latitude, longitude= longitude), variable= variable, mutable= False)
                trend = linreg.slope* x + linreg.intercept
                if temporary is None:
                    temporary = data[variable].copy(deep = True) * np.nan
                temporary = temporary.combine_first(data[variable].sel(latitude= latitude, longitude= longitude) - trend)#, concat_dim=["time"])
        if mutable:
            data[variable + "_detrend"] = temporary
            print(variable + "_detrend added to DataSet")
            return None
        return temporary


def deseason(data, variable = "swbgt", groupby = "week", rolling = False, mutable= True):
    '''
    Calculates deseasoned values of "varaible" using groupby time length
    INPUT: 
    - data : xarray DataSet containing the "variable"
    - variable : name of variable you want to detrend (default: "swbgt")
    - mutable = True : detrended data will be added to DataSet
    OUTPUT:
        mutable = True : detrended data will be added to DataSet , calculate_linreg(mutable= False) will be returned
        mutable 0 False : detrended data will be returned
    '''
    data[variable + "_deseason"] = data[variable].copy(deep= True )
    
    # with warnings.catch_warnings():
    warnings.filterwarnings("ignore")
    deseason = data.groupby("time." + groupby).mean('time')
    # print(deseason)
    # warnings.simplefilter("ignore")
    # print(deseason[groupby].min(),deseason[groupby].max())
    for step in np.arange(deseason[groupby].min(),deseason[groupby].max() + 1,1):
        mask = data.time["time." + groupby] == step
        data[variable + "_deseason"][mask] = (data[variable][mask] - deseason[variable].sel({groupby : step}))
    if rolling is not False:
        data[variable + "_deseason"] = data[variable + "_deseason"].rolling(time = rolling).mean()
    if mutable:
        print(variable + "_deseason added to DataSet")
        return None
    return data

def calc_daily_mean(data, values_per_day = 4):
    '''
    NOTE:
    Calculates daily DataSet if consistent values per day exist
    INPUT: 
    - data : xarray DataSet with dimension time
    - values_per_day : how values per day exist
    OUTPUT:
        xarray DataSet of daily mean data for each variable
        (hour of each day will be as hour of first time value in the DataSet)
    '''
    if data.time.size % values_per_day != 0:
        print("last day might miss as different values per day")
    DataSet_temporary = False
    for i in range(4):
        selection = np.arange(i,data.time.size, values_per_day)
        if not DataSet_temporary:
            DataSet_temporary = data.isel(time = selection).copy(deep = True)
            time_values = DataSet_temporary.time.values
        else:
            DataSet_temporary_old = DataSet_temporary.copy(deep = True)
            DataSet_temporary = data.isel(time = selection).copy(deep = True)
            DataSet_temporary.__setitem__("time", time_values)
            DataSet_temporary = DataSet_temporary + DataSet_temporary_old
    return DataSet_temporary/values_per_day

def calc_daily_max(data, variable = False, values_per_day = 4):
    '''
    NOTE:
    Calculates daily DataSet if consistent values per day exist
    INPUT: 
    - data : xarray DataSet with dimension time
    - variable : on which variable shall the maximum calculation be based on ?
    - values_per_day : how values per day exist
    OUTPUT:
        xarray DataSet of daily max data
        (hour of each day will be as hour of first time value in the DataSet)
    '''
    if data.time.size % values_per_day != 0:
        print("last day might miss as different values per day")
    DataSet_temporary = False
    for i in range(4):
        selection = np.arange(i,data.time.size, values_per_day)
        if not DataSet_temporary:
            DataSet_temporary = data.isel(time = selection).copy(deep = True)
            time_values = DataSet_temporary.time.values
        else:
            DataSet_temporary_old = DataSet_temporary.copy(deep = True)
            DataSet_temporary = DataSet_temporary* 0 # set to 0 and add the bigger value weherever it occures
            DataSet_temporary_new = data.isel(time = selection).copy(deep = True)
            DataSet_temporary_new.__setitem__("time", time_values)
            if not variable:
                DataSet_temporary = DataSet_temporary + DataSet_temporary_old.where(DataSet_temporary_new[variable] < DataSet_temporary_old[variable], other = 0)
                DataSet_temporary = DataSet_temporary + DataSet_temporary_new.where(DataSet_temporary_new[variable] >= DataSet_temporary_old[variable], other = 0)
            else :
                DataSet_temporary = DataSet_temporary + DataSet_temporary_old.where(DataSet_temporary_new < DataSet_temporary_old, other = 0)
                DataSet_temporary = DataSet_temporary + DataSet_temporary_new.where(DataSet_temporary_new >= DataSet_temporary_old, other = 0)
    return DataSet_temporary

def calc_hw_idx(data, window):
    threshold = 25
    mask_org = (data.t2m_mask_quantiles.isel(quantile= 0) + (data.t2m >= threshold))
    mask_new = (data.t2m_mask_quantiles.isel(quantile= 0) + (data.t2m >= threshold))
    mask_list = []

    for i in range(window):
        selection = np.arange(i,data.time.size, window)
        mask_list.append(mask_org.isel(time = selection).values)
    mask = None
    for mask in mask_list:
        if mask_res is None:
            mask_res = mask
        mask_res = mask & mask_res

def gufunc_linreg(first_samples,
                         second_samples,
                         dim='time'):
    '''
    xa.apply_ufunc linreg function gets 
    INupt:
    x,y data and the dimensino along which this linearregression shall be done
    Output
    dict with 
    slope, intercept, , r_value, p_value, std_err 
    '''
    
    def func(varx,vary) :
        mask = ~np.isnan(varx) & ~np.isnan(vary)
        if mask.sum() != 0:
            res = stats.linregress(varx[mask], vary[mask])
            return res
        else :
            return np.nan, np.nan, np.nan, np.nan, np.nan
    linreg_list = xr.apply_ufunc(func,
                       first_samples, second_samples,
                       input_core_dims=[[dim], [dim]],
                       output_core_dims=[[],[],[],[],[]],
                       dask = "parallelized",
                       vectorize=True,
                       output_dtypes = (float,float,float,float,float))
    return dict(slope = linreg_list[0], 
                intercept = linreg_list[1],
                r_value = linreg_list[2], 
                p_value = linreg_list[3], 
                std_err = linreg_list[4])


def gufunc_detrend(data, variable = "t2m" , dim = "time", mutable = True) :
    '''
    Calculates detrended values of "varaible" using gufunc_linreg
    INPUT: 
    - data : xarray DataSet containing the "variable"
    - variable : name of variable you want to detrend (default: "swbgt")
    - dim : dimension along which the trend shall be calculated
    - mutable = True : detrended data will be added to DataSet
    OUTPUT:
        1. dict with keywords slope, intercept, r_value, p_value, std_err 
        2. x-values ussed for calculation -> needed if change per time shall be analysed.
    '''
    y = data[variable]
    x = (y[dim] - y[dim][0]).astype(float)
    res = gufunc_linreg(x,y, dim = dim)
    if mutable :
        data[variable + '_detrend'] = data[variable] -(res['slope']*x + res['intercept'])
        return res, x
    else :
        return res, x
    
def quantile_calc(data, variable = "t2m", dim = "time", quantile_list = [0.9,0.95,0.98], mutable = True):
    '''
    Calculates mask quantile values of "varaible" along dim
    INPUT: 
    - data : xarray DataSet containing the "variable"
    - variable : name of variable quantile shall be calculated for
    - dim : dimension for calculatoion
    - mutable = True : quanitlemask data will be added to DataSet
    OUTPUT:
        data-array containing the quantile values
    '''
    if mutable :
        data[variable  + "_mask_quantiles"] = (data[variable] >= data[variable].quantile(quantile_list, dim = "time", keep_attrs = False)).astype(bool)
        data[variable  + "_mask_quantiles"].astype(bool)
        return None #data[variable].quantile(quantile_list, dim = "time")
    else :
        return data[variable].quantile(quantile_list, dim = "time")
    
    
def calc_mask_ydrunpctl(data, cdo_data, variable = "t2m"):
    '''
    This function uses the results from the cdo ydrunpctl.
    It creates a mask where the value for each calendar day is higher than the threshold given by cdo_data
    (data[d] >= cdo_data[d] , d = calender day)
    INPUT:
        data: DataSet containing the variable
        cdo_data : result if the cdo ydrunpctl function - Dataset containing the variable
        variable : variable which one wants to create the mask for
    OUTPUT :
        mask with: data[d] >= cdo_data[d] , d = calender day
    
    '''
    # get vector with corresponding months and days for the data
    months_data = data.time.dt.month
    days_data = data.time.dt.day
    
    # get vector with corresponding months and days for the cdo processed data
    time_of = cdo_data.time
    months_of = cdo_data.time.dt.month
    days_of = cdo_data.time.dt.day

    # create tuples with month and day for each timestep of the processed data
    tuples = []
    for of_step in time_of:
        tuples.append((of_step.dt.month.values, of_step.dt.day.values))
    
    # result will save the mask ->> create an fully False mask 
    result = (data[variable] * 0).astype(bool)
    for step in tuples:
        # for each tuple get the mask for data and cdo_data
        mask_data        = ((step[0] == months_data) & (step[1] == days_data))
        mask_cdo_data    = ((step[0] == months_of)   & (step[1] == days_of))
        # mask is the boolean array where the
        mask = data[variable].sel(time = mask_data).load() >= cdo_data[variable].sel(time = mask_cdo_data).values
        if mask.values.sum() > 0:
            result[mask_data] = result[mask_data]  + mask.values
        else :
            pass
    
    return result