
#%%
import numpy as np
import xarray as xr

def calc_heatwave_index(mask_original, duration = 1, as_DataArray = True):
        '''
        This function calculates the mask of heatwave days with a heat wave duration of multiple days (default 3)
        INPUT:
        - mask_original (bool) : boolean array as xarray dataarray
        - duration (int) : how many timesteps (not days or whatever) the mas needs to be True
        - as_DataArray : True returns a xarry datarray
        OUTPUT:
        - mask array with True for each timestep that is part of a heatwave
        
        NOTE:
        1. how does the calculation work here:
        (1 - True, 0 - False)
        example:
        mask_original = [1,0,1,1,1,0,0,1,1,1,1,0]
        duraiotn = 3
        1. for loop: mask_temporary
                it will create the following temporary mask_index (start index always moved to the right.)
                        [1,0,1,1,1,0,0,1,1,1]           [1,0,1,1,1,0,0,1,1,1]
                        [0,1,1,1,0,0,1,1,1,1]             [0,1,1,1,0,0,1,1,1,1]
                        [1,1,1,0,0,1,1,1,1,0]               [1,1,1,0,0,1,1,1,1,0]
                        - - - - - - - - - -     # with & operator the following will be done:
        mask_temporary: [0,0,1,0,0,0,0,1,1,0]   # if a heatwave occures, one can see that all values of the lists on the left side are 1 above each other
        
                but the mask_temporary is shorter ??
                -> for each value with 1 in mask_temporary the following 3 days (duration = 3) are part of the heatwave.
        
        2. for loop: =============
        mask_result     [0,0,0,0,0,0,0,0,0,0,0,0]       # mask_result starts with 0 ( need a array with the size to work on)
                         - - - - - - - - - - - -
        mask_temporary  [0,0,1,0,0,0,0,1,1,0]           # here the mask_temporary will be moved to the right for each itteration:
                ""        [0,0,1,0,0,0,0,1,1,0]         # now the + opertaor will be used (1+1=1, 1+0=1 , 0+0=0) :
                ""          [0,0,1,0,0,0,0,1,1,0]       # all 3 days after the 1 in mask_ temoprary will change to 1, as they are part of a heatwave.
                         - - - - - - - - - - - - 
        mask_result:    [0,0,1,1,1,0,0,1,1,1,1,0]       # you should be convinced that this works :D
        mask_original : [1,0,1,1,1,0,0,1,1,1,1,0]
        '''
        # duration = 3
        time_length = mask_original.time.size
        mask_temporary = None
        for index in range(duration):
                # the selection always chooses a slice of the mask along time.
                # the selection always hat the same length but changing starting position.
                selection = np.arange(index, time_length - duration + index, 1)
                # the mask_orginal will be sliced along the time axis based on the selection index vauels e.g. [0,1,2,3,4,...., time_length - duration]
                mask_index = mask_original.isel(time = selection).values
                if mask_temporary is None :
                        mask_temporary = mask_index
                # the & operator is used to get True for each day where 
                # afterwards at least for duration the mask_origianl is True (see Note) 
                mask_temporary = mask_temporary & mask_index

        mask_result = mask_original.values == -99 # this just creates a mask with False everywhere (see Note)
        for index in range(duration):
                # same slice as in 1. for loop
                selection = np.arange(index, time_length - duration + index, 1)
                # teh + operator will be used to identifiy each day ehich is part of the heat wave (see Note)
                mask_result[selection,:,:] = mask_result[selection,:,:] + mask_temporary

        if as_DataArray: # return a xarray DataArray with same coords and dimensions 
                return xr.DataArray(mask_result, dims = mask_original.dims, coords = mask_original.coords)
        else :  # else numpy.array will be returned
                return mask_result

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
    
