import xarray as xr
import matplotlib.pyplot as plt
import time
import numpy as np

def load_select_data(folderpath, timeslice) :
    #%% load data 
    folderpath = r"I:\Bachlor_thesis\Data"

    # load temperature
    st = time.time()
    temp = xr.open_dataset(folderpath + r"\t2m_era20c_1900-2010.nc")#, chunks={'time': 12*4*365})
    print("=======\nload temp: {:.2f}s".format(time.time()-st))
    # load dew point temperature
    st = time.time()
    temp_dew = xr.open_dataset(folderpath + r"\td2m_era20c_1900-2010.nc")#, chunks={'time': 12*4*365})
    print("=======\nload dew: {:.2f}s".format(time.time()-st))
    # merge both
    st = time.time()
    all_data = xr.merge([temp, temp_dew])
    print("=======\nmerge: {:.2f}s".format(time.time()-st))
    # print data
    print(all_data)

    # select by time slice
    st = time.time()
    data = all_data.sel(time= timeslice) - 273.15
    print("=======\nselect timeslice {} till {}: {:.2f}s".format(timestart, timeend, time.time()-st))
    print(data, '\n')
    return data
