#%% PERCEBTILE OF SWBGT
import os
# os.chdir(r"C:\Users\Nils Niebaum\Documents\Uni\Bachlor_thesis\Prog\Scripts")
import xarray as xr
import time
import numpy as np
import Calculation.humiditycalculation as humidcalc
import Calculation.trend_func as trends
from memory_profiler import profile

#%% =============== LOAD DATA ================== 

@profile
def save_calc_swbgt(save_dict, load_dict, save = True) :

        timeslice = load_dict["time_slice"]
        latitude_slice = load_dict["latitude_slice"]
        longitude_slice = load_dict["longitude_slice"]

        @profile
        def __main__(save = True):
                try : # if file exists load the data
                        data = load_DataSet_and_add(save_dict = save_dict, load_dict = load_dict)
                        data, added = add_missing_variables(data = data, save_dict = save_dict, load_dict = load_dict)
                        if save and added: save_data(data = data, save_dict = save_dict, load_dict = load_dict)
                except :
                        # if not found calc and save data 
                        print("file not complete\n => calc swbgt and create new file")
                        data = load_t2m_d2m_calc_swbgt(save_dict = save_dict, load_dict = load_dict)
                        data, added = add_missing_variables(data = data, save_dict = save_dict, load_dict = load_dict)
                        if save: 
                                save_data(data = data, save_dict = save_dict, load_dict = load_dict)
                return data
        # exceute __main__ and return data
        
        # if longitude_slice.stop > 180:
        #         longitude_slice = slice(longitude_slice.stop - 360, )

        # if the slice start longitude is west of 0° than use west and east data and combine later
        if longitude_slice.start < 0 :
                # east of 0° Longitude
                data_east = __main__(save = False)
                print("=============\nDATA EAST",data_east)
                input("ok?")
                print("calculated east side")
                # west of 0° Longitude
                load_dict["longitude_slice"] = slice(longitude_slice.start + 360, 360)
                data_west = __main__(save = False)
                print("=============\nDATA WEST",data_west)
                # turn longitude on west into west values (350° => -10°)
                data_west.__setitem__("longitude", data_west.longitude - 360)
                print("calculated east side")
                print("=============\nDATA WEST",data_west)
                input("ok?")
                data_all = xr.merge([data_east, data_west])
                print("=============\nDATA ALL",data_all)
                if save: 
                        save_data(data = data_all, save_dict = save_dict, load_dict = load_dict)
                return data_all
        else:
                return __main__()


def load_t2m_d2m_calc_swbgt(save_dict, load_dict) :
        ''' CALCULATION OF SWBGT '''

        timeslice = load_dict["time_slice"]
        latitude_slice = load_dict["latitude_slice"]
        longitude_slice = load_dict["longitude_slice"]
        quantiles = save_dict["quantiles"]
        deseason = save_dict["deseason"]
        detrend = save_dict["detrend"]
        print("=======\ncalculate and save data")
        # -------- load data from path --------
        try :
                st = time.time()
                t2m = xr.open_dataset(load_dict["t2m_path"])#, chunks={'time': 12*4*365})
                print("load temp took {:.2f}s".format(time.time()-st))
                st = time.time()
                d2m = xr.open_dataset(load_dict["d2m_path"])#, chunks={'time': 12*4*365})
                print("load dew took {:.2f}s".format(time.time()-st))
                all_data = xr.merge([t2m, d2m])
                print("-------\ninput data\n",all_data)
        except: 
                raise KeyError
        # -------- slice data based on input --------
        org_data = all_data.sel(time = timeslice, latitude= latitude_slice, longitude= longitude_slice)
        print("-------\nsliced data\n", org_data)
        # input("-------\nsliced data ok?\n")
        #     input("-------\nsliced data ok?\n")

        # -------- rechunk data if wanted --------
        if save_dict["rechunked"]:
                print("-------\nuse chunked dask array\n")
                rechunk_dim = {
                "time" : int(np.shape(org_data.t2m)[0] / 8),
                }
                # ======== USE  °C  ========
                org_data = org_data.chunk(rechunk_dim) -273.15
        else:
                # ======== USE  °C  ========
                org_data = org_data -273.15
        
        # -------- calc daily mean or max if wanted ------------
        if save_dict["daily"] == "mean":
                print("----\ncalculate daily mean")
                st = time.time()
                org_data = trends.calc_daily_mean(data= org_data)
                print("took {:.2f}s".format(time.time()-st))
        elif save_dict["daily"] == "max":
                print("----\ncalculate daily mean")
                st = time.time()
                org_data = trends.calc_daily_mean(data= org_data)
                print("took {:.2f}s".format(time.time()-st))
        print(org_data)
        # input("is it daily ???\n")
        # -------- set attrs --------
        org_data.t2m.attrs = {"full name" : "temperature",
                        "units" : "degC",
                        "place" : save_dict["location"]}
        org_data.d2m.attrs = {"full name" : "dew point temperature",
                        "units" : "degC",
                        "place" : save_dict["location"]}
        # ======== CALC SWBGT ========
        swbgt = humidcalc.swbgt_func(org_data.t2m, org_data.d2m)
        swbgt.name = "swbgt"
        swbgt.attrs = {"place" : save_dict["location"], 
                        "full name" : "simplified wet bulb globe temperature",
                        "units" : "degC",
                        "note" : "unitless but unit needed if xclim needs to be used"}
        print("------ \ncalculated swbgt")
        
        # -------- merge t2m, d2m and swbgt --------
        data = xr.merge([swbgt, org_data])
        return data

def load_DataSet_and_add(save_dict, load_dict) :
        data = xr.open_dataset(save_dict["savefolder"] + save_dict["savename"])
        print("=======\nload data done")
        print(data)
        try :
                data.t2m[0]
                print("t2m found")
                data.d2m[0]
                print("d2m found")
                data.swbgt[0]
                print("swbgt found")
                
                # -------- calc daily mean or max if wanted ------------
                if save_dict["daily"]:
                        print(data.time[0:5])
                        doit = input("is it daily data?\n y/n\n") is not "y"
                        if not doit:
                                print("sorry 6h data is required")
                if doit and save_dict["daily"] is "mean":
                        print("----\ncalculate daily MEAN")
                        st = time.time()
                        data = trends.calc_daily_mean(data= data)
                        print("took {:.2f}s".format(time.time()-st))
                elif doit and save_dict["daily"] is "max":
                        print("----\ncalculate daily MAX")
                        st = time.time()
                        data = trends.calc_daily_max(data= data)
                        print("took {:.2f}s".format(time.time()-st))

                # if t2m, d2m and swbgt found: calc wanted parameters
                print("-------\nadd variables")
                data, added = add_missing_variables(data = data, save_dict = save_dict, load_dict = load_dict)
                print("-------\nsave data")
                if save and added : save_data(data = data, save_dict = save_dict, load_dict = load_dict)
                return data
        except: 
                return data
def save_data(data, save_dict, load_dict) :
        
        timeslice = load_dict["time_slice"]
        latitude_slice = load_dict["latitude_slice"]
        longitude_slice = load_dict["longitude_slice"]
        quantiles = save_dict["quantiles"]
        deseason = save_dict["deseason"]
        detrend = save_dict["detrend"]

        st = time.time()
        try :
                data.to_netcdf(save_dict["savefolder"] + save_dict["savename"])
        except :
                # if not overwritable create NEW_....nc file
                data.to_netcdf(save_dict["savefolder"] + 'NEW_' + save_dict["savename"], mode = "w")
        print("all saved {:.2f}s".format(time.time()-st))

@profile
def add_missing_variables(data, save_dict, load_dict = None) :
        ''' create the wanted varaibles if wanted '''

        quantiles = save_dict["quantiles"]
        deseason = save_dict["deseason"]
        detrend = save_dict["detrend"]
        added = False

        # -------- create quantile mask if wanted --------
        if quantiles is not False:
                variable_list = quantiles["variables"]
                quantile_list = quantiles["quantiles"]
                for variable in variable_list :
                        try :
                                data[variable  + "_mask_quantiles"][0]
                                print("mask_swbgt_quantiles found")
                        except:
                                data[variable  + "_mask_quantiles"] = data[variable] >= data[variable].quantile(quantile_list)
                                print("quaniles mask added")
                                added = True
        else :
                print('no quantile')
        # -------- create deseason if wanted --------
        
        if deseason is not False:
                print(deseason)
                variable_list = deseason["variables"]
                groupby = deseason["groupby"]
                rolling = deseason.setdefault('rolling', False)
                for variable in variable_list :
                        try :
                                data[variable  + "_deseason"][0]
                                print(variable + "_deseason found")
                        except:
                                trends.deseason(data= data, variable= variable, groupby= groupby, rolling = rolling, mutable=True)
                                print(variable + "_deseason added")
                                added = True
        else :
                print('no deseason')
        # -------- create detrend if wanted --------+

        if detrend is not False:
                print(detrend)
                variable_list = detrend["variables"]
                for variable in variable_list :
                        try :
                                data.swbgt_deseason_detrend[0]
                                print("swbgt_deseason_detrend found")
                        except:
                                trends.detrend(data= data, variable= "swbgt_deseason")
                                print("deseasoned data added")
                                added = True
        else :
                print('no detrend')

        return data, added

        # try to open dataset and and missing values


#%%
# def save_calc_swbgt(save_dict, load_dict) :

#         timeslice = load_dict["time_slice"]
#         latitude_slice = load_dict["latitude_slice"]
#         longitude_slice = load_dict["longitude_slice"]
#         quantiles = save_dict["quantiles"]
#         deseason = save_dict["deseason"]
#         detrend = save_dict["detrend"]


#         def load_t2m_d2m_calc_swbgt() :
#                 ''' CALCULATION OF SWBGT '''

#                 print("=======\ncalculate and save data")
#                 # -------- load data from path --------
#                 try :
#                         st = time.time()
#                         t2m = xr.open_dataset(load_dict["t2m_path"])#, chunks={'time': 12*4*365})
#                         print("load temp took {:.2f}s".format(time.time()-st))
#                         st = time.time()
#                         d2m = xr.open_dataset(load_dict["d2m_path"])#, chunks={'time': 12*4*365})
#                         print("load dew took {:.2f}s".format(time.time()-st))
#                         all_data = xr.merge([t2m, d2m])
#                         print("-------\ninput data\n",all_data)
#                 except: 
#                         raise KeyError
#                 # -------- slice data based on input --------
#                 org_data = all_data.sel(time = timeslice, latitude= latitude_slice, longitude= longitude_slice)
#                 print("-------\nsliced data\n", org_data)
#                 # input("-------\nsliced data ok?\n")
#                 #     input("-------\nsliced data ok?\n")

#                 # -------- rechunk data if wanted --------
#                 if save_dict["rechunked"]:
#                         print("-------\nuse chunked dask array\n")
#                         rechunk_dim = {
#                         "time" : int(np.shape(org_data.t2m)[0] / 8),
#                         }
#                         # ======== USE  °C  ========
#                         org_data = org_data.chunk(rechunk_dim) -273.15
#                 else:
#                         # ======== USE  °C  ========
#                         org_data = org_data -273.15
                
#                 if save_dict["daily"] :
#                         print("----\ncalculate daily mean")
#                         st = time.time()
#                         org_data = trends.calc_daily_mean(data= org_data)
#                         print("took {:.2f}s".format(time.time()-st))
#                 print(org_data)
#                 # input("is it daily ???\n")
#                 # -------- set attrs --------
#                 org_data.t2m.attrs = {"full name" : "temperature",
#                                 "units" : "degC",
#                                 "place" : save_dict["location"]}
#                 org_data.d2m.attrs = {"full name" : "dew point temperature",
#                                 "units" : "degC",
#                                 "place" : save_dict["location"]}
#                 # ======== CALC SWBGT ========
#                 swbgt = humidcalc.swbgt_func(org_data.t2m, org_data.d2m)
#                 swbgt.name = "swbgt"
#                 swbgt.attrs = {"place" : save_dict["location"], 
#                                 "full name" : "simple wet bulb globe temperature",
#                                 "units" : "degC",}
#                 print("------ \ncalculated swbgt")
                
#                 # -------- merge t2m, d2m and swbgt --------
#                 data = xr.merge([swbgt, org_data])
#                 return data

#         def save_data(data) :
#                 st = time.time()
#                 try :
#                         data.to_netcdf(save_dict["savefolder"] + save_dict["savename"])
#                 except :
#                         # if not overwritable create NEW_....nc file
#                         data.to_netcdf(save_dict["savefolder"] + 'NEW_' + save_dict["savename"], mode = "w")
#                 print("all saved {:.2f}s".format(time.time()-st))

#         def add_missing_variables(data) :
#                 ''' create the wanted varaibles if wanted '''

#                 added = False

#                 # -------- create quantile mask if wanted --------
#                 if quantiles is not False:
#                         variable_list = quantiles["variables"]
#                         quantile_list = quantiles["quantiles"]
#                         for variable in variable_list :
#                                 try :
#                                         data[variable  + "_mask_quantiles"][0]
#                                         print("mask_swbgt_quantiles found")
#                                 except:
#                                         data[variable  + "_mask_quantiles"] = data[variable] >= data[variable].quantile(quantile_list)
#                                         print("quaniles mask added")
#                                         added = True
#                 else :
#                         print('no quantile')
#                 # -------- create deseason if wanted --------
                
#                 if deseason is not False:
#                         print(deseason)
#                         variable_list = deseason["variables"]
#                         groupby = deseason["groupby"]
#                         for variable in variable_list :
#                                 try :
#                                         data[variable  + "_deseason"][0]
#                                         print(variable + "_deseason found")
#                                 except:
#                                         trends.deseason(data= data, variable= variable, groupby= groupby,mutable=True)
#                                         print(variable + "_deseason added")
#                                         added = True
#                 else :
#                         print('no deseason')
#                 # -------- create detrend if wanted --------+

#                 if detrend is not False:
#                         print(detrend)
#                         variable_list = detrend["variables"]
#                         for variable in variable_list :
#                                 try :
#                                         data.swbgt_deseason_detrend[0]
#                                         print("swbgt_deseason_detrend found")
#                                 except:
#                                         trends.detrend(data= data, variable= "swbgt_deseason")
#                                         print("deseasoned data added")
#                                         added = True
#                 else :
#                         print('no detrend')

#                 return added

#         # try to open dataset and and missing values
#         def __main__(save = True) :
                
#                 try : # if file exists load the data
#                         data = xr.open_dataset(save_dict["savefolder"] + save_dict["savename"])
#                         print("=======\nload data done")
#                         print(data)
#                         try :
#                                 data.t2m[0]
#                                 print("t2m found")
#                                 data.d2m[0]
#                                 print("d2m found")
#                                 data.swbgt[0]
#                                 print("swbgt found")
#                                 # if t2m, d2m and swbgt found: calc wanted parameters
#                                 print("-------\nadd variables")
#                                 added = add_missing_variables(data = data)
#                                 print("-------\nsave data")
#                                 if save and added : save_data(data)
#                                 return data
#                         except: 
#                                 None
#                 except :
#                         # if not found calc and save data 
#                         print("file not complete\n => calc swbgt and create new file")
#                         data = load_t2m_d2m_calc_swbgt()
#                         added = add_missing_variables(data = data)
#                         if save and added : save_data(data = data)
#                 return data
#         # exceute __main__ and return data
        
#         # if longitude_slice.stop > 180:
#         #         longitude_slice = slice(longitude_slice.stop - 360, )

#         # if the slice start longitude is west of 0° than use west and east data and combine later
#         # if longitude_slice.start < 0 :
#         #         longitude_slice_temp = slice(longitude_slice.start + 360, 360)
#         #         data_east = __main__(save = False)
#         #         longitude_slice = longitude_slice_temp
#         #         data_west = __main__(save = False)
#         #         data_west.roll(longitude = 180).__setitem__("longitude", data.west.longitude - 180)
#         #         data = xr.merge
#         return __main__()





# %%
