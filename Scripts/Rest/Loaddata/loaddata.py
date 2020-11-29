
import xarray as xr
import time


def loaddata(folderpath= r"I:\Bachlor_thesis\Data", filepath):
        
        found = False
        while not found:
        try : # if file exists load the data
                data = xr.open_dataset(save_dict["savefolder"] + save_dict["savename"])
                print("=======\nload data done")
                try :
                        data.t2m[0]
                        print("t2m found" )
                        data.d2m[0]
                        print("d2m found" )
                        data.swbgt[0]
                        print("swbgt found" )

                        found = True
                except :
                        print("time_since_start missing\n => calc new")
                        found = False
        except: 
                found = False
        
        if not found :
                #%% ==================== CALCULATION OF SWBGT ================
                print("=======\nfile not found => calculate and save data")
                folderpath = r"I:\Bachlor_thesis\Data"

                st = time.time()
                t2m = xr.open_dataset(folderpath + r"\t2m_era20c_1900-2010.nc")#, chunks={'time': 12*4*365})
                print("load temp took {:.2f}s".format(time.time()-st))
                st = time.time()
                d2m = xr.open_dataset(folderpath + r"\td2m_era20c_1900-2010.nc")#, chunks={'time': 12*4*365})
                print("load dew took {:.2f}s".format(time.time()-st))
                st = time.time()
                all_data = xr.merge([t2m, d2m])
                print("---------------input data\n",all_data)

                data = all_data.sel(time = timeslice, latitude= latitude_slice, longitude= longitude_slice)

                # #%% ========== CHUNK DATA ==============
                # rechunk data for perfomance
                rechunk_dim = {
                        "time" : int(np.shape(data.t2m)[0] / 8),
                        }
                # !!!!!!!!!!!!! USE  °C !!!!!!!!!!!!!!!!!!!!
                data = data.chunk(rechunk_dim) -273.15
                # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                data.resample(time="1d").mean() # resample data
                data.t2m.attrs = {"full name" : "temperature", "unit" : "°C"}
                data.d2m.attrs = {"full name" : "dew point temperature", "unit" : "°C"}
                swbgt = humidcalc.swbgt_func(data.t2m, data.d2m)
                swbgt.name = "swbgt"
                swbgt.attrs = {"place" : city, "full name" : "simple wet bulb globe temperature", "unit" : "°C"}
                print("swbgt calc done")
                st = time.time()

                data = xr.merge([swbgt, data])
                print(data)
                input("ok?")
                data.to_netcdf(save_dict["savefolder"] + save_dict["savename"])
                print("swbgt saved {:.2f}s".format(time.time()-st))
