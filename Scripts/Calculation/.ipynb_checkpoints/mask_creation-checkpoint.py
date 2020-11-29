#%%
import regionmask
import matplotlib.pyplot as plt
import xarray as xr
import numpy as np

# land = regionmask.defined_regions.natural_earth.land_110
# regions = regionmask.defined_regions.ar6.land.mask(data.longitude, data.latitude, lon_name="longitude", lat_name="latitude")

# data= xr.open_dataset(r"C:\Users\Nils Niebaum\Documents\Uni\Bachlor_thesis\Data\europe_all_1900-01-01_1930-12-31.nc")

def mask_land(data) :
    return regionmask.defined_regions.natural_earth.land_110.mask(data.longitude, data.latitude, lon_name="longitude", lat_name="latitude") == 0


def mask_regions(data) :
    return regionmask.defined_regions.ar6.land.mask(data.longitude, data.latitude, lon_name="longitude", lat_name="latitude")

def mask_region(data, region) :
    if region == None:
        return (mask_regions(data) == 16) + (mask_regions(data) == 17) + (mask_regions(data) == 18) + (mask_regions(data) == 19) 
    else:
        return mask_regions(data) == region

def calc_mask_ydrunpctl(data, cdo_data, variable = "t2m"):
    '''
    INPUT : 
    data: data mask shall be created of based on cdo_data
    cdo_data: xarraydataset created by cdo ydrunpctl which is made of 366 timesteps
    variabel: variable the calculation shaöö be done for
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
#         for each tuple get the mask for data and cdo_data
        mask_data        = ((step[0] == months_data) & (step[1] == days_data))
        mask_cdo_data    = ((step[0] == months_of)   & (step[1] == days_of))

#         print('calc mask')
        mask = data[variable].sel(time = mask_data).load() >= cdo_data[variable].sel(time = mask_cdo_data).values
        if mask.values.sum() > 0:
            result[mask_data] = result[mask_data]  + mask.values
        else :
            pass
    
    return result

# import cartopy.crs as ccrs

# def area_plot(data = None, ax = None, 
#         levels = np.arange(20,38,1), contour_levels = np.arange(20,40,2),
#         add_colorbar=False, cmap = "Reds",
#         landcolor = "lightgrey", oceancolor = "lightskyblue",
#         zorder_ocean = 2, zorder_land = 2) :
#     #create colormap
#     cmap = plt.get_cmap(cmap, len(levels))

#     # plot pcolormesh of the data with according colormap
#     pmesh = data.plot.pcolormesh(ax = ax, cmap = cmap, \
#         x = "longitude", y = "latitude",\
#         yincrease = True, transform = ccrs.PlateCarree(), \
#         # norm = mplcolors.BoundaryNorm(levels, ncolors=len(levels) -1, clip=False), \
#         zorder=3, add_colorbar= False)

#     ct = contourplot = data.plot.contour(ax= ax, levels= contour_levels, colors= 'k',\
#         linestyles = '-', linewidths = 0.5, x="longitude", y="latitude", \
#         yincrease=True, transform=ccrs.PlateCarree(),\
#         add_colorbar= False,\
#         zorder=4)

#     # plot gridlines and coastlines
#     ax.coastlines(resolution = 'auto', color = 'k', )
#     gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,\
#                 linewidth=1, color='k', alpha=0.5, linestyle=':',\
#                 zorder=8)
#     gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False,\
#                 linewidth=1, color='w', alpha=0.5, linestyle='-',\
#                 zorder=8)
#     # gl.xlabels_top = False
#     # gl.ylabels_right = False
#     gl.xlines = True
#     # gl.xlocator = mticker.FixedLocator(range(-30,50,10))
#     # gl.xformatter = LONGITUDE_FORMATTER
#     # gl.yformatter = LATITUDE_FORMATTER

#     # gl.xlabel_style = {'size': 15, 'color': 'gray'}
#     # gl.xlabel_style = {'color': 'red', 'weight': 'bold'}
#     # ax.add_feature(cart.feature.OCEAN, facecolor= oceancolor, zorder=zorder_ocean)
#     # ax.add_feature(cart.feature.LAND, facecolor= landcolor, zorder=zorder_land)


#     return pmesh, ct
# #%%


# #%%
# fig = plt.figure()
# ax = fig.add_subplot(111,\
#         projection=ccrs.NearsidePerspective(
#                         central_latitude=50.72,
#                         central_longitude=10.53,
#                         satellite_height=10000000.0))

# # area_plot(data= data.t2m.where(mask == 1).where(regions == 16).isel(time = 0) *0, ax = ax,
# #         levels = np.arange(-2, 2, 0.1),
# #         contour_levels= np.arange(0,100), cmap = "Reds")
# # area_plot(data= data.t2m.where(mask == 1).where(regions == 17).isel(time = 0)*0, ax = ax,
# #         levels = np.arange(-2, 2, 0.1),
# #         contour_levels= np.arange(0,100), cmap = "Greens")
# # area_plot(data= data.t2m.where(mask == 1).where(regions == 19).isel(time = 0)*0, ax = ax,
# #         levels = np.arange(-2, 2, 0.1),
# #         contour_levels= np.arange(0,100), cmap = "Oranges")
# # area_plot(data= data.t2m.where(mask == 1).where(regions == 18).isel(time = 0)*0, ax = ax,
# #         levels = np.arange(-2, 2, 0.1),
# #         contour_levels= np.arange(0,100), cmap = "Blues")
# ax.coastlines(resolution='110m', zorder=20)

# ax.gridlines()
# area_plot(regions.where(mask == 1).where(regions <= 21), ax = ax,
#         levels = [15,16,17,18,19,20],
#         contour_levels= None, cmap = "twilight")
# plt.show()
# # %%

# # calculate russon rolling seasonal mean
# mean_data = data.isel(latitude =0, longitude = 10).rolling(dict(time = 30), center= True).mean() #.groupby("time.year").mean()



# timeslice = slice("1903-01-01", "1904-01-01")
# data.isel(latitude =0, longitude = 10).t2m.sel(time= timeslice).plot()
# mean_data.t2m.sel(time= timeslice).plot()

# # %%
# res =   ((data.time <= np.datetime64("1900-01-05")) & \
#         (data.time > np.datetime64("1900-01-04"))) + \
#         ((data.time <= np.datetime64("1901-01-05")) & \
#         (data.time > np.datetime64("1901-01-04"))) + \
#         ((data.time <= np.datetime64("1902-01-05")) & \
#         (data.time > np.datetime64("1902-01-04")))