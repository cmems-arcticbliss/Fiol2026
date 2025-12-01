"""
Diverse functions to work with Lagrangian trajectories.
"""

from libensdiv import fmean_mb
import numpy as np
import cartopy.crs as ccrs

print(f"Name: {__name__}")
print(f"Package: {__package__}")

print("This is a collection of diverse functions to work with Lagrangian trajectories.")

def projGeo2Cartesian(lat,lon):
    '''''
    This function transforms geographic coordinates (lat,lon) in degrees into Cartesian (x,y) in km with RGPS' 'NorthPolarStereo' projection.
    It is the same projection as the one used in sitrack (https://github.com/brodeau/sitrack).
    lat = latitude values
    lon = longitude values
    Warning: lat and lon should just contain values (not a DataArray or equivalent).
    '''''
    crs_src=ccrs.PlateCarree()
    crs_trg =ccrs.NorthPolarStereo(central_longitude=-45., true_scale_latitude=70.)

    zx,zy,_ = crs_trg.transform_points(crs_src, lon, lat).T
    
    return np.array([ zx/1000., zy/1000. ]).T #dividing by 1000 to convert into km

def fbuoys_alive(data):
    '''''
    This function computes the time series of the number of buoys alive 
    and the buoys that are no longer alive at the end of the tracking.
    ...
    data = a trajectory file of a member of an experiment 
    ...
    The function is designed to work with trajectory files generated with sitrack (https://github.com/brodeau/sitrack).
    '''''
    time_serie=data.mask.sum("buoy")
    nb_buoys=data.sizes["buoy"]
    if time_serie[-1]==nb_buoys: #all the buoys are alive at the end
        return (time_serie,[])
    else:
        buoys_disregarded=np.where(data.mask[-1,:]==0)[0] #storing of the buoys that are not alive until the end
        return (time_serie,buoys_disregarded)

def ffind_gaps(mask,time):
    '''''
    Given a 1D mask, this function identifies where there are gaps (mask=0) and the time lag between each gap. 
    The 1D mask should be a time series for this function to make sense.
    ...
    Inputs:
    mask = a 1D mask array equal to 1 if the data exists and to 0 otherwise (it should just contain values, not a DataArray or equivalent)
    time = 1D time array associated with the mask; it should have a dtype = datetime64[ns]
    ...
    Outputs:
    pos_end_gap = list containing the index of the data following gaps
    time_lag = associated list of the time lag induced by the temporal gaps
    '''''
    
    pos_end_gap=np.where((mask[1:]-mask[:-1])==1)[0]+1
    
    #test of if there is an open gap at the beginning (we do not consider open gaps)
    if mask[0]==0:
        pos_end_gap=pos_end_gap[1:]

    #computation of the time lag
    time_lag=[] #initialisation of the time lag list

    #finding where are the last data before the gaps
    for it in pos_end_gap: #iteration over the gaps
        cpt=1
        while mask[it-cpt]==0:
            cpt+=1

        time_lag.append((time[it]-time[it-cpt])/np.timedelta64(1,'h')) #computation of the time lag

    return (pos_end_gap,time_lag)

def fdistance(x1,y1,x2,y2):
    '''''
    This function computes the Euclidean distance between (x1,y1) and (x2,y2).
    x1, y1, x2, and y2 could be a singular position or an array.
    '''''
    return np.sqrt((x1-x2)**2+(y1-y2)**2)