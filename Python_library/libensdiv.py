"""
Diverse functions for ensembles.
"""

print(f"Name: {__name__}")
print(f"Package: {__package__}")

print("This is a collection of diverse functions for ensemble.")

import numpy as np

def fspatialmean(data,mask,e1,e2):
    '''''
    This function computes the spatial mean of a variable. 
    ...
    The arguments needed are the following:
    - data = data of a variable of one member of an experiment (DataArray with at least "x" and "y" as dimensions)
    - mask = the mask to use 
    /!\ should be 2D with the same dimensions as the variable of interest
    /!\ should only contain the 0 or 1 values
    - e1 and e2 = horizontal mesh sizes
    ...
    /!\ the mask, e1, and e2 should be given at the correct point according to the variable (possibilities for Arakawa C-grid: T, U, V, W, F)
    '''''
    mean=(data*mask*e1*e2).sum(("y","x"))/(mask*e1*e2).sum(("y","x"))
    return mean

def maskpt(masktot,maskdomain,pt):
    '''''
    This function gives the mask, e1 and e2 at the right point (Arakawa C-grid)
    ...
    The arguments needed are the following:
    - masktot = dataset containing e1 and e2 (horizontal mesh sizes), should have 3 dimensions: "time" (or equivalent), "y", and "x"
    - maskdomain = dataset containing masks of the domain (named tmask, umask, vmask)
    - pt = the Arakawa C-grid point of interest
    /!\ masktot and maskdomain should have the same dimensions
    ...
    The output is as follow: (mask, e1, e2)
    '''''

    if pt=="T": #if the variable is located at T-point
        #print(pt,True,"T")
        mask=maskdomain.tmask
        e1=masktot.e1t[0,:,:]
        e2=masktot.e2t[0,:,:]
    
    elif pt=="U": #if the variable is located at U-point
        #print(pt,True,"U")
        mask=maskdomain.umask
        e1=masktot.e1u[0,:,:]
        e2=masktot.e2u[0,:,:]

    elif pt=="V": #if the variable is located at V-point
        #print(pt,True,"V")
        mask=maskdomain.vmask
        e1=masktot.e1v[0,:,:]
        e2=masktot.e2v[0,:,:]
    
    return (mask,e1,e2)

def fmean_mb(data,var,nb_member,n):
    '''''
    This function computes the mean over all members (except one if asked) without ponderation.
    ...
    data = all the members of an experiment; it should be a list of members, with each member as a DataSet/DataArray
    var = variable on which to compute the mean; var=False if there is no variable defined
    nb_member = total number of members
    n = member to ignore (must be an integer); if n=-1 all members are considered
    '''''
    isum=0.0
    for i in range(nb_member):
        if i!=n:
            if var==False:
                isum+=data[i]
            else:
                isum+=data[i][var]
    ntot=nb_member-(n!=-1)*1.0
    #print(ntot)
    return isum/ntot

def fstd_mb(data,var,nb_member,n):
    '''''
    This function computes the standard deviation (unbiased definition) over all members (except one if asked).
    ...
    data = all the members of an experiment; it should be a list of members, with each member as a DataSet/DataArray
    var = variable on which to compute the std; var=False if there is no variable defined
    nb_member = total number of members
    n = member to ignore (must be an integer); if n=-1 all members are considered
    '''''
    mean=fmean_mb(data,var,nb_member,n)
    isum=0.0
    for i in range(nb_member):
        if i!=n:
            if var==False:
                isum+=(data[i]-mean)**2
            else:
                isum+=(data[i][var]-mean)**2
    ntot=nb_member-(n!=-1)*1.0
    #print(ntot)
    return np.sqrt(isum/(ntot-1))
