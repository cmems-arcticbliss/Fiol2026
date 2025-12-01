"""
Functions to look at the sea ice edge.
"""

print(f"Name: {__name__}")
print(f"Package: {__package__}")

print("This is a collection of functions to look at the sea ice edge.")

import numpy as np
import xarray as xr
import libdiv
    
def fmask_ce(data,ce):
    '''''
    This function computes a mask equal to 1 where data>ce and 0 elsewhere. 
    ...
    The arguments needed are the following:
    - data = data of a variable of one member of an experiment (a DataArray)
    - ce = the condition
    /!\ data should have 2 (or 3) dimensions in the following order: ("time_counter",) "y", "x".
    '''''
    size_y=data.sizes["y"] ; size_x=data.sizes["x"] #horizontal dimensions sizes
    nb_dim=len(np.shape(data))
    
    if nb_dim==2: #2D
        mask=xr.DataArray(data=np.zeros((size_y,size_x),dtype=int), dims=["y","x"])+1 #initialization
    elif nb_dim==3: #3D
        size_t=data.sizes["time_counter"] #time dimension size
        mask=xr.DataArray(data=np.zeros((size_t,size_y,size_x),dtype=int),\
                          dims=["time_counter","y","x"])+1 #initialization
    else:
        print("problem in the number of dimensions")
    
    mask=mask.where(data.values>ce,0) #the mask is equal to 0 where data<=ce
    return mask

def fproba_ce(data,var,ce,n=-1):
    '''''
    This function computes the frequency over all the members (except one if n!=-1) of the event data>ce.
    ...
    The arguments needed are the following:
    - data = data of an experiment /!\ to work should have 3 dimensions in that order: "time_counter", "y", "x".
    data should be a list of all the members, with members as a DataSet
    - var = name of the variable of interest 
    - ce = the condition
    - n = the member to exclude; if you don't want to exclude any member put n=-1
    By default it computes the frequency over all members.
    '''''
    size_y=data[0].sizes["y"] ; size_x=data[0].sizes["x"] #horizontal dimension sizes (all the members have the same sizes)
    size_t=data[0].sizes["time_counter"] #time dimension size (all the members have the same size)
        
    proba=xr.DataArray(data=np.zeros((size_t,size_y,size_x)), dims=["time_counter","y","x"]) #initialization

    nb_member=len(data) #number of members of data
    
    for i in range(nb_member): #for all members
        if i!=n: #except n
            mask=fmask_ce(data[i][var],ce) #computation of a mask equal to 1 where data[i][var]>ce and 0 elsewhere
            proba+=mask
    if n==-1:
        return proba/nb_member #return a probability (frequency) array of having the value of var strictly superior to ce, taking into account all members
    else:
        return proba/(nb_member-1) #return a probability (frequency) array of having the value of var strictly superior to ce, without taking into account member n

def fSPS_IIEE(data,type_ref,data_ref,var,ce,e1,e2,corners=False,width=False,height=False):
    '''''
    This function computes the SPS for the ice edge. 
    If specified, the SPS is computed over the domain defined with the indexes:
    [corners[1]:corners[1]+height+1,corners[0]:corners[0]+width+1]
    The reference chosen to compute the SPS is indicated by type_ref. There are 3 possible kinds of reference:
    - same_ens: each member of data is taken one after another as the reference
    - other_ens: each member of data_ref is taken one after another as the reference
    - masked_field: the reference is a 2D (or 3D if there is a time dimension) field with values only between 0 and 1.
    /!\ the reference given should use the same ce criterion to define the ice edge
    This last type of reference could be processed observations, a probability map of another ensemble, etc.
    ...
    The arguments needed are the following:
    - data = data of the experiment of which we want to compute the score
    data should be a list of all the members, with members as a DataSet
    - type_ref = character string indicating the kind of reference that is data_ref (3 possibilities; see above)
    - data_ref = data to take as the reference should be coherent with type_ref 
    if type_ref="same_ens" could be equal to False
    if type_ref="other_ens" data_ref should be a list of all the members, with members as a DataSet
    - var = name of the variable of interest over which to compute the score
    - ce = the condition/threshold to compute the score
    - e1, e2 = horizontal mesh sizes (should be 2D with the same dimensions "y" and "x" as data and data_ref!)
    - corners = position of the lower left corner. A tuple as follows: (x position, y position).
    If corners=False computation over all the domain (default).
    - width = length of the domain in the x-direction (default=False)
    - height = length of the domain in the y-direction (default=False)
    /!\ This function is designed to work with NEMO outputs!
    '''''
    #saving the sizes of the dimensions "x" and "y" of data
    size_y=data[0].sizes["y"] ; size_x=data[0].sizes["x"]
    
    #computing the mask of the domain over which to compute the SPS
    if corners: #computation over the domain defined by fmask_square_domain(size_x,size_y,corners,width,height)
        print("computation over the specified domain")
        mask_domain=libdiv.fmask_square_domain(size_y,size_x,corners,width,height)
    else: #computation over all the domain
        mask_domain=1
    
    SPS=[] #list that will stock all the SPS
    SPS_id=[type_ref,[]] #list that will indicate the kind of reference and more information about SPS

    if type_ref=="same_ens": #members of data are taken one after another as the reference
        nb_member=len(data) #number of members of data
        for n in range(nb_member): #iteration over the members of data
            SPS_id[1].append("ref mber: "+str(n+1))
            proba=fproba_ce(data,var,ce,n)*mask_domain #compute the probability (frequency) array of having the value of var strictly superior to ce, without taking into account member n
            proba_ref=fmask_ce(data[n][var],ce)*mask_domain #compute the reference probability (member n as reference) equal to 1 where data[n][var]>ce and 0 elsewhere
            print(np.shape(proba),np.shape(proba_ref),np.shape(e1),np.shape(e2))
            SPS.append((((proba-proba_ref)**2)*e1*e2).sum(("y","x"))) #computing of the SPS
    
    elif type_ref=="other_ens": #members of data_ref are taken one after another as the reference
        nb_member_ref=len(data_ref) #number of members of data_ref
        proba=fproba_ce(data,var,ce)*mask_domain #compute the probability (frequency) array of having the value of var strictly superior to ce for data
        for n in range(nb_member_ref):#iteration over the members of data_ref
            SPS_id[1].append("ref mber: "+str(n+1))
            proba_ref=fmask_ce(data_ref[n][var],ce)*mask_domain #compute the reference probability (member n of data_ref as reference) equal to 1 where data_ref[n][var]>ce and 0 elsewhere
            print(np.shape(proba),np.shape(proba_ref),np.shape(e1),np.shape(e2))
            SPS.append((((proba-proba_ref)**2)*e1*e2).sum(("y","x"))) #computing of the SPS
    
    elif type_ref=="masked_field": #data_ref is taken as the reference
        proba=fproba_ce(data,var,ce)*mask_domain #compute the probability (frequency) array of having the value of var strictly superior to ce for data
        proba_ref=data_ref*mask_domain
        print(np.shape(proba),np.shape(proba_ref),np.shape(e1),np.shape(e2))
        SPS.append((((proba-proba_ref)**2)*e1*e2).sum(("y","x"))) #computing of the SPS
    else:
        print("This kind of reference is not possible.")
        return ()
    
    return (SPS,SPS_id)

def fOU_IIEE(data,type_ref,data_ref,var,ce,e1,e2,corners=False,width=False,height=False):
    '''''
    This function computes the components O and U of the IIEE for the ensemble-median ice edge.
    If specified, the components O and U are computed over the domain defined with the indexes:
    [corners[1]:corners[1]+height+1,corners[0]:corners[0]+width+1]
    The reference chosen is indicated by type_ref. There are 3 possible kinds of reference:
    - same_ens: each member of data is taken one after another as the reference
    - other_ens: each member of data_ref is taken one after another as the reference
    - masked_field: the reference is a 2D (or 3D if there is a time dimension) field with values only between 0 and 1.
    /!\ the reference given should use the same ce criterion to define the ice edge
    This last type of reference could be processed observations, a probability map of another ensemble, etc.
    ...
    The arguments needed are the following:
    - data = data of the experiment of which we want to compute the score
    data should be a list of all the members, with members as a DataSet
    - type_ref = character string indicating the kind of reference that is data_ref (3 possibilities; see above)
    - data_ref = data to take as the reference should be coherent with type_ref
    if type_ref="same_ens" could be equal to False
    if type_ref="other_ens" data_ref should be a list of all the members, with members as a DataSet
    - var = name of the variable of interest over which to compute the score
    - ce = the condition/threshold to compute the score
    - e1, e2 = horizontal mesh sizes (should be 2D with the same dimensions "y" and "x" as data and data_ref!)
    - corners = position of the lower left corner. A tuple as follows: (x position, y position).
    If corners=False computation over all the domain (default).
    - width = length of the domain in the x-direction (default=False)
    - height = length of the domain in the y-direction (default=False)
    /!\ This function is designed to work with NEMO outputs!
    '''''
    #saving the sizes of the dimensions "x" and "y" of data
    size_y=data[0].sizes["y"] ; size_x=data[0].sizes["x"]
    
    #computing the mask of the domain over which to compute the SPS
    if corners: #computation over the domain defined by fmask_square_domain(size_x,size_y,corners,width,height)
        print("computation over the specified domain")
        mask_domain=libdiv.fmask_square_domain(size_y,size_x,corners,width,height)
    else: #computation over all the domain
        mask_domain=1
    
    OU=[] #list that will stock O and U components as (O,U)
    OU_id=[type_ref,[]] #list that will indicate the kind of reference and more information about OU

    if type_ref=="same_ens": #members of data are taken one after another as the reference
        nb_member=len(data) #number of members of data
        for n in range(nb_member): #iteration over all members of data
            OU_id[1].append("ref mber: "+str(n+1))

            #computation of a mask equal to 1 where the probability to have the value of var strictly superior to ce is strictly superior to 0.5
            proba=fproba_ce(data,var,ce,n) #computation of probability map of having the value of var strictly superior to ce
    
            mask_iemed=fmask_ce(proba,0.5)*mask_domain #computation of a mask equal to 1 where proba > 0.5 and 0 elsewhere

            #computation of a mask equal to 1 where data[n][var] > ce and 0 elsewhere with n the reference member
            mask_ref=fmask_ce(data[n][var],ce)*mask_domain

            #computing O component
            Otmp=np.maximum(mask_iemed-mask_ref,0.)
            O=(Otmp*e1*e2).sum(("y","x")) #spatial integral

            #computing U component
            Utmp=np.maximum(mask_ref-mask_iemed,0.)
            U=(Utmp*e1*e2).sum(("y","x")) #spatial integral
            
            OU.append((O,U))
            
    elif type_ref=="other_ens": #members of data_ref are taken one after another as the reference
        nb_member_ref=len(data_ref) #number of members of data_ref
        
        #computation of a mask equal to 1 where the probability to have the value of var strictly superior to ce is strictly superior to 0.5
        proba=fproba_ce(data,var,ce) #computation of probability map of having the value of var strictly superior to ce
    
        mask_iemed=fmask_ce(proba,0.5)*mask_domain #computation of a mask equal to where proba > 0.5 and 0 elsewhere
        
        for n in range(nb_member_ref):#iteration over all members of data_ref
            OU_id[1].append("ref mber: "+str(n+1))

            #computation of a mask equal to 1 where data_ref[n][var] > ce and 0 elsewhere with n the reference member
            mask_ref=fmask_ce(data_ref[n][var],ce)*mask_domain

            #computing O component
            Otmp=np.maximum(mask_iemed-mask_ref,0.)
            O=(Otmp*e1*e2).sum(("y","x")) #spatial integral

            #computing U component
            Utmp=np.maximum(mask_ref-mask_iemed,0.)
            U=(Utmp*e1*e2).sum(("y","x")) #spatial integral
            
            OU.append((O,U))
            
    elif type_ref=="masked_field": #data_ref is taken as the reference
        #computation of a mask equal to 1 where the probability to have the value of var strictly superior to ce is strictly superior to 0.5
        proba=fproba_ce(data,var,ce) #computation of probability map of having the value of var strictly superior to ce
    
        mask_iemed=fmask_ce(proba,0.5)*mask_domain #computation of a mask equal to where proba > 0.5 and 0 elsewhere

        #data_ref is the reference
        mask_ref=data_ref*mask_domain

        #computing O component
        Otmp=np.maximum(mask_iemed-mask_ref,0.)
        O=(Otmp*e1*e2).sum(("y","x")) #spatial integral

        #computing U component
        Utmp=np.maximum(mask_ref-mask_iemed,0.)
        U=(Utmp*e1*e2).sum(("y","x")) #spatial integral
 
        OU.append((O,U))
    else:
        print("This kind of reference is not possible.")
        return ()

    return (OU,OU_id)

def fproba_occurence_ce(data,var,ce,t0,tend): 
    '''''
    This function computes the frequency over all members of the event: the event data>ce occurs at least once over the period [t0,tend].
    ...
    The arguments needed are the following:
    - data = data of an experiment /!\ to work should have 3 dimensions in that order: "time_counter", "y", "x".
    data should be a list of all the members, with members as a DataSet
    - var = name of the variable of interest 
    - ce = the condition
    - t0 = time index of the first instant of the considered period
    - tend = time index of the last instant of the considered period
    '''''
    
    size_y=data[0].sizes["y"] ; size_x=data[0].sizes["x"] #horizontal dimension sizes (all the members have the same sizes)
    nb_mber=len(data) #number of members

    proba=xr.DataArray(data=np.zeros((size_y,size_x)), dims=["y","x"]) #initialisation
    
    for imb in range(nb_mber):#iteration over the members
        
        mask=fmask_ce(data[imb][var],ce) #mask in 3 dimensions equal to 1 if data[imb][var]>ce
        proba+=np.where(mask[t0:tend+1,:,:].sum("time_counter")>0,1,0) #+1 if the event data[imb][var]>ce occurs at least once
        
    return proba/nb_mber