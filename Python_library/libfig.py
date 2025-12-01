"""
Diverse functions for creating figures.
"""

print(f"Name: {__name__}")
print(f"Package: {__package__}")

print("This is a collection of diverse functions for creating figures.")

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
import numpy as np

def faxes(nsubfig,ncol,nrow,srow=-1):
    '''''
    This function returns axes of a figure with nsubfig subfigures, nrow rows and ncol columns.
    If the number of subfigures is not equal to ncol*nrow, you can choose which row will have fewer subfigures with the argument srow.
    '''''
    print("number of choosen subfigures = " +str(nsubfig))
    print("number of choosen columns = " +str(ncol))
    print("number of choosen rows = " +str(nrow))
    print("choosen row with unequal number of columns: row "+str(srow))
    
    if nsubfig>ncol*nrow:
        print("The number of rows and columns does not match with the number of subfigures")
        return ()
    
    axs=[] #initialisation of the list that will contain the axes
    axs_tmp=[] #initialisation of the list that will contain the axes (in the wrong order)
    pos_axs_tmp=np.zeros((nrow,ncol),dtype="int8")+999 #initialisation of the list that will contain the position number of the associated axis for a given row and column.  If no axis is associated with the row-column combination, the array is equal to 999
    gs = gridspec.GridSpec(nrow, (nsubfig==ncol*nrow)*nsubfig + (nsubfig!=ncol*nrow)*ncol*nrow)

    if nrow==1:#if only one row
        for i in range(nsubfig): #iteration over the subfigures
            axs.append(plt.subplot(gs[0, i]))
    else:
        ncolsolo= nsubfig -(nrow-1)*ncol #number of columns on srow
        for i in range(nsubfig): #iteration over the subfigures
            icol=i%ncol #column of interest
            irow_tmp=i//ncol #row of interest if all rows have the same number of columns
        
            ##some useful Booleans
            selectrow=(irow_tmp==(srow-1))*1 #row with fewer subfigures
            selectshift=selectrow*(icol<=(ncolsolo-1)) #if we should apply a shift or not
            test_row=(icol>((ncolsolo-1)))*(irow_tmp>=(srow-1))*1 #if we should go to the next row
        
            irow=irow_tmp+test_row #row of interest

            begincol=2 * icol + selectshift * (ncol - ncolsolo)
        
            pos_axs_tmp[irow,icol]=i
            axs_tmp.append(plt.subplot(gs[irow, begincol:begincol+2]))

        #putting the axes in the right order for the following
        for irow in range(nrow):
            for icol in range(ncol):
                if pos_axs_tmp[irow,icol]!=999:
                    axs.append(axs_tmp[pos_axs_tmp[irow,icol]])
    return axs

def fgraph_attributesv2(exp_name,type_color,type_label="full",focus_model="both"):
    '''''
    This function returns the color and label associated with this experiment according to our choices.
    —-
    exp_name = name of the experiment
    type_color = the kind of color palette to use. 2 possibilities: "atmos" and "rheo". The different initial amplitudes of perturbation correspond to different colors.
        "atmos" different color between "ABL" and "BLK" but the same color regardless of the rheology used
        "rheo" different color between "BBM" and "EVP" but the same color regardless of the atmosphere used
    type_label = the kind of label we want. 3 possibilities: "IC", "model" or "full". "full" is the default.
        "IC" gives the amplitude of perturbation as "IC STD X km"
        "model" gives the kind of model according to what has been asked in focus_model
        "full" gives a combination of the two above, for example: "ABL+BBM ; IC STD 50 km"
    focus_model = the part of the kind of the model we want to display. 3 possibilities: "atmos", "rheo" or "both". "both" is the default
        "atmos" gives the kind of atmos used ("BLK" or "ABL")
        "rheo" gives the kind of rheology used ("BBM" or "EVP")
        "both" gives more complete information about the model used, for example: "ABL+BBM"   
    '''''
    #The different color palettes that could be used
    colors1=["indigo","royalblue","coral"]
    colors2=["violet","darkgreen","gold"]

    #The different initial amplitudes of perturbation possible
    std_poss=["1","10","50"]

    #The kind of model
    label_model_atmos=exp_name[1:4]
    label_model_rheo=exp_name[-6:-3]
    if focus_model=="atmos": #focus on the atmosphere
        label_model=label_model_atmos
    elif focus_model=="rheo": #focus on the rheology
        label_model=label_model_rheo
    elif focus_model=="both": #no specific focus
        label_model=label_model_atmos+"+"+label_model_rheo
    else:
        print("Not in the 3 possibilities: atmos, rheo or both")
        return ()

    #The kind of initial amplitude of perturbation
    IC_perturbation=str(int(exp_name[-2:]))
    
    if IC_perturbation not in std_poss:
        print("Not a possible amplitude of perturbation")
        return()

    #The color associated with the experiment
    pos_IC=np.where(np.array(std_poss)==IC_perturbation)[0][0]
    if type_color=="atmos": #the colors depend on the kind of atmos
        if label_model_atmos=="ABL":
            color=colors1[pos_IC]
        elif label_model_atmos=="BLK":
            color=colors2[pos_IC]
        else:
            print("The name of the experiment given is not correct")
            return()
    
    elif type_color=="rheo": #the colors depend on the kind of rheology
        if label_model_rheo=="BBM":
            color=colors1[pos_IC]
        elif label_model_rheo=="EVP":
            color=colors2[pos_IC]
        else:
            print("The name of the experiment given is not correct")
            return()

    else:
        print("Not in the two possibilities: atmos or rheo")
        return()

    #The label to return
    if type_label=="IC": #label focus on the initial perturbations
        label="IC STD "+IC_perturbation+" km"
    elif type_label=="model": #label focus on the kind of model
        label=label_model
    elif type_label=="full": #full label
        label=label_model+" ; IC STD "+IC_perturbation+" km"
    else:
        print("Not in the 3 possibilities: IC, model or full")
        return ()

    return (color,label)
