"""
Functions relative to ellipses (fit, display, capture test)
"""

print(f"Name: {__name__}")
print(f"Package: {__package__}")

print("This is a collection of functions relative to ellipses.")

import numpy as np
from scipy.stats import chi2
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt

def fit_ellipse(points,proba):
    '''''
    This function returns properties of the proba x 100% confidence ellipse.
    The ellipse is obtained assuming that the positions (x and y) of the points follow a bivariate normal distribution.
    ...
    Inputs:
    points = points over which we fit the bivariate normal distribution 
        shape(points) should be (n,2), with n the number of observations 
        the first column corresponds to positions on the x-axis, and the second column to positions on the y-axis
    proba = probability of the confidence ellipse (between 0 and 1)
    ...
    Outputs:
    mu = mean of the distribution
        shape(mu)=(1,2)
        the first column corresponds to the mean of the x-positions, and the second column to the mean of the y-positions
    sigma = covariance matrix of the distribution (unbiased definition)
        shape(sigma)=(2,2)
        the first line is associated with the x-positions (in the first column: variance of the x-positions and in the second column: covariance with the y-positions)
        the second line is associated with the y-positions (in the first column: covariance with the x-positions and in the second column: variance of the y-positions)
        the matrix is symmetrical
    s = eigenvalues of sigma
        shape(s)=(1,2)
    U = change-of-basis matrix
        shape(U)=(2,2)
        first column: eigenvector associated with the first eigenvalue
        second column: eigenvector associated with the second eigenvalue
    axlength = length of the axes of the ellipse
        shape(axlength)=(1,2)
        the first column corresponds to the size (total) of the axis defined by the first eigenvector
        the second column corresponds to the size (total) of the axis defined by the second eigenvector
    angle = angle between the x-axis and the eigenvectors
        shape(angle)=(1,2)
        the first column corresponds to the angle between the first eigenvector and the x-axis
        the second column corresponds to the angle between the second eigenvector and the x-axis
    '''''

    ##Computation of the mean and covariance matrix
    
    mu=np.mean(points,axis=0)
    sigma=np.cov(points,rowvar=False,ddof=1)

    ##Computation of the eigenvalues and eigenvectors
    
    s, U = np.linalg.eig(sigma)
    
    ##Deducing properties of the ellipse
    
    if s[0]!=s[1]: #check that the two eigenvalues are different
        qchi2=chi2.ppf(proba,df=2) #; print("chi2 quantile proba "+str(proba)+":",qchi2)
        axlength=2*np.sqrt(qchi2*s) #total length of the axes defined by each of the eigenvectors
        angle=np.degrees(np.arctan2(U[1, :], U[0, :])) #angle between the x-axis and each eigenvector
        return (mu, sigma, s, U, axlength, angle)
    else:
        print("the two eigenvalues are the same!!! ")
        return (mu,sigma,s,U)


def draw_ellipse(center,axlength,angle,ax,\
                 legend=False,proba=False,loc_legend="upper right",size_legend=12,text_legend="",\
                 draw_axes=False,axes=False,axes_color="black",axes_style="-",axes_width=1,\
                 edge_color="blue",edge_width=1,edge_style="-",fill_color="None",alpha=1,zorder=2):
    '''''
    This function draws the ellipse whose properties are given.
    ...
    Inputs:
    center = (x_pos,y_pos) coordinates of the center of the ellipse
    axlength = length (total) of the axes of the ellipse 
        shape(axlength)=(1,2)
    angle = angles between each axis of the ellipse and the x-axis
        shape(angle)=(1,2)
    ax = subfigure on which the ellipse is drawn

    optional inputs:
    legend = boolean equal to True if we want to plot a legend and to False otherwise (default=False)
    proba = probability corresponding to the ellipse
        if legend=True, it should be a number between 0 and 1
        if legend=True, display of the probability as the legend of the ellipse (if text_legend=="")
    loc_legend = position of the legend, default="upper right"
    size_legend = size of the text of the legend, default = 12
    text_legend = chosen text of the legend (default = "")
    draw_axes = boolean equal to True if we want to draw the axes and to False otherwise (default=False)
    axes = axes of the ellipse (vectors), if draw_axes=True, we draw them
        shape(axes)=(2,2)
        the first column corresponds to the first axis (associated with the value in the first column of axlength)
        the second column corresponds to the second axis (associated with the value in the second column of axlength)
    axes_color = color of the axes if drawn, default="black"
    axes_style = style of the axes if drawn, default="-"
    axes_width = width of the lines of the drawn axes, default=1
    edge_color = color of the edge of the ellipse; if equal to "None", the edge is not drawn; default="blue"
    edge_style = style of the line of the edge if drawn, default="-"
    edge_width = width of the lines of the edge, default=1
    fill_color = color of the interior of the ellipse; if equal to "None", the interior is not filled; default="None"
    alpha = the opacity of the ellipse, between 0 and 1: 1 is total opacity, and 0 is no opacity (no ellipse); default=1
    zorder = order of the ellipse compared to the other elements plotted on the figure, default=2

    /!\ To work, you have to plot another thing (using plt.plot or plt.scatter for example) than just the ellipse
    '''''
    #plot of the ellipse
    label=text_legend+("proba "+str(proba*100)+" %")*(text_legend=="")
    ax.add_patch(Ellipse(center, axlength[0], axlength[1], angle=angle[0],\
                         facecolor=fill_color,edgecolor=edge_color,ls=edge_style,lw=edge_width,zorder=zorder,label=""+label*legend,alpha=alpha))
    
    #plot of the axes if asked
    if draw_axes:
        print("make sure to have axes with the same scale to not distort the angles!!")
        ax.plot([center[0]-0.5*axlength[0]*axes[0,0],center[0],center[0]+0.5*axlength[0]*axes[0,0]],\
            [center[1]-0.5*axlength[0]*axes[1,0],center[1],center[1]+0.5*axlength[0]*axes[1,0]],\
            color=axes_color,ls=axes_style,lw=axes_width,zorder=zorder)
        ax.plot([center[0]-0.5*axlength[1]*axes[0,1],center[0],center[0]+0.5*axlength[1]*axes[0,1]],\
            [center[1]-0.5*axlength[1]*axes[1,1],center[1],center[1]+0.5*axlength[1]*axes[1,1]],\
            color=axes_color,ls=axes_style,lw=axes_width,zorder=zorder)

    #plot of a legend if asked
    if legend:
        ax.legend(loc=loc_legend,fontsize=size_legend)


def isin_ellipse(points,mu,s,U,proba):
    '''''
    This function returns booleans: True for the points within the ellipse defined by mu, s, U, and proba, and False for the points outside.
    ...
    Inputs:
    points = points over which we do the test 
        shape(points) should be (n,2), with n the number of observations 
        the first column corresponds to positions on the x-axis and the second column to positions on the y-axis
    mu = mean of the distribution associated with the confidence ellipse
        shape(mu)=(1,2)
        the first column corresponds to the mean of the x-positions and the second column to the mean of the y-positions
    s = eigenvalues of sigma the covariance matrix of the distribution associated with the confidence ellipse
        shape(s)=(1,2)
    U = change-of-basis matrix associated with the confidence ellipse
        shape(U)=(2,2)
        first column: eigenvector associated with the first eigenvalue
        second column: eigenvector associated with the second eigenvalue
    proba = probability of the confidence ellipse (between 0 and 1)
    ...
    Outputs:
    isin = boolean matrix equal to True for the points within the ellipse and to False otherwise
    '''''
    #computation of the new coordinates (in the basis defined by the eigenvectors)
    new_coord=np.dot(points-mu,U)
    
    #computation of the test
    qchi2=chi2.ppf(proba,df=2)
    test=(new_coord[:,0]**2/(qchi2*s[0])+new_coord[:,1]**2/(qchi2*s[1])) #; print(test)
    isin=test<=1
    isin=np.where(test<0,np.nan,isin) #test should be positive; if not, it is that the inputs are not correct, and so the output value of the function doesn't exist 

    return isin