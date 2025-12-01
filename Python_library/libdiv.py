"""
Diverse functions.
"""

print(f"Name: {__name__}")
print(f"Package: {__package__}")

print("This is a collection of diverse functions.")

import xarray as xr
import numpy as np

def fmask_square_domain(size_y,size_x,corners,width,height):
    '''''
    This function computes a mask equal to 1 in the square which is defined with the indexes:
    [corners[1]:corners[1]+height+1,corners[0]:corners[0]+width+1]
    The mask returned is 2D with the dimensions in the following order: "y", "x".
    ...
    The arguments needed are the following:
    - size_x, size_y = sizes of the dimensions "x" and "y"
    - corners = position of the lower left corner. A tuple as follows: (x position, y position).
    - width = length of the domain in the x-direction
    - height = length of the domain in the y-direction
    '''''
    mask=xr.DataArray(data=np.zeros((size_y,size_x),dtype=int), dims=["y","x"])#initialization
    mask[corners[1]:corners[1]+height+1,corners[0]:corners[0]+width+1]=1
    return mask