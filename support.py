import numpy as np 
import pandas as pd

def Combine_Plates(red, blue, yellow, green):
    '''mimicks combining 4 plates into 1 on the rotor'''
    #makes sure we are combining same size arrays
    assert (red.shape == blue.shape) & (yellow.shape == green.shape) & (red.shape == green.shape)
   
    
    Combined= np.empty((red.shape[0]*2,red.shape[1]*2)).astype(object)

    Combined[::2,::2 ]= red
    Combined[::2,1::2 ]= blue
    Combined[1::2,::2 ]= yellow
    Combined[1::2,1::2 ]= green

    return(Combined)