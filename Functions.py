import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from InputData import *


#function to create voltage array to be used as axis and to obtain the field
#d= intrinsic region thickness
#N = number of voltage values

def volt(N, Vmin, Vmax, Vbi, d):
    
    voltage_values = np.linspace(Vmin, Vmax, N)
    field_values = - ( voltage_values - Vbi ) / d*10000  # field in kV/cm
    DelV = (Vmax-Vmin)/(N-1.)
    
    return voltage_values, field_values, DelV

#%% function to import and read dataset

def read_data(import_dir, data_filename):
    
    read_data = import_dir + data_filename
    data = np.genfromtxt(read_data)
    
    return data
#%% fitting functions

# Gaussian fit
def gaussian_gen(x,Amp,sig,mu,Bg):
    return(Bg+Amp*np.exp(-np.power(x-mu,2.)/(2*np.power(sig,2.))))

# parabolic fit
def polynom2(x,a0,a1,a2):
  return(a0+a1*x + a2*np.power(x,2.))

#%%
#to plot map
# takes X and Y as axis values, with extent Xmin-Xmax and Ymin-Ymax
# Imin and Imax can be used to set max and min value of colorbar
# X/Ymin and X/Ymax are the indexes of the arrays to set the boundaries of the image
#you may want to add the title

def MapPlot(Matrix,Fig,X,Y,Xmin,Xmax,Ymin,Ymax,Imin,Imax,Xlabel,Ylabel):
    fig, ax = plt.subplots(figsize=(10,10), num=Fig)
    pos = ax.imshow(Matrix, aspect='auto', origin='lower', interpolation='none',
                extent=[X[Xmin],X[Xmax],Y[Ymin],Y[Ymax]], 
                cmap='rainbow',vmin=Imin, vmax=Imax)
    
    cbar = plt.colorbar(pos, pad=0.1, shrink=0.7)
    cbar.set_label('PL intensity (arb. units)', rotation=270, labelpad=25)
    plt.xlabel(Xlabel)
    plt.ylabel(Ylabel)
    
    plt.show()

   
