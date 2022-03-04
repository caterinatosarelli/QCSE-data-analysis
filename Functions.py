import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Function to create voltage array to be used as axis
def volt(N_val, Vmin, Vmax):
    """
    This function creates the voltage values used during the voltage sweep.

    Parameters
    ----------
    N_val : number of voltage values used.
    Vmin : minimum voltage value.
    Vmax : maximum voltage value.

    Returns
    -------
    voltage_values : array of voltage values equally spaced.
    DelV : step size of the voltage array.

    """
    
    voltage_values = np.linspace(Vmin, Vmax, N_val)
    DelV = (Vmax-Vmin)/(N_val-1.)
    
    return voltage_values, DelV

#%% 

# Function to import and read the dataset
def read_data(import_dir, data_filename):
    """ 
    This function creates the data frames reading from txt files.
    
    Parameters
    ----------
    import_dir : path of the imported data.
    data_filename : name of the file to be imported.

    Returns
    -------
    data : matrix of pl intensity during voltage sweep.

    """
    
    read_data = import_dir + data_filename
    data = np.genfromtxt(read_data)
    
    return data


#%% Functions that will be used to perform fit

# Gaussian function
def gaussian_gen(x,Amp,sig,mu,Bg):
    """
    This function creates a Gaussian function.

    Parameters
    ----------
    x : array of energy values
    Amp : float, amplitude of curve
    sig : float, standard deviation
    mu : float, mean
    Bg : float, background
    
    Returns
    -------
    array of float following the formula : Bg + Amp * exp(-(x-mu)^2/sig^2)

    """
    return(Bg + Amp*np.exp(-np.power(x-mu,2.)/(2*np.power(sig,2.))))

# Quadratic function
def polynom2(x,a0,a1,a2):
    """
    This function creates a parabolic function.

    Parameters
    ----------
    x: array of float
    a0 : intercept
    a1 : coefficient of x
    a2 : coefficient of x^2

    Returns
    -------
    array of float following the formula : a0 + a1 * x + a2 * x^2

    """
    return(a0 + a1*x + a2*np.power(x,2.))

#%%
# Function to create the colormap: 
# takes X and Y as axis values, with extent Xmin-Xmax and Ymin-Ymax
# Imin and Imax can be used to set max and min value of colorbar
# X/Ymin and X/Ymax are the indexes of the arrays to set the boundaries of the image

def map_plot(Matrix,Fig,X,Y,Xmin,Xmax,Ymin,Ymax,Imin,Imax,Xlabel,Ylabel):
    """
    This function creates a colormap of the intensity of PL recorded during the 
    voltage sweep.

    Parameters
    ----------
    Matrix : data of the sweep
    Fig : str, name of the figure
    X : float, values of x axis (voltage)
    Y : float, values of y axis (energy)
    Xmin : index of min x value
    Xmax : index of max x value 
    Ymin : index of max y value
    Ymax : index of max y value
    Imin : int, min intensity level in the colorbar
    Imax : int, max intensity level in the colorbar
    Xlabel : str, name of x axis
    Ylabel : str, name of y axis

    Returns
    -------
    fig : colormap of photoluminescence during voltage sweep

    """

    fig, ax = plt.subplots(figsize=(10,10), num=Fig)
    pos = ax.imshow(Matrix, aspect='auto', origin='lower', interpolation='none',
                extent=[X[Xmin],X[Xmax],Y[Ymin],Y[Ymax]], 
                cmap='rainbow',vmin=Imin, vmax=Imax)
    
    cbar = plt.colorbar(pos, pad=0.1, shrink=0.7)
    cbar.set_label('PL intensity (arb. units)', rotation=270, labelpad=25)
    plt.xlabel(Xlabel)
    plt.ylabel(Ylabel)

    return fig

#%%
def peak_lim(wav, en, peak_pos, sigma) : 
    """
    This function finds the index at which the wavelength value corresponds to the 
    energy value given as input, following the conversion wav(nm) = 1239.8/en(eV). 
    It also finds two points distant 10*sigma from peak position.

    Parameters
    ----------
    wav : array of wavelength values 
    peak_pos : float, energy value taken by gauss. fit
    sigma : float, sigma value taken by gauss. fit, can be used to find index of mu+sigma
    
    Returns
    -------
    pos : int, index at which wav(nm)=1239.8/en(eV)
    lim1 : int, lower limit of range around peak pos.
    lim2 : int, upper limit of range aroung peak pos.

    """

    pos = np.where(wav>=1239.8/peak_pos)[0][0]
    midi = en[pos]
    delta = -10 * (np.where(wav>=1239.8/(peak_pos+sigma))[0][0] - pos)
    lim1 = pos - delta
    lim2 = pos + delta
    
    return midi, lim1, lim2
#%%
def direction(forward, Izero, voltage, vlim, i):
    """
    This function sets parameters to be used in fit_loop function.

    Parameters
    ----------
    forward : bool, set direction in which the fit is performed along voltage array
    Izero : int, index of voltage as starting point for the loop
    voltage : array of voltage values
    vlim : float, limit of voltage range selected
    i : int, value of the iterable variable of the loop

    Returns
    -------
    v_index : int, index of current voltage value
    vi : float, value of current voltage
    end : bool, true if out of the voltage range selected

    """
   
    if forward:
        v_index = Izero+i+1
        vi = voltage[v_index]
        end = vi>(np.round(voltage[vlim],2))

    else:
        v_index = Izero-i-1
        vi = voltage[v_index]
        end = vi<(np.round(voltage[vlim],2))
        
    return v_index, vi, end

#%%        
def fit_loop(Izero, voltage, VLim, params, spectrum, matrix, forward):
    """
    This function performes a series of gaussian fit in matrix of data over 
    the voltage range selected.

    Parameters
    ----------
    Izero : int, index of voltage starting point of the loop
    voltage : array of voltage values
    VLim : int, index of limit of voltage range
    params : array with starting parameters for gaussian fit
    spectrum : dataframe containing energy and wavelength arrays
    matrix : matrix of intensity values
    forward : bool, defines the direction of the loop from the starting point

    Returns
    -------
    table : dataframe containing parameters from fit for every voltage value in the range.

    """
    
    v0 = voltage[Izero]
    Wav = spectrum[:,0]
    En = spectrum[:,1]
    
    A, sigma, mu, bg = params

    # LimC0 = peak_lim(Wav, mu, 0)[0]
    midi = peak_lim(Wav, En, mu, 0)[0] 
    fit_vals=[[A, sigma, mu, bg, midi, v0]]
    columns_name = ('A','sigma','mu','Bg','mid','vfit')
    table = pd.DataFrame(fit_vals, columns=columns_name)

    for i in range(len(voltage)):
        
        v_index, vi, end = direction(forward, Izero, voltage, VLim, i)
        
        if end :
            print('End of selected voltage range')
            break
        else : 
            # to avoid division by 0
            if table['mu'][i]==0.:
                midi = table['mid'][i]
            else:
                midi = table['mu'][i]
                # redefines energy range of the peak to be fitted using sigma
                Lim1 = peak_lim(Wav,midi,table['sigma'][i])[1]
                Lim2 = peak_lim(Wav,midi,table['sigma'][i])[2]
            try:
                Gauss_fit = curve_fit(gaussian_gen,En[Lim1:Lim2],matrix[Lim1:Lim2,v_index],
                                      p0=(table['A'][i],table['sigma'][i],midi,0),
                                      )
                parsi = Gauss_fit[0]
                
            except RuntimeError: 
                print(vi,'V error in fit - just use peak')
                parsi = [0,0,0,0]
            
            # reinitialize midi in the new range
            midi = En[Lim1+np.argmax(matrix[Lim1:Lim2,v_index])]
            
            # append new parameters to dataframe
            parsi = np.append(parsi,[midi,vi]).reshape(1,6)
            dfi = pd.DataFrame(parsi, columns = columns_name)
            table = pd.concat([table, dfi]).reset_index(drop=True)
    
    return table
