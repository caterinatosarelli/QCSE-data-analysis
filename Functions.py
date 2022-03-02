import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Function to create voltage array to be used as axis
def volt(N_val, Vmin, Vmax):
    """
    This method creates the voltage values used during the voltage sweep.

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
    This method creates the data frames reading from txt files.
    
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
    This method creates a Gaussian function.

    Parameters
    ----------
    x, Amp, sig, mu, Bg
    
    Returns
    -------
    Bg + Amp * exp(-(x-mu)^2/sig^2)

    """
    return(Bg + Amp*np.exp(-np.power(x-mu,2.)/(2*np.power(sig,2.))))

# Quadratic function
def polynom2(x,a0,a1,a2):
    """
    This method creates a parabolic function.

    Parameters
    ----------
    x, a0, a1, a2

    Returns
    -------
    a0 + a1 * x + a2 * x^2

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
    Fig : name of the figure
    X : values of x axis
    Y : values of y axis
    Xmin : index of min x value
    Xmax : index of max x value 
    Ymin : index of max y value
    Ymax : index of max y value
    Imin : min intensity level in the colorbar
    Imax : max intensity level in the colorbar
    Xlabel : name of x axis
    Ylabel : name of y axis

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

#%% fitting loop

def fit_loop(Izero, voltage, VLim, params, spectrum, matrix, forward):
    """
    This method performs a gaussian fit on each spectrum present in the matrix 
    of the voltage sweep using a loop on the voltage range selected.

    Parameters
    ----------
    Izero : index of the voltage of the starting point 
    voltage : voltage values used
    VLim : limit of the voltage range on which the fit will be performed
    params : starting parameters for gaussian fit
    spectrum : dataframe from which wavelength and energy values are taken
    matrix : data of the spectrum during voltage sweep
    forward : boolean variable to set direction of the loop

    Returns
    -------
    table : dataframe with results of gaussian fit for each spectrum in matrix

    """
    
    v0 = voltage[Izero]
    Wav = spectrum[:,0]
    En = spectrum[:,1]
    
    LimC0 = np.where(Wav>=1239.8/params[2])[0][0] 
    # sets midi as the position of the mean from the first fit
    midi = En[LimC0] 
    fit_vals=[[params[0],params[1],params[2],params[3],midi,v0]]
    columns_name = ('A','sigma','mu','Bg','mid','vfit')
    table = pd.DataFrame(fit_vals, columns=columns_name)

    for i in range(len(voltage)):
        if forward:
            v_index = Izero+i+1
            vi = voltage[v_index]
            end = vi>(np.round(voltage[VLim],2))
        else:
            v_index = Izero-i-1
            vi = voltage[v_index]
            end = vi<(np.round(voltage[VLim],2))
        if end:
            print('End of selected voltage range')
            break
        else:
            try:
                if table['mu'][i]==0.:
                    midi = table['mid'][i]
                else:
                    midi = table['mu'][i]
                # at every cycle the position of the peak is updated with the mean from the previous fit
                LimC0 = np.where(Wav>=1239.8/midi)[0][0]
                LimDel = -10*(np.where(Wav>=1239.8/(midi+table['sigma'][i]))[0][0] - LimC0)
                Lim1 = LimC0 - LimDel
                Lim2 = LimC0 + LimDel
                Gauss_fit = curve_fit(gaussian_gen,En[Lim1:Lim2],matrix[Lim1:Lim2,v_index],
                            p0=(table['A'][i],table['sigma'][i],midi,0),
                            )
                parsi = Gauss_fit[0]
                midi = En[Lim1+np.argmax(matrix[Lim1:Lim2,v_index])]
                parsi = np.append(parsi,[midi,vi]).reshape(1,6)
                dfi = pd.DataFrame(parsi, columns = columns_name)
                table = pd.concat([table, dfi]).reset_index(drop=True)
               
            except RuntimeError:
                # if the fit is not performed, it takes the position of the highest intensity as peak
                print(vi,'V error in fit - just use peak')
                midi = En[Lim1+np.argmax(matrix[Lim1:Lim2,v_index])]
                dfi = pd.DataFrame([[0,0,0,0,midi,v0]], columns = columns_name)
                table = pd.concat([table, dfi]).reset_index(drop=True)
    return table
  
