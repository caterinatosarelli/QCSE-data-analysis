import numpy as np
import pandas as pd
import argparse
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import Functions as func
import Clickfind as click


plt.close('all')


parser = argparse.ArgumentParser(description=("input data"))

parser.add_argument('N_val', type = int, help = 'Number of voltage values', default=281)
parser.add_argument('Vmin', type = float, help = 'Minimum voltage value used', default=-2.2)
parser.add_argument('Vmax', type = float, help = 'Maximum voltage value used',default=0.6)
parser.add_argument('Vbi', type = float, help = 'Built-in voltage',default=1.76)
parser.add_argument('d', type = float, help = 'Intrinsic region thickness (nm)',default=118)
parser.add_argument('data_path', type = str, help = 'directory of input data',default='data/')
parser.add_argument('name_of_file', type = str, help = 'data filename', default='Dev1_C3_-2,2to0,6_970nW_960nm_1sec_281stp_pos3.txt')
parser.add_argument('grid_shape',type = str, help = 'grid filename', default = 'gr960nm.txt')
parser.add_argument('saving', type = bool, help = 'save plots', default = True)

args = parser.parse_args()

nval = args.N_val
vmin = args.Vmin
vmax = args.Vmax
vbi = args.Vbi
t = args.d
path = args.data_path
data = args.name_of_file
grid = args.grid_shape
save = args.saving

voltage, dV = func.volt(nval, vmin, vmax)
matrix = func.read_data(path, data)
spectrum = func.read_data(path, grid)

#from spectrum take wavelength and energy values
Wav = spectrum[:,0]
En = spectrum[:,1]


#%% plot the data in a colormap

colormap = func.MapPlot(matrix, 'Map', voltage, Wav, 0, -1, 0, -1, 0, 600, 
        'Voltage(V)', 'Wavelength(nm)')

if save:
    plt.savefig("Colormap.png",dpi=600)

#%% select peak to be fitted clicking on the plot using function that takes y coord

Zoom_map, Lim1, Lim2 = click.findY(colormap,Wav,matrix)
plt.close(colormap) 

func.MapPlot(Zoom_map, 'Zoom', voltage, Wav, 0, -1, Lim1, Lim2, 
        0, 800, 'Voltages(V)', 'Wavelength(nm)')

plt.close('Map')
#%% select the voltage range if needed

choice = input('Do you need to pick the range of voltages? y/n >> ')

if choice == 'y':
    zoom =plt.figure("Zoom")
    VLimLow, VLimHi = click.findX(zoom, voltage)
    plt.close('Zoom')
    
else:
    VLimLow = np.argmin(voltage)
    VLimHi = np.argmax(voltage)


#%%  fitting routine

# Chose V from which start the fitting routine. A good choice is approximately 
# the middle of the voltage range you picked
zeroPoint = float(input('Type approximate voltage that you want to use to get first \n approximation of fit >> '))

Izero = int((zeroPoint-vmin)/dV) #index of starting point
Vzero = voltage[Izero]

# position of peak
midi = En[Lim1+np.argmax(matrix[Lim1:Lim2,Izero])]

# perform gaussian fit in the initial point to use parameters as starting guess 
# to fit each spectrum
Gauss_fit = curve_fit(func.gaussian_gen,En[Lim1:Lim2],matrix[Lim1:Lim2,Izero],
                    p0=(100,0.00005,midi,0))
A0,sigma0,mu0,Bg0 = Gauss_fit[0]

#double check center of peak
LimC0 = np.where(Wav>=1239.8/mu0)[0][0]
midi = En[LimC0]
LimDel = -10*(np.where(Wav>=1239.8/(midi+sigma0))[0][0] - LimC0)

Gauss_fit = curve_fit(func.gaussian_gen,En[LimC0-LimDel:LimC0+LimDel],
                    matrix[LimC0-LimDel:LimC0+LimDel,Izero],
                    p0=(200,sigma0,midi,0))
A0,sigma0,mu0,Bg0 = Gauss_fit[0]

# fit two ranges of voltage separately starting from Izero, forward and backward
table_low = func.fit_loop(Izero, voltage = voltage, VLim = VLimLow, params = Gauss_fit[0],
                      spectrum = spectrum, matrix = matrix, forward = False)
table_high = func.fit_loop(Izero, voltage = voltage, VLim = VLimHi, params = Gauss_fit[0],
                      spectrum = spectrum, matrix = matrix, forward = True)

table = pd.concat([table_low, table_high], keys='vfit')
table = table.sort_values('vfit')


#%% plot peak position taken from gauss. mean and position of max intensity to check
#   that the fitting routine was performed on the correct peak

plt.figure()
plt.plot(table['vfit'],table['mu'], label = 'peak pos. from gauss. fit')
plt.plot(table['vfit'],table['mid'], label = 'pos. of maximum')
plt.xlabel('voltage (V)')
plt.ylabel('peak position (eV)')
plt.legend()

if save:
    plt.savefig("Gauss_fitroutine.png",dpi=600)
#%% fitting shif of peak position with parabolic function of the field

field = - ( table['vfit'] - vbi ) / t * 10000 # field values in kV/cm

starkshift_fit = curve_fit(func.polynom2, field, table['mu'])
parabolic_shift = func.polynom2(field, *starkshift_fit[0])

#plot peak position from gauss fit and parabolic shift
plt.figure()
fig, ax = plt.subplots(figsize=(11,9))
plt.plot(field,parabolic_shift, label = 'parabolic fit')
plt.plot(field,table['mu'], label = 'peak center from gauss.fit')
plt.xlabel('voltage (V)')
plt.ylabel('peak position (eV)')
plt.legend()

if save:
    plt.savefig("Starkshift_fit.png", dpi=600)

#%% plot sigma

plt.figure()
plt.plot(table['vfit'], table['sigma'])
plt.xlabel('voltage (V)')
plt.ylabel('sigma')

if save:
    plt.savefig("Sigma.png",dpi=600)
