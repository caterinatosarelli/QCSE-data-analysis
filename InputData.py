import os

# this file contains all the infos and the data to be analysed and plotted
# in my case all the data where in the same directory, add other import directories 
# if this is not you case
work_dir = 'C:/Users/catet_ch27s6d/OneDrive/Desktop/Software/'
import_dir = 'C:/Users/catet_ch27s6d/OneDrive/Desktop/Software/data/'
save_dir = 'C:/Users/catet_ch27s6d/OneDrive/Desktop/Software/data/plots/'

# The map file is the matrix created during the application of the voltage.
# To plot the matrix you need the voltage axis and the energy (or wavelength) axis,
# the voltages values are equispaced so you can create them using a np array.
# The energy axis (and corresponding wavelength axis) was taken from the file 
# of a single spectrum

map_filename = 'Dev1_C3_-2,2to0,6_970nW_960nm_1sec_281stp_pos3.txt'
spectrum_filename = 'gr960nm.txt'

#set the name that will be used to save all the plots generated for the dataset
plot_filename = 'Dev1_C3'

# to create the voltage axis you need to set number of values, Vmin and Vmax
# you may want to convert the voltage into field, for which you need to set the
# Vbi and the thickness of the intrinsic region of your diode (called d)

N = 281 #number of voltage values, must be equal to number of columnes of the map
Vmin = -2.2
Vmax = 0.6
Vbi = 1.76
d = 118 # thickness of intrins. region written in nm


