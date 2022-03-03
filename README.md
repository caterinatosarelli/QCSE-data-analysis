# General informations
This project is meant to be used by a person studying the application of field to a pin junction containing quantum dots. The code analyses the luminescence spectrum 
of the dots as the voltage is swept from negative to positive values, giving informations about the redshift that the emission lines undergo thanks to the 
Quantum Confined Stark Effect (QCSE). This phenomenon tells us that the emission energy of a confined structure is supposed to vary as a quadratic function of the field,
following the formula :   

<img src="https://latex.codecogs.com/svg.image?E&space;=&space;E_{0}&space;&plus;&space;p&space;\times&space;F&space;&plus;&space;\beta&space;\times&space;F^{2}" title="E = E_{0} + p \times F + \beta \times F^{2}" />    

where F is the applied field, obtained from the voltage using    

<img src="https://latex.codecogs.com/svg.image?F&space;=&space;-&space;\frac{\left&space;(&space;V&space;-&space;Vbi&space;\right&space;)}{d}" title="F = - \frac{\left ( V - Vbi \right )}{d}" />

Vbi being the built-in voltage of the junction and d the thickness of the intrinsic region.

The code is meant to receive two datasets:

1. The intensity map recorded while varying the voltage: each column is correspond to a PL spectrum taken at a certain voltage value.
Each row contains then the intensity at each voltage value for each energy/wavelength value. Note that this file contains only intensity values, not voltage or energy values. 
See "Dev1_C3_-2,2to0,6_970nW_960nm_1sec_281stp_pos3.txt" in the repository as example.

2. The grid of the spectrometer, which is not reproducible by using a simple formula so must be taken from a single spectrum file. 
This type of file contains 3 columns: one for energy values, one for wavelength values and one for counts. From this file the user will take the energy and wavelength values to be used. See "gr960nm.txt" in the repository as example.

# Code structure
The project is divided in 4 main parts:

* Functions.py contains all the functions used to create the matrix of the data, to plot it and to fit it.
* ClickFind.py contains two functions to be able to save position of the mouse clicked as index. The findY function 
is meant to select the peak the user wants to perform analysis on, clicking on either sides of the peak. The findX is meant 
to be used when the colormap shows problems in the dataset that do not allow to use the whole voltage range. To interact with the figure, the user needs to set in python Tools/Preferences/Graphics to be QT5.
* testing.py contains tests of the functions in Functions.py.
* main.py is the main part of the code. When runned, it needs arguments to create the voltage range used and to import the dataset. 
After the creation of the colormap and the selection of the ranges, a fitting loop is performed starting from a voltage choosen by the user. 
A first fit is performed on the spectrum at that voltage and the parameters obtained are then used as starting guess for the spectra in the next voltage value. 
This procedure is repeated in the function fit_loop, starting from the "zeroPoint" voltage, first in the "left side" voltage range (forward=False), 
then in the "right side" voltage range (forward=True). For this reason a good choice for zeroPoint is the middle of the voltage range selected. 
The parameters are then used to perform data analysis of most interest. The peak position is then plotted vs field and then fitted with a parabolic function. 
The parameters obtained are printed in a file called "Stark_shift_params". Also sigma vs voltage is plotted since it can be a useful information on how the peak evolves. 

# Configuration
To start the analysis the user needs to call main.py giving as input :
- number of voltage values used
- minimum voltage value used
- maximum voltage value used
- built-in voltage of the pin diode
- thickness of the intrinsic region of the diode (in nm)
- path to data
- data filename
- grid filename
- save (y/n) , to save plots and results
