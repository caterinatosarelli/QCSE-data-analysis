import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as curve_fit
import sys
sys.path.append('C:/Users/catet_ch27s6d/OneDrive/Desktop/Software/')
from Functions import *
from InputData import *
from Clickfind import *

plt.close('all')

runfile(work_dir+"/InputData.py")

#%%
#create voltage and field values used
voltage, field, dV = volt(N,Vmin,Vmax,Vbi,d)

#create data frames, spectrum contains an array of energy and the corresponding 
#array of wavelengths
matrix = read_data(import_dir, map_filename)
spectrum = read_data(import_dir, spectrum_filename)

#from spectrum take wavelength and energy values
Wav = spectrum[:,0]
En = spectrum[:,1]

#%% plot the data in a colormap
# Imin and Imax can be adjusted to regulate the contrast of the image
MapPlot(matrix, 'Raw', voltage, Wav, 0, -1, 0, -1, 0, 600, 
        'Voltage(V)', 'Wavelength(nm)')

#plt.savefig(save_dir+plot_filename+'_raw.png',dpi=600)

#%% select peak to be fitted clicking on the plot using function that takes y coord
f = plt.figure("Raw")
times = findY(f)
plt.close(f)

# use the coordinates to set the limits of the zoomed map
ZLimLow = np.min(times)
ZLimHi = np.max(times)
ZLim1=np.where(Wav<=ZLimLow)[0][-1] 
ZLim2=np.where(Wav>=ZLimHi)[0][0] 

Zoom_matrix = matrix[ZLim1:ZLim2,:]

MapPlot(Zoom_matrix, 'Zoom', voltage, Wav, 0, -1, ZLim1, ZLim2, 
        0, 800, 'Voltages(V)', 'Wavelength(nm)')

#%% select the voltage range if needed

f=plt.figure("Zoom")
print('Do you need to pick the range of voltages? y/n')

choice = input()
if choice == 'y':
    times = findX(f)

    print('limits', times)

    VLimLow = np.min(times)
    VLimHi = np.max(times)
else:
    VLimLow = voltage[0]
    VLimHi = voltage[-1]

#%% 
#DelV= (Vmax-Vmin)/(N-1.)

# Chose V from which start the fitting routine. A good choice is approximately 
# the middle of the voltage range you picked
print('Type approximate voltage that you want to use to get first ')
print('approximation of fit')
zeroPoint = input()
zeroPoint = float(zeroPoint)

Izero = np.int((zeroPoint-Vmin)/dV) #index of starting point
Vzero = voltage[Izero]

# position of peak
midi = En[ZLim1+np.argmax(matrix[ZLim1:ZLim2,Izero])]

# p0 contains the initial guess of the parameters
# perform gaussian fit in the initial point to use parameters as starting guess 
# to fit each spectrum
Gauss_fit = curve_fit(gaussian_gen,En[ZLim1:ZLim2],matrix[ZLim1:ZLim2,Izero],
                    p0=(100,0.00005,midi,0))
A0,sigma0,mu0,Bg0 = Gauss_fit[0]

#plot to check the fit on the first spectrum
Fit_Gauss=gaussian_gen(En,A0,sigma0, mu0,Bg0)
plt.figure()
plt.plot(En[ZLim1:ZLim2],matrix[ZLim1:ZLim2,Izero])
plt.plot(En[ZLim1:ZLim2],Fit_Gauss[ZLim1:ZLim2])
plt.show()

# sets limit for fit as wavelength corresponding to mu0 
LimC0=np.where(Wav>=1239.8/mu0)[0][0] # takes index

# double check that everything is centred
midi = En[LimC0] 
meno = np.where(Wav>=1239.8/(midi+sigma0))
LimDel = -10*(np.where(Wav>=1239.8/(midi+sigma0))[0][0] - LimC0) #perchè negativo?

Gauss_fit = curve_fit(gaussian_gen,En[LimC0-LimDel:LimC0+LimDel],
                    matrix[LimC0-LimDel:LimC0+LimDel,Izero],
                    p0=(200,sigma0,midi,0))
A0,sigma0,mu0,Bg0 = Gauss_fit[0]

Fit_Gauss=gaussian_gen(En,A0,sigma0, mu0,Bg0)
# plt.figure()
# plt.plot(En[LimC0-LimDel:LimC0+LimDel],matrix[LimC0-LimDel:LimC0+LimDel, Izero])
# plt.plot(En[LimC0-LimDel:LimC0+LimDel],Fit_Gauss[LimC0-LimDel:LimC0+LimDel])
# plt.show()

#%% 

A = []
sigma = []
mu = []
Bg = []
mid=[]
vfit= []
v0 = voltage[Izero]
A.append(A0)
sigma.append(sigma0)
mu.append(mu0)
Bg.append(Bg0)
mid.append(midi)
vfit.append(v0)
iBreakLow=0
for i in range(Izero):
    vi = voltage[Izero-i]
    if vi<(np.round(VLimLow,2) ):
#pad with zeros
        mu.append(0)
        Bg.append(0)
        mid.append(0)
        sigma.append(0)
        vfit.append(vi)
        H.append(0)   
    else:
        try:
            if mu[i]==0.:
                midi = mid[i]
            else:
                midi = mu[i]
            LimC0=np.where(Wav>=1239.8/midi)[0][0]
    #        LimDel=20
    #        midi = En[LimC0]
            LimDel = -10*(np.where(Wav>=1239.8/(midi+sigma[i]))[0][0] - LimC0)
            Lim1 = LimC0-LimDel
            Lim2 = LimC0+LimDel
            Gauss_fit=curve_fit(gaussian_gen,En[Lim1:Lim2],matrix[Lim1:Lim2,Izero-i],
                        p0=(A[i],sigma[i],midi,0))
            Ai,sigmai,mui,Bgi = Gauss_fit[0]
            midi =En[Lim1+np.argmax(matrix[Lim1:Lim2,Izero-i])]
            print('i, amp, sigma, mu, background', Ai,sigmai, mui, Bgi)
            
            # Fit_Gauss=gaussian_gen(En,Hi,sigmai, mui,Bgi)
            # plt.figure()
            # plt.plot(En[Lim1:Lim2],matrix_bg[Lim1:Lim2, Izero-i])
            # plt.plot(En[Lim1:Lim2],Fit_Gauss[Lim1:Lim2])
            # plt.title(str(i)+ '  ' + str(Izero-i) + ' ' + str(vi) + ' volts')
            # plt.show()
            # plt.pause(1)
            # plt.close()
            
            A.append(Ai)
            sigma.append(sigmai)
            mu.append(mui)
            Bg.append(Bgi)
            mid.append(midi)
            vfit.append(vi)
           
        except RuntimeError:
            print(i,'error in fit - just use peak')
            midi = En[Lim1+np.argmax(matrix[Lim1:Lim2,Izero-i])]
            mui = 0.
            Ai=0.
            sigmai=0.
            Bgi=0.
            A.append(Ai)
            sigma.append(sigmai)
            mu.append(mui)
            Bg.append(Bgi)
            mid.append(midi)
            vfit.append(vi)



plt.figure()
plt.plot(mu)
plt.plot(mid)
plt.show()