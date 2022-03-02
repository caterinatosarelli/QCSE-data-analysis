import numpy as np
import matplotlib.pyplot as plt

# Functions that allow to save coordinates of the point 
# one clicks on with the mouse and print them.
# To interact with the figure one needs to set:
# Tools/Preferences/Graphics to be QT5 not inline

def findY(f,ax,matrix):
    """
    This function allows to choose the peak to be fitted by clicking on its 
    sides.
    
    Parameters
    ----------
    f : figure to be clicked on.
    ax : axis along which the points will be picked.
    matrix : matrix of the data 

    Returns
    -------
    Zoom_matrix : matrix reduced in the limits chosen by the user.
    Zlim1 : left side of the peak.
    Zlim2 : right side of the peak.

    """
   
    print('Open the figure and click widely on either side of peak of interest.')
    
    def pickY(event):
        ymouse = event.ydata # takes y coordinate of the point clicked
        print ('y of mouse: {:.2f}'.format(ymouse))
        times.append(ymouse)
        if len(times)==2:
            f.canvas.mpl_disconnect(cid)
        return times
    
    times=[]
    
    cid=f.canvas.mpl_connect('button_press_event',pickY)
    while len(times)<2:
        plt.pause(5)

    print('limits', times)
    
    ZLimLow = np.min(times)
    ZLimHi = np.max(times)
    ZLim1=np.where(ax<=ZLimLow)[0][-1] 
    ZLim2=np.where(ax>=ZLimHi)[0][0] 
   
    Zoom_matrix = matrix[ZLim1:ZLim2,:]
    
    return Zoom_matrix, ZLim1, ZLim2

def findX(f,ax):
    """
    This function allows to select the voltage range on which we want to 
        perform the fit by picking two points on the figure.
    
    Parameters
    ----------
    f : figure on which point are picked.
    ax : array of voltage values.

    Returns
    -------
    VLimLow : index of lowest voltage value selected
    VLimHi : index of highest voltage value selected.

    """

    def pickX(event):
    
        xmouse = event.xdata #takes x coordinate of the point clicked
        print ('x of mouse: {:.2f}'.format(xmouse))
        times.append(xmouse)
        if len(times)==2:
            f.canvas.mpl_disconnect(cid)
        return times
    
    times=[]     
    print('Open the figure and click on two points to define the voltage range you want to fit.')
    cid=f.canvas.mpl_connect('button_press_event',pickX)
    while len(times)<2:
        plt.pause(5)

    print('limits', times)

    VLim1 = np.min(times)
    VLim2 = np.max(times)
    VLimLow = np.where(ax<=VLim1)[0][-1] 
    VLimHi = np.where(ax>=VLim2)[0][0] 

    return VLimLow, VLimHi

