
import numpy as np
import matplotlib.pyplot as plt

# functions that allow to save coordinates of the point 
# you click on with the mouse and print them
# to interact with the figure you need to set:
# tools/preferences/graphics to be QT5 not inline

def findY(f,ax,matrix):
    
    """This function allows to choose the peak to be fitted by clicking on its 
    sides.
    
    Parameters:
        f: figure to be clicked on
        ax: axis along which the points will be picked
        matrix: 
    """
    print('open figure and click WIDELY on either side of peak of interest')
    
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

    def pickX(event):
        
        """This function allows to select the voltage range on which we want to 
        perform the fit by picking two points on the figure.
        
        Parameters:
            f : figure on which point are picked
            
        Returns: 
            voltage range limits
        """
        xmouse = event.xdata #takes x coordinate of the point clicked
        print ('x of mouse: {:.2f}'.format(xmouse))
        times.append(xmouse)
        if len(times)==2:
            f.canvas.mpl_disconnect(cid)
        return times
    
    times=[]     
    print('open figure and click on two points to define the voltage range you want to fit')
    cid=f.canvas.mpl_connect('button_press_event',pickX)
    while len(times)<2:
        plt.pause(5)

    print('limits', times)

    VLim1 = np.min(times)
    VLim2 = np.max(times)
    VLimLow = np.where(ax<=VLim1)[0][-1] 
    VLimHi = np.where(ax>=VLim2)[0][0] 

    return VLimLow, VLimHi
