
import numpy as np
import matplotlib.pyplot as plt

# functions that allow to save coordinates of the point 
# you click on with the mouse and print them
# to interact with the figure you need to set:
# tools/preferences/graphics to be QT5 not inline

def findY(f):
    
    def pickY(event):
        ymouse = event.ydata # takes y coordinate of the point clicked
        print ('y of mouse: {:.2f}'.format(ymouse))
        times.append(ymouse)
        if len(times)==2:
            f.canvas.mpl_disconnect(cid)
        return times
    
    times=[]
    print('open figure and click WIDELY on either side of peak of interest')
    cid=f.canvas.mpl_connect('button_press_event',pickY)
    while len(times)<2:
        plt.pause(5)

    print('limits', times)
   
    return times

def findX(f):

    def pickX(event):
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
   
    return times
