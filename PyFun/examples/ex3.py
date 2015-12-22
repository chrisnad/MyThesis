import matplotlib.pyplot as plt
import random
import Tkinter
#import sys
 
 
def quitApp():
    plt.close()
    root.quit()
    sys.exit()
 
def interrupt():
    counter.set(counter.get() + 1)
    root.after(300, interrupt)
 
def addPlot():
    ax = plt.subplot(111)
    xdata = range(256)
    ydata = [random.randint(0,255) for counter in xdata]
    ax.plot(xdata, ydata)
    plt.show()
 
root = Tkinter.Tk()
Tkinter.Button(root, text='  ADD A PLOT  ', command=addPlot).grid()
counter = Tkinter.IntVar()
Tkinter.Label(root, textvariable=counter).grid()
root.protocol('WM_DELETE_WINDOW', quitApp)
root.after(0, interrupt)
root.mainloop()
