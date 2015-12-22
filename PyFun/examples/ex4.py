from Tkinter import *
import os, numpy, pylab, glob
 
def addPlot(x, y):
    pylab.plot(x, y)
    pylab.show()
 
root = Tk()

b = Button(root, text='  ADD A PLOT  ', command = addPlot([1,2,3],[1,4,9])).pack()

Label(root, textvariable = counter).pack()


mainloop()
