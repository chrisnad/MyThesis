from Tkinter import *
import tkMessageBox
import os, numpy, pylab, glob, matplotlib.pyplot, random

#def addPlot():
    


app = Tk()
app.title("Curve analyser")
app.geometry("450x300+200+200")

w1 = Label(text = "name of the file").pack()
b = Button(text = "Plot", command = addPlot).pack(side = LEFT)

w2 = Label(text = "change something").pack(side = LEFT)
e = Entry(app, bd = 1).pack(side = LEFT)

#c = Canvas(app, bg = "white", height = 250, width = 300)
#coord = 10, 50, 240, 210
#arc = c.create_arc(coord, start=0, extent=150, fill="red")
#c.pack()


mainloop()
