from Tkinter import *
master = Tk()
Label(master, text="Input:").grid(row=0)

e1 = Entry(master, width = 100)
e1.grid(row=0, column=1, columnspan=30)

Button(master, text='Q').grid(row=3, column=2)

Button(master, text='C').grid(row=3, column=1)

Button(master, text='Confirm').grid(row=3, column=0)
mainloop()
