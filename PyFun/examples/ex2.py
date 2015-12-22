from Tkinter import *
 
root=Tk()
e=StringVar()
ent=Entry(root,textvariable=e)
def action():
	print e.get(),type(e.get())
 
b=Button(root,text='PRINT',command=action)
b.pack()
ent.pack()
 
root.mainloop()
