from Tkinter import *
 
root=Tk()
e=StringVar()
ent=Entry(root,textvariable=e)

print e.get()
 
b=Button(root,text='PRINT',command="exit")
b.pack()
ent.pack()
 
root.mainloop()
