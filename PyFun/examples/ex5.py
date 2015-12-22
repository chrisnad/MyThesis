from Tkinter import *

app = Tk()
app.title("Curve analyser")
app.geometry("450x300+200+200")

scrollbar = Scrollbar(app)
scrollbar.pack( side = RIGHT, fill=Y )

mylist = Listbox(app, yscrollcommand = scrollbar.set )
for i in range(6):
   mylist.insert(END, str(i))

mylist.pack( side = RIGHT )
scrollbar.config( command = mylist.yview )

mainloop()
