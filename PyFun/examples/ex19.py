from Tkinter import *  # AFAIK Tkinter is always capitalized
#import os
#import easygui as eg

class App:
    characterPrefix = "character_"
    def __init__(self, master):
        self.master = master  # You'll want to keep a reference to your root window
        frame = Frame(master)
        frame.pack()

        # character box
        Label(frame, text = "Characters Editor").grid(row = 0, column = 0, rowspan = 1, columnspan = 2)
        self.charbox = Listbox(frame)  # You'll want to keep this reference as an attribute of the class too.
        for chars in []:
            self.charbox.insert(END, chars)
        self.charbox.grid(row = 1, column = 0, rowspan = 5)
        charadd = Button(frame, text = "   Add   ", command = self.addchar).grid(row = 1, column = 1)
        charremove = Button(frame, text = "Remove", command = self.removechar).grid(row = 2, column = 1)
        charedit = Button(frame, text = "    Edit    ", command = self.editchar).grid(row = 3, column = 1)

    def addchar(self, initialCharacter='', initialInfo=''):
        t = Toplevel(root)  # Creates a new window
        t.title("Add character")
        characterLabel = Label(t, text="Character name:")
        characterEntry = Entry(t)
        characterEntry.insert(0, initialCharacter)
        infoLabel = Label(t, text="Info:")
        infoEntry = Entry(t)
        infoEntry.insert(0, initialInfo)
        def create():
            characterName = characterEntry.get()
            self.charbox.insert(END, characterName)
            with open(app.characterPrefix + characterName, 'w') as f:
                    f.write(infoEntry.get())
            t.destroy()
        createButton = Button(t, text="Create", command=create)
        cancelButton = Button(t, text="Cancel", command=t.destroy)

        characterLabel.grid(row=0, column=0)
        infoLabel.grid(row=0, column=1)
        characterEntry.grid(row=1, column=0)
        infoEntry.grid(row=1, column=1)
        createButton.grid(row=2, column=0)
        cancelButton.grid(row=2, column=1)

    def removechar(self):
        for index in self.charbox.curselection():
            item = self.charbox.get(int(index))
            self.charbox.delete(int(index))
            try:
                os.remove(characterPrefix + item)
            except IOError:
                print "Could not delete file", characterPrefix + item
    #def editchar(self):
        # You can implement this one ;)

root = Tk()
root.wm_title("IA Development Kit")
app = App(root)
root.mainloop()
