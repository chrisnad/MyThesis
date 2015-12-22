from Tkinter import*
 
class MenuBar(Frame):
    """Barre de menus deroulants"""
    def __init__(self, boss =None):
        Frame.__init__(self, borderwidth =2, relief =GROOVE)
 
    #---------------------------------------------------------------------------------#
 
        ##### Menu <Client> #####
 
        fileMenu_client = Menubutton(self, text ="Client")
        fileMenu_client.grid(row=0,column=0)
 
        # Partie "deroulante" :
        menu_client = Menu(fileMenu_client,tearoff=0)
        menu_client.add_command(label ="Ajout")#,command = boss.effacer)
        menu_client.add_command(label="Supression")#, command=self.donothing)
        menu_client.add_command(label="Modification")#, command=self.donothing)
 
        # Integration du menu :
        fileMenu_client.configure(menu = menu_client )
 
    #---------------------------------------------------------------------------------#
 
        ##### Menu <Devis> #####
 
        fileMenu_devis = Menubutton(self, text ="Devis")
        fileMenu_devis.grid(row=0,column=1)
 
        # Partie "deroulante" :
        menu_devis = Menu(fileMenu_devis,tearoff=0)
        menu_devis.add_command(label="Creation")#,command = boss.effacer)
        menu_devis.add_command(label="Modification")#, command=self.donothing)
        menu_devis.add_command(label="Supression")#, command=self.donothing)
 
        # Integration du menu :
        fileMenu_devis.configure(menu = menu_devis)
 
    #---------------------------------------------------------------------------------#
 
        ##### Menu <facture> #####
 
        fileMenu_facture = Menubutton(self, text ="Facture")
        fileMenu_facture.grid(row=0,column=2)
 
        # Partie "deroulante" :
        menu_facture = Menu(fileMenu_facture,tearoff=0)
        menu_facture.add_command(label="Recherche")#,command = boss.effacer)
        menu_facture.add_command(label="Creation")#, command=self.donothing)
 
 
        # Integration du menu :
        fileMenu_facture.configure(menu = menu_facture)
 
    #---------------------------------------------------------------------------------#
 
        ##### Menu <stock> #####
 
        fileMenu_stock = Menubutton(self, text ="Stock")
        fileMenu_stock.grid(row=0,column=3)
 
        # Partie "deroulante" :
        menu_stock = Menu(fileMenu_stock, tearoff=0)
        menu_stock .add_command(label="Stock")#, command=donothing)
        menu_stock .add_command(label="Ajout")#, command=donothing)
        menu_stock .add_command(label="Modification")#, command=donothing)
        menu_stock .add_command(label="Supression")#, command=donothing)
 
 
        # Integration du menu :
        fileMenu_stock.configure(menu = menu_stock)
 
    #---------------------------------------------------------------------------------#
 
        ##### Menu <aide> #####
 
        fileMenu_aide = Menubutton(self, text ="Aide")
        fileMenu_aide.grid(row=0,column=4)
 
        # Partie "deroulante" :
        menu_aide = Menu(fileMenu_aide, tearoff=0)
        menu_aide.add_command(label="Client")#, command=donothing)
        menu_aide.add_command(label="Devis")#, command=donothing)
        menu_aide.add_command(label="Facture")#, command=donothing)
        menu_aide.add_command(label="Stock")#, command=donothing)
 
 
        # Integration du menu :
        fileMenu_aide.configure(menu = menu_aide)

class Test(MenuBar,Frame):
 
 
    def __init__(self, boss =None):
        Frame.__init__(self)
        self.grid()
 
        self.master.title('Fenetre avec menus')
        mBar = MenuBar(self)
        mBar.grid(row=0,column=0)
 
        self.master.geometry("400x300+200+200")
        #self.master.config(bg ="cadet blue")
 
        bouton= n() # pour signaler que la command concernant le bouton est dans la class n
        buton=Button(text="Quitter",fg="red",command=bouton.nou)
        buton.grid(row=1,column=0)
        #self.master.destroy()
 
 
class n(MenuBar):
 
    def __init__(self, boss =None):
        Frame.__init__(self)
        self.grid()
 
 
    def nou(self):
 
        #self.f=Tk()
        #self.f.title('Fenetre  menus')
        fenetre=Tk()
        mBar = MenuBar(self)
        mBar.grid(row=0,column=0)
 
        buton=Button(text="Quitter",fg="red")
        buton.grid(row=1,column=0)
 
        print("couou")
        fenetre.mainloop()
 
 
 
if __name__ =="__main__":              # --- Programme de test ---
 
    Test().mainloop()
    
