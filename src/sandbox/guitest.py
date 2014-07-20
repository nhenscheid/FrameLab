"""
This is a test of Tkinter for gui
"""

import Tkinter as tk
import tkMessageBox

top = tk.Tk()

def helloCallBack():
    tkMessageBox.showinfo("Hello, python", "Hello, world.")
    
    
B = tk.Button(top,text = "Hello", command = helloCallBack)

B.pack()
top.mainloop()