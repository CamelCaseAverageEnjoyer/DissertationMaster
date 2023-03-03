from tkinter import *
from tkinter import ttk
from all_objects import *
from PIL import Image, ImageTk
from tkinter_files.assembly0 import *


class Window(Tk):
    def __init__(self):
        super().__init__()
        self.title("Проект Г.И.Б.О.Н.")
        self.geometry("1000x685+200+100")


def click_button_assembly0():
    from tkinter_files.assembly0 import click_button_assembly
    global root
    root.destroy()
    click_button_assembly()


def click_button_test0():
    from tkinter_files.testing0 import click_button_test
    global root
    root.destroy()
    click_button_test()


def click_button_talk():
    global root
    talk()


def click_button_plot():
    from tkinter_files.plotting0 import click_button_plot
    global root
    root.destroy()
    click_button_plot()


def click_button_info0():
    from tkinter_files.information import click_button_info
    global root
    root.destroy()
    click_button_info()


def return_home():
    global root
    root = Tk()
    root.title("Проект Г.И.Б.О.Н.")
    root.geometry("1000x690+200+100")
    root.minsize(1000, 690)
    root.maxsize(1000, 690)
    photo_icon = PhotoImage(file="icons/satellite.png")
    root.iconphoto(True, photo_icon)

    img = ImageTk.PhotoImage(Image.open("icons/robot_talk1.png"))
    b = Label(image=img)
    b.grid(row=0, column=0, columnspan=5)

    photo_talk = PhotoImage(file="icons/discussing.png").subsample(10, 10)
    photo_test = PhotoImage(file="icons/processing.png").subsample(10, 10)
    photo_plot = PhotoImage(file="icons/vision.png").subsample(10, 10)
    photo_home = PhotoImage(file="icons/home.png").subsample(10, 10)
    photo_assembly = PhotoImage(file="icons/solution.png").subsample(10, 10)
    photo_idea = PhotoImage(file="icons/idea.png").subsample(10, 10)

    btn0 = Button(text="Поболтать", command=click_button_talk, image=photo_talk, compound=LEFT)
    btn0.grid(row=1, column=0, padx='7', pady='7')
    btn1 = Button(text="Начать сборку", command=click_button_assembly0, image=photo_assembly, compound=LEFT)
    btn1.grid(row=1, column=1, padx='7', pady='7')
    btn2 = Button(text="Тестировка", command=click_button_test0, image=photo_test, compound=LEFT)
    btn2.grid(row=1, column=2, padx='7', pady='7')
    btn3 = Button(text="Графики", command=click_button_plot, image=photo_plot, compound=LEFT)
    btn3.grid(row=1, column=3, padx='7', pady='7')
    btn4 = Button(text="Что я такое?", command=click_button_info0, image=photo_idea, compound=LEFT)
    btn4.grid(row=1, column=4, padx='7', pady='7')

    root.mainloop()


if __name__ == '__main__':
    return_home()
