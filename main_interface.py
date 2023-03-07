from tkinter_files.assembly0 import *
from tkinter_files.tk_functions import *


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
    icons = Icons()

    root.title("Проект Г.И.Б.О.Н.")
    root.geometry("1000x690+200+100")
    root.minsize(1000, 690)
    root.maxsize(1000, 690)
    root.iconphoto(True, icons.icon)

    img = ImageTk.PhotoImage(Image.open("icons/robot_talk1.png"))
    b = Label(image=img)
    b.grid(row=0, column=0, columnspan=5)

    btn0 = Button(text="Поболтать", command=click_button_talk, image=icons.talk, compound=LEFT)
    btn1 = Button(text="Начать сборку", command=click_button_assembly0, image=icons.assembly, compound=LEFT)
    btn2 = Button(text="Тестировка", command=click_button_test0, image=icons.test, compound=LEFT)
    btn3 = Button(text="Графики", command=click_button_plot, image=icons.plot, compound=LEFT)
    btn4 = Button(text="Что я такое?", command=click_button_info0, image=icons.idea, compound=LEFT)
    btn0.grid(row=1, column=0, padx='7', pady='7')
    btn1.grid(row=1, column=1, padx='7', pady='7')
    btn2.grid(row=1, column=2, padx='7', pady='7')
    btn3.grid(row=1, column=3, padx='7', pady='7')
    btn4.grid(row=1, column=4, padx='7', pady='7')

    root.mainloop()


if __name__ == '__main__':
    return_home()
