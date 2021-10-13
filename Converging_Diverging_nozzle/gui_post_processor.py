import numpy as np
import matplotlib.pyplot as plt
from tkinter import *
import time
import pandas as pd
from tkinter import filedialog




class PostProcessor:
    def buttonClose(self):
        self.window.destroy()

    def get_vals(self):
        return np.array(self.selecteditems1), np.array(self.selecteditems2)

    def buttonPrint(self):
        self.selecteditems1 = self.box1.curselection()
        self.selecteditems2 = self.box2.curselection()

    def __init__(self, data1, data2):
        self.window = Tk()
        self.selecteditems1 = ()
        self.selecteditems2 = ()
        self.box1 = Listbox(self.window, exportselection=0, height=10, width=50, selectmode='multiple')
        for num in data1:
            self.box1.insert(END, num)
        self.box1.place(x=10, y=150)

        self.box2 = Listbox(self.window, exportselection=0, height=10, width=50, selectmode='multiple')
        for num in data2:
            self.box2.insert(END, num)
        self.box2.place(x=300, y=150)

        self.label = Label(self.window, text="Available plots", fg='red', font=("Helvetica", 16))
        self.label.place(x=60, y=100)

        self.label2 = Label(self.window, text="Available tables", fg='red', font=("Helvetica", 16))
        self.label2.place(x=350, y=100)

        self.btn = Button(self.window, text="Confirm", fg='blue', font=('helvetica', 12, 'bold'),
                          command=lambda: [self.buttonPrint(), self.buttonClose()])
        self.btn.place(x=260, y=350)

        items = map(int, self.box1.curselection())
        print(items)

        self.window.title('Post Processor')
        self.window.geometry("600x400+10+10")

    def run(self):
        self.window.mainloop()
#
#
# data1 = ["direct quantities", "residuals", "m.flow vs distance", "dens vs Mach vs distance"]
# data2 = ["one", "two", "three", "four"]
# a = ['pressure', 'velocity', 'density', 'temperature']
#
# obj = PostProcessor(data1, data2, a)
# obj.run()
# selection1, selection2 = obj.get_vals()



def plotter(x, y, title, x_name, y_name, color, labels):
    if len(y) <= 6:
        for i in range(len(y)):
            plt.plot(x, y[i], color=color[i], label=labels[i])
    else:
        plt.plot(x, y[:], color=color[0], label=labels)
    plt.title(title)
    plt.xlabel(x_name)
    plt.ylabel(y_name)
    plt.legend()



def tabloter(titles, quantities, index1, flag):

    if index1 == 0:
        table1 = {titles[0]: quantities[0], titles[1]: quantities[1],
                titles[2]: quantities[2],titles[3]: quantities[3],
                titles[4]: quantities[4],titles[5]: quantities[5],
                titles[6]: quantities[6],titles[7]: quantities[7]}
        df = pd.DataFrame(table1, columns=[titles[0], titles[1],
                                           titles[2], titles[3],
                                           titles[4], titles[5],
                                           titles[6], titles[7]])
        print(df)
        df.to_excel(r'C:\Users\Dimitris Terzis\Desktop\Πανεπιστήμιο\Μαθήματα\10 εξάμηνο\cfd\3d_table1.xlsx', index = False, header=True)
    elif index1 == 1:
        if flag == 41:
            table1 = {titles[0]: quantities[0], titles[1]: quantities[1],
                      titles[2]: quantities[2], titles[3]: quantities[3],
                      titles[4]: quantities[4], titles[5]: quantities[5],
                      titles[6]: quantities[6],  titles[7]: quantities[7]}
            df = pd.DataFrame(table1, columns=[titles[0], titles[1],
                                               titles[2], titles[3],
                                               titles[4], titles[5],
                                               titles[6], titles[7]])
            print(df)
            df.to_excel(r'C:\Users\Dimitris Terzis\Desktop\Πανεπιστήμιο\Μαθήματα\10 εξάμηνο\cfd\3d_table2.xlsx',
                        index=False, header=True)
        else:
            pass
    else:
        table2 = {titles[0]: quantities[0], titles[1]: quantities[1],
                  titles[2]: quantities[2], titles[3]: quantities[3],
                  titles[4]: quantities[4], titles[5]: quantities[5],
                  titles[6]: quantities[6], titles[7]: quantities[7]}
        df = pd.DataFrame(table2, columns=[titles[0], titles[1],
                                           titles[2], titles[3],
                                           titles[4], titles[5],
                                           titles[6], titles[7]])
        print(df)
        df.to_excel(r'C:\Users\Dimitris Terzis\Desktop\Πανεπιστήμιο\Μαθήματα\10 εξάμηνο\cfd\3d_table3.xlsx', index = False, header=True)

