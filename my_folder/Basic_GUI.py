
import numpy as np
import matplotlib.pyplot as plt
from tkinter import Tk, Label, Entry, Button, Frame
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from tkinter import ttk

class SineWavePlotter:
    def __init__(self, root):
        self.root = root
        self.root.title("Sine Wave Plotter")
        self.root.configure(bg="#282C35")

        self.plot_frame = Frame(root, bg="#282C35")
        self.plot_frame.pack(side="left", padx=20, pady=20)

        self.figure = Figure(figsize=(10,6), dpi=100, facecolor='#282C35') 
        self.plot_area = self.figure.add_subplot(111)
        self.plot_area.set_facecolor('#1E1E1E')
        self.plot_area.grid(color='gray', linestyle='dashed')

        self.canvas = FigureCanvasTkAgg(self.figure, master=self.plot_frame)
        self.canvas.get_tk_widget().pack()

        self.params_frame = Frame(root, bg="#1E1E1E", bd=2, relief="solid")
        self.params_frame.pack(side="right", padx=20, pady=20)

        self.a_label = Label(self.params_frame, text="a:", bg="#1E1E1E", fg="white")
        self.a_label.pack()

        self.a_entry = Entry(self.params_frame)
        self.a_entry.pack()

        self.b_label = Label(self.params_frame, text="b:", bg="#1E1E1E", fg="white")
        self.b_label.pack()

        self.b_entry = Entry(self.params_frame)
        self.b_entry.pack()
        
        
        self.w_label = Label(self.params_frame, text="w:", bg="#1E1E1E", fg="white")
        self.w_label.pack()

        self.w_entry = Entry(self.params_frame)
        self.w_entry.pack()
        
        
        
        
        
        

        self.generate_button = Button(self.params_frame, text="Generate", command=self.generate_plot, bg="#1E90FF", fg="white")
        self.generate_button.pack()

    def generate_plot(self):
        try:
            a = float(self.a_entry.get())
            b = float(self.b_entry.get())
            w = float(self.w_entry.get())
            
            if a >= b:
                raise ValueError("a should be less than b")

            x = np.linspace(a, b, 1000)
            #w = 1.0  # You can adjust this value as needed
            y = np.sin(w * x)

            
            
            self.plot_area.clear()
            self.plot_area.plot(x, y, color='orange')
            
            self.plot_area.set_xlabel('x', color='white')  # Set x-axis label color
            self.plot_area.set_ylabel('sin(wx)', color='white')  # Set y-axis label color
            self.plot_area.set_title('Sine Wave Plot', color='white')  # Set title color
            
            self.plot_area.grid(color='gray', linestyle='dashed')
            self.plot_area.spines['bottom'].set_color('white')  # Set x-axis color
            self.plot_area.spines['left'].set_color('white')    # Set y-axis color
            
            # Customize tick label colors
            self.plot_area.tick_params(axis='x', colors='white')
            self.plot_area.tick_params(axis='y', colors='white')
            
            self.canvas.draw()


        except ValueError as e:
            self.show_error(str(e))

    def show_error(self, message):
        error_label = Label(self.root, text=message, fg="red")
        error_label.pack()

if __name__ == "__main__":
    root = Tk()
    root.configure(bg="#282C35")
    app = SineWavePlotter(root)
    root.mainloop()
