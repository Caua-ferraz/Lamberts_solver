import tkinter as tk
from tkinter import ttk, messagebox
import sys
import math
from main import LambertSolver, vector_norm, vector_subtract, orbital_energy, propagate_orbit

class RedirectText:
    def __init__(self, text_widget):
        self.text_widget = text_widget

    def write(self, string):
        self.text_widget.insert(tk.END, string)
        self.text_widget.see(tk.END)
        print(string, end='')  # Also print to terminal

    def flush(self):
        pass

class LambertSolverGUI:
    def __init__(self, master):
        self.master = master
        master.title("Lambert Solver")
        master.geometry("600x500")

        self.mu = 398600.4418  # Earth's gravitational parameter (km^3/s^2)
        self.solver = LambertSolver(self.mu)

        self.create_widgets()

    def create_widgets(self):
        # Input frame
        input_frame = ttk.Frame(self.master, padding="10")
        input_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        # r1 inputs
        ttk.Label(input_frame, text="r1 (km):").grid(row=0, column=0, sticky=tk.W)
        self.r1_x = ttk.Entry(input_frame, width=10)
        self.r1_x.grid(row=0, column=1)
        self.r1_y = ttk.Entry(input_frame, width=10)
        self.r1_y.grid(row=0, column=2)
        self.r1_z = ttk.Entry(input_frame, width=10)
        self.r1_z.grid(row=0, column=3)

        # r2 inputs
        ttk.Label(input_frame, text="r2 (km):").grid(row=1, column=0, sticky=tk.W)
        self.r2_x = ttk.Entry(input_frame, width=10)
        self.r2_x.grid(row=1, column=1)
        self.r2_y = ttk.Entry(input_frame, width=10)
        self.r2_y.grid(row=1, column=2)
        self.r2_z = ttk.Entry(input_frame, width=10)
        self.r2_z.grid(row=1, column=3)

        # Time of flight input
        ttk.Label(input_frame, text="Time of flight (s):").grid(row=2, column=0, sticky=tk.W)
        self.dt = ttk.Entry(input_frame, width=10)
        self.dt.grid(row=2, column=1)

        # Solve button
        self.solve_button = ttk.Button(input_frame, text="Solve", command=self.solve)
        self.solve_button.grid(row=3, column=0, columnspan=4, pady=10)

        # Output text widget
        self.output_text = tk.Text(self.master, wrap=tk.WORD, width=70, height=20)
        self.output_text.grid(row=1, column=0, padx=10, pady=10, sticky=(tk.W, tk.E, tk.N, tk.S))

        # Redirect stdout to both terminal and text widget
        sys.stdout = RedirectText(self.output_text)

        # Configure grid weights
        self.master.columnconfigure(0, weight=1)
        self.master.rowconfigure(1, weight=1)

    def solve(self):
        try:
            r1 = [float(self.r1_x.get()), float(self.r1_y.get()), float(self.r1_z.get())]
            r2 = [float(self.r2_x.get()), float(self.r2_y.get()), float(self.r2_z.get())]
            dt = float(self.dt.get())

            if dt <= 0:
                raise ValueError("Time of flight must be positive")

            v1, v2 = self.solver.solve(r1, r2, dt)

            self.output_text.delete(1.0, tk.END)  # Clear previous output
            print("Initial velocity vector (km/s):", v1)
            print("Final velocity vector (km/s):", v2)

            # Additional calculations
            e1 = orbital_energy(r1, v1, self.mu)
            e2 = orbital_energy(r2, v2, self.mu)
            print(f"Initial orbital energy: {e1:.2f} km^2/s^2")
            print(f"Final orbital energy: {e2:.2f} km^2/s^2")
            print(f"Energy difference: {abs(e1 - e2):.2e} km^2/s^2")

            r2_prop, v2_prop = propagate_orbit(r1, v1, dt, self.mu)
            r2_error = vector_norm(vector_subtract(r2, r2_prop))
            v2_error = vector_norm(vector_subtract(v2, v2_prop))
            print(f"Position error after propagation: {r2_error:.2f} km")
            print(f"Velocity error after propagation: {v2_error:.2f} km/s")

            transfer_dist = vector_norm(vector_subtract(r2, r1))
            t_transfer = (transfer_dist**1.5 / math.sqrt(8 * self.mu)) * math.pi
            print(f"Estimated minimum transfer time: {t_transfer:.2f} s")
            print(f"Input transfer time: {dt:.2f} s")

        except ValueError as e:
            messagebox.showerror("Error", str(e))
        except Exception as e:
            messagebox.showerror("Unexpected Error", f"An unexpected error occurred: {str(e)}")

if __name__ == "__main__":
    root = tk.Tk()
    app = LambertSolverGUI(root)
    root.mainloop()