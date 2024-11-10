import numpy as np
from matplotlib import pyplot as plt
c = 0.45 # chord length
# density calculator https://www.engineeringtoolbox.com/air-density-specific-weight-d_600.html
rho = 1.237
# kinematic viscosity calculator https://www.engineeringtoolbox.com/air-absolute-kinematic-viscosity-d_601.html
v = 14.34

# define file path to data set
file_path = ('group_8/group_8_test_4.dat')

def extract_data(file_path):
    with open(file_path, 'r') as f:
        # Read the first row
        first_line = f.readline().strip().split()

        # Extract AoA from the third column (index 2) and Uinf from the seventh column (index 6)
        AoA = float(first_line[2])
        Uinf = float(first_line[6])

        # Load the remaining data, skipping the first two rows
        data = np.loadtxt(file_path, delimiter=' ', skiprows=2)

        # Split the data into x position, y position, and pressure difference arrays
        x = np.array(data[:, 0])
        y = np.array(data[:, 1])
        p = np.array(data[:, 2])
    return AoA, Uinf, x, y, p

def compute_cp(x, p, Uinf):
    cp = np.zeros_like(p)
    for i in range(len(p)):
        cp[i] = (p[i]) / (0.5 * rho * Uinf)
        print(cp[i])
    return cp

def reynolds(Uinf):
    return Uinf * c / v

def plot_cp(x, cp, AoA, Re):
    plt.figure(figsize=(10, 6))
    plt.plot(x/c, cp, label='NASA exp', color='black', linestyle='-', marker='o', markersize=4)
    plt.xlabel('x/c')
    plt.ylabel('Cp')
    plt.title(f"NACA018 Pressure Coefficient, AoA = {AoA}°, Re = {Re}")
    plt.gca().invert_yaxis()
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.show()

def main():
    AoA, Uinf, x, y, p = extract_data(file_path)
    Re = reynolds(Uinf)
    cp = compute_cp(x, p, Uinf)
    plot_cp(x, cp, AoA, Re)

main()





