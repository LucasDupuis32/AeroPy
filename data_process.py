import numpy as np
from matplotlib import pyplot as plt
import math
from scipy.integrate import quad

c = 0.45 # chord length
s = 1 # wing span
# density calculator https://www.engineeringtoolbox.com/air-density-specific-weight-d_600.html
rho = 1.237
# kinematic viscosity calculator https://www.engineeringtoolbox.com/air-absolute-kinematic-viscosity-d_601.html
v = 14.34e-6
S = 2.5 * 1.8 # section area of the wind tunnel
naca0018area = 0.123289 # surface area of a unit chord NACA0018, given by XFoil 6.99
naca0018area = naca0018area * c ** 2




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

def compute_cp(p, Uinf):
    cp = np.zeros_like(p)
    for i in range(len(p)):
        cp[i] = (p[i]) / (0.5 * rho * Uinf ** 2)
    return cp

def reynolds(Uinf):
    return Uinf * c / v

def plot_cp(x, cp, AoA, Re, y):
    plt.figure(figsize=(10, 6))
    plt.plot(x/c, cp, color='black', linestyle='-', marker='o', markersize=4)
    plt.xlabel('x/c')
    plt.ylabel('Cp')
    plt.title(f"NACA018 Pressure Coefficient, AoA = {AoA}°, Re = {round(Re)}")
    plt.gca().invert_yaxis()
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.show()

def compute_cl(x, cp, y):
    cp_u = []
    cp_l = []
    x_u = []
    x_l = []
    # separate lower and upper data points
    for i in range((len(y))):
        if y[i] > 0:
            cp_u.append(cp[i])
            x_u.append(x[i])
        elif y[i] < 0:
            cp_l.append(cp[i])
            x_l.append(x[i])
    cp_u = np.flip(np.array(cp_u))
    cp_l = np.array(cp_l)
    x_u = np.flip(np.array(x_u))
    x_l = np.array(x_l)

    # There are 20 data points on the lower side, and 19 on the upper side.
    # -> separate the integral in two
    cl = np.trapezoid(cp_l, x_l) - np.trapezoid(cp_u, x_u)

    return cl

def blockage_correction(U, Re, V, S, K):
    e = (K * V) / S ** (3 / 2)
    U = U * (1 + e)
    Re = Re * (1 + e)
    return U, Re

def main():
    # print the AoA and Re for every data set
    for i in range(0, 8):
        file_path = (f'group_8/group_8_test_{i}.dat')
        AoA, Uinf, x, y, p = extract_data(file_path)
        Re = reynolds(Uinf)
        V = naca0018area * s
        _, Re = blockage_correction(Uinf, Re, V, S, K=0.52)
        print(f"data set {i} : AoA = {AoA}, Uinf = {Uinf}, Re = {Re}")

    # compute c_l for specific data set
    file_path = ('group_8/group_8_test_7.dat')
    AoA, Uinf, x, y, p = extract_data(file_path)
    Re = reynolds(Uinf)
    V = naca0018area * s # volume of the work piece
    Uinf, Re = blockage_correction(Uinf, Re, V, S, K=0.52) # apply solid blockage correction

    print()
    print(f"Reynolds number: Re = {Re:.0f}")
    print(f"Angle of attack: alpha = {AoA}°")
    cp = compute_cp(p, Uinf)
    cl = compute_cl(x, cp, y)
    print(f"Lift coefficient: c_l = {cl:.4f}")
    plot_cp(x, cp, AoA, Re, y)

if __name__ == '__main__':
    main()





