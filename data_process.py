import numpy as np
from matplotlib import pyplot as plt

c = 0.45 # chord length
# density calculator https://www.engineeringtoolbox.com/air-density-specific-weight-d_600.html
rho = 1.237
# kinematic viscosity calculator https://www.engineeringtoolbox.com/air-absolute-kinematic-viscosity-d_601.html
v = 14.34e-6

# define file path to data set
file_path = ('group_8/group_8_test_4.dat')
#file_path = ('common_data/common_test_3.dat')

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


def main():
    AoA, Uinf, x, y, p = extract_data(file_path)
    Re = reynolds(Uinf)
    print(f"Reynolds number: Re = {Re:.0f}")
    print(f"Angle of attack: alpha = {AoA}°")
    cp = compute_cp(p, Uinf)
    cl = compute_cl(x, cp, y)
    print(f"Lift coefficient: c_l = {cl:.4f}")
    plot_cp(x, cp, AoA, Re, y)

if __name__ == '__main__':
    main()





