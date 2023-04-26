# Import the required libraries
import tkinter as tk
from tkinter import *

import matplotlib
import numpy as np
import scipy.constants as C

matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# gui part
root = tk.Tk()  # Displays a window
root.geometry("450x250")  # Sets the initial window size
root.title('Rutherford scattering simulation')  # Sets the window title

# Alpha particle mass and atomic number
m_alpha = C.physical_constants["alpha particle mass energy equivalent in MeV"][0]
m_alpha_kg = C.physical_constants["alpha particle mass"][0]
Z_alpha = 2
s_target = StringVar()
s_target.set(79)
E_target = Entry(root, width=20, textvariable=s_target)

# Constant coefficient:
# k1=e^2/4\pi\varepsilon_0
# k2=Z_\alpha Z_Au/m_\alpha
MeV = 1e6 * C.e
Z_Au = 79
k1 = 1.44  # fm.Mev

# The differential equation used to solve the four columns are the x-positive velocity, x-positive acceleration, y-positive velocity, and y-positive acceleration of the particle
def differential_equations(y, t, k2):
    return np.array([y[1], k1 * k2 * y[0] / (np.sqrt((y[0]) ** 2 + (y[2]) ** 2)) ** 3, y[3],
                     k1 * k2 * y[2] / (np.sqrt((y[0]) ** 2 + (y[2]) ** 2)) ** 3])


# A function of coordinates (x,y) at different times t solved by initial conditions and differential equations
def solving(val, b, k2):
    initial_value = np.array([-3e3, val, b, 0])
    t = np.arange(0, 1e-18, (1e-20))
    res = odeint(differential_equations, initial_value, t, args=(k2,))
    x, y = res[:-1, 0], res[:-1, 2]
    return x, y


# Create a list of colors so that you can choose from different colors for subsequent drawings
colorlist1 = ['skyblue', 'cyan', 'palegreen', 'red', 'hotpink', 'gold', 'brown']

# Create a matrix of the required size, the number of rows is the number of data points solved for the same incident energy, and the number of columns is the number of groups set by different initial values
x_m = np.ones((100, 6), dtype='float')  # m is the matrix required for manual input of the drawing
y_m = np.ones((100, 6), dtype='float')
x_r = np.ones((100, 10000), dtype='float')  # r is the matrix required for randomly generated plotting
y_r = np.ones((100, 10000), dtype='float')

v_x = np.ones((100, 10000), dtype='float')  # r is the matrix required for randomly generated plotting
v_y = np.ones((100, 10000), dtype='float')
tan_theta = np.ones((1, 10000), dtype='float')

# Create three empty dictionaries, set the composition text, and set the initial values of incidence energy and aiming distance to 0
energy = {}
aimdistance = {}
particlerank = {}
for i in range(0, 6):
    j = 1 + i
    s1 = StringVar()
    s1.set(0)
    s2 = StringVar()
    s2.set(0)
    energy['E' + str(i)] = Entry(root, width=20, textvariable=s1)
    aimdistance['A' + str(i)] = Entry(root, width=20, textvariable=s2)
    particlerank['R' + str(i)] = Label(root, text='alpha particle%d' % j)

S1 = Scale(root, from_=10, to=60, tickinterval=6, resolution=10, length=120)


# Manually enter the drawing button binding function
def draw(event):
    Z_target = float(s_target.get())
    k2 = Z_alpha * Z_target / m_alpha * (C.c * 1e15) ** 2  # fm^2/(s^2 Mev)
    a = [0, 0, 0, 0, 0, 0]
    Energycompare = np.array([0])
    SUM = 0
    for i in range(0, 6): # Cycle to different incidence initial conditions
        Energy = float(energy['E' + str(i)].get()) * MeV  # Gets the user's manually entered incident energy value
        v0 = np.sqrt(2 * Energy / m_alpha_kg) * 1e15  # The initial velocity of incidence is calculated from the incident energy
        b = float(aimdistance['A' + str(i)].get())  # Gets the user manually entering the aiming distance value
        x_m[:, i], y_m[:, i] = solving(v0, b, k2)  # Use the solution function defined earlier to solve the coordinates corresponding to the initial conditions of group i, and import the data into the i-th column of the matrix for backup
        a[i] = b
        SUM = SUM + Energy

    for j in range(
            len(x_m)):  # The number of rows of the matrix loops, and the loop order of the first row and then the column makes the image appear on the picture with different initial conditions of incidence at the same time
        plt.cla()  # clear the chart
        plt.plot(0, 0, 'ko', ms=10)  # target particle position
        # plt.text(2,-2,'target particle',fontsize=20)#target particle position
        if max(a) < 2000:
            plt.xlim((-3e3, 3e3))  # The abscissa range setting
            plt.ylim((-3e3, 3e3))  # Ordinate axis range setting
        else:
            plt.xlim((-3e3, 3e3))  # The abscissa range setting
            plt.ylim((-1.5 * max(a), 1.5 * max(a)))  # Ordinate axis range setting
        plt.xticks(fontsize=7)
        plt.yticks(fontsize=7)
        plt.title('alpha particle scattering\n time: %d0zs' % (j), fontsize=15)  # Image caption
        plt.xlabel("x axis/fm", fontsize=15)  # The abscissa title
        plt.ylabel("y axis/fm")  # The vertical axis title
        for k in range(0, 6):  # Loop through the number of columns of the matrix
            plt.plot(x_m[:j + 1, k], y_m[:j + 1, k], color=colorlist1[k],
                     linewidth=1)  # Draw a curve, and each time you draw from the first row of data to the j row, the trajectory of the particle is preserved
            plt.scatter(x_m[j, k], y_m[j, k], color='red', s=15)  # Plot the latest point
        # plt.pause(1/S1.get())
        plt.pause(0.01)
        # plt.pause(20000*SUM)
    plt.show()


# Randomly generate a drawing button binding function
def random_draw(event):
    theta_hist = []
    missing_list = []
    Z_target = float(s_target.get())
    k2 = Z_alpha * Z_target / m_alpha * (C.c * 1e15) ** 2  # fm^2/(s^2 Mev)
    cycle_index = 10000
    for i in range(0, cycle_index):
        Energy = float(np.random.rand() * 10 + 1) * MeV  # Randomly generate different incident energies
        v0 = np.sqrt(2 * Energy / m_alpha_kg) * 1e15
        b = float(np.random.rand() * 2000 - 1000)  # Different aiming distances are randomly generated
        x_r[:, i], y_r[:, i] = solving(v0, b,
                                       k2)  # Use the solution function defined earlier to solve the coordinates corresponding to the initial conditions of group i, and import the data into the i-th column of the matrix for backup
        # v_x[:,i],v_y[:,i]=solving_speed(v0,b)#Use the solution function defined earlier to solve the velocity corresponding to the initial conditions of group i, and import the data into the ith column of the matrix for later use
        # tan_theta[:,i]=(v_y[:,i]/v_x[:,i])
    # Random Plot
    figure, axes = plt.subplots(nrows=1, ncols=2, figsize=(18, 9))
    for j in range(len(x_m)):
        # The same principle as manually entering a drawing
        plt.subplot(1, 2, 1)
        plt.cla()  # Clear the chart
        plt.plot(0, 0, 'ko', ms=10)  # Target particle position
        # plt.text(2,-2,'target particle',fontsize=20)
        plt.xlim((-3e3, 3e3))
        plt.ylim((-3e3, 3e3))
        plt.title('alpha particle scattering\n time: %d0zs' % (j))
        plt.xlabel("x axis/fm")
        plt.ylabel("y axis/fm")
        # plt.plot(x_r[:j+1,:], y_r[:j+1,:],'--', color='orange',alpha=0.1,linewidth=1)# Draw a curve, and each time you draw from the first row of data to the j row, the trajectory of the particle is preserved
        plt.scatter(x_r[j, :], y_r[j, :], color='red', s=15)  # Plot the latest point
        plt.subplot(1, 2, 2)
        plt.cla()
        plt.xlim((0, 180))  # The abscissa range setting
        plt.ylim((0, 1000))  # Ordinate axis range settings
        for k in range(0, cycle_index):
            if abs(x_r[j, k]) > 3e3 or abs(y_r[j, k]) > 3e3:
                if k not in missing_list:
                    tan_theta[0, k] = ((y_r[j + 1, k] - y_r[j, k]) / (x_r[j + 1, k] - x_r[j, k]))
                    if np.arctan(tan_theta[0, k]) < 0 and x_r[j, k] < 0:  # Second quadrant
                        theta_hist.append((np.arctan(tan_theta[0, k]) + np.pi) / np.pi * 180)
                    elif np.arctan(tan_theta[0, k]) < 0 and x_r[j, k] > 0:  # Forth quadrant
                        theta_hist.append(abs(np.arctan(tan_theta[0, k])) / np.pi * 180)
                    elif np.arctan(tan_theta[0, k]) > 0 and x_r[j, k] > 0:  # First quadrant
                        theta_hist.append((np.arctan(tan_theta[0, k])) / np.pi * 180)
                    elif np.arctan(tan_theta[0, k]) > 0 and x_r[j, k] < 0:  # Third quadrant
                        theta_hist.append((np.pi - np.arctan(tan_theta[0, k])) / np.pi * 180)
                    missing_list.append(k)
        plt.hist(theta_hist, bins=180)
        n, bins, pach = plt.hist(theta_hist, bins=180)
        k = n[0] * np.sin(bins[1] * np.pi / 180) ** 2
        plt.xlabel(r"$\theta$/°")
        plt.ylabel(r"n($\theta$)")
        plt.plot(bins[1::], k * (np.sin(bins[1::] * np.pi / 360)) ** -2)
        left, bottom, width, height = 0.6, 0.3, 0.28, 0.5
        plt.axes([left, bottom, width, height])
        bigtheta = np.array(theta_hist)[np.array(theta_hist) > 45]
        plt.hist(bigtheta, bins=180 - 45)
        plt.xlabel(r"$\theta$/°")
        plt.ylabel(r"n($\theta$)")
        plt.pause(0.001)
    plt.show()


# Mini Program interface layout
L1 = Label(root, text='Kinetic energies(MeV)')
L2 = Label(root, text='Impact parameters(fm)')

L4 = Label(root, text='nuclear charge Z')
L1.grid(row=0, column=1)
L2.grid(row=0, column=2)

L4.grid(row=7, column=0)
B1 = Button(root, text='Run with parameters\n in the input box', width=20, fg='ivory', bg='saddlebrown')
B1.bind('<Button-1>', draw)  # Bind manual input drawing functions
B2 = Button(root, text='Run with\n random parameters', width=20, fg='ivory', bg='saddlebrown')
B2.bind('<Button-1>', random_draw)  # Bind a randomly generated drawing function
B1.grid(row=8, column=1)
B2.grid(row=8, column=2)

E_target.grid(row=7, column=1, columnspan=2)
# Use loops to lay out the input box
for i in range(0, 6):
    particlerank['R' + str(i)].grid(row=i + 1, column=0)
    energy['E' + str(i)].grid(row=i + 1, column=1)
    aimdistance['A' + str(i)].grid(row=i + 1, column=2)
root.mainloop()  # Display the window
