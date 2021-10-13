import numpy as np
import matplotlib.pyplot as plt
import time as ti
from tkinter import messagebox


def riti(y, Re, sfalma, omega, tropos):
    dy = 1 / y
    dt = 0.5 * Re * dy ** 2
    y_tot = np.linspace(0, 1, y + 2)
    n = 1
    u = np.zeros(len(y_tot))
    u[-1] = 1
    logic = True
    velocity_profiles = u.copy()
    plt.figure(1)
    u_avg = np.zeros(len(y_tot))
    u_prev = np.zeros(len(y_tot))
    while logic:
        for j in range(1, len(y_tot) - 1):
            if tropos == 1:
                u_avg[j] = (1 - (2 * dt / (Re * dy ** 2))) * u[j] + (dt / (Re * (dy ** 2))) * (u_prev[j + 1] + u[j - 1])
                u[j] = u_prev[j] + omega * (u_avg[j] - u_prev[j])
            else:
                u[n, j] = u[n - 1, j] + (dt / (Re * (dy ** 2))) * (u[n - 1, j + 1] - 2 * u[n - 1, j] + u[n - 1, j - 1])
        sigklisi = np.abs(u_prev[:] - u[:]).copy
        u_prev = u.copy()
        if all(sigklisi <= sfalma):
            logic = False
        n = n + 1
        velocity_profiles = np.vstack([velocity_profiles, u])
    n = n - 1
    time = dt * n
    plt.plot(velocity_profiles[0], y_tot, color='r', label='Δt = %f ' % dt)
    plt.plot(velocity_profiles[int(0.25 * n)], y_tot, color='b', label='Δt = 25% της μόνιμης κατάστασης')
    plt.plot(velocity_profiles[int(0.5 * n)], y_tot, color='g', label='Δt = 50% της μόνιμης κατάστασης')
    plt.plot(velocity_profiles[int(0.75 * n)], y_tot, color='y', label='Δt = 75% της μόνιμης κατάστασης')
    plt.plot(velocity_profiles[int(n)], y_tot, color='m', label=' Mόνιμη κατάσταση')
    plt.title('Προφίλ ταχυτήτων ροής Couette με πεπλεγμένη μέθοδο')
    plt.xlabel('u/Ue', fontsize=13)
    plt.ylabel('y/D', fontsize=13)
    plt.text(0.5, 0.32, 'Συνολικός Χρόνος για μόνιμη ροή: %.3f' % time + ' sec', fontsize=8)
    plt.text(0.5, 0.36, 'Re = %.0f' % Re, fontsize=8)
    plt.text(0.5, 0.42, 'Υπολογιστικά χρονικά βήματα: %.0f' % n, fontsize=8)
    plt.text(0.5, 0.48, 'Υπολογιστικά χωρικά βήματα: %.0f' % y, fontsize=8)
    plt.legend()
    plt.show()


def peplegmeni(y, Re, E, sfalma, omega, tropos):
    dy = 1 / (y)
    dt = E * Re * (dy) ** 2
    y_tot = np.linspace(0, 1, y + 2)
    u = np.zeros(len(y_tot))
    u[-1] = 1
    logic = True
    n = 0
    velocity_profiles = u.copy()
    while logic:
        u_prev = u.copy()
        K = np.zeros(y)
        n = n + 1
        A = -(E / 2) * np.ones(y)
        B = (1 + E) * np.ones(y)
        C = -(E / 2) * np.ones(y)
        for j in range(1, len(y_tot) - 1):
            K[j - 1] = (1 - E) * u_prev[j] + E/2 * (u_prev[j + 1] + u_prev[j - 1])
        K[-1] = K[-1] - C[-1] * u[-1]
#Thomas
        for i in range(1, len(A)):
            mc = A[i] / B[i - 1]
            B[i] = B[i] - mc * C[i - 1]
            K[i] = K[i] - mc * K[i - 1]
        Thomas = A.copy()
        Thomas[-1] = K[-1] / B[-1]
        for j in range(len(A)-2, -1, -1):
            Thomas[j] = (K[j] - C[j] * Thomas[j + 1]) / B[j]
        for j in range(1, len(y_tot) - 1):
            if tropos == 1:
                u[j] = (omega * Thomas[j - 1]) + (1 - omega) * u_prev[j]
            else:
                u[j] = Thomas[j - 1]
        velocity_profiles = np.vstack([velocity_profiles, u])
        sigkrisi = abs(u - u_prev)
        if all(sigkrisi <= sfalma):
            logic = False
    plt.figure(2)
    time = n * dt
    plt.plot(velocity_profiles[0], y_tot, color='r', label='Δt = %f ' % dt)
    plt.plot(velocity_profiles[int(0.25 * n)], y_tot, color='b', label='Δt = 25% της μόνιμης κατάστασης')
    plt.plot(velocity_profiles[int(0.5 * n)], y_tot, color='g', label='Δt = 50% της μόνιμης κατάστασης')
    plt.plot(velocity_profiles[int(0.75 * n)], y_tot, color='y', label='Δt = 75% της μόνιμης κατάστασης')
    plt.plot(velocity_profiles[int(n)], y_tot, color='m', label=' Mόνιμη κατάσταση')
    plt.title('Προφίλ ταχυτήτων ροής Couette με πεπλεγμένη μέθοδο')
    plt.xlabel('u/Ue', fontsize=13)
    plt.ylabel('y/D', fontsize=13)
    plt.text(0.5, 0.32, 'Συνολικός Χρόνος για μόνιμη ροή: %.3f' % time + ' sec', fontsize=8)
    plt.text(0.5, 0.36, 'Re = %.0f' % Re, fontsize=8)
    plt.text(0.5, 0.42, 'Υπολογιστικά χρονικά βήματα: %.0f' % n, fontsize=8)
    plt.text(0.5, 0.48, 'Υπολογιστικά χωρικά βήματα: %.0f' % y, fontsize=8)
    plt.legend()
    plt.show()


control = 0
while control == 0:
    option = float(input("Define solving method [1 for explicid, 2 for implicid]: "))
    if option == 1 or option == 2:
        control = 1
    else:
        messagebox.showinfo("Wrong input", "Please define a proper solving method (1 or 2)")
if option == 1:
    control = 0
    while control == 0:
        y = int(input("Define number of grid points: [>1]"))
        if y > 1:
            control = 1
        else:
            messagebox.showinfo("Wrong input", "Please define a proper number of grid points: ")
    control = 0
    while control == 0:
        Re = float(input("Define Reynolds number [>0]"))
        if Re > 0 :
            control = 1
        else:
            messagebox.showinfo("Wrong input", "Please define a proper number : ")
    control = 0
    while control == 0:
        error = float(input("Define desired error threshold: [<1]"))
        if 1 > error > 0:
            control = 1
        else:
            messagebox.showinfo("Wrong input", "Please define a proper error threshold: [<1]")
    control = 0
    while control == 0:
        relax = input("Do u wish to add relaxation method?[y, n]")
        if relax == 'y':
            tropos = 1
            control = 1
            control2 = 0
            while control2 == 0:
                omega = float(input("Define desired omega factor: [0 < ω < 2]"))
                if 2 > omega > 0:
                    control2 = 1
                else:
                    messagebox.showinfo("Wrong input", "Please define a proper value: [0 < ω < 2]")
        elif relax == 'n':
            tropos = 0
            control = 1
        else:
            messagebox.showinfo("Wrong input", "Please give a proper answer: [y, n]")
    riti(y, Re, error, omega, tropos)
else:
    control = 0
    while control == 0:
        y = int(input("Define number of grid points: [>1]"))
        if y > 1:
            control = 1
        else:
            messagebox.showinfo("Wrong input", "Please define a proper number of grid points: ")
    control = 0
    while control == 0:
        Re = float(input("Define Reynolds number [>0]"))
        if Re > 0:
            control = 1
        else:
            messagebox.showinfo("Wrong input", "Please define a proper number : ")
    control = 0
    while control == 0:
        E = float(input("Define desired E coefficient: "))
        if 4000 > E > 1 :
            control = 1
        else:
            messagebox.showinfo("Wrong input", "Please define a proper error threshold: [<1]")
    control = 0
    while control == 0:
        error = float(input("Define desired error threshold: [<1]"))
        if 1 > error > 0:
            control = 1
        else:
            messagebox.showinfo("Wrong input", "Please define a proper error threshold: [<1]")
    control = 0
    while control == 0:
        relax = input("Do u wish to add relaxation method?[y, n]")
        if relax == 'y':
            tropos = 1
            control = 1
            control2 = 0
            while control2 == 0:
                omega = float(input("Define desired omega factor: [0 < ω < 2]"))
                if 2 > omega > 0:
                    control2 = 1
                else:
                    messagebox.showinfo("Wrong input", "Please define a proper value: [0 < ω < 2]")
        elif relax == 'n':
            tropos = 0
            control = 1
        else:
            messagebox.showinfo("Wrong input", "Please give a proper answer: [y, n]")
    peplegmeni(y, Re, E, error, omega, tropos)





