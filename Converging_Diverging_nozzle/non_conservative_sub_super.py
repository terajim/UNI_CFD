import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math
import time
import pandas as pd
from tkinter import *
import gui_post_processor
plt.style.use('ggplot')

L = 3
gamma = 1.4

# use ggplot style for more sophisticated visuals
plt.style.use('ggplot')


def live_plotter(x_vec, y1_data, line1, identifier='', pause_time=0.1):
    if line1 == []:
        # this is the call to matplotlib that allows dynamic plotting
        plt.ion()
        fig = plt.figure(figsize=(13, 6))
        ax = fig.add_subplot(111)
        # create a variable for the line so we can later update it
        line1, = ax.plot(x_vec, y1_data, '-o', alpha=0.8)
        # update plot label/title
        plt.ylabel('Y Label')
        plt.title('Title: {}'.format(identifier))
        plt.show()

    # after the figure, axis, and line are created, we only need to update the y-data
    line1.set_ydata(y1_data)
    # adjust limits if new data goes beyond bounds
    if np.min(y1_data) <= line1.axes.get_ylim()[0] or np.max(y1_data) >= line1.axes.get_ylim()[1]:
        plt.ylim([np.min(y1_data) - np.std(y1_data), np.max(y1_data) + np.std(y1_data)])
    # this pauses the data so the figure/axis can catch up - the amount of pause can be altered above
    plt.pause(pause_time)

    # return line so we can update it again in the next iteration
    return line1

def nonconservative(n, C, interations, sfalma):
    x_40 = np.linspace(0, L, 41)
    A_40 = 1 + 2.2 * (x_40 - 1.5)**2

    dx = L / n
    x_total = np.linspace(0, L, n)
    A_x = 1 + 2.2 * (x_total - 1.5) ** 2
    dt = np.zeros([interations, len(x_total)])
    time_table = np.zeros([2, interations])

    dens = np.zeros([interations, len(x_total)])
    V = np.zeros([interations, len(x_total)])
    T = np.zeros([interations, len(x_total)])

    ddens_dt = np.zeros([interations, len(x_total)])
    dV_dt = np.zeros([interations, len(x_total)])
    dT_dt = np.zeros([interations, len(x_total)])

    dens_bar = np.zeros([interations, len(x_total)])
    V_bar = np.zeros([interations, len(x_total)])
    T_bar = np.zeros([interations, len(x_total)])

    ddens_bar_dt = np.zeros([interations, len(x_total)])
    dV_bar_dt = np.zeros([interations, len(x_total)])
    dT_bar_dt = np.zeros([interations, len(x_total)])

    ddens_avg = np.zeros([interations, len(x_total)])
    dV_avg = np.zeros([interations, len(x_total)])
    dT_avg = np.zeros([interations, len(x_total)])
    a = np.zeros([interations, len(x_total)])

    dens[0, :] = 1 - 0.3146 * x_total[:]
    T[0, :] = 1 - 0.2314 * x_total[:]
    V[0, :] = (0.1 + 1.09 * x_total[:]) * T[0, :] ** 0.5
    a[0, :] = np.sqrt(T[0, :])

    data = np.array([x_total[:], A_x[:], V[0, :], dens[0, :], T[0, :]])
    data = data.T
    df = pd.DataFrame(data)


    T[:, 0] = 1
    T_bar[:, 0] = 1
    dens[:, 0] = 1
    dens_bar[:, 0] = 1
    k = 0
    line1 = []
    y_vector = V[0, 16]
    dt[0, :] = C * dx / (a[0, :] + V[0, :])
    #dt = min(C * (dx / (T[k, :] ** 0.5 + V[k, :])))


    #analytical solution
    M_an = np.array([0.09782, 0.10656, 0.11646, 0.12769, 0.14049, 0.15513 ,0.17195,
                     0.19134, 0.21377, 0.23981, 0.27012, 0.30549, 0.3468, 0.39505,
                     0.45128, 0.51657, 0.59186, 0.67786, 0.77485, 0.88254, 1 , 1.12566,
                     1.25755, 1.39348, 1.53136, 1.66934, 1.80597, 1.94015, 2.07115,
                     2.17853, 2.32205, 2.44163, 2.5572, 2.66915, 2.77733, 2.88201, 2.98336,
                     3.08156, 3.17678, 3.26919, 3.35895])

    p_an = (1 + (gamma - 1) / 2 * M_an[:]**2) ** (-gamma / (gamma - 1))
    dens_an = (1 + (gamma - 1) / 2 * M_an[:]**2) **(-1/(gamma - 1))
    T_an = (1 + (gamma - 1) / 2 * M_an[:]**2) ** (-1)
    m_flow_an = dens_an[:] * A_40[:] * M_an[:] * np.sqrt(T_an[:])

    control = True

    while k <= interations - 2 and control:

        for i in range(1, len(x_total) - 1):
            dt = min(C * (dx / ((T[k, :]) ** 0.5 + V[k, :])))

            ddens_dt[k, i] = -dens[k, i] * ((V[k, i + 1] - V[k, i]) / dx) - dens[k, i] * V[k, i] * (
                    (np.log(A_x[i + 1]) - np.log(A_x[i])) / dx) - V[k, i] * ((dens[k, i + 1] - dens[k, i]) / dx)

            dV_dt[k, i] = -V[k, i] * (V[k, i + 1] - V[k, i]) / dx - (1 / gamma) * (
                    (T[k, i + 1] - T[k, i]) / dx + (T[k, i] / dens[k, i]) * (dens[k, i + 1] - dens[k, i]) / dx)

            dT_dt[k, i] = -V[k, i] * ((T[k, i + 1] - T[k, i]) / dx) - (gamma - 1) * T[k, i] * (
                    ((V[k, i + 1] - V[k, i]) / dx) + V[k, i] * ((np.log(A_x[i + 1]) - np.log(A_x[i])) / dx))

            dens_bar[k + 1, i] = dens[k, i] + ddens_dt[k, i] * dt
            V_bar[k + 1, i] = V[k, i] + dV_dt[k, i] * dt
            T_bar[k + 1, i] = T[k, i] + dT_dt[k, i] * dt

        # boundary Conditions
            V_bar[k + 1, 0] = V[k, 0]
            V_bar[k, 0] = V[k, 0]

        #V_bar[k, -1] = 2 * V_bar[k, -2] - V_bar[k, -3]
        #T_bar[k, -1] = 2 * T_bar[k, -2] - T_bar[k, -3]
        #dens_bar[k, -1] = 2 * dens_bar[k, -2] - dens_bar[k, -3]

            ddens_bar_dt[k + 1, i] = -dens_bar[k + 1, i] * ((V_bar[k + 1, i] - V_bar[k + 1, i - 1]) / dx) - \
                             dens_bar[k + 1, i] * V_bar[k + 1, i] * (np.log(A_x[i]) - np.log(A_x[i - 1])) / dx - \
                                     V_bar[k + 1, i] * ((dens_bar[k + 1, i] - dens_bar[k + 1, i - 1]) / dx)

            dV_bar_dt[k + 1, i] = -V_bar[k + 1, i] * ((V_bar[k + 1, i] - V_bar[k + 1, i - 1]) / dx) -\
                                  (1 / gamma) * ((T_bar[k + 1, i] - T_bar[k + 1, i - 1] ) / dx +
                                       T_bar[k + 1, i] / dens_bar[k + 1, i] * (dens_bar[k + 1, i] - dens_bar[k + 1, i - 1]) / dx)

            dT_bar_dt[k + 1, i] = -V_bar[k + 1, i] * ((T_bar[k + 1, i] - T_bar[k + 1, i - 1]) / dx) -\
                          (gamma - 1) * T_bar[k + 1, i] * (((V_bar[k + 1, i] - V_bar[k + 1, i - 1]) / dx) + V_bar[k + 1, i] *
                                                   ((np.log(A_x[i]) - np.log(A_x[i - 1])) / dx))

            ddens_avg[k + 1, i] = 0.5 * (ddens_dt[k, i] + ddens_bar_dt[k + 1, i])
            dV_avg[k + 1, i] = 0.5 * (dV_dt[k, i] + dV_bar_dt[k + 1, i])
            dT_avg[k + 1, i] = 0.5 * (dT_dt[k, i] + dT_bar_dt[k + 1, i])

            dens[k + 1, i] = dens[k, i] + ddens_avg[k + 1, i] * dt
            V[k + 1, i] = V[k, i] + dV_avg[k + 1, i] * dt
            T[k + 1, i] = T[k, i] + dT_avg[k + 1, i] * dt
            a[k + 1, i] = np.sqrt(abs(T[k + 1, i]))
            ######      dt[k + 1, i] = C * (dx / (a[k + 1, i] + V[k + 1, i]))

        # boundary conditions

        V[k + 1, 0] = 2 * V[k + 1, 1] - V[k + 1, 2]
        #####       dt[k + 1, 0] = C * (dx / (T[k + 1, 0] ** 0.5 + V[k + 1, 0]))
        a[k + 1, 0] = np.sqrt(T[k + 1, 0])

        dens[k + 1, -1] = 2 * dens[k + 1,  -2] - dens[k + 1,  -3]
        T[k + 1, -1] = 2 * T[k + 1,  -2] - T[k + 1, -3]
        V[k + 1, -1] = 2 * V[k + 1,  -2] - V[k + 1, -3]
        ######      dt[k + 1, -1] = C * (dx / (a[k + 1, -1] + V[k + 1,-1]))

        data = np.array([x_total[:], A_x[:], V[k + 1, :], dens[k + 1, :], T[k + 1, :]])
        data = data.T
        df = pd.DataFrame(data)
        # print(df)

        time_table[0, k] = k
        time_table[1, k] = dt
        if all(np.abs(V[k + 1,:] - V[k, :]) <=sfalma) and all(np.abs(T[k + 1,:] - T[k, :]) <=sfalma) and all(np.abs(dens[k + 1,:] - dens[k, :]) <=sfalma):
            control = False
        k += 1
    print(k)
    a = np.sqrt(T)
    M = V[0:k, :] / a[0:k, :]
    p = dens[0:k, :] * T[0:k, :]
    A = np.zeros([k, len(x_total)])
    A[:, :] = A_x[:]
    k = k -1
########################################## PLOTS ######################################################


    plot_titles = ['Flow Field quantities vs number of iterations', 'Residuals vs number of iterations',
                   'Mass flow different timesteps vs nozzle non-d. distance', 'Density over Mach vs non-d. distance']

    x_labels = ['Number of time steps', 'Number of timesteps', 'Non-dimensional distance [x/L]',
                'Non-dimensional distance [x/L]']

    y_labels = ['Non-Dimensional quantities', 'Residuals', 'Non-dimensional mass flow', 'Non-dimensional density']

    colorpallete = ['r', 'b', 'y', 'm', 'g']


    line_labels = ['nond. density', 'nond. temperature', 'nond. pressure', 'Mach number',
                   '(dρ/dt)avg', '(dT/dt)avg', '(dV/dt)avg', '0.05t_total', '0.1t_total', '0.2t_total', '0.7t_total', 'convergence',
                   'nond. density', 'Mach number']
    table_titles = ['Calculated converged numerical values ', 'Calculated and analytical values']

    values_titles = ['x/L', 'A/A*', 'V/a', 'ρ/ρο', 'T/Tο', 'p/po', 'M', 'm_dot']
    values2_titles = ['x/L', 'A/A*',  'ρ/ρο', 'ρ/ρο(analytic)', '%difference', 'M', 'M(analytical)', '%difference']

    tables = [x_total, A_x, V[k, :], dens[k, :], T[k, :], p[k, :], M[-1, :],  dens[k, :] * V[k, :] * A_x]
    if n == 41:
        tables2 = [x_total, A_x, dens[k, :], dens_an, np.abs(dens[k, :] - dens_an) / dens_an, M[-1, :], M_an,
               np.abs(M[-1, :] - M_an) / M_an]

    plots = [dens[0:k, int(n / 2)], T[0:k, int(n / 2)], p[0:k, int(n / 2)], M[0:k, int(n / 2)],
             np.abs(ddens_avg[0:k, int(n / 2)]), np.abs(dT_avg[0:k, int(n / 2)]), np.abs(dV_avg[0:k, int(n / 2)])
        , dens[int(.05 * k), :] * V[int(.05 * k), :] * A[int(.05 * k), :],
             dens[int(.1 * k), :] * V[int(.1 * k), :] * A[int(.1 * k), :],
             dens[int(.2 * k), :] * V[int(.2 * k), :] * A[int(.2 * k), :],
             dens[int(.7 * k), :] * V[int(.7 * k), :] * A[int(.7 * k), :], dens[k, :] * V[k, :] * A[-1, :], dens[k, :],
             M[k, :]]

    obj = gui_post_processor.PostProcessor(plot_titles, table_titles)
    obj.run()
    selection1, selection2 = obj.get_vals()

    for i in range (0, len(selection1)):
        plt.figure(i + 1)
        if selection1[i] == 0:
            gui_post_processor.plotter(time_table[0, 0:k], plots[0:4], plot_titles[0], x_labels[0], y_labels[0],
                                       colorpallete, line_labels[0:4])
        elif selection1[i] == 1:
            gui_post_processor.plotter(time_table[0, 0:k], plots[4:7], plot_titles[1], x_labels[1], y_labels[1],
                                       colorpallete, line_labels[4:7])
            #plt.ylim(0.000001, 0.01)
        elif selection1[i] == 2:
            gui_post_processor.plotter(x_total, plots[7:12], plot_titles[2], x_labels[2], y_labels[2],
                                       colorpallete, line_labels[7:12])
            plt.scatter(0, m_flow_an[2], color='red', s=20, label='Ακριβής τιμή')
            plt.legend()
        else:
            fig, ax1 = plt.subplots()
            ax2 = ax1.twinx()
            ax1.plot(x_total, plots[12], 'r')
            ax2.plot(x_total, plots[-1], 'b-')

            ax1.set_xlabel(x_labels[3])
            ax1.set_ylabel(y_labels[3], color='r')
            ax2.set_ylabel(y_labels[3], color='b')
            plt.title(plot_titles[3])
            ax1.scatter(x_40, M_an, color='black', s=10)
            ax2.scatter(x_40, dens_an, color='black', s=10)

    for i in range(0, len(selection2)):
        if selection2[i] == 0:
            gui_post_processor.tabloter(values_titles, tables, selection2[i], n)
        elif selection2[i] == 1 and n == 41:
            gui_post_processor.tabloter(values2_titles, tables2, selection2[i], n)
    total_time = sum(time_table[1, :])
    print(total_time, 'Συνολικός χρόνος')
    print(dens_an[20], T_an[20], p_an[20], M_an[20])
    print(dens[k, int(n/2)], T[k, int(n/2)], p[k, int(n/2)], M[k, int(n/2)])
    print('total time:', np.sum(time_table[1, :]))
    plt.show()

#nonconservative(41, 0.70, 4000, 10**(-5))
