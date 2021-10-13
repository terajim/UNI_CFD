import numpy as np
import matplotlib.pyplot as plt
import math
import time
import pandas as pd
import tkinter.messagebox as tk
from inspect import currentframe, getframeinfo
import gui_post_processor
gamma = 1.4
L = 3
plt.style.use('ggplot')
def conservative(n, C, interations, dp, Cx, sfalma):
    x_40 = np.linspace(0, L, 41)
    A_40 = 1 + 2.2 * (x_40 - 1.5) ** 2


    x_total = np.linspace(0, L, n)
    dx = L / n
    A_x = np.zeros(len(x_total))
    # dt = np.zeros(len(x_total))

    # Arxikopoiish
    time_table = np.zeros([2, interations])



    U1 = np.zeros([interations, len(x_total)])
    U2 = np.zeros([interations, len(x_total)])
    U3 = np.zeros([interations, len(x_total)])
    F1 = np.zeros([interations, len(x_total)])
    F2 = np.zeros([interations, len(x_total)])
    F3 = np.zeros([interations, len(x_total)])
    J2 = np.zeros([interations, len(x_total)])

    V = np.zeros([interations, len(x_total)])
    dens = np.zeros([interations, len(x_total)])
    T = np.zeros([interations, len(x_total)])

    dU1_dt = np.zeros([interations, len(x_total)])
    dU2_dt = np.zeros([interations, len(x_total)])
    dU3_dt = np.zeros([interations, len(x_total)])

    U1_bar = np.zeros([interations, len(x_total)])
    U2_bar = np.zeros([interations, len(x_total)])
    U3_bar = np.zeros([interations, len(x_total)])

    dens_bar = np.zeros([interations, len(x_total)])
    T_bar = np.zeros([interations, len(x_total)])
    V_bar = np.zeros([interations, len(x_total)])

    F1_bar = np.zeros([interations, len(x_total)])
    F2_bar = np.zeros([interations, len(x_total)])
    F3_bar = np.zeros([interations, len(x_total)])

    dU1_bar_dt = np.zeros([interations, len(x_total)])
    dU2_bar_dt = np.zeros([interations, len(x_total)])
    dU3_bar_dt = np.zeros([interations, len(x_total)])

    dU1_avg = np.zeros([interations, len(x_total)])
    dU2_avg = np.zeros([interations, len(x_total)])
    dU3_avg = np.zeros([interations, len(x_total)])

    S1 = np.zeros([interations, len(x_total)])
    S2 = np.zeros([interations, len(x_total)])
    S3 = np.zeros([interations, len(x_total)])

    S1_bar = np.zeros([interations, len(x_total)])
    S2_bar = np.zeros([interations, len(x_total)])
    S3_bar = np.zeros([interations, len(x_total)])

    # initial conditions
    for i in range(len(x_total)):
        A_x[i] = 1 + 2.2 * (x_total[i] - 1.5) ** 2
        if x_total[i] <= 0.5:
            dens[0, i] = 1
            T[0, i] = 1
        elif x_total[i] <= 1.5:
            dens[0, i] = 1 - 0.366 * (x_total[i] - 0.5)
            T[0, i] = 1 - 0.167 * (x_total[i] - 0.5)
        elif x_total[i] <= 2.1:
            dens[0, i] = 0.634 - 0.702 * (x_total[i] - 1.5)
            T[0, i] = 0.833 - 0.4908 * (x_total[i] - 1.5)
        else:
            dens[0, i] = 0.5892 + 0.10228 * (x_total[i] - 2.1)
            T[0, i] = 0.93968 + 0.0622 * (x_total[i] - 2.1)

    # Boundary conditions
    T[:, 0] = 1
    T_bar[:, 0] = 1
    dens[:, 0] = 1
    dens_bar[:, 0] = 1
    U1[:, 0] = A_x[0]
    U1_bar[:, 0] = A_x[0]

    # initial conditions
    V[0, :] = 0.59 / (dens[0, :] * A_x[:])
    U1[0, 1:] = dens[0, 1:] * A_x[1:]
    U2[0, :] = 0.59
    U3[0, :] = dens[0, :] * A_x[:] * ((T[0, :] / (gamma - 1)) + ((gamma / 2) * V[0, :] ** 2))

    # analytical solution
    M_an = np.array([0.09782, 0.10656, 0.11646, 0.12769, 0.14049, 0.15513,
                     0.17195, 0.19134, 0.21377, 0.23981, 0.27012, 0.30549,
                     0.3468, 0.39505, 0.45128, 0.51657, 0.59185, 0.67785,
                     0.77485, 0.88254, 1, 1.12566, 1.25755, 1.39348,
                     1.53136, 1.66934, 1.80597, 1.94015, 0.56498, 0.480978,
                     0.41626, 0.36423, 0.321384, 0.285538, 0.255199, 0.229285,
                     0.20698, 0.187656, 0.170816, 0.156061, 0.14307])

    p_an = (1 + (gamma - 1) / 2 * M_an[:] ** 2) ** (-(gamma) / (gamma - 1))
    p_an[x_40 >= 2.1] = p_an[x_40 >= 2.1] * 0.6882
    dens_an = (1 + (gamma - 1) / 2 * M_an[:] ** 2) ** (-1 / (gamma - 1))
    dens_an[x_40 >= 2.1] = dens_an[x_40 >= 2.1] * 0.6882
    T_an = (1 + (gamma - 1) / 2 * M_an[:] ** 2) ** (-1)
    V_an = M_an * np.sqrt(T_an)
    m_flow_an = dens_an[:] * A_x[:] * V_an[:]

    # Visualization
    data = np.array([x_total[:], A_x[:], V[0, :], dens[0, :], T[0, :], U1[0, :], U2[0, :], U3[0, :]])
    data = data.T
    df = pd.DataFrame(data)

    control = True
    k = 0
    while k <= interations -2 and control:
        dt = min(C * (dx / ((T[k, :]) ** 0.5 + V[k, :])))

        F1[k, :] = U2[k, :]

        for i in range(len(x_total)):
            F2[k, i] = ((U2[k, i] ** 2) / U1[k, i]) + (gamma - 1) / gamma * (
                    U3[k, i] - (gamma / 2) * (U2[k, i] ** 2 / U1[k, i]))
            F3[k, i] = ((gamma * U2[k, i] * U3[k, i] / U1[k, i]) - (
                    gamma * (gamma - 1) / 2) * (U2[k, i] ** 3 / U1[k, i] ** 2))

        for i in range(1, len(x_total) - 1):
            J2[k, i] = ((1/gamma) * dens[k, i] * T[k, i] * (A_x[i + 1] - A_x[i])/dx)


            dU1_dt[k, i] = -(F1[k, i + 1] - F1[k, i]) / dx
            dU2_dt[k, i] = -(F2[k, i + 1] - F2[k, i]) / dx + J2[k, i]
            dU3_dt[k, i] = -(F3[k, i + 1] - F3[k, i]) / dx

            S1[k, i] = Cx * (np.abs(dens[k, i + 1] * T[k, i + 1] - 2 * dens[k, i] * T[k, i] +
                                dens[k, i - 1] * T[k, i - 1]) /
                                (dens[k, i + 1] * T[k, i + 1] + 2 * dens[k, i] * T[k, i] +
                                  dens[k, i - 1] * T[k, i - 1])) * (U1[k, i + 1] -
                                                                   2 * U1[k, i] + U1[k, i - 1])
            S2[k, i] = Cx * (np.abs(dens[k, i + 1] * T[k, i + 1] - 2 * dens[k, i] * T[k, i] +
                                     dens[k, i - 1] * T[k, i - 1]) /
                                 (dens[k, i + 1] * T[k, i + 1] + 2 * dens[k, i] * T[k, i] +
                                  dens[k, i - 1] * T[k, i - 1])) * (U2[k, i + 1] -
                                                                   2 * U2[k, i] + U2[k, i - 1])
            S3[k, i] = Cx * (np.abs(dens[k, i + 1] * T[k, i + 1] - 2 * dens[k, i] * T[k, i] +
                                     dens[k, i - 1] * T[k, i - 1]) /
                                 (dens[k, i + 1] * T[k, i + 1] + 2 * dens[k, i] * T[k, i] +
                                  dens[k, i - 1] * T[k, i - 1])) * (U3[k, i + 1] -
                                                                   2 * U3[k, i] + U3[k, i - 1])

            U1_bar[k + 1, i] = U1[k, i] + dU1_dt[k, i] * dt + S1[k, i]
            U2_bar[k + 1, i] = U2[k, i] + dU2_dt[k, i] * dt + S2[k, i]
            U3_bar[k + 1, i] = U3[k, i] + dU3_dt[k, i] * dt + S3[k, i]

            # Boundaries
        U1_bar[k + 1, -1] = 2 * U1_bar[k + 1, -2] - U1_bar[k + 1, -3]
        U2_bar[k + 1, -1] = 2 * U2_bar[k + 1, -2] - U2_bar[k + 1, -3]
        U3_bar[k + 1, -1] = 2 * U3_bar[k + 1, -2] - U3_bar[k + 1, -3]

        U2_bar[k + 1, 0] = 2 * U2_bar[k + 1, 1] - U2_bar[k + 1, 2]
        V_bar[k + 1, :] = U2_bar[k + 1, :] / U1_bar[k + 1, :]
        U3_bar[k + 1, 0] = U1_bar[k + 1, 0] * (T[k + 1, 0] / (gamma - 1) + (gamma / 2) * V_bar[k + 1, 0] ** 2)


        for i in range(1, len(x_total)):
            dens_bar[k + 1, i] = U1_bar[k + 1, i] / A_x[i]

            T_bar[k + 1, i] = (gamma - 1) * ((U3_bar[k + 1, i] / U1_bar[k + 1, i]) -
                                             (gamma / 2) * (U2_bar[k + 1, i] / U1_bar[k + 1, i]) ** 2)

        F1_bar[k + 1, :] = U2_bar[k + 1, :]


        for i in range(len(x_total)):
            F2_bar[k + 1, i] = (U2_bar[k + 1, i] ** 2 / U1_bar[k + 1, i]) + (gamma - 1) / gamma * (
                    U3_bar[k + 1, i] - (gamma / 2) * (U2_bar[k + 1, i] ** 2 / U1_bar[k + 1, i]))

            F3_bar[k + 1, i] = (gamma * U2_bar[k + 1, i] * U3_bar[k + 1, i] / U1_bar[k + 1, i]) - (
                    gamma * (gamma - 1) / 2 * (U2_bar[k + 1, i] ** 3 / U1_bar[k + 1, i] ** 2))


        for j in range(1, len(x_total) - 1):
            dU1_bar_dt[k + 1, j] = - (F1_bar[k + 1, j] - F1_bar[k + 1, j - 1]) / dx
            dU2_bar_dt[k + 1, j] = - (F2_bar[k + 1, j] - F2_bar[k + 1, j - 1]) / dx + \
                                   (1 / gamma) * dens_bar[k + 1, j] * T_bar[k + 1, j] * (A_x[j] - A_x[j - 1]) / dx
            dU3_bar_dt[k + 1, j] = - (F3_bar[k + 1, j] - F3_bar[k + 1, j - 1]) / dx


        for i in range(1, len(x_total) - 1):
            dU1_avg[k + 1, i] = 0.5 * (dU1_dt[k, i] + dU1_bar_dt[k + 1, i])
            dU2_avg[k + 1, i] = 0.5 * (dU2_dt[k, i] + dU2_bar_dt[k + 1, i])
            dU3_avg[k + 1, i] = 0.5 * (dU3_dt[k, i] + dU3_bar_dt[k + 1, i])

            S1_bar[k + 1, i] = Cx * (np.abs(dens_bar[k + 1, i + 1] * T_bar[k + 1, i + 1] - 2 * dens_bar[k + 1, i] * T_bar[k + 1, i] +
                                dens_bar[k + 1, i - 1] * T_bar[k + 1, i - 1]) /
                         (dens_bar[k + 1, i + 1] * T_bar[k + 1, i + 1] + 2 * dens_bar[k + 1, i] * T_bar[k + 1, i] +
                          dens_bar[k + 1, i - 1] * T_bar[k + 1, i - 1])) * (U1_bar[k + 1, i + 1] -
                                                            2 * U1_bar[k + 1, i] + U1_bar[k + 1, i - 1])
            S2_bar[k + 1, i] = Cx * (np.abs(dens_bar[k + 1, i + 1] * T_bar[k + 1, i + 1] - 2 * dens_bar[k + 1, i] * T_bar[k + 1, i] +
                                dens_bar[k + 1, i - 1] * T_bar[k + 1, i - 1]) /
                         (dens_bar[k + 1, i + 1] * T_bar[k + 1, i + 1] + 2 * dens_bar[k + 1, i] * T_bar[k + 1, i] +
                          dens_bar[k + 1, i - 1] * T_bar[k + 1, i - 1])) * (U2_bar[k + 1, i + 1] -
                                                            2 * U2_bar[k + 1, i] + U2_bar[k + 1, i - 1])
            S3_bar[k + 1, i] = Cx * (np.abs(dens_bar[k + 1, i + 1] * T_bar[k + 1, i + 1] - 2 * dens_bar[k + 1, i] * T_bar[k + 1, i] +
                                dens_bar[k + 1, i - 1] * T_bar[k + 1, i - 1]) /
                         (dens_bar[k + 1, i + 1] * T_bar[k + 1, i + 1] + 2 * dens_bar[k + 1, i] * T_bar[k + 1, i] +
                          dens_bar[k + 1, i - 1] * T_bar[k + 1, i - 1])) * (U3_bar[k + 1, i + 1] -
                                                            2 * U3_bar[k + 1, i] + U3_bar[k + 1, i - 1])

        for i in range(1, len(x_total) - 1):
            U1[k + 1, i] = U1[k, i] + dU1_avg[k + 1, i] * dt + S1_bar[k + 1, i]
            U2[k + 1, i] = U2[k, i] + dU2_avg[k + 1, i] * dt + S2_bar[k + 1, i]
            U3[k + 1, i] = U3[k, i] + dU3_avg[k + 1, i] * dt + S3_bar[k + 1, i]

        U1[k + 1, -1] = 2 * U1[k + 1, -2] - U1[k + 1, -3]
        U2[k + 1, -1] = 2 * U2[k + 1, -2] - U2[k + 1, -3]
        U3[k + 1, -1] = dp * A_x[-1] / (gamma - 1) + dens[k-1, -1] * A_x[-1] * gamma * V[k-1, -1]**2 / 2

        U2[k + 1, 0] = 2 * U2[k + 1, 1] - U2[k + 1, 2]
        V[k + 1, :] = U2[k + 1, :] / U1[k + 1, :]

        U3[k + 1, 0] = U1[k + 1, 0] * (T[k + 1, 0] / (gamma - 1) + (gamma / 2) * V[k + 1, 0] ** 2)

        dens[k + 1, 1:] = U1[k + 1, 1:] / A_x[1:]
        T[k + 1, 1:] = (gamma - 1) * ((U3[k + 1, 1:] / U1[k + 1, 1:]) - (gamma / 2) * V[k + 1, 1:]**2)

        time_table[0, k] = k
        time_table[1, k] = dt
        if all(np.abs(U1[k + 1, :] - U1[k, :]) <= sfalma) and all(
                np.abs(U2[k + 1, :] - U2[k, :]) <= sfalma) and all(
            np.abs(U3[k + 1, :] - U3[k, :]) <= sfalma):
            control = False
        k += 1

        #print(k, 'k')
        data = np.array([x_total[:], A_x[:], V[k, :], dens[k, :], T[k, :], U1[k, :], U2[k, :], U3[k, :]])
        data = data.T
        df = pd.DataFrame(data)
        #print(df)
    a = np.sqrt(T)
    M = V[0:k, :] / a[0:k, :]
    p = dens[0:k, :] * T[0:k, :]
    A = np.zeros([k, len(x_total)])
    A[:, :] = A_x[:]
    k = k - 1

######################################### PLOTS ###########################################

    plot_titles = ['Flow Field quantities vs number of iterations', 'Residuals vs number of iterations',
                       'Mass flow different timesteps vs nozzle non-d. distance', 'Nond. pressure vs nond. distance',
                   'Mach number vs non-d distance']

    x_labels = ['Number of time steps', 'Number of timesteps', 'Non-dimensional distance [x/L]',
                    'Non-dimensional distance [x/L]', 'Non-dimensional distance [x/L]']

    y_labels = ['Non-Dimensional quantities', 'Residuals', 'Non-dimensional mass flow', 'Non-dimensional pressure',
                'Mach number']

    colorpallete = ['r', 'b', 'y', 'm', 'g']

    line_labels = ['nond. density', 'nond. temperature', 'nond. pressure', 'Mach number',
                       '(dρ/dt)avg', '(dT/dt)avg', '(dV/dt)avg', '0.05t_total', '0.1t_total', '0.2t_total', '0.7t_total', 'convergence',
                       '0.1t_total', '0.5t_total', '0.7t_total', 'convergence', 'Mach number']
    table_titles = ['Calculated converged numerical values ', 'Calculated and analytical values',
                    'Calculated first step analytical values']

    values_titles = ['x/L', 'A/A*', 'V/a', 'ρ/ρο', 'T/Tο', 'p/po', 'M', 'm_dot']
    values2_titles = ['x/L', 'A/A*','ρ/ρο', 'ρ/ρο(analytic)', '%difference', 'M', 'M(analytical)', '%difference']

    tables = [x_total, A_x, V[k, :], dens[k, :], T[k, :], p[k, :], M[-1, :], dens[k, :] * V[k, :] * A_x]
    if n == 41:
        tables2 = [x_total, A_x, dens[k, :], dens_an, np.abs(dens[k, :] - dens_an) / dens_an, M[-1, :], M_an,
                   np.abs(M[-1, :] - M_an) / M_an]
    tables3 = [x_total, A_x, V[1, :], dens[1, :], T[1, :], p[1, :], M[1, :], dens[1, :] * V[1, :] * A_x]

    plots = [dens[0:k, int(n / 2) + 1], T[0:k, int(n / 2) + 1], p[0:k, int(n / 2) + 1], M[0:k, int(n / 2) + 1],
                 np.abs(dU1_avg[0:k, int(n / 2)]), np.abs(dU2_avg[0:k, int(n / 2)]), np.abs(dU3_avg[0:k, int(n / 2)])
            , dens[int(.05 * k), :] * V[int(.05 * k), :] * A[int(.05 * k), :],
                 dens[int(.1 * k), :] * V[int(.1 * k), :] * A[int(.1 * k), :],
                 dens[int(.2 * k), :] * V[int(.2 * k), :] * A[int(.2 * k), :],
                 dens[int(.7 * k), :] * V[int(.7 * k), :] * A[int(.7 * k), :], dens[k, :] * V[k, :] * A[-1, :],
                 p[int(.1 * k), :], p[int(.5 * k), :], p[int(.7 * k), :], p[int(k), :], M[k, :]]
    print(plots[-1])
    print(x_total)

    obj = gui_post_processor.PostProcessor(plot_titles, table_titles)
    obj.run()
    selection1, selection2 = obj.get_vals()

    for i in range(0, len(selection1)):
        plt.figure(i + 1)
        if selection1[i] == 0:
            gui_post_processor.plotter(time_table[0, 0:k], plots[0:4], plot_titles[0], x_labels[0], y_labels[0],
                                           colorpallete, line_labels[0:4])
        elif selection1[i] == 1:
            gui_post_processor.plotter(time_table[0, 0:k], plots[4:7], plot_titles[1], x_labels[1], y_labels[1],
                                           colorpallete, line_labels[4:7])
        elif selection1[i] == 2:
            gui_post_processor.plotter(x_total, plots[7:12], plot_titles[2], x_labels[2], y_labels[2],
                                           colorpallete, line_labels[7:12])
            plt.scatter(0, m_flow_an[2], color='red', s=20, label='Ακριβής τιμή')
            plt.legend()
        elif selection1[i] == 3:
            gui_post_processor.plotter(x_total, plots[12:16], plot_titles[3], x_labels[3], y_labels[3],
                                           colorpallete, line_labels[12:16])
            plt.scatter(x_40, p_an, color='black', s=10)
        else:
            gui_post_processor.plotter(x_total, plots[-1], plot_titles[4], x_labels[4], y_labels[4],
                                       colorpallete, line_labels[-1])
            plt.scatter(x_40, M_an, color='black', s=10)

    for i in range(0, len(selection2)):
        if selection2[i] == 0:
            gui_post_processor.tabloter(values_titles, tables, selection2[i], n)
        elif selection2[i] == 1:
            gui_post_processor.tabloter(values2_titles, tables2, selection2[i], n)
        else:
            print((tables3))
            print((values_titles))
            gui_post_processor.tabloter(values_titles, tables3, selection2[i], n)
    print(k)
    total_time = sum(time_table[1, :])
    print(total_time, 'Συνολικός χρόνος')
    print(dens_an[-1], V_an[-1],  T_an[-1], p_an[-1], M_an[-1], m_flow_an[-1])
    print(dens[k, -1], V[k, -1], T[k, -1], p[k, -1],
          M[k, -1], dens[k, -1] * V[k,-1] * A_x[-1])
    print('total time:', np.sum(time_table[1, :]))
    plt.show()


#conservative(41, 0.5, 10000, 0.6784, 0.3, 10**(-4))