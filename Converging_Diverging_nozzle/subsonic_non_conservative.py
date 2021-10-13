import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math
import time
import pandas as pd
import tkinter.messagebox
import gui_post_processor
plt.style.use('ggplot')
gamma = 1.4
L = 3

def subsonic(C, dp, n, sfalma, interations):
    x_40 = np.linspace(0, L, 41)
    A_40 = np.zeros(len(x_40))
    for i in range(len(x_40)):
        if x_40[i] <= 1.5:
            A_40[i] = 1 + 2.2 * (x_40[i] - 1.5) ** 2
        else:
            A_40[i] = 1 + 0.2223 * (x_40[i] - 1.5) ** 2

    dx = L/n
    x_total = np.linspace(0, L, n)
    A_x = np.zeros(len(x_total))
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

    for i in range(len(x_total)):
        if x_total[i] <= 1.5:
            A_x[i] = 1 + 2.2 * (x_total[i] - 1.5) ** 2
            dens[0, i] = 1 - (0.023 * x_total[i] )
            T[0, i] = 1 - (0.009333 * x_total[i] )
            V[0, i] = 0.05 + (0.11 * x_total[i] )
        else:
            A_x[i] = 1 + 0.2223 * (x_total[i] - 1.5) ** 2
            dens[0, i] = 1 - (0.023 * x_total[i])
            T[0, i] = 1 - (0.009333 * x_total[i])
            V[0, i] = 0.05 + (0.11 * x_total[i])

    data = np.array([x_total[:], A_x[:], V[0, :], dens[0, :], T[0, :]])
    data = data.T
    df = pd.DataFrame(data)

    T[:, 0] = 1
    T_bar[:, 0] = 1
    dens[:, 0] = 1
    dens_bar[:, 0] = 1
    k = 0

    # analytical solution
    M_an = np.array([0.076955, 0.083804, 0.091538, 0.100302, 0.110268, 0.121638,
                     0.13465, 0.149582, 0.166757, 0.186537, 0.209324, 0.235529,
                     0.265534, 0.299606, 0.337743, 0.379433, 0.423279, 0.466514,
                     0.504598, 0.531457, 0.54125, 0.540239, 0.537236, 0.532325,
                     0.52564, 0.517353, 0.50766, 0.496771, 0.484898, 0.472246,
                     0.459008, 0.445358, 0.431454, 0.417432, 0.403409, 0.389482,
                     0.375734, 0.36223, 0.349025, 0.336157, 0.323658])

    p_an = (1 + (gamma - 1) / 2 * M_an[:] ** 2) ** (-(gamma) / (gamma - 1))
    dens_an = (1 + (gamma - 1) / 2 * M_an[:] ** 2) ** (-1 / (gamma - 1))
    T_an = (1 + (gamma - 1) / 2 * M_an[:] ** 2) ** (-1)
    m_flow_an = dens_an[:] * A_40[:] * M_an[:] * np.sqrt(T_an[:])



    ###dt[0, :] = C * dx / (a[0, :] + V[0, :])
    dt = min(C * (dx / (T[k, :] ** 0.5 + V[k, :])))
    control = True
    while k <= interations -2 and control:

        for i in range(1, len(x_total) - 1):  # oi paragogoi gemizoun apo 2 mexri n-1

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
            dens_bar[k + 1, 0] = dens[k, 0]
            V_bar[k + 1, 0] = V[k, 0]
            T_bar[k + 1, 0] = T[k, 0]
            V_bar[k, 0] = V[k, 0]


            ddens_bar_dt[k + 1, i] = -dens_bar[k + 1, i] * ((V_bar[k + 1, i] - V_bar[k + 1, i - 1]) / dx) - \
                                     dens_bar[k + 1, i] * V_bar[k + 1, i] * (np.log(A_x[i]) - np.log(A_x[i - 1])) / dx - \
                                     V_bar[k + 1, i] * ((dens_bar[k + 1, i] - dens_bar[k + 1, i - 1]) / dx)

            dV_bar_dt[k + 1, i] = -V_bar[k + 1, i] * ((V_bar[k + 1, i] - V_bar[k + 1, i - 1]) / dx) - \
                                  (1 / gamma) * ((T_bar[k + 1, i] - T_bar[k + 1, i - 1]) / dx +
                                                 T_bar[k + 1, i] / dens_bar[k + 1, i] * (
                                                             dens_bar[k + 1, i] - dens_bar[k + 1, i - 1]) / dx)

            dT_bar_dt[k + 1, i] = -V_bar[k + 1, i] * ((T_bar[k + 1, i] - T_bar[k + 1, i - 1]) / dx) - \
                                  (gamma - 1) * T_bar[k + 1, i] * (
                                              ((V_bar[k + 1, i] - V_bar[k + 1, i - 1]) / dx) + V_bar[k + 1, i] *
                                              ((np.log(A_x[i]) - np.log(A_x[i - 1])) / dx))

            ddens_avg[k + 1, i] = 0.5 * (ddens_dt[k, i] + ddens_bar_dt[k + 1, i])
            dV_avg[k + 1, i] = 0.5 * (dV_dt[k, i] + dV_bar_dt[k + 1, i])
            dT_avg[k + 1, i] = 0.5 * (dT_dt[k, i] + dT_bar_dt[k + 1, i])

            dens[k + 1, i] = dens[k, i] + ddens_avg[k + 1, i] * dt
            V[k + 1, i] = V[k, i] + dV_avg[k + 1, i] * dt
            T[k + 1, i] = T[k, i] + dT_avg[k + 1, i] * dt
            a[k + 1, i] = math.sqrt(abs(T[k + 1, i]))
            #### dt[k + 1, i] = C * (dx / (a[k + 1, i] + V[k + 1, i]))

        # boundary conditions
        V[k + 1, 0] = 2 * V[k + 1, 1] - V[k + 1, 2]
        #### dt [k + 1, 0] = C * (dx / (T[k + 1, 0] ** 0.5 + V[k + 1, 0]))
        a[k + 1, 0] = math.sqrt(T[k + 1, 0])

        #### dt [k + 1, 0] = C * (dx / (T[k + 1, 0] ** 0.5 + V[k + 1, 0]))

        dens[k + 1, -1] = 2 * dens[k + 1,  - 2] - dens[k + 1,  - 3]
        T[k + 1, -1] = dp / dens[k + 1, -1]
        V[k + 1, -1] = 2 * V[k + 1, - 2] - V[k + 1, - 3]

        #### dt [k + 1, n] = C * (dx / (a[k + 1, n] + V[k + 1, n]))

        dt = min(C * (dx / (T[k + 1, :] ** 0.5 + V[k + 1, :])))
        data = np.array([x_total[:], A_x[:], V[k + 1, :], dens[k + 1, :], T[k + 1, :]])
        data = data.T
        df = pd.DataFrame(data)
        # print(df)
        print(k)
        time_table[0, k] = k
        time_table[1, k] = dt
        if all(np.abs(dV_avg[k + 1, :] - dV_avg[k, :]) <= sfalma) and all(
                np.abs(dT_avg[k + 1, :] - dT_avg[k, :]) <= sfalma) and all(
                np.abs(ddens_avg[k + 1, :] - ddens_avg[k, :]) <= sfalma):
            control = False
        k += 1
    a = np.sqrt(T)
    M = V[0:k, :] / a[0:k, :]
    p = dens[0:k, :] * T[0:k, :]
    A = np.zeros([k, len(x_total)])
    A[:, :] = A_x[:]
    k = k - 1


######################################### PLOTS ###########################################

    plot_titles = ['Flow Field quantities vs number of iterations', 'Residuals vs number of iterations',
                   'Mass flow different timesteps vs nozzle non-d. distance', 'Nond. pressure vs nond. distance']

    x_labels = ['Number of time steps', 'Number of timesteps', 'Non-dimensional distance [x/L]',
                'Non-dimensional distance [x/L]']

    y_labels = ['Non-Dimensional quantities', 'Residuals', 'Non-dimensional mass flow', 'Non-dimensional pressure']

    colorpallete = ['r', 'b', 'y', 'm', 'g']


    line_labels = ['nond. density', 'nond. temperature', 'nond. pressure', 'Mach number',
                   '(dρ/dt)avg', '(dT/dt)avg', '(dV/dt)avg', '0.05t_total', '0.1t_total', '0.2t_total', '0.7t_total', 'convergence',
                   '0.1t_total', '0.5t_total', '0.7t_total', 'convergence']
    table_titles = ['Calculated converged numerical values ', 'Calculated and analytical values']

    values_titles = ['x/L', 'A/A*', 'V/a', 'ρ/ρο', 'T/Tο', 'p/po', 'M', 'm_dot']
    values2_titles = ['x/L', 'A/A*', 'ρ/ρο', 'ρ/ρο(analytic)', '%difference', 'M', 'M(analytical)', '%difference']

    tables = [x_total, A_x, V[k, :], dens[k, :], T[k, :], p[k, :], M[-1, :],  dens[k, :] * V[k, :] * A_x]
    tables2 = [x_total, A_x, dens[k, :], dens_an, np.abs(dens[k, :] - dens_an) / dens_an, M[-1, :], M_an,
               np.abs(M[-1, :] - M_an) / M_an]

    plots = [dens[0:k, int(n / 2) + 1], T[0:k, int(n / 2) + 1], p[0:k, int(n / 2) + 1], M[0:k, int(n / 2) + 1],
             np.abs(ddens_avg[0:k, int(n / 2)]), np.abs(dT_avg[0:k, int(n / 2)]), np.abs(dV_avg[0:k, int(n / 2)])
        , dens[int(.05 * k), :] * V[int(.05 * k), :] * A[int(.05 * k), :],
             dens[int(.1 * k), :] * V[int(.1 * k), :] * A[int(.1 * k), :],
             dens[int(.2 * k), :] * V[int(.2 * k), :] * A[int(.2 * k), :],
             dens[int(.7 * k), :] * V[int(.7 * k), :] * A[int(.7 * k), :], dens[k, :] * V[k, :] * A[-1, :],
             p[int(.3 * k), :], p[int(.5 * k), :], p[int(.7 * k), :],  p[int(k), :]]

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
            # plt.ylim(0.000001, 0.01)

        elif selection1[i] == 2:
            gui_post_processor.plotter(x_total, plots[7:12], plot_titles[2], x_labels[2], y_labels[2],
                                       colorpallete, line_labels[7:12])
            plt.scatter(0, m_flow_an[2], color='red', s=20, label='Ακριβής τιμή')
            plt.legend()
        else:
            gui_post_processor.plotter(x_total, plots[12:], plot_titles[3], x_labels[3], y_labels[3],
                                       colorpallete, line_labels[7:12])
            plt.scatter(x_40, p_an, color='black', s=20, label='ακριβής τιμή')
            plt.legend()

    for i in range(0, len(selection2)):
        if selection2[i] == 0:
            gui_post_processor.tabloter(values_titles, tables, selection2[i], n)
        else:
            gui_post_processor.tabloter(values2_titles, tables2, selection2[i], n)
    print('total time:', np.sum(time_table[1, :]))
    plt.show()


#subsonic(0.5, 0.88, 41,  10**(-4), 10000)