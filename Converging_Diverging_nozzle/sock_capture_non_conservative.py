import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math
import time
import pandas as pd
import tkinter.messagebox
plt.style.use('ggplot')
gamma = 1.4
L = 3

def noncons(C, dp, n, sfalma, interations, omega, Cx):
    dx = L / n
    x_total = np.arange(0, L + dx, dx)
    A_x = np.zeros(len(x_total))
    dt = np.zeros([interations, len(x_total)])

    dens = np.zeros([interations, len(x_total)])
    V = np.zeros([interations, len(x_total)])
    T = np.zeros([interations, len(x_total)])

    S1 = np.zeros([interations, len(x_total)])
    S2 = np.zeros([interations, len(x_total)])
    S3 = np.zeros([interations, len(x_total)])

    S1_bar = np.zeros([interations, len(x_total)])
    S2_bar = np.zeros([interations, len(x_total)])
    S3_bar = np.zeros([interations, len(x_total)])

    continuity = np.zeros([interations, len(x_total)])
    momentum = np.zeros([interations, len(x_total)])
    energy = np.zeros([interations, len(x_total)])

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

    # initial conditions
    for i in range(len(x_total)):
        A_x[i] = 1 + 2.2 * (x_total[i] - 1.5) ** 2
        if x_total[i] <= 0.5:
            dens[0, i] = 1
            T[0, i] = 1
            V[0, i] = (0.1 + 1.09 * x_total[i]) * T[0, i] ** 0.5
        elif x_total[i] <= 1.5:
            dens[0, i] = 1 - 0.366 * (x_total[i] - 0.5)
            T[0, i] = 1 - 0.167 * (x_total[i] - 0.5)
            V[0, i] = (0.1 + 1.09 * x_total[i]) * T[0, i] ** 0.5
        elif x_total[i] <= 2.1:
            dens[0, i] = 0.634 - 0.702 * (x_total[i] - 1.5)
            T[0, i] = 0.833 - 0.4908 * (x_total[i] - 1.5)
            V[0, i] = (0.1 + 1.09 * x_total[i]) * T[0, i] ** 0.5
        else:
            dens[0, i] = 0.5892 + 0.10228 * (x_total[i] - 2.1)
            T[0, i] = 0.93968 + 0.0622 * (x_total[i] - 2.1)
            V[0, i] = 0.05 + (0.11 * x_total[i])


    data = np.array([x_total[:], A_x[:], V[0, :], dens[0, :], T[0, :]])
    data = data.T
    df = pd.DataFrame(data)
    print(df)
    T[:, 0] = 1
    T_bar[:, 0] = 1
    dens[:, 0] = 1
    dens_bar[:, 0] = 1
    k = 0


    continuity[0, :] = dens[0, :] * A_x[:]
    momentum[0, :] = dens[0, :] * V[0, :] * A_x[:]
    energy[0, :] = dens[0, :] * (T[0, :] / (gamma - 1) + gamma / 2 * V[0, :]**2) * A_x[:]



    ###dt[0, :] = C * dx / (a[0, :] + V[0, :])
    dt = min(C * (dx / (T[k, :] ** 0.5 + V[k, :])))
    while k <= interations - 2:

        for i in range(1, len(x_total) - 1):

            ddens_dt[k, i] = -dens[k, i] * ((V[k, i + 1] - V[k, i]) / dx) - dens[k, i] * V[k, i] * (
                    (np.log(A_x[i + 1]) - np.log(A_x[i])) / dx) - V[k, i] * ((dens[k, i + 1] - dens[k, i]) / dx)

            dV_dt[k, i] = -V[k, i] * (V[k, i + 1] - V[k, i]) / dx - (1 / gamma) * (
                    (T[k, i + 1] - T[k, i]) / dx + (T[k, i] / dens[k, i]) * (dens[k, i + 1] - dens[k, i]) / dx)

            dT_dt[k, i] = -V[k, i] * ((T[k, i + 1] - T[k, i]) / dx) - (gamma - 1) * T[k, i] * (
                    ((V[k, i + 1] - V[k, i]) / dx) + V[k, i] * ((np.log(A_x[i + 1]) - np.log(A_x[i])) / dx))


            S1[k, i] = Cx * (np.abs(dens[k, i + 1] * T[k, i + 1] - 2 * dens[k, i] * T[k, i] +
                                    dens[k, i - 1] * T[k, i - 1]) /
                             (dens[k, i + 1] * T[k, i + 1] + 2 * dens[k, i] * T[k, i] +
                              dens[k, i - 1] * T[k, i - 1])) * (dens[k, i + 1] -
                                                                2 * dens[k, i] + dens[k, i - 1])
            S2[k, i] = Cx * (np.abs(dens[k, i + 1] * T[k, i + 1] - 2 * dens[k, i] * T[k, i] +
                                    dens[k, i - 1] * T[k, i - 1]) /
                             (dens[k, i + 1] * T[k, i + 1] + 2 * dens[k, i] * T[k, i] +
                              dens[k, i - 1] * T[k, i - 1])) * (V[k, i + 1] -
                                                                2 * V[k, i] + V[k, i - 1])
            S3[k, i] = Cx * (np.abs(dens[k, i + 1] * T[k, i + 1] - 2 * dens[k, i] * T[k, i] +
                                    dens[k, i - 1] * T[k, i - 1]) /
                             (dens[k, i + 1] * T[k, i + 1] + 2 * dens[k, i] * T[k, i] +
                              dens[k, i - 1] * T[k, i - 1])) * (T[k, i + 1] -
                                                                2 * T[k, i] + T[k, i - 1])


            dens_bar[k + 1, i] = dens[k, i] + ddens_dt[k, i] * dt + S1[k, i]
            V_bar[k + 1, i] = V[k, i] + dV_dt[k, i] * dt + S2[k, i]
            T_bar[k + 1, i] = T[k, i] + dT_dt[k, i] * dt + S3[k, i]

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

            S1_bar[k + 1, i] = Cx * (
                        np.abs(dens_bar[k + 1, i + 1] * T_bar[k + 1, i + 1] - 2 * dens_bar[k + 1, i] * T_bar[k + 1, i] +
                               dens_bar[k + 1, i - 1] * T_bar[k + 1, i - 1]) /
                        (dens_bar[k + 1, i + 1] * T_bar[k + 1, i + 1] + 2 * dens_bar[k + 1, i] * T_bar[k + 1, i] +
                         dens_bar[k + 1, i - 1] * T_bar[k + 1, i - 1])) * (dens_bar[k + 1, i + 1] -
                                                                           2 * dens_bar[k + 1, i] + dens_bar[k + 1, i - 1])
            S2_bar[k + 1, i] = Cx * (
                        np.abs(dens_bar[k + 1, i + 1] * T_bar[k + 1, i + 1] - 2 * dens_bar[k + 1, i] * T_bar[k + 1, i] +
                               dens_bar[k + 1, i - 1] * T_bar[k + 1, i - 1]) /
                        (dens_bar[k + 1, i + 1] * T_bar[k + 1, i + 1] + 2 * dens_bar[k + 1, i] * T_bar[k + 1, i] +
                         dens_bar[k + 1, i - 1] * T_bar[k + 1, i - 1])) * (V_bar[k + 1, i + 1] -
                                                                           2 * V_bar[k + 1, i] + V_bar[k + 1, i - 1])
            S3_bar[k + 1, i] = Cx * (
                        np.abs(dens_bar[k + 1, i + 1] * T_bar[k + 1, i + 1] - 2 * dens_bar[k + 1, i] * T_bar[k + 1, i] +
                               dens_bar[k + 1, i - 1] * T_bar[k + 1, i - 1]) /
                        (dens_bar[k + 1, i + 1] * T_bar[k + 1, i + 1] + 2 * dens_bar[k + 1, i] * T_bar[k + 1, i] +
                         dens_bar[k + 1, i - 1] * T_bar[k + 1, i - 1])) * (T_bar[k + 1, i + 1] -
                                                                           2 * T_bar[k + 1, i] + T_bar[k + 1, i - 1])

            dens[k + 1, i] = (1 - omega) * dens[k, i] + omega * ddens_avg[k + 1, i] * dt + S1_bar[k + 1, i]
            V[k + 1, i] = (1 - omega) * V[k, i] + omega *  dV_avg[k + 1, i] * dt + S2_bar[k + 1, i]
            T[k + 1, i] = (1 - omega) * T[k, i] + omega *  dT_avg[k + 1, i] * dt + S3_bar[k + 1, i]
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
        print(df)
        print(k)
        k += 1

# noncons(0.5, 0.6784, 61, 0.00001, 5000, 0.6, 0.1)
##################################################################################
