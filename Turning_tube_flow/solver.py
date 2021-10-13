import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
import time

from metrics import *

dens = 1.23
visc = 18.37e-6
omega = 0.5

Nx = 200
Ny = 60

u_old = np.zeros([len(xu[:, 0]), len(yu[0, :])])
u_new = np.zeros([len(xu[:, 0]), len(yu[0, :])])

v_old = np.zeros([len(xv[:, 0]), len(yv[0, :])])
v_new = np.zeros([len(xv[:, 0]), len(yv[0, :])])

v_u_p = np.zeros([len(xu[:, 0]), len(yu[0, :])])
u_v_p = np.zeros([len(xv[:, 0]), len(yv[0, :])])

p_old = np.zeros([len(xp[:, 0]), len(yp[0, :])])
p_new = np.zeros([len(xp[:, 0]), len(yp[0, :])])


# arxikes sun8ikes
for i in range(0, len(xu[:, 0])):
    if yy[i, -1] - yy[0, -1] < 0.00000001:
        u_old[i, 1:-1] = 0.1
        u_new[i, 1:-1] = 0.1
        v_old[i, 1:-1] = 0
        v_new[i, 1:-1] = 0
    elif np.abs(xx[i, 0] - xx[i-1, 0]) < 0.0000001:
        u_old[i, 1:-1] = 0
        u_new[i, 1:-1] = 0
        v_old[i, 1:-1] = 0.1
        v_new[i, 1:-1] = 0.1
    else:
        u_old[i, 1:-1] = np.sqrt(2)/2 * 0.1
        u_new[i, 1:-1] = np.sqrt(2)/2 * 0.1
        v_old[i, 1:-1] = np.sqrt(2)/2 * 0.1
        v_new[i, 1:-1] = np.sqrt(2)/2 * 0.1



q1 = b21hta + b22hta
q2 = b11hta + b12hta
q3 = b11ksi + b12ksi

q1w = b21htaw + b22htaw
q2w = b11htaw + b12htaw
q3w = b11ksiw + b12ksiw

q1s = b21htas + b22htas
q2s = b11htas + b12htas
q3s = b11ksis + b12ksis

De_u = []
Fe_u = []
ae_u = []
Dw_u = []
Fw_u = []
aw_u = []
Dn_u = []
Fn_u = []
an_u = []
Ds_u = []
Fs_u = []
as_u = []
dp_u_x = []
dp_u_y = []

Sp_u = []

ap_u_big = np.zeros([0, Ny-2])
ap_v_big = np.zeros([0, Ny-2])

for i in range(1, len(xu[:, 0]) - 2):
    ap_u = []
    ae_u = []
    aw_u = []
    an_u = []
    as_u = []
    dp_u = []
    Sp_u = []
    for j in range(1, len(yu[0, :]) - 1):
        # u
        v_u_p[i, j] = (v_old[i, j+1] + v_new[i-1, j+1] + v_new[i-1, j] + v_old[i, j]) / 4

        # west
        Dw_u = visc / Ju[i, j] * (q1w[i, j] +
                                  q2w[i, j] * ((u_new[i-1, j+1] + u_old[i, j+1] + u_old[i, j] + u_new[i-1, j]) / 4 -
                                               (u_new[i-1, j] + u_old[i, j] + u_new[i-1, j-1] + u_new[i, j-1]) / 4) / 2)
        Fw_u = 0.5 * dens * (
                    b22w[i, j] * (u_old[i, j] - u_new[i-1, j]) - b21w[i, j] * (v_u_p[i, j] + v_u_p[i+1, j]))
        aw_u = np.append(aw_u, Fw_u + Dw_u)

        # east
        De_u = visc / Ju[i, j] * (q1w[i, j] -
                            q2w[i, j] * ((u_old[i, j+1] + u_old[i+1, j+1] + u_old[i+1, j] + u_old[i, j]) / 4 -
                                               (u_old[i+1, j] + u_old[i, j] + u_new[i, j-1] + u_old[i+1, j-1]) / 4) / 2)
        Fe_u = - 0.5 * dens * (b22w[i, j] * (- u_old[i, j] + u_old[i+1, j]) - b21w[i, j] * (-v_u_p[i, j] + v_u_p[i+1, j]))
        ae_u = np.append(ae_u, Fe_u + De_u)


        # north
        Dn_u = visc / Ju[i, j] * (q3w[i, j]-
                                      q2w[i, j] * (-(u_new[i-1, j+1] + u_old[i, j] + u_old[i, j+1] + u_new[i-1, j]) / 4 +
                                            (u_old[i, j] + u_old[i, j+1] + u_old[i+1, j+1] + u_old[i+1, j]) / 4) / 2)
        Fn_u = - 0.5 * dens * (b11w[i, j] * (- v_u_p[i, j] + v_u_p[i, j+1]) - b12w[i, j] * (- u_old[i, j] + u_old[i, j+1]))
        an_u = np.append(an_u, Fn_u + Dn_u)

        # south
        Ds_u = visc / Ju[i, j] * (q3w[i, j] +
                                      q2w[i, j] * ((u_new[i-1, j] + u_old[i, j] + u_new[i-1, j-1] + u_new[i, j-1]) / 4 -
                                            (u_old[i, j] + u_old[i+1, j] + u_old[i+1, j-1] + u_new[i, j-1]) / 4) / 2)
        Fs_u = 0.5 * dens * (b11w[i, j] * (v_u_p[i, j] - v_u_p[i, j-1]) - b12w[i, j] * (u_old[i, j] - u_new[i, j-1]))
        as_u = np.append(as_u, Fs_u + Ds_u)

        dp_u = np.append(dp_u_x, b22w[i, j] * p_old[i+1, j] - b22w[i, j] * p_old[i-1, j] +
                         b12w[i, j] * p_old[i, j+1] - b12w[i, j] * p_old[i, j-1])
        if j == 1 or j == Ny - 1:
            Sp_u = np.append(Sp_u, - b12w[i, j] * visc)
        else:
            Sp_u = np.append(Sp_u, 0)

    ap_u = aw_u + ae_u + an_u + as_u
    ap_u_big = np.vstack([ap_u_big, ap_u])
    # Thomas for
    A = -an_u
    B = ap_u
    C = -as_u
    K = dp_u + ae_u + aw_u + Sp_u
    # K[-1] = K[-1] - C[-1] * u[-1] kapoiou eidous or. sun8iki
    #Thomas
    for k in range(1, len(A)):
        mc = A[k] / B[k - 1]
        B[k] = B[k] - mc * C[k - 1]
        K[k] = K[k] - mc * K[k - 1]
    Thomas = A.copy()
    Thomas[-1] = K[-1] / B[-1]
    for l in range(len(A)-2, -1, -1):
        Thomas[l] = (K[l] - C[l] * Thomas[l + 1]) / B[l]
    for l in range(1, Ny - 2):
        u_new[i, l] = Thomas[l - 1]


aa = pd.DataFrame(u_old)
bb = pd.DataFrame(u_new)
aa.to_excel('aa.xlsx')
bb.to_excel('bb.xlsx')
print(u_new)



for i in range(1, len(xv[:, 0]) - 2):
    ap_v = []
    ae_v = []
    aw_v = []
    an_v = []
    as_v = []
    dp_v = []
    Sp_v = []
    for j in range(1, len(yv[0, :]) - 1):
        # v
        u_v_p[i, j] = (u_old[i, j + 1] + u_new[i-1, j + 1] + u_new[i-1, j] + u_old[i, j]) / 4
        # west
        Dw_v = visc / Jv[i, j] * (q1s[i, j] +
                                      q2s[i, j] * (-(v_new[i-1, j] + v_new[i-1, j-1] + v_new[i, j-1] + v_old[i, j]) / 4 +
                                            (v_new[i-1, j+1] + v_old[i, j] + v_old[i, j+1] + v_new[i-1, j]) / 4) / 2)
        Fw_v = 0.5 * dens * (b22s[i, j] * (u_v_p[i, j] - u_v_p[i-1, j]) - b21s[i, j] * (v_old[i, j] - v_new[i-1, j]))
        aw_v = np.append(aw_v, Dw_v + Fw_v)

        # east
        De_v = De_u = visc / Jv[i, j] * (q1s[i, j] -
                                         q2s[i, j] * ((v_new[i, j-1] + v_old[i+1, j-1] + v_old[i, j] + v_old[
                    i+1, j]) / 4 -
                                                      (v_old[i, j+1] + v_old[i, j] + v_old[i+1, j+1] + v_old[
                                                          i+1, j]) / 4) / 2)
        Fe_v = - 0.5 * dens * (
                    b22s[i, j] * (-u_v_p[i, j] + u_v_p[i + 1, j]) - b21s[i, j] * (-v_old[i, j] + v_old[i + 1, j]))
        ae_v = np.append(ae_v, De_v + Fe_v)

        # north
        Dn_v = visc / Jv[i, j] * (q3s[i, j] -
                                      q2s[i, j] * ((v_new[i-1, j+1] + v_old[i, j] + v_old[i, j+1] + v_new[i-1, j]) / 4 -
                                            (v_old[i, j] + v_old[i, j+1] + v_old[i+1, j+1] + v_old[i+1, j]) / 4) / 2)
        Fn_v = - 0.5 * dens * (b11s[i, j] * (-v_old[i, j] + v_old[i, j+1]) - b12s[i, j] * (-u_v_p[i, j] + u_v_p[i, j+1]))
        an_v = np.append(an_v, Fn_v + Dn_v)

        # south
        Ds_v = visc / Jv[i, j] * (q3s[i, j]  +
                                      q2s[i, j] * (-(v_new[i-1, j] + v_old[i, j] + v_new[i-1, j-1] + v_new[i, j-1]) / 4 +
                                            (v_old[i, j] + v_old[i+1, j] + v_old[i+1, j-1] + v_new[i, j-1]) / 4) / 2)
        Fs_v = 0.5 * dens * (b11s[i, j] * (v_old[i, j] - v_new[i, j-1]) - b12s[i, j] * (v_u_p[i, j] - v_u_p[i, j-1]))
        as_v = np.append(as_v, Ds_v + Fs_v)

        dp_v = np.append(dp_v, b22s[i, j] * p_old[i+1, j] - b22s[i, j] * p_old[i-1, j]
                         + b12s[i, j] * p_old[i, j+1] - b12s[i, j] * p_old[i, j-1])
        if j == 1 or j == Ny - 1:
            Sp_v = np.append(Sp_v, - b12s[i, j] * visc)
        else:
            Sp_v = np.append(Sp_v, 0)
    ap_v = aw_v + ae_v + an_v + as_v
    ap_v_big = np.vstack([ap_v_big, ap_v])

    A = -an_v
    B = ap_v
    C = -as_v
    K = dp_v + ae_v + aw_v + Sp_v
    # K[-1] = K[-1] - C[-1] * u[-1]
    #Thomas
    for k in range(1, len(A)):
        mc = A[k] / B[k - 1]
        B[k] = B[k] - mc * C[k - 1]
        K[k] = K[k] - mc * K[k - 1]
    Thomas = A.copy()
    Thomas[-1] = K[-1] / B[-1]
    for l in range(len(A)-2, -1, -1):
        Thomas[l] = (K[l] - C[l] * Thomas[l + 1]) / B[l]
    for l in range(1, Ny - 2):
        v_new[i, l] = Thomas[l - 1]


ae_p = []
aw_p = []
an_p = []
as_p = []




print((ap_u_big.shape))

# Pressure Correction
for i in range(1, len(xp[:, 0])-4):
    temp_e = []
    temp_w = []
    temp_n = []
    temp_s = []
    Source = []
    ae_p = []
    aw_p = []
    an_p = []
    as_p = []


    for j in range(1, len(yp[0, :]) - 2):
        if j == 1:
            if i == 1:
                temp_s = 0
                temp_n = b21[i, j] ** 2 / ap_v_big[i, j] + \
                         b11[i, j] ** 2 / (
                                 ap_u_big[i, j] + ap_u_big[i, j+1]) / 4
                temp_w = 0
                temp_e = b12[i, j] ** 2 / ap_u_big[i + 1, j] + \
                         b22[i, j] ** 2 / (0 + 0 + ap_v_big[i + 1, j] + ap_v_big[i, j]) / 4

            elif i == len(xp[:, 0])-4:
                temp_s = 0
                temp_n = b21[i, j] ** 2 / ap_v_big[i, j] + \
                         b11[i, j] ** 2 / (
                                 ap_u_big[i, j] + ap_u_big[i, j + 1]) / 2
                temp_w = b12[i, j] ** 2 / ap_u_big[i, j] + \
                         b22[i, j] ** 2 / (
                                 ap_v_big[i - 1, j] + ap_v_big[i, j] + 0 + 0) / 4
                temp_e = 0

            else:
                temp_s = 0
                temp_n = b21[i, j] ** 2 / ap_v_big[i, j] + \
                         b11[i, j] ** 2 / (
                                     ap_u_big[i, j+1] + ap_u_big[i + 1, j+1] + ap_u_big[i + 1, j ] + ap_u_big[i, j]) / 4
                temp_w = b12[i, j] ** 2 / ap_u_big[i - 1, j] + \
                         b22[i, j] ** 2 / (
                                     ap_v_big[i-1, j] + ap_v_big[i, j] + 0 + 0) / 4
                temp_e = b12[i, j] ** 2 / ap_u_big[i+1, j] + \
                         b22[i, j] ** 2 / (0 + 0 + ap_v_big[i + 1, j] + ap_v_big[i, j]) / 4
        elif j >= len(yp[0, :]) - 4:
            if i == 1:
                print(j)
                temp_s = b21[i, j] ** 2 / ap_v_big[i, j] + \
                         b11[i, j] ** 2 / (
                                     0 + ap_u_big[i + 1, j] + ap_u_big[i + 1, j - 1] + 0) / 4
                temp_n = 0
                temp_w = 0
                temp_e = b12[i, j] ** 2 / ap_u_big[i + 1, j] + \
                         b22[i, j] ** 2 / (0 + 0 + ap_v_big[i, j] + ap_v_big[i+1, j]) / 4

            elif i == len(yp[0, :]) - 4:
                temp_s = b21[i, j] ** 2 / ap_v_big[i, j] + \
                         b11[i, j] ** 2 / (
                                     ap_u_big[i, j] + 0 + 0 + ap_u_big[i, j - 1]) / 4
                temp_n = 0
                temp_w = b12[i, j] ** 2 / ap_u_big[i, j] + \
                         b22[i, j] ** 2 / (
                                 ap_v_big[i-1, j] + ap_v_big[i, j] + 0 + 0) / 4
                temp_e = 0

        else:
            print(j)

            temp_e = b12[i, j]**2 / ap_u_big[i+1, j] +\
                     b22[i, j]**2 / (ap_v_big[i, j+1] + ap_v_big[i+1, j+1] + ap_v_big[i+1, j] + ap_v_big[i, j]) / 4

            temp_w = b12[i, j] ** 2 / ap_u_big[i - 1, j] + \
                     b22[i, j] ** 2 / (
                                 ap_v_big[i, j + 1] + ap_v_big[i - 1, j + 1] + ap_v_big[i - 1, j] + ap_v_big[i, j]) / 4
            temp_n = b21[i, j] ** 2 / ap_v_big[i, j + 1] + \
                     b11[i, j] ** 2 / (
                                 ap_u_big[i, j] + ap_u_big[i + 1, j] + ap_u_big[i + 1, j - 1] + ap_u_big[i, j - 1]) / 4

            temp_s = b21[i, j] ** 2 / ap_v_big[i, j - 1] + \
                     b11[i, j] ** 2 / (
                             ap_u_big[i, j] + ap_u_big[i + 1, j] + ap_u_big[i + 1, j + 1] + ap_u_big[i, j + 1]) / 4
        Source = np.append(Source, dens * 0.5 * (b22[i, j] * (u_new[i+1, j] - u_new[i-1, j]) -
                           b21[i, j] *0.25 *  ((u_new[i, j] + u_new[i+1, j] + u_new[i+1, j+1] + u_new[i, j+1]) -
                                               u_new[i, j] + u_new[i+1, j] + u_new[i+1, j-1] + u_new[i, j-1]) +
                           b11[i, j] * 0.25 * ((v_new[i, j] + v_new[i+1, j] + v_new[i+1, j+1] + v_new[i, j+1]) -
                                               (v_new[i, j] + v_new[i-1, j] + v_new[i-1, j+1] + v_new[i, j+1])) -
                           b12[i, j] * (v_new[i, j+1] - v_new[i, j])))

        ae_p = np.append(ae_p, temp_e)
        aw_p = np.append(aw_p, temp_w)
        an_p = np.append(an_p, temp_n)
        as_p = np.append(as_p, temp_s)

    ap_p = aw_p + ae_p + an_p + as_p
    A = -an_p
    B = ap_p
    C = -as_p
    K = Source + ae_p + aw_p
    # K[-1] = K[-1] - C[-1] * u[-1]
    # Thomas
    print(len(A))
    for k in range(1, len(A)):
        mc = A[k] / B[k - 1]
        B[k] = B[k] - mc * C[k - 1]
        K[k] = K[k] - mc * K[k - 1]
    Thomas = A.copy()
    Thomas[-1] = K[-1] / B[-1]
    for l in range(len(A) - 2, -1, -1):
        Thomas[l] = (K[l] - C[l] * Thomas[l + 1]) / B[l]
    print(len(Thomas))
    for l in range(1, len(A)-1):
        p_new[i, l] = Thomas[l - 1]


p_final = p_old + omega * p_new
cc = pd.DataFrame(v_new)
cc.to_excel('cc.xlsx')
dd = pd.DataFrame(p_new)
dd.to_excel('cc.xlsx')
ee = pd.DataFrame(p_final)
ee.to_excel('dd.xlsx')








