import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
from pylab import meshgrid

# Mesh reading
mesh=np.loadtxt('mesh.dat')
x=np.zeros(0)
y=np.zeros(0)
x=np.append(x,mesh[:,0])
y=np.append(y,mesh[:,1])

#Array initialization

xx=np.zeros((200, 60))
yy=np.zeros((200, 60))
xp=np.zeros((200, 60))
yp=np.zeros((200, 60))
xu=np.zeros((200, 60))
yu=np.zeros((200, 60))
xv=np.zeros((200, 60))
yv=np.zeros((200, 60))
x_y=np.zeros((60, 200))
y_y=np.zeros((60, 200))


b11=np.zeros((200, 60))
b12=np.zeros((200, 60))
b21=np.zeros((200, 60))
b22=np.zeros((200, 60))
b11w=np.zeros((200, 60))
b12w=np.zeros((200, 60))
b21w=np.zeros((200, 60))
b22w=np.zeros((200, 60))
b11s=np.zeros((200, 60))
b12s=np.zeros((200, 60))
b21s=np.zeros((200, 60))
b22s=np.zeros((200, 60))
Jp=np.zeros((200, 60))
Ju=np.zeros((200, 60))
Jv=np.zeros((200, 60))
volp=np.zeros((200, 60))
volu=np.zeros((200, 60))
volv=np.zeros((200, 60))
pcheck=np.zeros((200, 60))
ucheck=np.zeros((200, 60))
vcheck=np.zeros((200, 60))

#second order central
b11ksi=np.zeros((200, 60))
b11hta=np.zeros((200, 60))
b12ksi=np.zeros((200, 60))
b12hta=np.zeros((200, 60))
b21ksi=np.zeros((200, 60))
b21hta=np.zeros((200, 60))
b22ksi=np.zeros((200, 60))
b22hta=np.zeros((200, 60))

#second order west
b11ksiw=np.zeros((200, 60))
b11htaw=np.zeros((200, 60))
b12ksiw=np.zeros((200, 60))
b12htaw=np.zeros((200, 60))
b21ksiw=np.zeros((200, 60))
b21htaw=np.zeros((200, 60))
b22ksiw=np.zeros((200, 60))
b22htaw=np.zeros((200, 60))

#Second order south
b11ksis=np.zeros((200, 60))
b11htas=np.zeros((200, 60))
b12ksis=np.zeros((200, 60))
b12htas=np.zeros((200, 60))
b21ksis=np.zeros((200, 60))
b21htas=np.zeros((200, 60))
b22ksis=np.zeros((200, 60))
b22htas=np.zeros((200, 60))



# Grid points
a=0
while a <len(x):
    for j in np.arange(0, 60):
        for i in np.arange(0, 200):
            xx[i, j] = x[a]
            yy[i, j] = y[a]
            a += 1

# Reading by y
b=0
while b <len(y):
    for i in np.arange(0, 60):
        for j in np.arange(0, 200):
            x_y[i, j] = x[b]
            y_y[i, j] = y[b]
            b += 1


# Center of cell - Pressure Points
for i in range(1, len(xx[:, 0])):
    for j in range(1, len(xx[0, :])):
        xp[i-1, j-1] = 0.25*(xx[i-1, j-1]+xx[i, j-1]+xx[i, j]+xx[i-1, j])
        yp[i-1, j-1] = 0.25*(yy[i-1, j-1]+yy[i, j-1]+yy[i, j]+yy[i-1, j])
# Correction for the last row
xp[-1, :] = xp[-2, :]
yp[-1, :] = yp[-2, :]

# u-cell ( West face)

for i in range(len(xx[:, 0])):
    for j in range(1, len(xx[0, :])):
        xu[i, j-1] = 0.5*(xx[i, j-1]+xx[i, j])
        yu[i, j-1] = 0.5*(yy[i, j-1]+yy[i, j])

# v-cell ( South face)

for i in range(1, len(xx[:, 0])):
    for j in range(len(xx[0, :])):
        xv[i-1, j]=0.5*(xx[i-1, j]+xx[i, j])
        yv[i-1, j]=0.5*(yy[i-1, j]+yy[i, j])

# Correction for the last row
xv[-1, :] = xv[-2, :]
yv[-1, :] = yv[-2, :]

plt.figure()
plt.title("όγκοι ελέγχου physical plane")
# for i in range(0, len(xp[:, 0])):
#     plt.plot(xv[i, :-1], yv[i, :-1], c='red')
# plt.plot(xv, yv, c='red')

for i in range(0, len(xv[:, 0])):
    plt.plot(xu[i, :-1], yu[i, :-1], c='green')
plt.plot(xu, yu, c='green')
plt.scatter(xp, yp, c='black', s=1)

# Staggered grids
# plt.plot(xx, yy, c='black')
# plt.plot(x_y, y_y, c='black')
# plt.figure()
for i in range(0, len(xv[:, 0])):
    plt.plot(xu[i, :-1], yu[i, :-1], c='black')
plt.plot(xv, yv, c='black')
plt.xlabel('x')
plt.ylabel('y')

plt.show()
# Metrics and Jacobians

# Center of c.v.
for i in np.arange(1,199):
    for j in np.arange(1,59):
        b11[i, j] = xu[i+1,j] - xu[i, j]
        b22[i, j] = yv[i, j+1]-yv[i, j]
        b12[i, j] = yu[i+1, j]-yu[i, j]
        b21[i, j] = xv[i, j+1]-xv[i, j]

        # West face of c.v.
        b11w[i, j] = xp[i, j]-xp[i-1, j]
        b22w[i, j] = yy[i, j+1]-yy[i, j]
        b12w[i, j] = yp[i, j]-yp[i-1, j]
        b21w[i, j] = xx[i, j+1]-xx[i, j]
        if i==1:
            b11w[i, j] = -3*xu[i,j]+4*xp[i,j]-xu[i+1,j]
            b12w[i, j] = -3*yu[i,j]+4*yp[i,j]-yu[i+1,j]

        if i==-1:
            b11w[i+1,j] = 3*xu[i,j]-4*xp[i-1,j]+xu[i-1,j]
            b12w[i+1,j] = 3*yu[i,j]-4*yp[i-1,j]+yu[i-1,j]
            b22w[i+1,j] = yy[i+1,j+1]-yy[i+1,j]
            b21w[i+1,j] = xx[i+1,j+1]-xx[i+1,j]

        # South face of c.v
        b11s[i,j] = xx[i+1,j]-xx[i,j]
        b22s[i,j] = yp[i,j]-yp[i,j-1]
        b12s[i,j] = yy[i+1,j]-yy[i,j]
        b21s[i,j] = xp[i,j]-xp[i,j-1]
        if i==1:
            b22s[i,j] = -3*yv[i,j]+4*yp[i,j]-yv[i,j+1]
            b21s[i,j] = -3*xv[i,j]+4*xp[i,j]-xv[i,j+1]
        if i==-1:
            b11s[i,j+1] = xx[i+1,j+1]-xx[i,j+1]
            b12s[i,j+1] = yy[i+1,j+1]-yy[i,j+1]
            b21s[i,j+1] = 3*xv[i,j]-4*xp[i,j-1]+xv[i,j-1]
            b22s[i,j+1] = 3*yv[i,j]-4*yp[i,j-1]+yv[i,j-1]

        # Special treatment on Boundaries / south
        if i==199:
            b11s[i+1,j] = 0
            b12s[i+1,j] = 0
            b21s[i+1,j] = xu[i+1,j]-xu[i,j]
            b22[i+1,j] = yu[i+1,j]-yu[i,j]

        # Second order metrics
        #Pressure points

        # d2x/dksi2
        b11ksi[i,j] = b11[i+1,j]-b11[i,j]
        #d2x/dksidhta
        b11hta[i,j] = b11[i,j+1]-b11[i,j]
        #d2y/dksi2
        b12ksi[i,j] = b12[i+1,j]-b12[i,j]
        #d2y/dksidhta
        b12hta[i,j] = b12[i,j+1]-b12[i,j]
        #d2x/dksidhta
        b21ksi[i,j] = b21[i+1,j]-b21[i,j]
        #d2x/dhta2
        b21hta[i,j] = b21[i,j+1]-b21[i,j]
        #d2y/dksidhta
        b22ksi[i,j] = b22[i+1,j]-b22[i,j]
        #d2y/dhta2
        b22hta[i,j] = b22[i,j+1]-b22[i,j]


        # u face
        # d2x/dksi2
        b11ksiw[i, j] = b11w[i + 1, j] - b11w[i, j]
        # d2x/dksidhta
        b11htaw[i, j] = b11w[i, j + 1] - b11w[i, j]
        # d2y/dksi2
        b12ksiw[i, j] = b12w[i + 1, j] - b12w[i, j]
        # d2y/dksidhta
        b12htaw[i, j] = b12w[i, j + 1] - b12w[i, j]
        # d2x/dksidhta
        b21ksiw[i, j] = b21w[i + 1, j] - b21w[i, j]
        # d2x/dhta2
        b21htaw[i, j] = b21w[i, j + 1] - b21w[i, j]
        # d2y/dksidhta
        b22ksiw[i, j] = b22w[i + 1, j] - b22w[i, j]
        # d2y/dhta2
        b22htaw[i, j] = b22w[i, j + 1] - b22w[i, j]


        #v-face
        # d2x/dksi2
        b11ksis[i, j] = b11s[i + 1, j] - b11s[i, j]
        # d2x/dksidhta
        b11htas[i, j] = b11s[i, j + 1] - b11s[i, j]
        # d2y/dksi2
        b12ksis[i, j] = b12s[i + 1, j] - b12s[i, j]
        # d2y/dksidhta
        b12htas[i, j] = b12s[i, j + 1] - b12s[i, j]
        # d2x/dksidhta
        b21ksis[i, j] = b21s[i + 1, j] - b21s[i, j]
        # d2x/dhta2
        b21htas[i, j] = b21s[i, j + 1] - b21s[i, j]
        # d2y/dksidhta
        b22ksis[i, j] = b22s[i + 1, j] - b22s[i, j]
        # d2y/dhta2
        b22htas[i, j] = b22s[i, j + 1] - b22s[i, j]

        # Jacobians

        Jp[i,j]=b11[i,j]*b22[i,j]-b12[i,j]*b21[i,j]
        Ju[i, j] = b11w[i, j] * b22w[i, j] - b12w[i, j] * b21w[i, j]
        Jv[i, j] = b11s[i, j] * b22s[i, j] - b12s[i, j] * b21s[i, j]

        #Volumes of c.v

        #Pressure points
        volp[i,j] = abs(0.5*((xx[i+1,j+1]-xx[i,j])*(yy[i,j+1]-yy[i+1,j])-(xx[i,j+1]-xx[i+1,j])*(yy[i+1,j+1]-yy[i,j])))

        #u points
        volu[i,j] = abs(0.5*((xv[i+1,j+1]-xv[i,j])*(yv[i,j+1]-yv[i+1,j])-(xv[i,j+1]-xv[i+1,j])*(yv[i+1,j+1]-yv[i,j])))

        #v points
        volv[i,j] = abs(0.5*((xu[i+1,j+1]-xu[i,j])*(yu[i,j+1]-yu[i+1,j])-(xu[i,j+1]-xu[i+1,j])*(yu[i+1,j+1]-yu[i,j])))


for i in np.arange(1, 199):
    for j in np.arange(1, 59):
        pcheck[i,j] = Jp[i,j]-volp[i,j]
        ucheck[i,j] = Ju[i,j]-volu[i,j]
        vcheck[i,j] = Jv[i,j]-volv[i,j]





# Dataframes

# df1=pd.DataFrame(Jp)
# df1.to_excel('Jp.xlsx')
# df2=pd.DataFrame(Ju)
# df2.to_excel('Ju.xlsx')
# df3=pd.DataFrame(Jv)
# df3.to_excel('Jv.xlsx')
# df4=pd.DataFrame(pcheck)
# df4.to_excel('pcheck.xlsx')
# df5=pd.DataFrame(ucheck)
# df5.to_excel('ucheck.xlsx')
# df6=pd.DataFrame(vcheck)
# df6.to_excel('vcheck.xlsx')


