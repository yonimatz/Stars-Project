"""
Created on Fri Mar 20 22:59:31 2020
@author: matthewsmith
"""
import numpy as np
import matplotlib.pyplot as plt
import time

start_time = time.time()

# =============================================================================
# Defining constants
# =============================================================================

G = 6.674e-11  # Gravitational Constant
c = 2.998e8  # Speed of light
hbar = 6.626e-34 / (2 * np.pi)
k = 1.381e-23  # Boltzmann Constant
sigma = 5.670e-8  # Stefan-Boltzmann Constant
m_e = 9.109e-31  # Electron mass
M_s = 1.989e30  # Mass of Sun
R_s = 6.963e8  # Radius of Sun
L_s = 3.828e26  # Luminosity of Sun
m_p = 1.6726219e-27  # Mass of proton
a = 4 * sigma / c
pi = np.pi
ep_pp_coeff = 1.07e-7 * 1e-5 * (1e-6) ** 4
ep_cno_coeff = 8.24e-26 * 1e-5 * (1e-6) ** 19.9
nonrelgenpress = (3 * pi ** 2) ** (2 / 3) / 5 * hbar ** 2 / m_e * m_p ** (-5 / 3)
mach_ep = np.finfo(np.float64).eps
tiny_float = 1e-20
gamma = 5 / 3  # ideal gas constant

# =================Defining mu, XYZ=====================================
# X = 0.738
# Y = 0.249

X = 0.8
Y = 0.2
Z = 1.0 - (X + Y) + 1e-15

mu = 1.0 / (2.0 * X + 0.75 * Y + 0.5 * Z)

# mu = mean molecular weight for fully ionized gas

# =============================================================================
# Pressure equations
# =============================================================================
def Pressure(rho, T):
    const1 = (((3.0 * pi ** 2.0) ** (2.0 / 3.0)) / 5.0) * ((hbar ** 2.0) / m_e)
    term1 = const1 * (rho / m_p) ** (5.0 / 3.0)

    const2 = k / (mu * m_p)
    term2 = const2 * rho * T

    const3 = a / 3.0
    term3 = const3 * (T ** 4.0)
    return term1 + term2 + term3


# Partial derivative of pressure w.r.t density
def DelPressureDelrho(rho, T):
    const1 = (((3.0 * pi ** 2.0) ** (2.0 / 3.0)) / 3.0) * ((hbar ** 2.0) / (m_e * m_p))
    term1 = const1 * (rho / m_p) ** (2.0 / 3.0)

    const2 = k / (mu * m_p)
    term2 = T * const2

    return term1 + term2


# Partial derivitave of pressure w.r.t temperature
def DelPressuredelT(rho, T):
    const1 = k / (mu * m_p)
    term1 = rho * const1

    const2 = 4.0 * a / 3.0
    term2 = const2 * (T ** 3.0)

    return term1 + term2


# =============================================================================
# kappa definitions
# =============================================================================

# Kappa H-
def kappaHminus(rho, T):
    rho3 = rho / 1.0e3

    return 2.5e-32 * (Z / 0.02) * (rho3 ** 0.5) * (T ** 9.0)


# Kappa ff
def Kappaff(rho, T):
    rho3 = rho / 1.0e3
    return 1.0e24 * (Z + 0.0001) * (rho3 ** 0.7) * (T ** -3.5)


# kappa es
def Kappaes():
    return 0.02 * (1.0 + X)


# the real slim kappa
def Kappa(rho, T):
    term1 = 1.0 / kappaHminus(rho, T)
    # if Z == 0:
    #     term1 = 0
    if T > 1.0e4:
        term2 = 1.0 / max(
            [Kappaff(rho, T), Kappaes()])  # max fn take itterale and returns largest, so i put the kappas in a list
    else:
        term2 = 1.0 / max([Kappaff(rho, T), Kappaes()])  # what andrew did, changed to max
    return 1.0 / (term1 + term2)


# =============================================================================
# Energy generation epsilon
# =============================================================================

def EpsilonPP(rho, T):
    rho5 = rho / 1.0e5
    T6 = T / 1e6
    return 1.07e-7 * rho5 * (X ** 2.0) * T6 ** 4.0


def EpsilonCNO(rho, T):
    rho5 = rho / 1.0e5
    T6 = T / 1.0e6
    Xcno = 0.03 * X
    return 8.24e-26 * rho5 * X * Xcno * (T6 ** 19.9)


def Epsilon(rho, T):
    return EpsilonPP(rho, T) + EpsilonCNO(rho, T)


# =============================================================================
# The Big 5 Equations
# =============================================================================
# 1) Partial derivative of rho w.r.t radius

def DRhoDr(r, rho, M, kappa, DelPDelrho, DelPressuredelT, DTempDr):
    term1 = -1.0 * (G * M * rho) / (r ** 2.0)

    term2 = DelPressuredelT * DTempDr

    return (term1 + term2) / DelPDelrho


# 2) Derivitave of Temperature w.r.t radius
def DTDr(r, rho, T, L, M, kappa, P):
    const1 = 3.0 / (16.0 * pi * a * c)

    term1 = (const1 * rho * kappa * L) / ((T ** 3.0) * (r ** 2.0))
    # print(kappa)
    const2 = (1.0 - (1.0 / gamma)) * G
    term2 = const2 * (T * M * rho) / (P * (r ** 2.0))
    return -1.0 * min([term1, term2])


# 3) Derivative of Mass w.r.t radius
def DMDr(r, rho):
    return 4.0 * pi * (r ** 2.0) * rho


# 4) Derivative of Luminosity w.r.t radius
def DLDr(r, rho, epsilon):
    return 4.0 * pi * (r ** 2.0) * rho * epsilon


# 5) Derivative of optical depth w.r.t radius
def DTauDr(kappa, rho):
    return kappa * rho


# =============================================================================
# Runge Kutta
# =============================================================================
# For runge kutta Solving for K values


# Runge-Kutta-K 2's-RVal
RKK2R = 2.500000000000000e-01  # 1/4
# ...
RKK3R = 3.750000000000000e-01  # 3/8
RKK4R = 9.230769230769231e-01  # 12/13
RKK5R = 1.000000000000000e+00  # 1
RKK6R = 5.000000000000000e-01  # 1/2

# Runge-Kukka-K 2's- K1val
RK21 = 2.500000000000000e-01  # 1/4

# Runge-Kukka-K 3's- K1val
RK31 = 9.375000000000000e-02  # 3/32
# Runge-Kukka-K 3's- K2val
RK32 = 2.812500000000000e-01  # 9/32

RK41 = 8.793809740555303e-01  # 1932/2197
RK42 = -3.277196176604461e+00  # -7200/2197
RK43 = 3.320892125625853e+00  # 7296/2197

RK51 = 2.032407407407407e+00  # 439/216
RK52 = -8.000000000000000e+00  # -8
RK53 = 7.173489278752436e+00  # 3680/513
RK54 = -2.058966861598441e-01  # -845/4104

RK61 = -2.962962962962963e-01  # -8/27
RK62 = 2.000000000000000e+00  # 2
RK63 = -1.381676413255361e+00  # -3544/2565
RK64 = 4.529727095516569e-01  # 1859/4104
RK65 = -2.750000000000000e-01  # -11/40

# constants for summing all x together

# Runge-Kukka-4thorder-X- K1val
RK4X1 = 1.157407407407407e-01  # 25/216
RK4X3 = 5.489278752436647e-01  # 1408/2565
RK4X4 = 5.353313840155945e-01  # 2197/4104
RK4X5 = -2.00000000000000e-01  # -1/5

RK5X1 = 1.185185185185185e-01  # 16 / 135
RK5X3 = 5.189863547758284e-01  # 6656 / 12825
RK5X4 = 5.061314903420167e-01  # 28561 / 56430
RK5X5 = -1.80000000000000e-01  # -9 / 50
RK5X6 = 3.636363636363636e-02  # 2 / 55

r1 = 2.777777777777778e-03  # 1/360
r3 = -2.994152046783626e-02  # -128/4275
r4 = -2.919989367357789e-02  # -2197/75240
r5 = 2.000000000000000e-02  # 1/50
r6 = 3.636363636363636e-02  # 2/55


def MakeDerivitaves(vals, currentR, i, listofK):
    KTermsrho = np.sum(listofK[:, 0])
    KTermsT = np.sum(listofK[:, 1])
    KTermsM = np.sum(listofK[:, 2])
    KTermsL = np.sum(listofK[:, 3])
    KTermstau = np.sum(listofK[:, 4])

    rho = vals[0, i] + KTermsrho
    T = vals[1, i] + KTermsT
    M = vals[2, i] + KTermsM
    L = vals[3, i] + KTermsL
    tau = vals[4, i] + KTermstau
    kappa = Kappa(rho, T)
    epsilon = Epsilon(rho, T)
    pressure = Pressure(rho, T)
    dpdt = DelPressuredelT(rho, T)
    dpdrho = DelPressureDelrho(rho, T)

    dTau = DTauDr(kappa, rho)
    dL = DLDr(currentR, rho, epsilon)
    dM = DMDr(currentR, rho)
    dT = DTDr(currentR, rho, T, L, M, kappa, pressure)
    dRho = DRhoDr(currentR, rho, M, kappa, dpdrho, dpdt, dT)
    # print(dRho)
    derivitaves = np.array([dRho, dT, dM, dL, dTau])
    return derivitaves


def CheckToStop(rho, T, M, drho):
    kappa = Kappa(rho, T)
    deltaTau = kappa * rho ** 2.0 / abs(drho)
    # print(deltaTau)
    if deltaTau < 0.00001:  # @@@@@@@@this is just a placeholder!!!!@@@@@@@
        return False, deltaTau

    if M > 1e3 * M_s:
        return False, deltaTau

    return True, deltaTau


def extendArray(arrayToExtend, dimensions):
    extender = np.zeros(dimensions)
    arrayToExtend = np.column_stack((arrayToExtend, extender))
    return arrayToExtend


def RK45(r_0, h_0, vals, minError, maxError, hmin, hmax, abort):
    # merge the functions
    i = 0  # itter counter
    dL = [0]
    r = [r_0]  # initial r   next r is r+h
    h = h_0
    errorList = []
    keepGoing = True
    zero = np.array([0., 0., 0., 0., 0.])
    vals = extendArray(vals, (5, 1000))
    kappa = [Kappa(vals[0, 0], vals[1, 0])]
    delP = [DelPressuredelT(vals[0, 0], vals[1, 0])]

    while keepGoing and i < abort:

        if h < hmin:
            h = hmin

        elif h > hmax:
            h = hmax

        #        try:
        derivitaves = MakeDerivitaves(vals, r[i], i,
                                      np.array([zero, zero, zero, zero, zero]))  # keep this becuase its usefull later
        K1 = h * derivitaves

        K2 = h * MakeDerivitaves(vals, r[i] + RKK2R * h, i,
                                 np.array([RK21 * K1, zero, zero, zero, zero]))

        K3 = h * MakeDerivitaves(vals, r[i] + RKK3R * h, i,
                                 np.array([RK31 * K1, RK32 * K2, zero, zero, zero]))

        K4 = h * MakeDerivitaves(vals, r[i] + RKK4R * h, i,
                                 np.array([RK41 * K1, RK42 * K2, RK43 * K3, zero, zero]))

        K5 = h * MakeDerivitaves(vals, r[i] + RKK5R * h, i,
                                 np.array([RK51 * K1, RK52 * K2, RK53 * K3, RK54 * K4, zero]))

        K6 = h * MakeDerivitaves(vals, r[i] + RKK6R * h, i,
                                 np.array([RK61 * K1, RK62 * K2, RK63 * K3, RK64 * K4, RK65 * K5]))

        RK4X = vals[:, i] + (RK4X1 * K1 + RK4X3 * K3 + RK4X4 * K4 + RK4X5 * K5)

        RK5X = vals[:, i] + (RK5X1 * K1 + RK5X3 * K3 + RK5X4 * K4 + RK5X5 * K5 + RK5X6 * K6)

        # R = abs(r1*K1 + r3*K3 + r4*K4 + r5*K5 + r6*K6) / (RK5X)

        errorlist = abs((RK4X - RK5X) / RK5X)

        # print(abs(RK4X - RK5X)/RK5X)
        # print(R)

        error = np.max(errorlist)

        if error > maxError and h > hmin:

            # h = h * (maxError / error)**.2
            h = h / 2
            continue

        else:
            vals[:, i + 1] = RK4X
            # print(RK4X)

            kappa.append(Kappa(RK4X[0], RK4X[1]))
            delP.append(DelPressuredelT(RK4X[0], RK4X[1]))
            drhodr = derivitaves[0]

            keepGoing, deltaTau = CheckToStop(RK4X[0], RK4X[1], RK4X[2], drhodr)

            dL.append(derivitaves[3])

            r.append(r[i] + h)

            i += 1
            errorList.append(error)

            if i % 1000 == 0:
                if i == abort:
                    pass
                else:
                    try:
                        vals = extendArray(vals, (5, 1000))
                    except:
                        print('Error: Array could not be extended')

                # try:
                #     print('Keep Going Criteria: ' + str(deltaTau))
                #     print('Current Step (x1e3): ' + str(i/1000))
                #     print('-------------------------------------')
                # except:
                #     pass
            if error < minError:
                # print(h)
                # h = h * (maxError / error)**.2
                h = h * 2
                # print(h)
    return vals, r, dL, kappa, errorList, delP


# make maximum error 5% to the previous value IE if T=1e6, the max allowed error is a 1% difference
# make minimum error 0.1%

maxError = 1E-9
minError = 3E-11
hmin = 100
hmax = 10000
h_0 = 5000



## --------------------------------------------------
# Plotting the main sequence:
## --------------------------------------------------
# Full lists of T_c and rho_c values we need to plot MS, if you want more data points simply add
# corresponding T_c and rho_c values to the end of the list.

# Example data set:
T_c = [1e4,2e4,3e4,4e4,5e4,6e4,7e4,8e4,9e4,1e5,1e5, 2e5, 3e5, 4e5, 5e5,6e5,7e5,8e5,9e5,1e6,2e6,3e6,5e6,7e6,9e6]

rho_c_list = [3.00000000e+09, 2.98000000e+09, 2.50000000e+09, 1.41469594e+09, 9.00000000e+08,
              5.83200000e+08, 4.37400000e+08, 3.24000000e+08, 2.34812358e+08, 1.77156100e+08,
              2e8,3.8e7,1.5e7,7e6,5.5e6,3.5e6,2.5e6,1.9e6,1.7e6,1.4e6,9e5,8e5,6.3e5,5.4e5,4.5e5]

# Additional data points function:
T_c_extension = []
T_c_extension_2 = []
T_c_extension_3 = []

rho_c_list_extension = []
rho_c_list_extension_2 = []
rho_c_list_extension_3 = []

for t in T_c:
    T_c_extension.append(t * 1.35)
    T_c_extension_2.append(t/1.75)
    T_c_extension_3.append(t/1.35)

for r in rho_c_list:
    rho_c_list_extension.append(r / 1.35)
    rho_c_list_extension_2.append(r * 1.75)
    rho_c_list_extension_3.append(r * 1.3)

T_c = T_c + T_c_extension + T_c_extension_2 + T_c_extension_3
rho_c_list = rho_c_list + rho_c_list_extension + rho_c_list_extension_2 + rho_c_list_extension_3

# T_c = [1e4, 2e4, 3e4]
#
# rho_c_list = [3e09, 2e9, 1e9]

filename = 'Stars_2_main_extended_100stars.csv'
title = 'L (row 1), T (row 2), M (row 3), R (row 4), T_c (row 5), rho_c (row 6), tau_0 = 0.231739, X = 0.738, Y = 0.249'
X = 0.738
Y = 0.249

## --------------------------------------------------
## --------------------------------------------------

values = np.zeros((5, 1))
numstars = len(T_c) # Number of stars we use
Stars = np.zeros((6, numstars)) # Make stars array for the ability to do several stars

i = 0
while i <= numstars - 1:
    # rho_c = 777460 # 1 magnitude larger then stated
    # T_c = 2e+07 # idx1: T
    rho_c = rho_c_list[i]

    r_0 = 1000
    M_0 = 4.0 * np.pi * r_0 ** 3.0 * rho_c / 3.0  # idx2: M
    L_0 = 4.0 * np.pi * r_0 ** 3.0 * rho_c * Epsilon(rho_c, T_c[i]) / 3.0  # idx3: L
    tau_0 = 0.231739  # possible value 1 from OPAL TABLE

    # initialiize the vals array
    values[0] = rho_c  # idx0: rho
    values[1] = T_c[i]  # idx1: T
    values[2] = M_0  # idx2: M
    values[3] = L_0  # idx3: L
    values[4] = tau_0  # idx4: tau

    match = False

    while (match == False): # While the error in L is bigger than we allow this statement runs
        vals, r, dL, kappa, errors, delP = RK45(r_0, h_0, values, minError, maxError, hmin, hmax, 100000)
        tau_surface = np.max(vals[4, :]) - (2 / 3) # Defining the surface as tau(infinity) - tau(2/3)
        idx = np.abs(vals[4] - tau_surface).argmin() # index of the smallest difference between tau values and tau at the surface
        R_surface = r[idx] # Location of the surface

        T_surface = vals[1, idx] # Temperature at the surface
        M_surface = vals[2, idx]  # Temperature at the surface
        L_surface = vals[3, idx]  # Luminosity at the surface

        L_error = (L_surface - (4 * pi * sigma * (R_surface ** 2) * (T_surface ** 4))) / (
            np.sqrt(4 * pi * sigma * (R_surface ** 2) * (T_surface ** 4) * L_surface)) # error in luminosity using the surface conditions
        print(L_error)

        if L_error < 0.1 and L_error > -0.1: # If the error is small enough:
            match = True
            Stars[0, i] = L_surface  # Luminosity of particular star added to list of stars
            Stars[1, i] = T_surface  # Temperature of particular star added to list of stars
            Stars[2, i] = M_surface  # Mass of particular star added to list of stars
            Stars[3, i] = R_surface  # Radius of particular star added to list of stars
            Stars[4, i] = T_c[i]  # Radius of particular star added to list of stars
            Stars[5, i] = values[0]  # Radius of particular star added to list of stars

            print('i = ', i)
            print("R Surface: " + str(R_surface)) # Prints the radius of the star's surface
            i += 1
        else:
            if L_error > 0:
                values[0] = values[0] * 0.9 # If error is too large, make rho_c smaller
            else:
                values[0] = values[0] * 1.1 # If error is too large, make rho_c smaller

vals = vals[:, :len(r)] # vals get cut off at the star's radius
errors.append(errors[-1]) # Append the error onto the error list

# @@@@@@@ Abort variable must be greater then and divisible by 10000 @@@@@@
stepsize = [h_0]
for idx, step in enumerate(r):
    if idx != 0:
        stepsize.append(r[idx] - r[idx - 1])

maxr = np.max(r)

xaxis = r / maxr

# =======================================================================================================
## Optional code to plot individual star data:

# fig, axs = plt.subplots(2, 3, sharex=True, frameon=False, figsize=(10, 5))
# axs[0, 0].plot(xaxis, vals[0, :]/rho_c, 'k-')
# axs[0, 0].set_title('Density (rho/rho_c)', fontsize=8)
# axs[0, 1].plot(xaxis, vals[1, :]/T_c, 'k-')
# axs[0, 1].set_title('Temperature (T/T_c)', fontsize=8)
# axs[0, 2].plot(xaxis, vals[2, :]/np.max(vals[2,:]), 'k-')
# axs[0, 2].set_title('Mass (M/M_max)', fontsize=8)
# axs[1, 0].plot(xaxis, dL, 'k-')
# axs[1, 0].set_title('Derivitave of Luminosity', fontsize=8)
# axs[1, 1].plot(xaxis, kappa, 'k-')
# axs[1, 1].set_title('Kappa', fontsize=8)
# axs[1, 1].set_yscale('log')
# axs[1, 2].plot(xaxis, stepsize, 'k-')
# axs[1, 2].set_title('stepsize', fontsize=8)
# axs[1, 2].set_yscale('log')

# plt.xscale('log')

# =======================================================================================================
# Saving data to CSV file with given names

data = Stars
np.savetxt(filename, data, delimiter=",", header=title)

# =======================================================================================================



# =======================================================================================================
# Plotting the MS plots for the given data

Stars_L = Stars[0, :]
Stars_T = Stars[1, :]
Stars_M = Stars[2, :]
Stars_R = Stars[3, :]
Stars_T_c = Stars[4, :]
Stars_rho_c = Stars[5, :]
print(Stars[:, :])

print("--- %s Minutes ---" % str(int(time.time() - start_time) / 60))


# Plotting L-T
plt.title('L-T plot', fontsize=8)
plt.rcParams.update({'font.size': 8})
plt.yscale('log')
plt.ylim((min(Stars_L), max(Stars_L)))
plt.xscale('log')
plt.xlim((max(Stars_T), min(Stars_T)))
plt.scatter(Stars_T, Stars_L, c='r')
plt.tight_layout()
plt.show()

# --------------------------------------
# Additional MS plots
# Plotting L-M
plt.title('L-M plot', fontsize=8)
plt.rcParams.update({'font.size': 8})
plt.yscale('log')
plt.ylim((min(Stars_L), max(Stars_L)))
plt.xscale('log')
plt.xlim((max(Stars_M), min(Stars_M)))
plt.scatter(Stars_M, Stars_L, c='r')
plt.tight_layout()
plt.show()

# ---------------------------------------
# Plotting M-R
plt.title('M-R plot', fontsize=8)
plt.rcParams.update({'font.size': 8})
plt.yscale('log')
plt.ylim((min(Stars_M), max(Stars_M)))
plt.xscale('log')
plt.xlim((max(Stars_R), min(Stars_R)))
plt.scatter(Stars_R, Stars_M, c='r')
plt.tight_layout()
plt.show()
