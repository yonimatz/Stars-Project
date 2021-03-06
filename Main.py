

"""
Created on Fri Mar 20 22:59:31 2020
@author: matthewsmith
Edited by susherov
"""
import numpy as np
import matplotlib.pyplot as plt
import time
start_time = time.time()


# Questions
# do i need ot add a minimum  to the 'kappa' function? or is that implied in his given eqn im not sure
# =============================================================================
# Defining constants
# =============================================================================

G = 6.674e-11  # Gravitational Constant
c = 2.998e8  # Speed of light
h = 6.626e-34  # Planck's constant
hbar = h / (2 * np.pi)
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
X = 0.7381
Y = 0.2485
Z = 1 - (X + Y)
mu = 1 / (2 * X + 0.75 * Y + 0.5 * Z)
# mu = mean molecular weight for fully ionized gas


# =============================================================================
# Pressure equations
# =============================================================================
def Pressure(rho, T):
    const1 = (((3 * pi ** 2) ** (2 / 3)) / 5) * ((hbar ** 2) / m_e)
    term1 = const1 * (rho / m_p) ** (5 / 3)

    const2 = k / (mu * m_p)
    term2 = const2 * rho * T

    const3 = a / 3
    term3 = const3 * (T ** 4)
    return term1 + term2 + term3


# Partial derivative of pressure w.r.t density
def DelPressureDelrho(rho, T):
    const1 = (((3 * pi ** 2) ** (2 / 3)) / 3) * ((hbar ** 2) / (m_e * m_p))
    term1 = const1 * (rho / m_p) ** (2 / 3)

    const2 = k / (mu * m_p)
    term2 = T * const2

    return term1 + term2


# Partial derivitave of pressure w.r.t temperature
def DelPressuredelT(rho, T):
    const1 = k / (mu * m_p)
    term1 = rho * const1

    const2 = 4 * a / 3
    term2 = const2 * (T ** 3)

    return term1 + term2


# =============================================================================
# kappa definitions
# =============================================================================
kappaH = []
kappaes = []
kappaf =[]

# Kappa H-
def kappaHminus(rho, T):
    rho3 = rho / 1e3
    kappa = 2.5e-32 * (Z / 0.02) * (rho3 ** 0.5) * (T ** 9)

    return kappa


# Kappa ff
def Kappaff(rho, T):
    rho3 = rho / 1e3
    kappa = 1e24 * (Z + 0.0001) * (rho3 ** 0.7) * (T ** -3.5)
    
    return kappa


# kappa es
def Kappaes():
    kappa = 0.02 * (1 + X)
    
    return kappa


# the real slim kappa
def Kappa(rho, T):
    term1 = 1 / kappaHminus(rho, T)
    if T > 1e4:
        term2 = 1 / max(Kappaff(rho, T), Kappaes())  # max fn take itterale and returns largest, so i put the kappas in a list
    else:
        term2 = 1 / min(Kappaff(rho, T), Kappaes())
    

    return 1 / (term1 + term2)


# =============================================================================
# Energy generation epsilon
# =============================================================================

def EpsilonPP(rho, T):
    rho5 = rho / 1e5
    T6 = T / 1e6
    return 1.07e-7 * rho5 * (X ** 2) * T6 ** 4


def EpsilonCNO(rho, T):
    rho5 = rho / 1e5
    T6 = T / 1e6
    Xcno = 0.03 * X
    return 8.24e-26 * rho5 * X * Xcno * (T6 ** 19.9)


def Epsilon(rho, T):
    return EpsilonPP(rho, T) + EpsilonCNO(rho, T)


# =============================================================================
# The Big 5 Equations
# =============================================================================
# 1) Partial derivative of rho w.r.t radius

def DRhoDr(r, rho, M, kappa, DelPDelrho, DelPressuredelT, DTempDr):
    term1 = -1 * (G * M * rho) / (r ** 2)

    term2 = DelPressuredelT * DTempDr

    return (term1 + term2) / DelPDelrho


# 2) Derivitave of Temperature w.r.t radius
def DTDr(r, rho, T, L, M, kappa, P):
    const1 = 3 / (16 * pi * a * c)
    term1 = (const1 * rho * kappa * L) / ((T ** 3) * (r ** 2))

    const2 = (1 - (1 / gamma)) * G
    term2 = const2 * (T * M * rho) / (P * (r ** 2))

    return -1 * min([term1, term2])


# 3) Derivative of Mass w.r.t radius
def DMDr(r, rho):
    return 4 * pi * (r ** 2) * rho


# 4) Derivative of Luminosity w.r.t radius
def DLDr(r, rho, epsilon):
    return 4 * pi * (r ** 2) * rho * epsilon


# 5) Derivative of optical depth w.r.t radius
def DTauDr(kappa, rho):
    return kappa * rho


# for a  function x(t)
# we are computing f(t,x)
# Where t and x are the function inputs


# =============================================================================
# Runge Kutta
# =============================================================================
# For runge kutta Solving for K values


# Runge-Kutta-K 2's-RVal
RKK2R = 1 / 4
# ...
RKK3R = 3 / 8
RKK4R = 12 / 13
RKK5R = 1
RKK6R = 1 / 2

# Runge-Kukka-K 2's- K1val
RK21 = 1 / 4

# Runge-Kukka-K 3's- K1val
RK31 = 3 / 32
# Runge-Kukka-K 3's- K2val
RK32 = 9 / 32

RK41 = 1932 / 2197
RK42 = -7200 / 2197
RK43 = 7296 / 2197

RK51 = 439 / 216
RK52 = -8
RK53 = 3680 / 513
RK54 = -845 / 4104

RK61 = -8 / 27
RK62 = 2
RK63 = -3544 / 2565
RK64 = 1859 / 4104
RK65 = -11 / 40

# constants for summing all x together

# Runge-Kukka-4thorder-X- K1val
RK4X1 = 25 / 216
RK4X3 = 1408 / 2565
RK4X4 = 2197 / 4104
RK4X5 = -1 / 5

RK5X1 = 16 / 135
RK5X3 = 6656 / 12825
RK5X4 = 28561 / 56430
RK5X5 = -9 / 50
RK5X6 = 2 / 55

    # list of K, has indices corresponding to what k i am using
    # so if is solving k2, k1 is known and index 0 = 0.25*K1, all other indices=0
    # if solving K3, K1 and K2 knows so listofK[0] = 3/32 * K1 and listofK[1] = 9/32* K2, all others 0


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

    derivitaves = np.array([dRho, dT, dM, dL, dTau])
    return derivitaves


def CheckToStop(rho, T, M, drho):
    kappa = Kappa(rho, T)
    deltaTau = kappa * rho ** 2 / abs(drho)
    # print(deltaTau)
    if deltaTau < 0.01:  # @@@@@@@@this is just a placeholder!!!!@@@@@@@
        return False, deltaTau

    if M > 1e3 * M_s:
        return False, deltaTau

    return True, deltaTau


def extendArray(arrayToExtend,dimensions):
    extender=np.zeros(dimensions)
    arrayToExtend=np.column_stack((arrayToExtend,extender))
    return arrayToExtend



def RK45(r_0, h_0, vals, minError, maxError, hmin, hmax, abort):
    # merge the functions
    i = 0  # itter counter
    dL=[0]
    r = [r_0]  # initial r   next r is r+h
    hlist = [h_0]
    errorList = []
    keepGoing = True
    zero = np.array([0, 0, 0, 0, 0])
    vals = extendArray(vals,(5,10000))
    kappa = [Kappa(vals[0,0],vals[1, 0])]
    
    while keepGoing and i < abort:
        
        h = hlist[i]

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

        RK4X = vals[:, i] + RK4X1 * K1 + RK4X3 * K3 + RK4X4 * K4 + RK4X5 * K5

        RK5X = vals[:, i] + RK5X1 * K1 + RK5X3 * K3 + RK5X4 * K4 + RK5X5 * K5 + RK5X6 * K6

        error = abs(RK4X - RK5X)
        maxerror = np.max(error)

        #        except ValueError:
        #            # Finding negative values causing errors
        #            hlist[i] = h/2
        #            continue

        if maxerror > maxError and h > hmin:
            hlist[i] = h / 2
            continue

        else:
            vals[0, i + 1] = RK4X[0]  # idx0: rho
            vals[1, i + 1] = RK4X[1]  # idx1: T
            vals[2, i + 1] = RK4X[2]  # idx2: M
            vals[3, i + 1] = RK4X[3]  # idx3: L
            vals[4, i + 1] = RK4X[4]  # idx4: tau
            #print(vals)

            kappa.append(Kappa(RK4X[0],RK4X[1]))
            drhodr = derivitaves[0]
            keepGoing, deltaTau = CheckToStop(RK4X[0], RK4X[1], RK4X[2], drhodr)
            
            dL.append(derivitaves[3])
            
            r.append(r[i] + h)

            #kappaH.append(kappaHminus(rho, T))
            #kappaes.append(Kappaea())
            #kappaf.append(kappaff(rho, T))

            i += 1

            errorList.append(error)
            
            if i % 10000 == 0:
                if i==abort:
                    pass
                else:
                    try:
                        vals=extendArray(vals,(5,10000))
                    except:
                        print('Error: Array could not be extended')
                        
                try:
                    print('Keep Going Criteria: ' + str(deltaTau))
                    print('Current Step (x1e4): ' + str(i/10000))
                    print('-------------------------------------')
                except:
                    pass
            if maxerror < minError:
                hlist.append(2 * h)

            else:
                hlist.append(h)

    return vals, r, hlist, errorList, dL, kappa


# we  need to define our initial parameters, 
rs = 0.865*R_s 
ms = 0.673*M_s
ls = 5.86e-2*L_s  
ts = 3056

r_0 = 1
rho_c = 585600 # Correct rho for sun
T_c = 8.23e6  # idx1: T
M_0 = 4 * np.pi * r_0 ** 3 * rho_c / 3  # idx2: M
L_0 = 4 * np.pi * r_0 ** 3 * rho_c * Epsilon(rho_c, T_c) / 3  # idx3: L
tau_0 = 0.231739  # possible value 1 from OPAL TABLE

# initialiize the vals array
values = np.zeros((5, 1))   
values[0] = rho_c           # idx0: rho
values[1] = T_c             # idx1: T
values[2] = M_0             # idx2: M
values[3] = L_0             # idx3: L
values[4] = tau_0           # idx4: tau


vals, r, stepsize, errors, dL, kappa = RK45(1, 100, values, 0.02, 0.5, 100000, 500000, 1e5)
#@@@@@@@ Abort variable must be greater then and divisible by 10000 @@@@@@

vals = vals[:,:len(r)]
#errors.append(errors[-1])

r = np.array(r)
print(r)
R =r/np.max(r)

fig, axs = plt.subplots(3, 3, sharex=True, frameon=False, figsize=(10, 5))

plt.xlim(0, 1)

axs[0, 0].plot(R, vals[0, :]/rho_c, 'k-')
axs[0, 0].set_title('Density (rho/rho_c)', fontsize=8)

axs[0, 1].plot(R, vals[1, :]/T_c, 'k-')
axs[0, 1].set_title('Temperature (T/T_c)', fontsize=8)

axs[0, 2].plot(R, vals[2, :]/ms, 'k-')
axs[0, 2].set_title('Mass (M/M_0)', fontsize=8)

axs[1, 0].plot(R, vals[3, :]/ls, 'k-')
axs[1, 0].set_title('Luminosity (L/L_0)', fontsize=8)

axs[1, 1].plot(R, vals[4, :], 'k-')
axs[1, 1].set_title('Tau', fontsize=8)

axs[1, 2].plot(R, dL, 'k-')
axs[1, 2].set_title('Derivitave of Luminosity', fontsize=8)

axs[2, 0].plot(R, kappa, 'k-')
axs[2, 0].set_title('Kappa', fontsize=8)

axs[2, 1].set_yscale('log')
axs[2, 1].plot(R, stepsize, 'k-')
axs[2, 1].set_title('Stepsize', fontsize=8)

#axs[2, 2].plot(R, kappaH, 'k-')
#axs[2, 2].plot(R, kappaes, 'b-')
#axs[2, 2].plot(R, kappaf, 'r')
#axs[2, 2].set_title('Kappa Comps', fontsize=8)

#plt.xscale('log')

plt.rcParams.update({'font.size': 8})
plt.tight_layout()
plt.show()


print("--- %s Minutes ---" % str(int(time.time() - start_time)/60))
