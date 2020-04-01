import numpy as np
import matplotlib.pyplot as plt
import time
start_time = time.time()


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
a = 4 * sigma / (3*c)
pi = np.pi
ep_pp_coeff = 1.07e-7 * 1e-5 * (1e-6) ** 4
ep_cno_coeff = 8.24e-26 * 1e-5 * (1e-6) ** 19.9
nonrelgenpress = (3 * pi ** 2) ** (2 / 3) / 5 * hbar ** 2 / m_e * m_p ** (-5 / 3)
mach_ep = np.finfo(np.float64).eps
tiny_float = 1e-20
gamma = 5 / 3  # ideal gas constant

# =================Defining mu, XYZ=====================================
X = 0.70
Y = 0.265
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

    const2 =  a
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
    if T > 1.0e4:
        term2 = 1.0 / max([Kappaff(rho, T), Kappaes()])  # max fn take itterale and returns largest, so i put the kappas in a list
    else:
        term2 = 1.0 / max([Kappaff(rho, T), Kappaes()]) # what andrew did, changed to max
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
    #print(kappa)
    const2 = (1.0 - (1.0 / gamma)) * G * mu
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

# TestVals = [ 1.338e+33,    6.01533e+10,    8.14146e-06,        3161.18,      1.338e+33,    2.24985e+32,     0.00157801,     0.00157801,    6.44724e-63,            2.5,       0.999976,    1.39037e-07,    0.000383851,              1,              1,      0.0224871,      0.0224871,    5.26823e+07,           0.34,     3.4781e+06,        32780.6,    3.44531e+06,        0.25184,   -2.82018e-13,   -7.30015e-05,    3.70195e+17,     0.00157801]
# TestVals2 = [1.338e+33,    6.01541e+10,    7.93209e-06,        3106.75,      1.338e+33,    2.24985e+32,     0.00139739,     0.00139739,    4.33158e-63,            2.5,       0.999988,    1.35462e-07,    0.000377241,              1,              1,       0.018984,       0.018984,    5.49725e+07,           0.34,     3.3303e+06,        31387.6,    3.29891e+06,       0.234937,   -2.79572e-13,   -7.29997e-05,    3.60684e+17,     0.00139739]
TestVals = [1.338e+33 ,    4.62509e+10,      0.0874953,    1.13347e+06,    1.32636e+33,    2.24985e+32,     1.7809e+15,     1.7809e+15,    0.000301589 ,        4.5815,       0.768866,     0.00149422,       0.137633,       0.991303,              1,        40.0499,    2.28359e+23,        40.0499,           0.34,    1.34519e+13,    1.71562e+11,    1.32762e+13,    4.16261e+09,   -1.83096e-11,   -6.65503e-05,    2.35199e+21,     1.7809e+15]
TestVals2 = [1.338e+33,    4.63015e+10,      0.0865743,    1.13011e+06,    1.32648e+33,    2.24985e+32,    1.72677e+15,    1.72677e+15,    0.000278928,        4.57609,       0.769706,     0.00147849,       0.137225,       0.991392,              1,        40.1699,     2.2116e+23,        40.1699 ,          0.34,    1.32701e+13,    1.68562e+11,    1.30974e+13,    4.11343e+09,   -1.81281e-11,    -6.6493e-05,    2.33232e+21,    1.72677e+15] 

#               M		         r	             rho(r)	        T(r)	       M(r)	            L(r)	       dL/dr             dLpp/dr       dLcno/dr	        dlogP/dlogT	        r/R	  	     rho(r)/rhoc        T(r)/Tc         M(r)/M         L(r)/L    	   kappa        kappaHm        kappaff            kappaes          P           Pdeg           Pgas               Ppho        drho/dr          dT/dr          dM/dr          dL/dr
#               0                1                 2            3               4                 5              6               7               8                  9                10          11                  12              13             14              15              16              17              18             19           20              21                22          23              24              25              26                                                   

print(np.shape(TestVals))
Delatr = (TestVals[3]-TestVals2[3])

deltacalc = ((TestVals[19]-TestVals2[19])/10)/Delatr

# deltaim = TestVals[19]/10

print(TestVals[15]/10)
print(Kappa(TestVals[2]*1000,TestVals[3]))

# deltatest = DRhoDr(TestVals[1]/100,TestVals[2]*1000,TestVals[4]/1000,TestVals[15]/10,DelPressureDelrho(TestVals[2]*1000,TestVals[3]),DelPressuredelT(TestVals[2]*1000,TestVals[3]),TestVals[24]*100)
# deltatest = DLDr(TestVals[1]/100,TestVals[2]*1000,Epsilon(TestVals[2]*1000,TestVals[3]))
# deltatest = DMDr(TestVals[1]/100,TestVals[2]*1000)
# deltatest = DTDr(TestVals[1]/100,TestVals[2]*1000,TestVals[3],TestVals[5]*1e7,TestVals[4]/1000,TestVals[15]/10,TestVals[19]/10)
deltatest = DelPressuredelT(TestVals[2]*1000, TestVals[3])
# deltatest = Pressure(TestVals[2]*1000,TestVals[3])
# deltatest = DelPressureDelrho (TestVals[2]*1000,TestVals[3])

# print ("Sheet " + str(deltaim))
print ("Calc " + str(deltacalc))
print ("Test " + str(deltatest))

