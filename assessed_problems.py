from __future__ import division

import numpy as np
import pylab as plt

Omega = 7.292e-5
a = 6.37e6
Rd = 287
cpd = 1004
L = 2.5e6
Rv = 461.5

def f(phi):
    return 2 * Omega * np.sin(np.pi /180 * phi)


def B(phi):
    return 2 * Omega * np.cos(np.pi / 180 * phi)/a

def f_B(phi, phi0=0):
    f0 = f(phi0)
    B0 = B(phi0)
    return f0 + B0 * np.pi / 180 * phi * a

def PE(phi):
    return np.abs((f(phi) - f_B(phi))/f(phi)) * 100

def Q1():
    plt.clf()
    test_lats = [-5, 5, 10]

    for test_lat in test_lats:
        print('{}deg, f={}'.format(test_lat, f(test_lat)))
    for test_lat in test_lats:
        print('{}deg, PE={}'.format(test_lat, PE(test_lat)))

    phi = np.linspace(0, 90, 9001)
    e = PE(phi)
    e[0] = 0
    plt.plot(phi, e, label='$\\frac{f - f_\\beta}{f}$')
    plt.plot(phi, np.ones_like(phi) * 5, 'r', label='critical % error')
    plt.xlabel('$phi$')
    plt.ylabel('% error')
    plt.legend(loc='upper left')
    for ev, phiv in zip(e, phi):
        if ev > 5:
            break
    print('phi where PE>5%: {}'.format(phiv))
    plt.show()

def tetens(T):
    T = T - 273.15
    return 6.112 * np.exp(17.67 * T / (T + 243.5))

def theta_e(theta_d, rvs, T):
    return theta_d * np.exp(L * rvs / (cpd * T) ) * (1 - tetens(T) / 1000)**(-Rd / cpd)

def theta_e_approx(theta_d, rvs, T):
    return theta_d * np.exp(L * rvs / (cpd * T))

def theta_e_approx(theta_d, rvs, T):
    return theta_d * np.exp(L * rvs / (cpd * T))

def theta_e_naive(rvs, T):
    return T + L * rvs / cpd

def Q2():
    dt = 2.5
    T = 300

    percentage_change = L /(Rv * T**2) * dt * 100
    print('% change: {}'.format(percentage_change))
    pc2 = (tetens(T + dt) - tetens(T)) / tetens(T) * 100
    print('% change 2: {}'.format(pc2))

    # See Ambaum p114 and https://en.wikipedia.org/wiki/Equivalent_potential_temperature
    # for formulae.
    for T, rvs in [(273.15 + 25, 0.017),
                   (273.15 + 40, 0.007)]:
        print('T={}, rvs={}'.format(T, rvs))
        print('theta_e={}'.format(theta_e(T, rvs, T)))
        print('theta_e_approx={}'.format(theta_e_approx(T, rvs, T)))
        print('theta_e_naive={}'.format(theta_e_naive(rvs, T)))
        print('')


if __name__ == '__main__':
    plt.ion()
    #Q1()
    Q2()
