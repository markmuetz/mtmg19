from __future__ import division

import numpy as np
import pylab as plt

Omega = 7.292e-5
a = 6.37e6
Rd = 287
cpd = 1004
L = 2.5e6
Rv = 461.5
g = 9.81

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
    print('Q1')
    plt.clf()
    test_lats = [-5, 5, 10]

    print('(a)')
    for test_lat in test_lats:
        print('{}deg, f={}'.format(test_lat, f(test_lat)))
    print('(b)')
    for test_lat in test_lats:
        print('{}deg, PE={}'.format(test_lat, PE(test_lat)))

    print('(c)')
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
    print('')
    print('(d)')
    # e-foldinging dist:
    c = 80 # m/s, taken from Q3.
    ye = np.sqrt(2*c/B(0))
    print('ye={}'.format(ye))
    print('ye={} deg'.format(ye/a * 180 / np.pi))

def tetens(T):
    T = T - 273.15
    return 6.112 * np.exp(17.67 * T / (T + 243.5))

def theta_e(rvs, T):
    return T + L * rvs / cpd

def Q2():
    print('Q2')
    print('(a)')
    dt = 2.5
    T = 300

    percentage_change = L /(Rv * T**2) * dt * 100
    print('% change: {}'.format(percentage_change))
    pc2 = (tetens(T + dt) - tetens(T)) / tetens(T) * 100
    print('% change 2: {}'.format(pc2))

    print('(b)')
    # See https://en.wikipedia.org/wiki/Equivalent_potential_temperature
    # for formulae.
    for T, rvs in [(273.15 + 25, 0.017),
                   (273.15 + 40, 0.007)]:
        print('T={}, rvs={}'.format(T, rvs))
        print('theta_e={}'.format(theta_e(rvs, T)))
        print('')


def Q3():
    print('Q3')
    print('(a)/(b)')
    c = 80 # m/s
    u1 = 10 # m/s

    lat = 12 * np.pi / 180
    y = lat * a
    u2 = u1 * np.exp(-B(lat) * y**2 / (2 * c))
     
    for u in [u1, u2]:
        phi = u*c
        print('phi={} m2s-2'.format(phi))
        Z = phi/g
        print('Z={} m'.format(Z))
        rho = 1.22
        dp = -rho*g*2*Z
        print('Delta p={} Pa, {} hPa'.format(dp, dp/100))
        print('')


def Q4():
    print('Q4')
    print('(b)')
    for lat in [5, 10, 25]:
        print('lat={}'.format(lat))
        y = lat * np.pi / 180 * a
        u = Omega * y**2 / a
        print('u={}'.format(u))
        print('')
    print('(c)')
    yc = 12 * np.pi / 180 * a
    y = 0
    u = Omega * (y**2 - yc**2) / a
    print('u={}'.format(u))
    print('')



if __name__ == '__main__':
    plt.ion()
    Q1()
    Q2()
    Q3()
    Q4()
