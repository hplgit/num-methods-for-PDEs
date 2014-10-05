import numpy as np
import matplotlib.pyplot as plt
import scitools.std as plt

def minmax(t, u):
    """
    Compute all local minima and maxima of the function u(t),
    represented by discrete points in the arrays u and t.
    Return lists minima and maxima of (t[i],u[i]) extreme points.
    """
    minima = []; maxima = []
    for n in range(1, len(u)-1, 1):
        if u[n-1] > u[n] < u[n+1]:
            minima.append((t[n], u[n]))
        if u[n-1] < u[n] > u[n+1]:
            maxima.append((t[n], u[n]))
    return minima, maxima

def periods(extrema):
    """
    Given a list of (t,u) points of the maxima or minima,
    return an array of the corresponding local periods.
    """
    p = [extrema[n][0] - extrema[n-1][0]
         for n in range(1, len(extrema))]
    return np.array(p)

def amplitudes(minima, maxima):
    """
    Given a list of (t,u) points of the minima and maxima of
    u, return an array of the corresponding local amplitudes.
    """
    # Compare first maxima with first minima and so on
    a = [(abs(maxima[n][1] - minima[n][1]))/2.0
         for n in range(min(len(minima),len(maxima)))]
    return np.array(a)

import nose.tools as nt

def test_empirical_analysis():
    t = np.linspace(0, 6*np.pi, 1181)
    # Modulated amplitude and period
    u = np.exp(-(t-3*np.pi)**2/12.0)*np.cos(np.pi*(t + 0.6*np.sin(0.25*np.pi*t)))
    plt.plot(t, u, label='signal')
    plt.hold('on')
    minima, maxima = minmax(t, u)
    t_min = [ti for ti, ui in minima]
    t_max = [ti for ti, ui in maxima]
    u_min = [ui for ui, ui in minima]
    u_max = [ui for ui, ui in maxima]
    plt.plot(t_min, u_min, 'bo', label='minima')
    plt.plot(t_max, u_max, 'ro', label='maxima')
    plt.legend()

    plt.figure()
    p = periods(maxima)
    a = amplitudes(minima, maxima)
    plt.plot(range(len(p)), p, 'g--', label='periods')
    plt.hold('on')
    plt.plot(range(len(a)), a, 'y-', label='amplitudes')
    plt.legend()

    p_ref = np.array([
        1.48560059,  2.73158819,  2.30028479,  1.42170379,  1.45365219,
        2.39612999,  2.63574299,  1.45365219,  1.42170379])
    a_ref = np.array([
        0.00123696,  0.01207413,  0.19769443,  0.59800044,  0.90044961,
        0.96007725,  0.42076411,  0.08626735,  0.0203696 ,  0.00312785])
    p_diff = np.abs(p - p_ref).max()
    a_diff = np.abs(a - a_ref).max()
    nt.assert_almost_equal(p_diff, 0, places=7)
    nt.assert_almost_equal(a_diff, 0, places=7)

if __name__ == '__main__':
    test_empirical_analysis()
    plt.show()

