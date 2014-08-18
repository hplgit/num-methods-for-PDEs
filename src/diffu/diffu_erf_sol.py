"""
Demonstrate exact analytical solution for diffusion of a step
into a straight line (and further on into u=1/2).
"""
from scipy.special import erfc
from numpy import sqrt, linspace
import scitools.std as plt
import time

def u(x, t):
    eta = (x - c)/sqrt(4*a*t)
    return 0.5*erfc(eta)

c = 0.5
a = 1

x = linspace(0, 1, 1001)
T = 5*a*1.6E-2
t_values = linspace(0, T, 1001)[1:]  # skip t=0

axis = [0, 1, -0.1, 1.1]
plt.plot([0, 0.5, 0.5, 1], [1, 1, 0, 0], 'r-',
         axis=axis, title='t=0')
time.sleep(1)
for t in t_values:
    y = u(x, t)
    plt.plot(x, y, 'r-',
             axis=axis, title='t=%f' % t)

