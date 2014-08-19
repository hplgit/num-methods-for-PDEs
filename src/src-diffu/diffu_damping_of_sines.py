from matplotlib.pyplot import *
from numpy import *

def component(Q, x, t, k, a=1):
    return Q*exp(-a*k**2*t)*sin(k*x)

def u(x, t):
    return component(1, x, t, pi, 1) + component(0.1, x, t, 100*pi, 1)

x = linspace(0, 1, 2001)
a = 1
amplitudes = array([0.1, 0.01])
k = 100*pi
times1 = log(amplitudes)/(-a*k**2)
k = pi
times2 = log(amplitudes)/(-a*k**2)
times = [0] + times1.tolist() + times2.tolist()

for t in times:
    figure()
    plot(x, u(x,t))
    title('t=%.2E' % t)
    xlabel('x')
    ylabel('u')
    axis([0, 1, -0.1, 1.1])
    savefig('tmp_%.2E.pdf' % t)
    savefig('tmp_%.2E.png' % t)

import os
times = times[:1] + times[2:]
os.system('doconce combine_images tmp_%.2E.pdf tmp_%.2E.pdf tmp_%.2E.pdf tmp_%.2E.pdf diffusion_damping.pdf' % tuple(times))
os.system('doconce combine_images tmp_%.2E.png tmp_%.2E.png tmp_%.2E.png tmp_%.2E.png diffusion_damping.png' % tuple(times))
show()



