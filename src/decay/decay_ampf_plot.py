"""Plot of amplification factors for the theta scheme."""

from numpy import linspace, exp
#from matplotlib.pyplot import *
from scitools.std import *

def A_exact(p):
    return exp(-p)

def A(p, theta):
    return (1-(1-theta)*p)/(1+theta*p)

def amplification_factor(names):
    curves = {}
    p = linspace(0, 3, 15)
    curves['exact'] = A_exact(p)
    plot(p, curves['exact'])
    hold('on')
    name2theta = dict(FE=0, BE=1, CN=0.5)
    for name in names:
        curves[name] = A(p, name2theta[name])
        plot(p, curves[name])
    plot([p[0], p[-1]], [0, 0], '--')  # A=0 line
    title('Amplification factors')
    grid('on')
    legend(['exact'] + names, loc='lower left', fancybox=True)
    xlabel('a*dt')
    ylabel('A')
    savefig('A_factors.png')
    savefig('A_factors.pdf')
    savefig('A_factors.eps')
    show()

if __name__ == '__main__':
    if len(sys.argv) > 1:
        names = sys.argv[1:]
    else: # default
        names = 'FE BE CN'.split()
    amplification_factor(names)
