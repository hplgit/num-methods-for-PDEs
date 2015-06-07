from scitools.easyviz.matplotlib_ import *
from numpy import *

def Gaussian(x, x0, s):
    return 1/(sqrt(2*pi)*s)*exp(-(x-x0)**2/(2*s**2))

x = linspace(0, 1, 10001)
s = 0.01
for x0 in (0.2, 0.4, 0.6, 0.8):
    y = Gaussian(x, x0, s)
    max_y = y.max()
    plot(x, y, '-')
    hold('on')
xlabel('x')
ylabel('w')
axis([x[0], x[-1], 0, 1.2*max_y])
savefig('tmp.pdf')
savefig('tmp.png')
show()
input()
