from diffusion1D_vc import *

# Test problem: start with u=u_L in left part and u=u_R in right part,
# let diffusion work and make a linear function from u_L to u_R as
# time goes to infinity. For a=1, u=1-x when u_L=1, u_R=0.

def I(x):
    return u_L if x <= L/2. else u_R

theta = 1
L = 1
Nx = 20
#Nx = 400
a = zeros(Nx+1) + 1
u_L = 1
u_R = 0
dx = L/float(Nx)
D = 500
dt = dx**2*D
dt = 1.25
D = dt/dx**2
T = 2.5
umin = u_R
umax = u_L

a_consts = [[0, 1]]
a_consts = [[0, 1], [0.5, 8]]
a_consts = [[0, 1], [0.5, 8], [0.75, 0.1]]
a = fill_a(a_consts, L, Nx)
#a = random.uniform(0, 10, Nx+1)

from scitools.std import plot, hold, subplot, figure, show

figure()
subplot(2,1,1)
u, x, cpu = viz(I, a, L, Nx, D, T, umin, umax, theta, u_L, u_R)

v = u_exact_stationary(x, a, u_L, u_R)
print 'v', v
print 'u', u
hold('on')
symbol = 'bo' if Nx < 32 else 'b-'
plot(x, v, symbol, legend='exact stationary')

subplot(2,1,2)
plot(x, a, legend='a')
show()
