from diffusion1D_vc import *

# Test problem: start with u=u_L in left part and u=u_R in right part,
# let diffusion work and make a linear function from u_L to u_R as
# time goes to infinity. For a=1, u=1-x when u_L=1, u_R=0.

theta = 1
L = 1
Nx = 41
a = zeros(Nx+1) + 1
u_L = 1
u_R = 0
D = 50
T = 0.4    # "infinite time"
umin = u_R
umax = u_L

def I(x):
    return u_L if x <= L/2. else u_R

def test2(p):
    """
    Given the values ``p=(a1,a2)`` of the diffusion coefficient
    in two subdomains [0, 0.5] and [0.5, 1], return u(0.5, inf).
    """
    assert len(p) == 2
    a1, a2 = p
    a_consts = [[0, a1], [0.5, a2]]
    a = fill_a(a_consts, L, Nx)

    u, x, t, cpu = solver_theta(I, a, L, Nx, D, T,
                                theta=theta, u_L=u_L, u_R=u_R)
    return u[Nx/2]

def test2_fast(p):
    """Fast version of test2 using the analytical solution directly."""
    a1, a2 = p
    a_consts = [[0, a1], [0.5, a2]]
    a = fill_a(a_consts, L, Nx)
    x = linspace(0, L, Nx+1)
    v = u_exact_stationary(x, a, u_L, u_R)

    return v[Nx/2]


def visualize(p):
    a1, a2 = p
    a_consts = [[0, a1], [0.5, a2]]
    a = fill_a(a_consts, L, Nx)
    # Choose smaller time step to better see the evolution
    global D
    D = 50.

    u, x, cpu = viz(I, a, L, Nx, D, T, umin, umax, theta, u_L, u_R)

if __name__ == '__main__':
    #visualize((8, 1))
    print test2((8, 1))
