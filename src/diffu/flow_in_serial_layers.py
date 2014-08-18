from Heaviside import PiecewiseConstant, IntegratedPiecewiseConstant
import numpy as np
from sets import Set

class SerialLayers:
    """
    b: coordinates of boundaries of layers, b[0] is left boundary
    and b[-1] is right boundary of the domain [0,L].
    a: values of the functions in each layer (len(a) = len(b)-1).
    U_0: u(x) value at left boundary x=0=b[0].
    U_L: u(x) value at right boundary x=L=b[0].
    """

    def __init__(self, a, b, U_0, U_L, eps=0):
        self.a, self.b = np.asarray(a), np.asarray(b)
        assert len(a) == len(b)-1, 'a and b do not have compatible lengths'
        self.eps = eps  # smoothing parameter for smoothed a
        self.U_0, self.U_L = U_0, U_L

        # Smoothing parameter must be less than half the smallest material
        if eps > 0:
            assert eps < 0.5*(self.b[1:] - self.b[:-1]).min(), 'too large eps'

        a_data = [[bi, ai] for bi, ai in zip(self.b, self.a)]
        domain = [b[0], b[-1]]
        self.a_func = PiecewiseConstant(domain, a_data, eps)

        # inv_a = 1/a is needed in formulas
        inv_a_data = [[bi, 1./ai] for bi, ai in zip(self.b, self.a)]
        self.inv_a_func = PiecewiseConstant(domain, inv_a_data, eps)
        self.integral_of_inv_a_func = \
             IntegratedPiecewiseConstant(domain, inv_a_data, eps)
        # Denominator in the exact formula is constant
        self.inv_a_0L = self.integral_of_inv_a_func(b[-1])

    def u_exact(self, x):
        solution = self.U_0 + (self.U_L-self.U_0)*\
                   self.integral_of_inv_a_func(x)/self.inv_a_0L
        return solution

    __call__ = u_exact

    def plot(self):
        x, y_a = self.a_func.plot()
        x = np.asarray(x); y_a = np.asarray(y_a)
        y_u = self.u_exact(x)
        import matplotlib.pyplot as plt
        plt.figure()
        plt.plot(x, y_u, 'b')
        plt.hold('on')  # Matlab style
        plt.plot(x, y_a, 'r')
        ymin = -0.1
        ymax = 1.2*max(y_u.max(), y_a.max())
        plt.axis([x[0], x[-1], ymin, ymax])
        plt.legend(['solution $u$', 'coefficient $a$'], loc='upper left')
        if self.eps > 0:
            plt.title('Smoothing eps: %s' % self.eps)
        plt.savefig('tmp.pdf')
        plt.savefig('tmp.png')
        plt.show()

def simulate_case1(
    mesh_resolution_factor=1.,
    uncertainty_type=Set(['material parameters',
                          'internal boundary 1',
                          'u at 0',
                          'u at L']),
    plot=True):  # include plot in simulator (ok when playing)
    """
    Return simulator that takes just the uncertain parameters
    as "unknowns" (i.e., works in parameter space), but uses
    other prescribed values for all remaining input.
    """

    def quantity_of_interest(solution, mesh):
        """Return gradient in last material."""
        return (solution[-1] - solution[-2])/(mesh[-1] - mesh[-2])

    def simulator(*args):
        """
        Function to be called from UQ libraries.
        `args` holds the uncertain parameters only.
        """
        # Note: args has different storing conventions for different
        # types of uncertainty parameters, see below.

        # *Initial all data needed for the simulation*

        # Should use a numerical mesh that has the boundaries of the
        # layers as mesh points. Therefore, the mesh is made first.
        L = 5
        n = 10 # no of mesh intervals with mesh_resolution_factor=1
        i = mesh_resolution_factor*n    # total no of intervals
        mesh = np.linspace(0, L, i+1)
        # Material boundaries matches mesh (as long as n=10)
        b = [0, 0.25, 0.5, 1]
        # Material values
        a = [0.2, 0.4, 4]
        # Boundary conditions
        U_0 = 0.5
        U_L = 5 # should approx match a for a nice plot...

        # Override with data from *args

        # Use Set to avoid a particular sequence of the elements
        # in the uncertainty_type set
        if uncertainty_type == Set(['material parameters']):
            a = np.asarray(*args)
        elif uncertainty_type == Set(['internal boundary 1']):
            b[1] = args[0]
        elif uncertainty_type == Set(['material parameters',
                                      'internal boundaries']):
            a = np.asarray(*args)[:-1]
            b[1] = args[-1]
        elif uncertainty_type == Set(['internal boundary 1', 'u at L']):
            b[1] = args[0]
            U_L = args[1]
        # else: no redefinition of data because of args

        # Create exact solution
        #eps=0.05
        eps = 0
        u = SerialLayers(a, b, U_0, U_L, eps=eps)
        if plot:
            u.plot()
        solution = u(mesh)
        q = quantity_of_interest(solution, mesh)
        return q
    return simulator

simulator = simulate_case1(
    mesh_resolution_factor=1.,
    uncertainty_type = Set(['internal boundary 1']),
    plot=True)

# Play around with moving b[1] and see the effect :-)
b1 = 0.25
while 0 < b1 < 0.5:
    b1 = float(raw_input('location of boundary 1 in (0, 0.5): '))
    q = simulator(b1)
    print q

