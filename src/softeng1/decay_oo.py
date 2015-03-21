
class Parameters:
    def __init__(self):
        """
        Subclasses must initialize self.prms with
        parameters and default values, self.types with
        the corresponding types, and self.help with
        the corresponding descriptions of parameters.
        self.types and self.help are optional, but
        self.prms must be complete and contain all parameters.
        """
        pass

    def ok(self):
        """Check if attr. prms, types, and help are defined."""
        if hasattr(self, 'prms') and \
           isinstance(self.prms, dict) and \
           hasattr(self, 'types') and \
           isinstance(self.types, dict) and \
           hasattr(self, 'help') and \
           isinstance(self.help, dict):
            return True
        else:
            raise ValueError(
                'The constructor in class %s does not '\
                'initialize the dictionaries\n'\
                'self.prms, self.types, self.help.' %
                self.__class__.__name__)

    def _illegal_parameter(self, name):
        """Raise exception about illegal parameter name."""
        raise ValueError(
            'parameter "%s" is not registered.\nLegal '\
            'parameters are\n%s' %
            (name, ' '.join(list(self.prms.keys()))))

    def set(self, **parameters):
        for name in parameters:
            if name in self.prms:
                self.prms[name] = parameters[name]
            else:
                self._illegal_parameter(name)

    def get(self, name):
        if name in self.prms:
            return self.prms[name]
        else:
            self._illegal_parameter(name)

    def define_command_line_options(self, parser=None):
        self.ok()
        if parser is None:
            import argparse
            parser = argparse.ArgumentParser()

        for name in self.prms:
            tp = self.types[name] if name in self.types else str
            help = self.help[name] if name in self.help else None
            parser.add_argument(
                '--' + name, default=self.get(name), metavar=name,
                type=tp, help=help)

        return parser

    def init_from_command_line(self, args):
        for name in self.prms:
            self.prms[name] = getattr(args, name)


class Problem(Parameters):
    """
    Physical parameters for the problem u'=-a*u, u(0)=I,
    with t in [0,T].
    """
    def __init__(self):
        self.prm  = dict(I=1, a=1, T=10)
        self.type = dict(I=float, a=float, T=float)
        self.help = dict(I='initial condition, u(0)',
                         a='coefficient in ODE',
                         T='end time of simulation')

    def u_exact(self, t):
        I, a = self.get('I'), self.get('a')
        return I*np.exp(-a*t)

class Solver(Parameters):
    def __init__(self, problem):
        self.problem = problem
        self.prm  = dict(dt=0.5, theta=0.5)
        self.type = dict(dt=float, theta=float)
        self.help = dict(dt='time step value',
                         theta='time discretization parameter')

    def solve(self):
        from decay_mod import solver
        self.u, self.t = solver(
            self.problem.get('I'),
            self.problem.get('a'),
            self.problem.get('T'),
            self.get('dt'),
            self.get('theta'))

    def error(self):
        try:
            u_e = self.problem.u_exact(self.t)
            e = u_e - self.u
            E = np.sqrt(self.get('dt')*np.sum(e**2))
        except AttributeError:
            E = None
        return E

def experiment_classes():
    problem = Problem()
    solver = Solver(problem)

    # Read input from the command line
    parser = problem.define_command_line_options()
    parser = solver. define_command_line_options(parser)
    args = parser.parse_args()
    problem.init_from_command_line(args)
    solver. init_from_command_line(args)

    # Solve and plot
    solver.solve()
    import matplotlib.pyplot as plt
    t_e = np.linspace(0, T, 1001)    # very fine mesh for u_e
    u_e = problem.u_exact(t_e)

    plt.plot(t,   u,   'r--o')       # dashed red line with circles
    plt.plot(t_e, u_e, 'b-')         # blue line for u_e
    plt.legend(['numerical, theta=%g' % theta, 'exact'])
    plt.xlabel('t')
    plt.ylabel('u')
    plotfile = 'tmp'
    plt.savefig(plotfile + '.png');  plt.savefig(plotfile + '.pdf')

    error = problem.u_exact(t) - u
    E = np.sqrt(dt*np.sum(error**2))
    print 'Error norm:', E
    plt.show()

if __name__ == '__main__':
    experiment_classes()
