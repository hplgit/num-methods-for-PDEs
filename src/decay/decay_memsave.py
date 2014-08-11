import numpy as np
import matplotlib.pyplot as plt

def solver_memsave(I, a, T, dt, theta, filename='sol.dat'):
    """
    Solve u'=-a*u, u(0)=I, for t in (0,T] with steps of dt.
    Minimum use of memory. The solution is stored in a file
    (with name filename) for later plotting.
    """
    dt = float(dt)         # avoid integer division
    Nt = int(round(T/dt))  # no of intervals

    outfile = open(filename, 'w')
    # u: time level n+1, u_1: time level n
    t = 0
    u_1 = I
    outfile.write('%.16E  %.16E\n' % (t, u_1))
    for n in range(1, Nt+1):
        u = (1 - (1-theta)*a*dt)/(1 + theta*dt*a)*u_1
        u_1 = u
        t += dt
        outfile.write('%.16E  %.16E\n' % (t, u))
    outfile.close()
    return u, t

def read_file(filename='sol.dat'):
    infile = open(filename, 'r')
    u = [];  t = []
    for line in infile:
        words = line.split()
        if len(words) != 2:
            print 'Found more than two numbers on a line!', words
            sys.exit(1)  # abort
        t.append(float(words[0]))
        u.append(float(words[1]))
    return np.array(t), np.array(u)

def read_file_numpy(filename='sol.dat'):
    data = np.loadtxt(filename)
    t = data[:,0]
    u = data[:,1]
    return t, u

def exact_solution(t, I, a):
    return I*np.exp(-a*t)

def explore(I, a, T, dt, theta=0.5, makeplot=True):
    """
    Run a case with the solver_memsave, load t and u data
    from file into arrays, compute error measure,
    and plot the numerical and exact solutions (if makeplot=True).
    """
    filename = 'u.dat'
    u, t = solver_memsave(I, a, T, dt, theta, filename)

    t, u = read_file(filename)
    u_e = exact_solution(t, I, a)
    e = u_e - u
    E = np.sqrt(dt*np.sum(e**2))
    if makeplot:
        plt.figure()                     # create new plot
        t_e = np.linspace(0, T, 1001)    # very fine mesh for u_e
        u_e = exact_solution(t_e, I, a)
        plt.plot(t,   u,   'r--o')       # red dashes w/circles
        plt.plot(t_e, u_e, 'b-')         # blue line for u_e
        plt.legend(['numerical', 'exact'])
        plt.xlabel('t')
        plt.ylabel('u')
        plt.title('Method: theta-rule, theta=%g, dt=%g' %
                  (theta, dt))
        theta2name = {0: 'FE', 1: 'BE', 0.5: 'CN'}
        plt.savefig('%s_%g.png' % (theta2name[theta], dt))
        plt.show()
    return E

def verify():

    def exact_discrete_solution(n, I, a, theta, dt):
        A = (1 - (1-theta)*a*dt)/(1 + theta*dt*a)
        return I*A**n

    theta = 0.8; a = 2; I = 0.1; dt = 0.8
    Nt = int(8/dt)  # no of steps
    filename = 'sol.dat'
    u, t = solver_minmem(I=I, a=a, T=Nt*dt, dt=dt, theta=theta,
                             filename=filename)

    # Test both methods of reading the file
    t, u = read_file(filename)
    t_, u_ = read_file_numpy(filename)
    file_reading = (t == t_).all() and (u == u_).all()

    # Compare numerical values
    u_de = np.array([exact_discrete_solution(n, I, a, theta, dt)
                     for n in range(Nt+1)])
    difference = abs(u_de - u).max()
    tol = 1E-15  # tolerance for comparing floats
    success = difference <= tol
    return success and file_reading

def define_command_line_options():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--I', '--initial_condition', type=float,
                        default=1.0, help='initial condition, u(0)',
                        metavar='I')
    parser.add_argument('--a', type=float,
                        default=1.0, help='coefficient in ODE',
                        metavar='a')
    parser.add_argument('--T', '--stop_time', type=float,
                        default=1.0, help='end time of simulation',
                        metavar='T')
    parser.add_argument('--makeplot', action='store_true',
                        help='display plot or not')
    parser.add_argument('--theta', type=float,
                        default=0.5, help='parameter in scheme',
                        metavar='theta')
    parser.add_argument('--dt', '--time_step', type=float,
                        default=1.0, help='length of time step',
                        metavar='dt')
    return parser

def read_command_line():
    parser = define_command_line_options()
    args = parser.parse_args()
    print 'I={}, a={}, T={}, makeplot={}, theta={}, dt={}'.format(
        args.I, args.a, args.T, args.makeplot, args.theta, args.dt)
    return args.I, args.a, args.T, args.makeplot, \
           args.theta, args.dt

def main():
    I, a, T, makeplot, theta, dt = read_command_line()
    E = explore(I, a, T, dt, theta, makeplot)
    print 'theta=%.1f dt=%g Error=%.3E' % (theta, dt, E)

if __name__ == '__main__':
    if verify():
        main()
    else:
        print 'Bug in the implementation!'
