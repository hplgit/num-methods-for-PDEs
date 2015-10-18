# The equation solved is the parabolic equaiton
#
#       du     d  du
#       -- = k -- --
#       dt     dx dx
#
# along with boundary conditions

import matplotlib.pyplot as plt
import matplotlib
# Change some default values to make plots more readable on the screen
LNWDT=2; FNT=15
matplotlib.rcParams['lines.linewidth'] = LNWDT; matplotlib.rcParams['font.size'] = FNT

import numpy as np
import scipy as sc
import scipy.sparse
import scipy.sparse.linalg


class Grid1d:
    """A simple grid class for grid information and the solution."""
    def __init__(self, nx=10, xmin=0.0, xmax=1.0):
        self.xmin, self.xmax = xmin, xmax
        self.dx = float(xmax-xmin)/(nx)
        self.nx = nx                           # Number of dx  
        self.u = np.zeros((nx+1, 1), 'd')      # Number of x-values is nx+1
        self.x = np.linspace(xmin,xmax,nx+1)


class HeatSolver1d:
    """A simple 1dheat equation solver that can use different schemes to solve the problem."""
    def __init__(self, grid, scheme='explicit',k=1.0,r=0.5,theta=1.0):
        self.grid = grid
        self.setSolverScheme(scheme)
        self.k = k
        self.r = r
        self.theta = theta #Used for implicit solver only

    def setSolverScheme(self, scheme='explicit'):
        """Sets the scheme to be which should be one of ['slow', 'explicit', 'implicit']."""
        if scheme == 'slow':
            self.solver = self.pythonExplicit
            self.name = 'python'
        elif scheme == 'explicit':
            self.solver = self.numpyExplicit
            self.name = 'explicit'
        elif scheme == 'implicit':
            self.solver = self.numpyImplicit
            self.name   = 'implicit'
        else:
            self.solver = self.numpyImplicit
            self.name   = 'implicit'

    def numpyExplicit(self, tmin, tmax, nPlotInc):
        """Solve equation for all t in time step using a NumPy expression."""
        g = self.grid
        k = self.k       #Diffusivity
        r = self.r       #Numerical Fourier number
        u, x, dx  = g.u, g.x, g.dx
        xmin, xmax = g.xmin, g.xmax
        
        dt=r*dx**2/k     #Compute timestep based on Fourier number, dx and diffusivity

        m=round((tmax-tmin)/dt) # Number of temporal intervals
        time=np.linspace(tmin,tmax,m)
        
        for t in time:
            u[1:-1] =  r*(u[0:-2]+ u[2:]) + (1.0-2.0*r)*u[1:-1]
    
        g.u=u

    def pythonExplicit(self, tmin, tmax, nPlotInc):
        """Solve equation for all t in time step using a NumPy expression."""
        g = self.grid
        k = self.k       #Diffusivity
        r = self.r       #Numerical Fourier number
        u, x, dx, n  = g.u, g.x, g.dx, g.nx
        xmin, xmax = g.xmin, g.xmax
        
        dt=r*dx**2/k     #Compute timestep based on Fourier number, dx and diffusivity

        m=round((tmax-tmin)/dt) # Number of temporal intervals
        time=np.linspace(tmin,tmax,m)
        
        for t in time:
            u0 = u
#            u[1:-1] =  r*(u[0:-2]+ u[2:]) + (1.0-2.0*r)*u[1:-1]
            for i in range(1,n):
                u[i] =  r*(u[i-1]+ u[i+1]) + (1.0-2.0*r)*u[i]


    def numpyImplicit(self, tmin, tmax,nPlotInc):
        g = self.grid
        k = self.k         #Diffusivity
        r = self.r         #Numerical Fourier number
        theta =self.theta  #Parameter for implicitness: theta=0.5 Crank-Nicholson, theta=1.0 fully implicit
        u, x, dx  = g.u, g.x, g.dx
        xmin, xmax = g.xmin, g.xmax
        
        dt=r*dx**2/k     #Compute timestep based on Fourier number, dx and diffusivity
 
        m=round((tmax-tmin)/dt) # Number of temporal intervals
        time=np.linspace(tmin,tmax,m)

        #Create matrix for sparse solver. Solve for interior values only (nx-1)
        diagonals=np.zeros((3,g.nx-1))   
        diagonals[0,:] = -r*theta                       #all elts in first row is set to 1
        diagonals[1,:] = 1+2.0*r*theta  
        diagonals[2,:] = -r*theta 
        As = sc.sparse.spdiags(diagonals, [-1,0,1], g.nx-1, g.nx-1,format='csc') #sparse matrix instance

        #Crete rhs array
        d=np.zeros((g.nx-1,1),'d')
        
        #Advance in time an solve tridiagonal system for each t in time
        for t in time:
            d[:] = u[1:-1]+r*(1-theta)*(u[0:-2]-2*u[1:-1]+u[2:])  
            d[0] += r*theta*u[0]
            w = sc.sparse.linalg.spsolve(As,d) #theta=sc.linalg.solve_triangular(A,d)
            u[1:-1] = w[:,None]
           
        g.u=u
        
           
    def solve(self, tmin, tmax,nPlotInc=5):
        return self.solver(tmin,tmax,nPlotInc)

    def initialize(self,U0=1.0):
        self.grid.u[0] = U0
        

## Main program

## Make grids for the solvers
nx = 120
L  = 1.0
# mg  = Grid1d(nx,0,L)
# mg2 = Grid1d(nx,0,L)
# mg3 = Grid1d(nx,0,L)
# mg4 = Grid1d(nx,0,L)
# mg5 = Grid1d(nx,0,L)

## Make various solvers.
solvers=[]
#solvers.append(HeatSolver1d(mg,  scheme = 'slow',     k=1.0, r=0.5))
#solvers.append(HeatSolver1d(Grid1d(nx,0,L),  scheme = 'slow',     k=1.0, r=0.5))
solvers.append(HeatSolver1d(Grid1d(nx,0,L), scheme = 'explicit', k=1.0, r=0.5))
solvers.append(HeatSolver1d(Grid1d(nx,0,L), scheme = 'explicit', k=1.0, r=0.5))
solvers.append(HeatSolver1d(Grid1d(nx,0,L), scheme = 'implicit', k=1.0,  r=3.0, theta=0.5))


U0=1.0
(tmin, tmax)=(0,0.025)

## Compute a solution for all solvers
for solver in solvers:
    solver.initialize(U0=U0)
    solver.solve(tmin,tmax,nPlotInc=2)
lstyle=['r-',':','.','-.','--']
mylegends=[]
i=0
for solver in solvers:
    plt.plot(solver.grid.x,solver.grid.u,lstyle[i])
    mylegends.append(str('%s r = %3.1f' % (solver.name, solver.r)))
    i+=1

plt.legend(mylegends)
plt.show()
#plt.pause(5)
#plt.close()
