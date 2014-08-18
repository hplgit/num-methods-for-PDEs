import numpy as np
import scipy.linalg
import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import scipy as sc
import scipy.linalg
import scipy.sparse
import scipy.sparse.linalg

class Grid:
    """A simple grid class that stores the details and solution of the
    computational grid."""
    def __init__(self, nx=10, ny=10, xmin=0.0, xmax=1.0,
                 ymin=0.0, ymax=1.0):
        self.xmin, self.xmax, self.ymin, self.ymax = xmin, xmax, ymin, ymax
        self.dx = float(xmax-xmin)/(nx-1)
        self.dy = float(ymax-ymin)/(ny-1)
        self.T = np.zeros((nx+2, ny+2), 'd')
        self.nx = nx
        self.ny = ny
    
    def setBC(self, top=30, bottom=10, left=10, right = 10):
        self.T[0 ,:] = bottom; self.bottom=bottom   
        self.T[-1,:] = top;    self.top=top
        self.T[:, 0] = left;   self.left=left
        self.T[:,-1] = right;  self.right=right

class LaplaceSolver:
    """Solvers for the laplacian equation. Scheme can be one of
        ['direct','slow', 'numeric', 'blitz', 'inline', 'fastinline','fortran'].  """
    def __init__(self,grid,scheme='direct'):
        self.grid = grid

        if scheme == 'slow':
            self.solver = self.slowTimeStep
        elif scheme == 'direct':
            self.solver =self.directSolver
        elif scheme == 'blitz':
             self.solver = self.blitzTimeStep
        else:
            self.solver = self.numericTimeStep

    def directSolver(self,dt=0):
        g = self.grid
        T = self.grid.T

        bottom, top, left, right = g.bottom, g.top, g.left, g.right
        nx, ny = g.T.shape #number of grid pts excluding bndry
        n = nx*ny
        
        dgs=np.zeros((5,N))
        dgs[0,:] = -1.0
        dgs[1,:] = -1.0
        dgs[2,:] =  4.0 
        dgs[3,:] = -1.0
        dgs[4,:] = -1,0

        #impose BCs for sparse solver
        dgs[1,nx-1::nx] = 0.0
        dgs[3,nx::nx]   = 0.0
        
        # assemble the right hand side vector b
        b = np.zeros(n)
        b[:nx]  = b[:nx]  + bottom
        b[-nx:] = b[-nx:] + top
        b[::nx] = b[::nx] + left
        b[nx-1::nx] = b[nx-1::nx] + right

        Asp  = sc.sparse.spdiags(dgs, [-nx,-1,0,1,nx], n, n,format='csc') #sparse matrix instance, 
        Tv = sc.sparse.linalg.spsolve(Asp,b) # Compute the solution with sparse solver
        Tmatrix = np.reshape(Tv,(nx,ny))

        # Assign the computed values to the field of the T
        T[1:nx+1,1:ny+1]=Tmatrix[0:nx-1,0:ny-1]
        T[5:8,5:8] = 30
#        T[1:nx,1:ny]=Tmatrix[:,:]

        g.T = T
        return 'Ok'

    def solve(self):
        return self.solver()


# discretizing geometry
width = 1.0
height = 1.0
Nx = 20 # number of unknowns in the x-direction 
Ny = 20 # number of unknowns in the y-direction 
N = Nx*Ny



# BCs
bottom = 10.0
left = 10.0
top = 30.0
right = 10.0


diagonals=np.zeros((5,N))
diagonals[0,:] = -1.0                       #all elts in first row is set to 1
diagonals[1,:] = -1.0
diagonals[2,:] =  4.0
diagonals[3,:] = -1.0
diagonals[4,:] = -1.0

#impose BCs for sparse solver
diagonals[1,Nx-1::Nx] = 0.0
diagonals[3,Nx::Nx]   = 0.0

# assemble the right hand side vector b
b = np.zeros(N)
b[:Nx]  = b[:Nx] + bottom
b[-Nx:] = b[-Nx:] + top# assemble coefficient matrix A
b[::Nx] = b[::Nx] + left
b[Nx-1::Nx] = b[Nx-1::Nx] + right


# solve the equation system
As = sc.sparse.spdiags(diagonals, [-Nx,-1,0,1,Nx], N, N,format='csc') #sparse matrix instance, 
Tvector = sc.sparse.linalg.spsolve(As,b) # Compute the solution with a sparse solver.

Tmatrix = np.reshape(Tvector,(Nx,Ny))

# adding the boundary points (for plotting)
T = np.zeros((Nx+2,Ny+2))

myGrid=Grid(nx=Nx,ny=Nx)
myGrid.setBC(top=top,bottom=bottom,left=left,right=right)
#mySolver=LaplaceSolver(myGrid,scheme='direct')
mySolver=LaplaceSolver(Grid(nx=Nx,ny=Nx),scheme='direct')
mySolver.grid.setBC(top=top,bottom=bottom,left=left,right=right)
mySolver.solve()

# Set the boundary values 
T[:,0]    = left
T[:,Ny+1] = right
T[0,:]    = bottom
T[Nx+1,:] = top

# Assign the computed values to the field of the T
T[1:Nx+1,1:Ny+1]=Tmatrix[0:Nx,0:Ny]


mySolver.grid.T[5:5,5:5]=30.0
T=mySolver.grid.T[:,:]
# plotting
x = np.linspace(0,width,Nx+2)
y = np.linspace(0,height,Ny+2)
Tmax = np.max(T)
Tmin = np.min(T)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
X,Y = np.meshgrid(x,y)
surf = ax.plot_surface(X, Y, T, rstride=1, cstride=1, cmap=cm.jet)
ax.set_zlim3d(Tmin,Tmax)
fig.colorbar(surf)
plt.title('Temperature field in beam cross section')
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Temperature')
ax.view_init(elev=10., azim=-140)
plt.show()
# #plt.savefig('oppg1.pdf')
