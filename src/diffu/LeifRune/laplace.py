import numpy as np
import scipy.linalg
import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import scipy as sc
import scipy.linalg
import scipy.sparse
import scipy.sparse.linalg


# discretizing geometry
width = 1.0
height = 1.0
Nx = 20 # number of points in x-direction
Ny = 20 # number of points in y-direction
N = Nx*Ny

# BCs
bottom = 10.0
left = 10.0
top = 30.0
right = 10.0

# diagonals
diag1, diag5 = np.zeros(N-3), np.zeros(N-3)
diag2, diag4 = np.zeros(N-1), np.zeros(N-1)
diag3 = np.zeros(N)
diag1[:], diag2[:], diag3[:], diag4[:], diag5[:] = -1.0, -1.0, 4.0, -1.0, -1.0

diagonals=np.zeros((5,N))
diagonals[0,:] = -1.0                       #all elts in first row is set to 1
diagonals[1,:] = -1.0
diagonals[2,:] =  4.0
diagonals[3,:] = -1.0
diagonals[4,:] = -1.0

# impose BCs on diag2 and diag4
diag2[Nx-1::Nx] = 0.0
diag4[Nx-1::Nx] = 0.0

#impose BCs for sparse solver
diagonals[1,Nx-1::Nx] = 0.0
diagonals[3,Nx::Nx]   = 0.0

# assemble coefficient matrix A
A = np.zeros((N,N))
for i in range(0,N):
    A[i,i] = diag3[i]
for i in range(0,N-1):
    A[i+1,i] = diag2[i]
    A[i,i+1] = diag4[i]
for i in range(0,N-Nx):
    A[i,i+Nx] = diag5[i]
    A[i+Nx,i] = diag1[i]

# assemble the right hand side vector b
b = np.zeros(N)
b[:Nx]  = b[:Nx] + bottom
b[-Nx:] = b[-Nx:] + top
b[::Nx] = b[::Nx] + left
b[Nx-1::Nx] = b[Nx-1::Nx] + right


# solve the equation system
As = sc.sparse.spdiags(diagonals, [-Nx,-1,0,1,Nx], N, N,format='csc') #sparse matrix instance, 
Tvector = sc.sparse.linalg.spsolve(As,b) # Compute the solution with a sparse solver.
Tvector2 = scipy.linalg.solve(A,b)

print 'max diff = ', np.max(Tvector - Tvector2)
print 'min diff = ', np.min(Tvector - Tvector2)

Tmatrix = np.reshape(Tvector,(Nx,Ny))

# adding the boundary points (for plotting)
T = np.zeros((Nx+2,Ny+2))

# Set the boundary values 
T[:,0]    = left
T[:,Ny+1] = right
T[0,:]    = bottom
T[Nx+1,:] = top

# Assign the computed values to the field of the T
T[1:Nx+1,1:Ny+1]=Tmatrix[0:Nx,0:Ny]

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
#plt.savefig('oppg1.pdf')
