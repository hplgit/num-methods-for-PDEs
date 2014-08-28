def dispersion_relation_2D(kh, theta, C):
    arg = C*sqrt(sin(0.5*kh*cos(theta))**2 +
                 sin(0.5*kh*sin(theta))**2)
    c_frac = 2./(C*kh)*arcsin(arg)

    return c_frac

from numpy import exp, sin, cos, linspace, \
     pi, meshgrid, arcsin, sqrt
r = kh = linspace(0.001, pi, 101)
theta = linspace(0, 2*pi, 51)
r, theta = meshgrid(r, theta)

# Make 2x2 filled contour plots for 4 values of C
import matplotlib.pyplot as plt
C_max = 1/sqrt(2)
C = [[C_max, 0.9*C_max], [0.5*C_max, 0.2*C_max]]
fix, axes = plt.subplots(2, 2, subplot_kw=dict(polar=True))
for row in range(2):
    for column in range(2):
        error = 1 - dispersion_relation_2D(
            kh, theta, C[row][column])
        print error.min(), error.max()
        cax = axes[row][column].contourf(
            theta, r, error, 50, vmin=0, vmax=0.36)
        axes[row][column].set_xticks([])
        axes[row][column].set_yticks([])

# Add colorbar to the last plot
cbar = plt.colorbar(cax)
cbar.ax.set_ylabel('error in wave velocity')
plt.savefig('disprel2D.png')
plt.savefig('disprel2D.pdf')
plt.show()

# See
# http://blog.rtwilson.com/producing-polar-contour-plots-with-matplotlib/
# for polar plotting in matplotlib
