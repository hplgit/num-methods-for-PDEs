from pysketcher import *
from numpy import exp, linspace

L = 10.
H = 6.
cylinder_center = (L/4, H/2)
radius = 0.8
wall_thickness = 0.5
margin = 1

drawing_tool.set_coordinate_system(xmin=-margin, xmax=L+margin,
                                   ymin=-margin, ymax=H+margin,
                                   axis=False)

drawing_tool.set_linecolor('blue')
outer = Rectangle((0,0), L, H).set_linestyle('dashed')
cylinder = Circle(cylinder_center, radius).set_linewidth(3)


def velprofile(y):
    return [0.75, 0]

inlet_profile = VelocityProfile((0,0), H, velprofile, 10)

fig = Composition({
    'outer': outer,
    'inlet': inlet_profile,
    'cylinder': cylinder,
    })

fig.draw()  # send all figures to plotting backend

symbols = {
    'A': Text('A', (-0.1,0), alignment='right'),
    'B': Text('B', (L+0.1,0), alignment='left'),
    'D': Text('D', (-0.1,H), alignment='right'),
    'C': Text('C', (L+0.1,H), alignment='left'),
}

symbols = Composition(symbols)
symbols.draw()

drawing_tool.display()
drawing_tool.savefig(os.path.splitext(__file__)[0])


raw_input()
