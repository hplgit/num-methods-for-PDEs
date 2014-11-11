from pysketcher import *
from numpy import exp, linspace

L = 10.
H = 3.
margin = 2
wall_thickness = 0.5

drawing_tool.set_coordinate_system(xmin=-margin, xmax=L+margin,
                                   ymin=-margin, ymax=H+margin,
                                   axis=False)

drawing_tool.set_linecolor('blue')
drawing_tool.set_fontsize(20)

bottom1 = Rectangle((0,H/2-wall_thickness),
                    L/3-wall_thickness, wall_thickness).\
                    set_filled_curves(pattern='/')
bottom2 = Rectangle((L/3-wall_thickness,0), wall_thickness, H/2).\
          set_filled_curves(pattern='/')
bottom3 = Rectangle((L/3-wall_thickness,-wall_thickness),
                    2*L/3+wall_thickness, wall_thickness).\
                    set_filled_curves(pattern='/')
top    = Line((0, H), (L,H)).set_linestyle('dashed')
inlet  = Line((0,H/2), (0,H)).set_linestyle('dashed')
outlet = Line((L,0), (L,H)).set_linestyle('dashed')

def velprofile(y):
    return [1, 0]

inlet_profile = VelocityProfile((0,H/2), H/2, velprofile, 5)

fig = Composition({
    'bottom': Composition(
        dict(bottom1=bottom1, bottom2=bottom2, bottom3=bottom3)),
    'top': top,
    'inlet': inlet_profile,
    'outlet': outlet,
    })

fig.draw()  # send all figures to plotting backend

symbols = {
    'u_x inlet':
    Text('$u_x=1$', (-0.2, 1.5*H/2), alignment='right'),
    'u_y inlet':
    Text('$u_y=0$', (-0.2, 1.5*H/2-0.6), alignment='right'),
    'bottom':
    Text('$u_x=u_y=0$', (L/2, -1), alignment='center'),
    'top':
    Text(r'$u_y=0,\ \frac{\partial u_x}{\partial n}=0$', (L/2, H+0.4), alignment='center'),
    'outlet1':
    Text(r'$\frac{\partial u_y}{\partial n}=0$', (L+0.1,H/2-1), alignment='left'),
    'outlet2':
    Text(r'$\frac{\partial u_x}{\partial n}=0$', (L+0.1,H/2), alignment='left'),
    'outlet3':
    Text('$p=p_0$', (L+0.1, H/2+0.75), alignment='left')}

symbols = Composition(symbols)
symbols.draw()

drawing_tool.display()
drawing_tool.savefig('os.path.splitext(__file__)[0]')


raw_input()
