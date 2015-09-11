"""
Illustrate a spring with point masses and constant tension.
For derivation of the wave equation for a spring.
"""

from pysketcher import *

#test_test()
xaxis = 1
drawing_tool.set_coordinate_system(0.8, 5.2, 0, 4.2, axis=False)

x = 3
xb = 2
xf = 4
xf2 = 5
xb2 = 1

p = point(x, 3)
pb = point(xb, 2.5)
pf = point(xf, 2.7)
pf2 = point(xf2, 2.2)
pb2 = point(xb2, 2.6)

force_L = Arrow1(p, p + 0.7*unit_vec(pb-p)).\
          set_linecolor('black').set_linewidth(1)
force_R = Arrow1(p, p + 0.7*unit_vec(pf-p)).\
          set_linecolor('black').set_linewidth(1)
text_L = Text('$T$', force_L.geometric_features()['end'] + point(0,0.1))
text_R = Text('$T$', force_R.geometric_features()['end'] + point(0,0.1))
forces = Composition([force_L, force_R, text_L, text_R])

# Masses
r = 0.1
c   = Circle(p, r).  set_linecolor('black').set_filled_curves('black')
cf  = Circle(pf, r). set_linecolor('black').set_filled_curves('black')
cb  = Circle(pb, r). set_linecolor('black').set_filled_curves('black')
cf2 = Circle(pf2, r).set_linecolor('black').set_filled_curves('black')
cb2 = Circle(pb2, r).set_linecolor('black').set_filled_curves('black')

# u values
v = point(0,0.25)  # vertical displacement of text over masses
fnts = 18
u_values = Composition([Text('$u_{i-1}$', pb + v, fontsize=fnts),
                        Text('$u_{i+1}$', pf + v, fontsize=fnts),
                        Text('$u_{i}$',   p  + v, fontsize=fnts)
                        ])

# Dashed lines along the string
line_L1 = Line(p, pb).set_linestyle('dashed').\
          set_linecolor('black').set_linewidth(1)
line_R1 = Line(p, pf).set_linestyle('dashed').\
          set_linecolor('black').set_linewidth(1)
line_L2 = Line(pb, pb2).set_linestyle('dashed').\
          set_linecolor('black').set_linewidth(1)
line_R2 = Line(pf, pf2).set_linestyle('dashed').\
          set_linecolor('black').set_linewidth(1)

masses = Composition([c, cf, cb, cf2, cb2, line_L1, line_R1, line_L2, line_R2])

# Points in the mesh
p0 = point(x, xaxis)
pf0 = point(xf, xaxis)
pb0 = point(xb, xaxis)
tick = 0.05

# 1D mesh with three points
mesh = Composition({
    'xim1': Text('$x_{i-1}$', pb0 - point(0, 0.3)),
    'xi': Text('$x_{i}$', p0 - point(0, 0.3)),
    'xip1': Text('$x_{i+1}$', pf0 - point(0, 0.3)),
    'axis': Composition({
        'hline': Line(pf0-point(3,0), pb0+point(3,0)).\
        set_linecolor('black').set_linewidth(1),
        'tick_m1': Line(pf0+point(0,tick), pf0-point(0,tick)).\
        set_linecolor('black').set_linewidth(1),
        'tick_i':  Line(p0+point(0,tick), p0-point(0,tick)).\
        set_linecolor('black').set_linewidth(1),
        'tick_p1': Line(pb0+point(0,tick), pb0-point(0,tick)).\
        set_linecolor('black').set_linewidth(1)}),
    })

fig = Composition(dict(masses=masses, mesh=mesh,
                       forces=forces, u_values=u_values))

drawing_tool.erase()
fig.draw()
drawing_tool.display()
drawing_tool.savefig('wave_on_string')
raw_input()

