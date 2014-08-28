"""
Illustrate forward, backward and centered finite differences
in four figures.
"""

from pysketcher import *

#test_test()
xaxis = 2
drawing_tool.set_coordinate_system(0, 7, 2, 6, axis=True)

f = SketchyFunc1('$B(x)$')
x = 4                 # point where we want the slide
p = (x, f(x))         # curve point
r = 0.1               # radius of circles placed at key points
c = Circle(p, r).set_linecolor('blue')


# Tangents illustrating the derivative
domain = [2.5, 5]
h = 1E-3  # h in finite difference approx used to compute the exact tangent
tangent = Line((x+h, f(x+h)), (x-h, f(x-h))).\
          new_interval(x=domain).\
          set_linestyle('dotted').set_linecolor('black')

exact = Composition(dict(graph=f, tangent=tangent))

def g(x, g_m, g_s=1, g_0=0, g_a=1):
    return g_0 + g_a*exp(-((x-g_m)/float(g_s))**2)

from numpy import linspace, exp, arctan, pi
g_s = 0.3
g_m = 1.5
g_0 = 2.5
xc = linspace(g_m - 2.5*g_s, g_m + 2.5*g_s, 21)
slide = Curve(xc, g(xc, g_m, g_s, g_0, 0.4))

drawing_tool.erase()
slide.draw()
drawing_tool.display()
drawing_tool.savefig('slide1')

L1 = tangent.geometric_features()['start']
L0 = tangent.geometric_features()['end']

angle = -arctan((L0[1]-L1[1])/(L1[0]-L0[0]))*180/pi
print (x-g_m, p[1]-g_0)
print xc
print g(xc, g_m, g_s, g_0, 0.4)
print xc + x-g_m
print g(xc, g_m, g_s, g_0, 0.4) + p[1]-g_0

slide.rotate(angle, (g_m, g_0))
slide.translate((x-g_m, p[1]-g_0))
fig = Composition(dict(curve=f, slide=slide, tangent=tangent, circle=c))

drawing_tool.erase()
fig.draw()
drawing_tool.display()
drawing_tool.savefig('curve_with_slide')
raw_input()

