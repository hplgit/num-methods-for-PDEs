#!/usr/bin/env python
import os, sys, glob
import decay_exper1

I = 1; a = 2; T = 5
decay_exper1.run_experiments(I=I, a=a, T=T)
dt_values_str = ', '.join(sys.argv[1:])  # needed in report

# Write HTML report
html = open('tmp_report.html', 'w')
# Need raw strings because of latex math with backslashes
html.write(r'''
<html>
<head>

<!-- Use MathJax to render mathematics -->
<script type="text/x-mathjax-config">
MathJax.Hub.Config({
  TeX: {
     equationNumbers: {  autoNumber: "AMS"  },
     extensions: ["AMSmath.js", "AMSsymbols.js", "autobold.js"]
  }
});
</script>
<script type="text/javascript"
 src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML">
</script>
</head>

<body bgcolor="white">

<title>Numerical investigations</title>

<center><h1>Experiments with Schemes for Exponential Decay</h1></center>

<center><b>Hans Petter Langtangen</b></center>
<center><b>Simula Research Laboratory</b></center>
<center><b>University of Oslo</b></center>
<center><h4>August 20, 2012</h4></center>

<b>Summary.</b> This report investigates the accuracy of three
finite difference schemes for the ordinary differential equation
\( u'=-au \) with the aid of numerical experiments.  Numerical
artifacts are in particular demonstrated.

<h2>Mathematical problem</h2>

We address the initial-value problem

$$
\begin{align}
u'(t) &= -au(t), \quad t \in (0,T], \label{ode}\\
u(0)  &= I,                         \label{initial:value}
\end{align}
$$
where \( a \), \( I \), and \( T \) are prescribed parameters,
and \( u(t) \) is the unknown function to be estimated.
This mathematical model is relevant for physical phenomena
featuring exponential decay in time, e.g., vertical pressure
variation in the atmosphere, cooling of an object, and
radioactive decay.


<h2>Numerical solution method</h2>

We introduce a mesh in time with points
\( 0= t_0< t_1 \cdots < t_{N_t}=T \).
For simplicity, we assume constant spacing \( \Delta t \)
between the mesh points: \( \Delta t = t_{n}-t_{n-1} \),
\( n=1,\ldots,N_t \). Let \( u^n \) be the numerical approximation
to the exact solution at \( t_n \).

The \( \theta \)-rule <a href="#Iserles_2009">[1]</a>
is used to solve \eqref{ode} numerically:

$$
u^{n+1} = \frac{1 - (1-\theta) a\Delta t}{1 + \theta a\Delta t}u^n,
$$
for \( n=0,1,\ldots,N_t-1 \). This scheme corresponds to

<ul>
  <li> The Forward Euler scheme when \( \theta=0 \)
  <li> The Backward Euler scheme when \( \theta=1 \)
  <li> The Crank-Nicolson scheme when \( \theta=1/2 \)
</ul>

<h2>Implementation</h2>

The numerical method is implemented in a Python function
<a href="#Langtangen_2012">[2]</a>:

<pre>
def solver(I, a, T, dt, theta):
    """Solve u'=-a*u, u(0)=I, for t in (0,T] with steps of dt."""
    Nt = int(round(T/float(dt)))  # no of intervals
    u = zeros(Nt+1)
    t = linspace(0, T, Nt+1)

    u[0] = I
    for n in range(0, Nt):
        u[n+1] = (1 - (1-theta)*a*dt)/(1 + theta*dt*a)*u[n]
    return u, t
</pre>

<h2>Numerical experiments</h2>

We define a set of numerical experiments where
\( I \), \( a \), and \( T \) are fixed, while
\( \Delta t \) and \( \theta \) are varied.
In particular, \( I=%(I)g \), \( a=%(a)g \),
\( \Delta t= %(dt_values_str)s \).

''' % vars())

short2long = dict(FE='The Forward Euler method',
                  BE='The Backward Euler method',
                  CN='The Crank-Nicolson method')
for method in 'BE', 'CN', 'FE':
    html.write('\n<h3>%s</h3>\n' % short2long[method])
    html.write('<img src="%s.png" width="800">\n\n' % method)
# Remember raw strings for latex math with backslashes
html.write(r"""

<h3>Error versus time discretization</h3>

How \( E \) varies with \( \Delta t \) for
\( \theta = 0, 0.5, 1 \) is shown below.

<p><b>Observe:</b>
The data points for the three largest \( \Delta t \) values in the
Forward Euler method are not relevant as the solution behaves
non-physically.

<p>
<img="error.png", width="400">

<h3>Summary</h3>

<ol>
  <li> \( \theta =1 \): \( E\sim \Delta t \) (first-order convergence).
  <li> \( \theta =0.5 \): \( E\sim \Delta t^2 \) (second-order convergence).
  <li> \( \theta =1 \) is always stable and gives qualitatively corrects
       results.
  <li> \( \theta =0.5 \) never blows up, but may give oscillating solutions
       if \( \Delta t \) is not sufficiently small.
  <li> \( \theta =0 \) suffers from fast-growing solution if \( \Delta t \) is
       not small enough, but even below this limit one can have oscillating
       solutions that disappear if \( \Delta t \) is sufficiently small.
</ol>

<h2>Bibliography</h2>

<ol>
 <li> <a name="Iserles_2009"></a> <b>A. Iserles</b>.
    <em>A First Course in the Numerical Analysis of Differential Equations</em>,
    Cambridge University Press, 2009.
 <li> <a name="Langtangen_2012"></a> <b>H. P. Langtangen</b>.
    <em>A Primer on Scientific Programming With Python</em>,
    Springer, 2012.
</ol>


</body>
</html>
""")
