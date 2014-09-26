#!/usr/bin/env python
import os, sys, glob
import decay_exper1

I = 1; a = 2; T = 5
decay_exper1.run_experiments(I=I, a=a, T=T)

# Write Doconce report
do = open('tmp_report.do.txt', 'w')
title = 'Experiments with Schemes for Exponential Decay'
author1 = 'Hans Petter Langtangen Email:hpl@simula.no at '\
          'Center for Biomedical Computing, '\
          'Simula Research Laboratory & '\
          'Department of Informatics, University of Oslo.'
date = 'today'

dt_values_str = ', '.join(sys.argv[1:])  # needed in report

# Remember to use raw string because of latex commands for math
do.write(r"""
TITLE: %(title)s
AUTHOR: %(author1)s
DATE: %(date)s

__Summary.__
This report investigates the accuracy of three finite difference
schemes for the ordinary differential equation $u'=-au$ with the
aid of numerical experiments. Numerical artifacts are in particular
demonstrated.

## Include table of contents (latex and html; sphinx always has a toc).
## (Lines starting with ## are not propagated to the output file,
## otherwise comments lines starting with # are visible in the
## output file.)

TOC: on


======= Mathematical problem =======
label{math:problem}

## Purpose: section with multi-line equation.
## idx{keyword} adds keyword to the index
## (to be place before the actual paragraph).

idx{model problem} idx{exponential decay}

We address the initial-value problem

!bt
\begin{align}
u'(t) &= -au(t), \quad t \in (0,T], label{ode}\\
u(0)  &= I,                         label{initial:value}
\end{align}
!et
where $a$, $I$, and $T$ are prescribed parameters, and $u(t)$ is
the unknown function to be estimated. This mathematical model
is relevant for physical phenomena featuring exponential decay
in time, e.g., vertical pressure variation in the atmosphere,
cooling of an object, and radioactive decay.

======= Numerical solution method =======
label{numerical:problem}

## Purpose: section with single-line equation, equation reference,
## a bullet list, and links to URLs ("link text": "url").

idx{mesh in time} idx{$\theta$-rule} idx{numerical scheme}
idx{finite difference scheme}

We introduce a mesh in time with points $0= t_0< t_1 \cdots < t_{N_t}=T$.
For simplicity, we assume constant spacing $\Delta t$ between the
mesh points: $\Delta t = t_{n}-t_{n-1}$, $n=1,\ldots,N_t$. Let
$u^n$ be the numerical approximation to the exact solution at $t_n$.

The $\theta$-rule cite{Iserles_2009}
is used to solve (ref{ode}) numerically:

!bt
\[
u^{n+1} = \frac{1 - (1-\theta) a\Delta t}{1 + \theta a\Delta t}u^n,
\]
!et
for $n=0,1,\ldots,N_t-1$. This scheme corresponds to

  * The "Forward Euler":
    "http://en.wikipedia.org/wiki/Forward_Euler_method"
    scheme when $\theta=0$
  * The "Backward Euler":
    "http://en.wikipedia.org/wiki/Backward_Euler_method"
    scheme when $\theta=1$
  * The "Crank-Nicolson":
    "http://en.wikipedia.org/wiki/Crank-Nicolson"
    scheme when $\theta=1/2$


======= Implementation =======

## Purpose: section with computer code taken from a part of
## a file. The fromto: f@t syntax copies from the regular
## expression f up to the line, but not including, the regular
## expression t.

The numerical method is implemented in a Python function
cite{Langtangen_2012} `solver` (found in the "`decay_mod`":
"https://github.com/hplgit/INF5620/blob/gh-pages/src/decay/experiments/dc_mod.py" module):

@@@CODE ../decay_mod.py  fromto: def solver@def verify_three


======= Numerical experiments =======

## Purpose: section with figures.

idx{numerical experiments}

We define a set of numerical experiments where $I$, $a$, and $T$ are
fixed, while $\Delta t$ and $\theta$ are varied.
In particular, $I=%(I)g$, $a=%(a)g$, $\Delta t = %(dt_values_str)s$.

""" % vars())

short2long = dict(FE='The Forward Euler method',
                  BE='The Backward Euler method',
                  CN='The Crank-Nicolson method')
methods = 'BE', 'CN', 'FE'
inline_figures = True  # True: subsections with inline graphics
                       # False: no subsecs, figures with captions
if inline_figures:
    for method in methods:
        do.write("""

===== %s =====

## Purpose: subsection with inline figure (figure without caption).


idx{%s}

FIGURE: [%s, width=800]

""" % (short2long[method], method, method))
else:
    do.write("""

===== Time series =====

Figures ref{fig:BE}-ref{fig:FE} display the results.

""")
    # Full figures with captions
    for method in methods:
        fig = 'FIGURE: [%s, width=800] %s. '\
              'label{fig:%s}\n\n' % \
              (method, short2long[method], method)
        do.write(fig)

# Remember raw string for latex math with backslashes
do.write(r"""

===== Error vs $\Delta t$ =====

## Purpose: exemplify referring to a figure with label and caption.

idx{error vs time step}

How $E$ varies with $\Delta t$ for $\theta=0,0.5,1$
is shown in Figure ref{fig:error}.

## Here is an admonition box for warnings

!bwarning Observe:
The data points for the three largest $\Delta t$ values in the
Forward Euler method are not relevant as the solution behaves
non-physically.
!ewarning

FIGURE: [error, width=400] Variation of the error with the time step. label{fig:error}

## A corresponding table

The numbers corresponding to the figure above are given in the table below.

|------c--------------c--------------c--------------c-------|
| $\Delta t$   | $\theta=0$   | $\theta=0.5$ | $\theta=1$   |
|------r--------------r--------------r--------------r-------|
| 1.25         | 7.4630       | 0.2161       | 0.2440       |
| 0.75         | 0.6632       | 0.0744       | 0.1875       |
| 0.50         | 0.2797       | 0.0315       | 0.1397       |
| 0.10         | 0.0377       | 0.0012       | 0.0335       |
|-----------------------------------------------------------|

## Here is an admonition box for summaries

!bsummary
 o $\theta =1$: $E\sim \Delta t$ (first-order convergence).
 o $\theta =0.5$: $E\sim \Delta t^2$ (second-order convergence).
 o $\theta =1$ is always stable and gives qualitatively corrects results.
 o $\theta =0.5$ never blows up, but may give oscillating solutions
   if $\Delta t$ is not sufficiently small.
 o $\theta =0$ suffers from fast-growing solution if $\Delta t$ is
   not small enough, but even below this limit one can have oscillating
   solutions that disappear if $\Delta t$ is sufficiently small.
!esummary

======= Bibliography =======

## Publish (https://bitbucket.org/logg/publish is used to
## handle references. The line below specifies the name of
## the Publish database file (see the doconce manual for details).

BIBFILE: .publish_references.pub
""")

# Good habits when writing for latex, sphinx and mathjax-html
# at once:
#
# Minimize math in captions as the caption becomes the figure
# name in sphinx, when referring to figures, and any math
# is deleted in the name.
#
# Use \[, equation or align enviroments in latex math.
# Sphinx cannot handle labels in align environments.
#
# Figures float around in latex, but are placed at where
# they are defined in sphinx and html. Figures without captions
# are placed inline in latex and may be convenient.
#
# Remember raw strings for any text with latex math with
# backslashes.
