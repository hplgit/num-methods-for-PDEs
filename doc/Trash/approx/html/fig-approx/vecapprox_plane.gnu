#unset tics

a=1.
b=0.6
set xrange [0:6]
set yrange [0:6]
plot (b/a)*x notitle with lines lt 2
set arrow from 0,0 to a,b linewidth 1
set label "(a,b)" at a+0.025,b-0.05
set arrow from 0,0 to 3,5 linewidth 1
set label "(3,5)" at 3+0.025,5+0.05
c0 = (3*a + 5*b)/(a*a + b*b)
set arrow from 0,0 to c0*a,c0*b linewidth 1
set label "c_0(a,b)" at c0*a+0.025,c0*b-0.05
set arrow from c0*a,c0*b to 3,5 ls 3
replot
# the size command is important for keeping aspect ratio correct
set term postscript eps enhanced size 10cm,10cm
set output "vecapprox_plane.eps"
replot
set term png fontscale 1.0 size 480, 480 background "#ffffff"
set output "vecapprox_plane.png"
replot
