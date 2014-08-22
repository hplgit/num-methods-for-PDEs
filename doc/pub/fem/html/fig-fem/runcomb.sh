#!/bin/sh
sh comb.sh parabola_ls_sines4.pdf parabola_ls_sines12.pdf
mv -f tmp.pdf parabola_ls_sines4_12.pdf

sh comb.sh parabola_interp1_linear.pdf parabola_interp2_linear.pdf
mv -f tmp.pdf parabola_interp12_linear.pdf

sh comb.sh Lagrange_ls_sin_4.pdf Lagrange_interp_sin_4.pdf
mv -f tmp.pdf Lagrange_ls_interp_sin_4.pdf

sh comb.sh Lagrange_interp_abs_8.pdf Lagrange_interp_abs_15.pdf
mv -f tmp.pdf Lagrange_interp_abs_8_15.pdf

sh comb.sh Lagrange_interp_abs_Cheb_8.pdf Lagrange_interp_abs_Cheb_15.pdf
mv -f tmp.pdf Lagrange_interp_abs_Cheb_8_15.pdf

sh comb.sh fe_p1_x2_2e.pdf fe_p1_x2_4e.pdf
mv -f tmp.pdf fe_p1_x2_2e_4e.pdf

sh comb.sh fe_p1_x9_4e.pdf fe_p2_x9_2e.pdf fe_p1_x9_8e.pdf fe_p2_x9_4e.pdf
mv -f tmp.pdf fe_p1_p2_x9_248e.pdf

sh comb.sh ElmT6n2Da.map.pdf ElmT6n2D.map.pdf
mv -f tmp.pdf ElmT6n2Da_ElmT6n2D.pdf

sh comb.sh ElmT3n2D.map.fig.pdf ElmB4n2D.map.fig.pdf
mv -f tmp.pdf ElmT3n2D_ElmB4n2D.pdf

sh comb.sh ElmB8n2D.map.fig.pdf ElmB9n2D.map.fig.pdf
mv -f tmp.pdf ElmB8n2D_ElmB9n2D.pdf

