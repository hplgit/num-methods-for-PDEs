from sympy import *
u_np1, u_n, u_nm1, u_nm2, dt, w = \
       symbols('u_np1 u_n u_nm1 u_nm2 dt w')
half = Rational(1,2)

def energy(u_np1, u_n, u_nm1):
    DtDtu = (u_np1 - u_nm1)/(2*dt)
    return half*DtDtu**2 + half*w**2*(u_n)**2

def energy2(u_np1, u_n, u_nm1):
    Dtu = (u_np1 - u_nm1)/(2*dt)
    return half*Dtu**2 + half*w**2*(u_n)**2 - half*dt*w**2*Dtu*u_n

# Energy at t_n
E_n = energy2(u_np1, u_n, u_nm1)
E_n = E_n.subs(u_n, 2*u_nm1 - u_nm2 - dt**2*w**2*u_nm1).\
           subs(u_np1, 2*u_n - u_nm1 - dt**2*w**2*u_n)
E_n = expand(E_n)
print 'E_n:', E_n

# Energy at t_n-1
E_nm1 = energy2(u_n, u_nm1, u_nm2)
print 'E_nm1:', E_nm1

E_diff = E_n - E_nm1
E_diff = simplify(E_diff)
print 'Difference:', E_diff
print 'factor:', factor(E_diff, deep=True)
print 'series:', E_diff.series(dt, 0, 2)
