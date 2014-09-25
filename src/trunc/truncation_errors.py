import sympy as sp

class TaylorSeries:
    """Class for symbolic Taylor series."""
    def __init__(self, f, num_terms=4):
        self.f = f
        self.N = num_terms
        # Introduce symbols for the derivatives
        self.df = [f]
        for i in range(1, self.N+1):
            self.df.append(sp.Symbol('D%d%s' % (i, f.name)))

    def __call__(self, h):
        """Return the truncated Taylor series at x+h."""
        terms = self.f
        for i in range(1, self.N+1):
            terms += sp.Rational(1, sp.factorial(i))*self.df[i]*h**i
        return terms


class DiffOp:
    """Class for discrete difference operators."""
    def __init__(self, f, independent_variable='x',
                 num_terms_Taylor_series=4):
        self.Taylor = TaylorSeries(f, num_terms_Taylor_series)
        self.f = self.Taylor.f
        self.h = sp.Symbol('d%s' % independent_variable)

        # Finite difference operators
        h, f, f_T = self.h, self.f, self.Taylor # short names
        theta = sp.Symbol('theta')
        self.diffops = {
            'Dtp':  (f_T(h) - f)/h,
            'Dtm':  (f - f_T(-h))/h,
            'Dt':   (f_T(h/2) - f_T(-h/2))/h,
            'D2t':  (f_T(h) - f_T(-h))/(2*h),
            'DtDt': (f_T(h) - 2*f + f_T(-h))/h**2,
            'barDt': (f_T((1-theta)*h) - f_T(-theta*h))/h,
            }
        self.diffops = {diffop: sp.simplify(self.diffops[diffop])
                        for diffop in self.diffops}

        self.diffops['weighted_arithmetic_mean'] = \
             self._weighted_arithmetic_mean()
        self.diffops['geometric_mean'] = self._geometric_mean()
        self.diffops['harmonic_mean'] = self._harmonic_mean()

    def _weighted_arithmetic_mean(self):
        # The expansion is around n*h + theta*h
        h, f, f_T = self.h, self.f, self.Taylor
        theta = sp.Symbol('theta')
        f_n = f_T(-h*theta)
        f_np1 = f_T((1-theta)*h)
        a_mean = theta*f_np1 + (1-theta)*f_n
        return sp.expand(a_mean)

    def _geometric_mean(self):
        h, f, f_T = self.h, self.f, self.Taylor
        f_nmhalf = f_T(-h/2)
        f_nphalf = f_T(h/2)
        g_mean = f_nmhalf*f_nphalf
        return sp.expand(g_mean)

    def _harmonic_mean(self):
        h, f, f_T = self.h, self.f, self.Taylor
        f_nmhalf = f_T(-h/2)
        f_nphalf = f_T(h/2)
        h_mean = 2/(1/f_nmhalf + 1/f_nphalf)
        return sp.expand(h_mean)

    def D(self, i):
        """Return the symbol for the i-th derivative."""
        return self.Taylor.df[i]

    def __getitem__(self, operator_name):
        return self.diffops.get(operator_name, None)

    def operator_names(self):
        """Return all names for the operators."""
        return list(self.diffops.keys())

def truncation_errors():
    # Make a table
    u, theta = sp.symbols('u theta')
    diffop = DiffOp(u, independent_variable='t',
                    num_terms_Taylor_series=5)
    D1u = diffop.D(1)   # symbol for du/dt
    D2u = diffop.D(2)   # symbol for d^2u/dt^2
    print 'R Dt:', diffop['Dt'] - D1u
    print 'R Dtm:', diffop['Dtm'] - D1u
    print 'R Dtp:', diffop['Dtp'] - D1u
    print 'R barDt:', diffop['barDt'] - D1u
    print 'R DtDt:', diffop['DtDt'] - D2u
    print 'R weighted arithmetic mean:', diffop['weighted_arithmetic_mean'] - u
    print 'R arithmetic mean:', diffop['weighted_arithmetic_mean'].subs(theta, sp.Rational(1,2)) - u
    print 'R geometric mean:', diffop['geometric_mean'] - u
    dt = diffop.h
    print 'R harmonic mean:', (diffop['harmonic_mean'] - u).\
          series(dt, 0, 3).as_leading_term(dt)

if __name__ == '__main__':
    truncation_errors()
