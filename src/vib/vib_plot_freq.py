from numpy import linspace, arcsin as asin, pi

def tilde_w(w, dt):
    return (2./dt)*asin(w*dt/2.)

def plot_frequency_approximations():
    w = 1  # relevant value in a scaled problem
    stability_limit = 2./w
    dt = linspace(0.2, stability_limit, 111)  # time steps
    series_approx = w + (1./24)*dt**2*w**3
    P = 2*pi/w  # one period
    num_timesteps_per_period = P/dt  # more instructive
    from scitools.std import plot, savefig
    plot(num_timesteps_per_period, tilde_w(w, dt), 'r-',
         num_timesteps_per_period, series_approx, 'b--',
         legend=('exact discrete frequency', '2nd-order expansion'),
         xlabel='no of time steps per period',
         ylabel='numerical frequency')
    savefig('discrete_freq.png')
    savefig('discrete_freq.eps')
    savefig('discrete_freq.pdf')

plot_frequency_approximations()
