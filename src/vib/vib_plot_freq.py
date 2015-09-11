import numpy as np

def tilde_w(w, dt):
    return (2./dt)*np.arcsin(w*dt/2.)

def plot_frequency_approximations():
    w = 1  # relevant value in a scaled problem
    stability_limit = 2./w
    dt = np.linspace(0.2, stability_limit, 111)  # time steps
    series_approx = w + (1./24)*dt**2*w**3
    P = 2*np.pi/w  # one period
    num_timesteps_per_period = P/dt  # more instructive
    import scitools.std as plt
    plt.plot(num_timesteps_per_period, tilde_w(w, dt), 'r-',
             num_timesteps_per_period, series_approx, 'b--',
             legend=('exact discrete frequency', '2nd-order expansion'),
             xlabel='no of time steps per period',
             ylabel='numerical frequency')
    plt.savefig('discrete_freq.png')
    plt.savefig('discrete_freq.pdf')

plot_frequency_approximations()
