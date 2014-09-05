def find_k(u0, u1, Ts, dt):
    return 1./dt*((u0-u1)/(u1-Ts)) #(1-(u1-Ts)/float(u0-Ts))/dt

k = find_k(26.7, 25.8, 20, 3600)
print 'k=%g' % k

t = 0
dt = 1
if abs(k*dt) > 1:
    print("To large time step")

T = 37
Ts = 20
while T > 25.8:
    T = Ts+(T-Ts)*(1-k*dt)
    t+= dt

minutes, seconds = divmod(t, 60)
hours, minutes = divmod(minutes, 60)
print("The death occurred %d hours, %d minutes, and %g seconds before 3am" % \
        (hours, minutes, seconds))
