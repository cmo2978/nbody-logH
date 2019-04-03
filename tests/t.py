import nb
import numpy as np
from numpy import pi, sqrt
import matplotlib.pyplot as plt

def rand(): return 2 * pi * np.random.rand()

m0 = 1.
m1, a1, e1, omega1, theta1 = 1e-5, 1., .5, 0., 0.
m2, a2, e2, omega2, theta2 = 1e-5, 1., .5, pi, 0.

sim = nb.Simulation()
p0 = sim.add(m=m0)
p1 = sim.add(m=m1, a=a1, e=e1, omega=omega1, theta=theta1, primary=p0)
p2 = sim.add(m=m2, a=a2, e=e2, omega=omega2, theta=theta2, primary=p0)

p0.label = r"$m=%s$" % (m0)
p1.label = r"$m=%s$, $a=%.2f$, $e=%.2f$, $\omega=%.2f$, $\theta=%.2f$" % (m1,a1,e1,omega1,theta1)
p2.label = r"$m=%s$, $a=%.2f$, $e=%.2f$, $\omega=%.2f$, $\theta=%.2f$" % (m2,a2,e2,omega2,theta2)

P  = 2 * pi * sqrt(p1.a**3 / (sim.G * (p0.m + p1.m)))
P2 = 2 * pi * sqrt(p2.a**3 / (sim.G * (p0.m + p2.m)))

t = 100 * P2
h = P / (40 * p0.m / p1.m)

sim.run(t, h)

fig = nb.plot(sim)
fig.savefig("3b4_3.png")
