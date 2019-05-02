import nb
import numpy as np
from numpy import pi, sqrt
import matplotlib.pyplot as plt

hdiv = 150
file = "kozai_h"+str(hdiv)+"_gbs"

m0 = 1.
m1, a1, e1, i1, omega1, theta1 = 1e-6,   1.,  0.,   0., 0.,   0.
m2, a2, e2, i2, omega2, theta2 = 1e-7, .002, .05, pi/3, 0., pi/2

sim = nb.Simulation()
p0 = sim.add(m=m0)
p1 = sim.add(m=m1, a=a1, e=e1, i=i1, omega=omega1, theta=theta1, primary=p0)
p2 = sim.add(m=m2, a=a2, e=e2, i=i2, omega=omega2, theta=theta2, primary=p1)

P1 = 2 * pi * sqrt(p1.a**3 / (sim.G * (p0.m + p1.m)))
P2 = 2 * pi * sqrt(p2.a**3 / (sim.G * (p1.m + p2.m)))
E = sim.T - sim.U

t = 50 * P1
h = P2 / hdiv * abs(E) * 2 * 4

# sim.run(t, h, "../dat/"+file+".dat")
sim = nb.load("../dat/"+file+".dat")
p = sim.particles

fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, dpi=300, figsize=(7,5))
ax1.plot(sim.data.t, p[2].data.e, '-o', lw=.1, ms=.2)
ax2.plot(sim.data.t, p[2].data.i, '-o', lw=.1, ms=.2)
ax3.plot(sim.data.t, (sim.data.T-sim.data.U-sim.E0)/sim.E0, '-o', lw=.1, ms=.2)
ax1.set_ylabel('e')
ax2.set_ylabel('i')
ax3.set_ylabel('rel.err.E')
ax3.set_xlabel('t')
ax3.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
fig.subplots_adjust(wspace=0)
fig.suptitle(file)
fig.savefig("tests/"+file+".png", bbox_inches='tight')
