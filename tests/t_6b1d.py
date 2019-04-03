import nb
import numpy as np
import matplotlib.pyplot as plt

sim = nb.Simulation(com=False)
p0 = sim.add(m=1., x=-3.)
p1 = sim.add(m=2., x=-2.)
p2 = sim.add(m=3., x=-1.)
p3 = sim.add(m=4., x= 0.)
p4 = sim.add(m=5., x= 2.)
p5 = sim.add(m=6., x= 4.)

h = 1e-4 *sim.m**2.5 / (sim.N * np.sqrt(abs(sim.T-sim.U)))

sim.run(5.,h)

fig = nb.plot(sim)
fig.savefig("6b1d.png")
