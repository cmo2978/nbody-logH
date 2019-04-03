import nb
import numpy as np
import matplotlib.pyplot as plt

sim = nb.Simulation(G=1.)
p0 = sim.add(m=3., x= 1., y= 3.)
p1 = sim.add(m=4., x=-2., y=-1.)
p2 = sim.add(m=5., x= 1., y=-1.)

h = 1e-4 *sim.m**2.5 / (sim.N * np.sqrt(abs(sim.T-sim.U)))

sim.run(10.,h)

fig = nb.plot(sim)
fig.savefig("t345.png")
