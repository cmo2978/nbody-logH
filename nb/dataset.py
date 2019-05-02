import pickle
import numpy as np

class Dataset:
    pass

def new_file(sim, filename):
    f = open(filename, 'w+b')
    pickle.dump(sim, f)
    return open(filename, 'a+b')

def save_data(sim, f):
    s = "%s %s %s" % (sim.t, sim.T, sim.U)
    for p in sim.particles:
        s += " %s %s %s %s %s %s %s %s %s %s %s %s" % (p.r[0], p.r[1], p.r[2],
                                                       p.p[0], p.p[1], p.r[2],
                                                       p.a, p.e, p.i,
                                                       p.Omega, p.omega, p.theta)
    pickle.dump(s, f)

def load(filename):
    f = open(filename, 'r+b')
    sim = pickle.load(f)
    data = []
    try:
        while True:
            data.append(pickle.load(f).split())
    except EOFError:
        pass
    data = np.array(data, dtype=float)

    sim.data = Dataset()
    sim.data.t = data[:,0]
    sim.data.T = data[:,1]
    sim.data.U = data[:,2]
    i = 0
    for p in sim.particles:
        p.data = Dataset()
        p.data.r = np.vstack([data[:,12*i+3], data[:,12*i+4], data[:,12*i+5]]).T
        p.data.p = np.vstack([data[:,12*i+6], data[:,12*i+7], data[:,12*i+8]]).T
        p.data.a = data[:,12*i+9]
        p.data.e = data[:,12*i+10]
        p.data.i = data[:,12*i+11]
        p.data.Omega = data[:,12*i+12]
        p.data.omega = data[:,12*i+13]
        p.data.theta = data[:,12*i+14]
        i += 1

    f.close()
    return sim
