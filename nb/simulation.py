from .particle import Particle
from .coords import move_to_com, convert_cart_to_orb
import numpy as np
from numpy.linalg import norm
import sys

class Simulation:
    def __init__(self, G=1., com=True):
        """
        Initialize a simulation with options
        G   : Gravitational constant
        com : Use center of mass frame (bool)
        """
        self.G = G
        self.com = com
        self.particles = []
        self.N = 0  # Number of particles
        self.m = 0. # Total mass
        self.t = 0. # Current time
        self.data = None

    def add(self, **kwargs):
        """
        Create a particle object and add it to the simulation.
        Returns the particle.
        """
        p = Particle(sim=self, **kwargs)
        self.particles.append(p)
        self.N += 1
        self.m += p.m
        return p

    @property
    def T(self):
        """
        Returns the kinetic energy of the system
        """
        return sum([0.5 * norm(p.p)**2 / p.m for p in self.particles])

    @property
    def U(self):
        """
        Returns the potential energy of the system
        """
        U = 0.
        for i in range(self.N-1):
            for j in range(i+1, self.N):
                M = self.particles[i].m * self.particles[j].m
                R = norm(self.particles[i].r - self.particles[j].r)
                U += M / R
        return self.G * U

    def run(self, t, h):
        """
        Start the simulation, integrating until time t, time-step h
        """
        if self.data is None:
            if self.com:
                move_to_com(self)

            for p in self.particles:
                convert_cart_to_orb(p)

            self.data = Dataset(self)
            self.B = self.U - self.T    # Binding energy = - Initial total energy

        # Leapfrog:
        while self.t < t:

            dt = h / (self.T + self.B)
            self.t += 0.5 * dt
            for p in self.particles:
                p.r += 0.5 * dt * p.p / p.m

            for p in self.particles:
                p.p += h / self.U * p.F

            dt = h / (self.T + self.B)
            self.t += 0.5 * dt
            for p in self.particles:
                p.r += 0.5 * dt * p.p / p.m

            for p in self.particles:
                convert_cart_to_orb(p)

            self.data.save()

            sys.stdout.write("\r%.2f%%, %d" %(self.t/t*100, len(self.data.t)))
            sys.stdout.flush()

class Dataset:
    def __init__(self, parent):
        if isinstance(parent, Simulation):
            self._arrays = []
            self._lists = ['t', 'T', 'U']
            for p in parent.particles:
                p.data = Dataset(p)
        elif isinstance(parent, Particle):
            self._arrays = ['r', 'p']
            self._lists = ['a','e','i','Omega','omega','theta']
        for key in self._arrays + self._lists:
            self.__dict__[key] = np.array([parent.__getattribute__(key)])

        self.parent = parent

    def save(self):
        if isinstance(self.parent, Simulation):
            for p in self.parent.particles:
                p.data.save()
        for key in self._arrays:
            val = self.parent.__getattribute__(key)
            self.__dict__[key] = np.vstack((self.__dict__[key], val))
        for key in self._lists:
            val = self.parent.__getattribute__(key)
            self.__dict__[key] = np.hstack((self.__dict__[key], val))

# class Dataset:
#     def __init__(self, parent):
#         if isinstance(parent, Simulation):
#             self.t = np.array([ parent.t ])
#             self.T = np.array([ parent.T ])
#             self.U = np.array([ parent.U ])
#             for p in parent.particles:
#                 p.data = Dataset(p)
#         elif isinstance(parent, Particle):
#             self.r = np.array([ parent.r ])
#             self.p = np.array([ parent.p ])
#             self.a = np.array([ parent.a ])
#             self.e = np.array([ parent.e ])
#             self.i = np.array([ parent.i ])
#             self.Omega = np.array([ parent.Omega ])
#             self.omega = np.array([ parent.omega ])
#             self.theta = np.array([ parent.theta ])
#
#         self.parent = parent
#
#     def save(self):
#         if isinstance(self.parent, Simulation):
#             self.t = np.hstack((self.t, self.parent.t))
#             self.T = np.hstack((self.T, self.parent.T))
#             self.U = np.hstack((self.U, self.parent.U))
#             for p in self.parent.particles:
#                 p.data.save()
#         elif isinstance(self.parent, Particle):
#             self.r = np.vstack((self.r, self.parent.r))
#             self.p = np.vstack((self.p, self.parent.p))
#             self.a = np.hstack((self.a, self.parent.a))
#             self.e = np.hstack((self.e, self.parent.e))
#             self.i = np.hstack((self.i, self.parent.i))
#             self.Omega = np.hstack((self.Omega, self.parent.Omega))
#             self.omega = np.hstack((self.omega, self.parent.omega))
#             self.theta = np.hstack((self.theta, self.parent.theta))
