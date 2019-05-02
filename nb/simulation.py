from .particle import Particle
from .coords import move_to_com, convert_cart_to_orb
from .integrator import GBS_step
from .dataset import new_file, save_data
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

    def run(self, t, h, filename=None, gbs=True):
        """
        Start the simulation, integrating until time t, using time-step H
        """
        if self.com:
            move_to_com(self)
        for p in self.particles:
            convert_cart_to_orb(p)
        self.E0 = self.T - self.U

        if filename: f = new_file(self, filename)
        while self.t < t:
            GBS_step(self, h)
            # Leapfrog
            # self.t += 0.5 * h / (self.T - self.E0)
            # for p in self.particles:
            #     p.r += 0.5 * h / (self.T - self.E0) * p.p / p.m
            #
            # for p in self.particles:
            #     p.p += h / self.U * p.F
            #
            # self.t += 0.5 * h / (self.T - self.E0)
            # for p in self.particles:
            #     p.r += 0.5 * h / (self.T - self.E0) * p.p / p.m

            for p in self.particles:
                convert_cart_to_orb(p)

            if filename: save_data(self, f)
            sys.stdout.write("\r%.2f%%" %(self.t/t*100))
            sys.stdout.flush()
        if filename: f.close()
