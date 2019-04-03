import numpy as np
from numpy.linalg import norm
from .coords import convert_orb_to_cart, get_com

class Particle:
    def __init__(self, m, sim=None, label=None,
                 x=None, y=None, z=None, vx=None, vy=None, vz=None,
                 primary=None, a=None, e=None, i=None, Omega=None, omega=None, theta=None):
        """
        Initialize a particle object using either cartesian coordinates or orbital elements.
        Use .simulation.Simulation.add() to create a particle and add it directly to a simulation.

        m           : Mass
        sim         : Reference to parent simulation
        label       : Label string for plots
        x, y, z     : Positions
        vx, vy, vz  : Velocities
        primary     : Primary object the particle is bound to
        a, e, i, Omega, omega, theta    : Semi-major axis, eccentricity, inclination, long. of asc. node, arg. of periapsis, true anomaly
        """
        self.sim = sim
        self.label = label
        self.m = m

        cartesian = [x,y,z,vx,vy,vz] != 6*[None]
        orbital = [primary,a,e,i,Omega,omega,theta] != 7*[None]
        if cartesian and orbital:
            raise ValueError("Cannot pass both cartesian and orbital coordinates.")

        self.primary = primary
        if orbital:
            if a is None:
                raise ValueError("Missing parameter a.")
            e = e or 0.
            if e == 1.:
                raise ValueError("e cannot be exactly 1.")
            if e > 1. and a > 0.:
                raise ValueError("e > 1 requires a < 0.")
            if e < 1. and a < 0.:
                raise ValueError("e < 1 requires a > 0.")
            if primary and primary.sim != self.sim:
                raise ValueError("Primary object must be part of the simulation.")
            self.a = a
            self.e = e
            self.i = i or 0.
            self.Omega = Omega or 0.
            self.omega = omega or 0.
            self.theta = theta or 0.
            convert_orb_to_cart(self, primary)
        else:
            self.r = np.array([x or 0., y or 0., z or 0.])
            self.p = self.m * np.array([vx or 0., vy or 0., vz or 0.])

        self.data = None

    @property
    def F(self):
        """
        Total gravitational force acting on the particle.
        """
        F = 0.
        for p in self.sim.particles:
            if p != self:
                F += self.m * p.m * (p.r - self.r) / norm(p.r - self.r)**3
        return self.sim.G * F
