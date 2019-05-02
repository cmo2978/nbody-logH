import numpy as np
from numpy import sin, cos, sqrt, cross, dot, pi, arccos
from numpy.linalg import norm

def convert_orb_to_cart(P):
    """
    Convert orbital elements to cartesian state vectors
    """
    primary = P.primary or get_com([p for p in P.sim.particles if p != P])
    mu = P.sim.G * (P.m + primary.m)
    p = P.a * (1. - P.e**2)
    r = p / (1. + P.e * cos(P.theta))
    v = sqrt(mu * (2. / r - 1. / P.a))
    sing = sqrt(mu / p) * P.e / v * sin(P.theta)
    cosg = sqrt(1. - sing**2)
    u = P.omega + P.theta

    x = cos(P.Omega) * cos(u) - sin(P.Omega) * cos(P.i) * sin(u)
    y = sin(P.Omega) * cos(u) + cos(P.Omega) * cos(P.i) * sin(u)
    z = sin(P.i) * sin(u)

    vx = v * (x * sing - cosg * (cos(P.Omega) * sin(u) + sin(P.Omega) * cos(P.i) * cos(u)))
    vy = v * (y * sing - cosg * (sin(P.Omega) * sin(u) - cos(P.Omega) * cos(P.i) * cos(u)))
    vz = v * (z * sing + cosg * sin(P.i) * cos(u))

    P.r = r * np.array([x, y, z]) + primary.r
    P.p = P.m * (np.array([vx, vy, vz]) + primary.p/primary.m)

def convert_cart_to_orb(P):   # TODO: singularity checks
    """
    Convert cartesian state vectors to orbital elements
    """
    primary = P.primary or get_com([p for p in P.sim.particles if p != P])
    mu = P.sim.G * (P.m + primary.m)
    r = P.r - primary.r
    v = P.p / P.m - primary.p / primary.m

    h = cross(r, v)
    n = cross([0,0,1], h)

    ev = 1 / mu * ((norm(v)**2 - mu / norm(r)) * r - dot(r, v) * v)
    e = norm(ev)

    E = norm(v)**2 / 2 - mu / norm(r)
    a = - mu / (2 * E)

    i = arccos(h[2] / norm(h))

    # Singularity checks
    # Non-inclined orbit
    if abs(i) < 1e-15:
        Omega = 0.
        # Circular orbit
        if abs(e) < 1e-15:
            omega = 0.
        else:
            omega = arccos(ev[0]/e)
    else:
        Omega = arccos(n[0] / norm(n))
        if n[1] < 0: Omega = 2 * pi - Omega
        omega = arccos(dot(n, ev) / (norm(n) * norm(ev)))

    # Circular orbit
    if abs(e) < 1e-15:
        # Non-inclined orbit
        if abs(i) < 1e-15:
            theta = arccos(r[0] / norm(r))
            if v[0] > 0: theta = 2 * pi - theta
        else:
            theta = arccos(dot(n, r) / (norm(n) * norm(r)))
            if dot(n, v) > 0: theta = 2 * pi - theta
    else:
        theta = arccos(dot(ev, r) / (e * norm(r)))
        if dot(r, v) < 0: theta = 2 * pi - theta

        if ev[2] < 0:
            omega = 2 * pi - omega

    P.a, P.e, P.i, P.Omega, P.omega, P.theta = a, e, i, Omega, omega, theta


def get_com(particles):
    """
    Return the center of mass of the list of particles
    """
    class Com:
        m = 0.
        r = np.zeros(3)
        p = np.zeros(3)

    com = Com()
    for p in particles:
        com.m += p.m
        com.r += p.m * p.r
        com.p += p.p
    com.r /= com.m
    return com

def move_to_com(sim):
    """
    Move all particles in sim to center of mass frame
    """
    com = get_com(sim.particles)
    for p in sim.particles:
        p.r -= com.r
        p.p -= com.p / com.m * p.m
