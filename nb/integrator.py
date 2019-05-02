import numpy as np
from numpy.linalg import norm

n = 2 * (np.arange(8) + 1)
kmax = 8

def GBS_step(sim, H, tol=1):
    T = []
    W_old, t_old = save_state(sim)  # Save the state vector at start

    # Set T00 using leapfrog with n[0]=2 substeps
    leapfrog(sim, H, n[0])
    W, t = save_state(sim)
    set_Tk(sim, T, 0, W)
    set_state(sim, W_old, t_old)

    k = 1
    while k < kmax:
        # Set kth row of T using leapfrog with n[k] substeps
        leapfrog(sim, H, n[k])
        W, t = save_state(sim)
        set_Tk(sim, T, k, W)

        err = norm( T[k][k] - T[k][k-1] ) / norm(T[k][k])   # Relative error check. ??
        if err < tol:
            break
        set_state(sim, W_old, t_old)
        k += 1

    if k == kmax:
        print("\nk = kmax. err = %e. H too large?" % err)
        GBS_step(sim, H/2, tol)
        GBS_step(sim, H/2, tol)
    else:
        set_state(sim, T[k][k], t)

def leapfrog(sim, H, n):
    '''
    Leapfrog with step-size H, # of substeps n.
    '''
    h = H / n

    sim.t += 0.5 * h / (sim.T - sim.E0)
    for p in sim.particles:
        p.r += 0.5 * h / (sim.T - sim.E0) * p.p / p.m
    for p in sim.particles:
        p.p += h / sim.U * p.F

    for _ in range(n-1):
        sim.t += h / (sim.T - sim.E0)
        for p in sim.particles:
            p.r += h / (sim.T - sim.E0) * p.p / p.m
        for p in sim.particles:
            p.p += h / sim.U * p.F

    sim.t += 0.5 * h / (sim.T - sim.E0)
    for p in sim.particles:
        p.r += 0.5 * h / (sim.T - sim.E0) * p.p / p.m

def set_Tk(sim, T, k, W):
    T.append([W])
    for j in range(k):
        T[k].append( T[k][j] + (T[k][j] - T[k-1][j]) / ( (n[k]/n[k-j-1])**2 - 1 ) )

def save_state(sim):
    W = []
    for p in sim.particles:
        W.append(p.r)
        W.append(p.p)
    t = sim.t
    return np.array(W), t

def set_state(sim, W, t):
    i = 0
    for p in sim.particles:
        p.r = W[2*i]
        p.p = W[2*i+1]
        i += 1
    sim.t = t
