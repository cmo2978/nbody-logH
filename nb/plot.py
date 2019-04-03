import matplotlib.pyplot as plt

def plot(sim):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14,5), dpi=300)

    colors = ax1._get_lines.prop_cycler
    for p in sim.particles:
        c = next(colors)['color']
        ax1.plot(p.data.r[:,0], p.data.r[:,1], 'o', ms=.5, c=c)
        ax1.plot(p.r[0], p.r[1], 'o', ms=3., c=c, label=p.label)
    ax1.axis('equal')
    ax1.legend(loc='upper right', fontsize='x-small')

    ax2.axhline(y=0., linestyle='-', c='gray', lw=.5)
    ax2.plot(sim.data.t, (sim.data.T - sim.data.U + sim.B) / sim.B, '-o', lw=.5, ms=.7)

    fig.tight_layout()

    return fig
