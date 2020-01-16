"""
Visualize ebtel++ example results
"""
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import MaxNLocator
import seaborn


def make_figure(results, filename):
    seaborn.set_context('notebook', font_scale=1.1)
    seaborn.set(font='serif')
    seaborn.set_style("white", {
        "font.family": "serif",
        "font.serif": ["Times", "Palatino", "serif"]
    })

    fig = plt.figure(figsize=(9, 6))
    gs = gridspec.GridSpec(3, 2)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)
    ax3 = fig.add_subplot(gs[2, 0], sharex=ax1)
    ax4 = fig.add_subplot(gs[:, 1])

    ax1.plot(results['time'], results['heat'])
    ax2.plot(results['time'], results['electron_temperature']/1e+6,
             color=seaborn.color_palette()[0],
             label=r'$T_e$')
    ax2.plot(results['time'], results['ion_temperature']/1e+6,
             color=seaborn.color_palette()[2],
             label=r'$T_i$')
    ax3.plot(results['time'], results['density']/1e+8)
    ax4.plot(results['dem_temperature'],
             results['dem_tr'],
             color=seaborn.color_palette()[0],
             label=r'$\mathrm{DEM}_{\mathrm{TR}}$')
    ax4.plot(results['dem_temperature'],
             results['dem_corona'],
             color=seaborn.color_palette()[2],
             label=r'$\mathrm{DEM}_{\mathrm{C}}$')
    ax4.plot(results['dem_temperature'],
             results['dem_total'],
             color=seaborn.color_palette()[1],
             label=r'$\mathrm{DEM}_{\mathrm{total}}$')

    seaborn.despine()
    seaborn.despine(ax=ax1, bottom=True)
    seaborn.despine(ax=ax2, bottom=True)
    ax1.set_ylabel(r'$H$ (erg cm$^{-3}$ s$^{-1}$)')
    ax1.yaxis.set_major_locator(MaxNLocator(prune='lower', nbins=5))
    ax1.tick_params(axis='x', labelbottom=False)
    ax2.set_ylabel(r'$T$ (MK)')
    ax2.yaxis.set_major_locator(MaxNLocator(prune='lower', nbins=5))
    ax2.tick_params(axis='x', labelbottom=False)
    ax3.set_ylabel(r'$n$ ($10^8$ cm$^{-3}$)')
    ax3.set_xlabel(r'$t$ (s)')
    ax3.yaxis.set_major_locator(MaxNLocator(nbins=5))
    ax3.set_xlim(results['time'][[0, -1]])
    ax2.legend(loc='best')
    ax4.set_xlabel(r'$T$ (K)')
    ax4.set_ylabel(r'$\mathrm{DEM}$ (cm$^{-5}$ K$^{-1}$)')
    ax4.set_xlim([10**(4.5), 10**(7.5)])
    ax4.set_ylim([10**(20.0), 10**(23.5)])
    ax4.set_xscale('log')
    ax4.set_yscale('log')
    ax4.legend(loc='best')

    plt.subplots_adjust(hspace=0.0, wspace=0.3)
    plt.savefig(filename)
