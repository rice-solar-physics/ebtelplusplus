"""
Second example, ion heating, multiple pulses
"""
import os

import numpy as np

from plot import make_figure
from util import run_ebtel, read_xml


top_dir = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))


if __name__ == '__main__':
    # Configure and run ebtel++
    base = read_xml(os.path.join(top_dir, 'config', 'ebtel.example.cfg.xml'))
    base['calculate_dem'] = True
    base['total_time'] = 10000.0
    base['heating']['partition'] = 0.0
    events = []
    for i in range(5):
        events.append({'event': {
            'rise_start': i*2000.0,
            'rise_end': i*2000.0+100.0,
            'decay_start': i*2000.0+100.0,
            'decay_end': i*2000.0+200.0,
            'magnitude': np.random.uniform(0.001, 0.1)
        }})
    base['heating']['events'] = events
    results = run_ebtel(base, top_dir)
    # Time average DEM
    results['dem_total'] = np.average(
        results['dem_tr']+results['dem_corona'],
        axis=0,
        weights=np.gradient(results['time']))
    results['dem_tr'] = np.average(
        results['dem_tr'],
        axis=0,
        weights=np.gradient(results['time']))
    results['dem_corona'] = np.average(
        results['dem_corona'],
        axis=0,
        weights=np.gradient(results['time']))
    # Plot
    figure_filename = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), 'ex2.png')
    make_figure(results, figure_filename)
