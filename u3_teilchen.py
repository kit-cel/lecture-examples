#!/usr/bin/env python
#  -*- coding: utf-8 -*-
"""
Wurfweite bei einer festen initialen Geschwindigkeit v und einem zufälligen
Winkel phi gegenüber der horizontalen Achse.
"""

from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as pp

# Parameter ###################################################################

num_samples = 1e5
num_bins = 30

v = 20  # m²/s (Beispielwert)
g = 9.81  # m²/s

# Wurfweite
c = v ** 2 / g
h = lambda phi: c * np.sin(2 * phi)

# Daten für Teilaufgaben
data = [[
    'a',

    #  analytische Lösung
    lambda d: 2 / np.pi / np.sqrt(c**2 - d**2),

    # N Realisierungen der ZV Phi
    np.random.uniform(0, np.pi/2, num_samples), 'b'
],
[
    "b",

    # analytische Lösung
    lambda d: d / c / np.sqrt(c**2 - d**2),

    # N Realisierungen der ZV Phi
    np.arccos(np.random.uniform(-1, 1, num_samples)) / 2, 'r',
]]

# Plot vorbereiten ############################################################
d_plot = np.linspace(0, c, 251)[:-1]
fig = pp.figure(figsize=(8,4.1))
ax = fig.add_axes((0.10, 0.12, 0.90, 0.88))

ax.set_xlim(0, c)
ax.set_ylim(0, data[0][1](d_plot[-4]))
ax.set_xticks((0, c/2, c), (u'0', u'v²/2g', u'v²/g'))
ax.set_xlabel("Wurfweite d")
ax.set_ylabel("rel. Haeufigkeit / Dichte")

# Berechnung / Plot ###########################################################
for section, f_D, phi_samples, color in data:

    # Bestimme für jede Realisierung die resultierende Weite
    d_samples = h(phi_samples)

    # relative Häufigkeit der diskretisierten Wurfweite
    h_D, bins = np.histogram(d_samples, bins=num_bins, range=(0, c), normed=True)

    # Plot
    ax.plot(d_plot, f_D(d_plot), color, label=section + ") analytisch")
    ax.plot((bins[:-1] + bins[1:]) / 2, h_D, color + 'o', ms=5,
            label=section + ") simulation")


# Ausgabe #####################################################################
ax.legend(loc="upper left")
pp.draw()
#pp.savefig("../figures/u3_teilchen.pdf")
pp.show()
