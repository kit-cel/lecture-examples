#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
WT-Übung 1: Aufgabe 4

Elementarereignisse für Ereignis A 'Strom fließt von E nach A'
"""

from __future__ import division, print_function
import time
from itertools import product


## Aufgabenstellung
# #############################################################################

print(__doc__)

# Indikatorfunktionen für Ereignis A
check_event_A = lambda: ((r[0] or r[1]) and r[2]) or (r[3] and r[4])

# Funktion zur Ausgabe des Schaltwerks
network = lambda relais, pre='', post='': """\
{i:}        ┌─{r[0]:}─┐
{s:}  E ──┬─┤    ├──────{r[2]:}─────┬── A {e:}
{i:}      │ └─{r[1]:}─┘             │
{i:}      └─────────{r[3]:}─────{r[4]:}──┘
""".format(s=pre, r=relais, e=post, i=' '*len(pre))


## Simulation
# #############################################################################

print("Fall             Schaltwerk           Geschlossen  Anzahl")
print("─────────────────────────────────────────────────────────")

event_A_counter = 0
for n, r in enumerate(product(range(2), repeat=5)):  # alle möglichen Ergebnisse

    # Prüfe ob Ereignis A für Ereignis r eintritt
    event_A = check_event_A()
    if event_A:
        event_A_counter += 1  # Zähle Elementarereignisse in Ereignis A

    # Ausgabe
    print(network(
        relais=['──' if ri else '┘└' for ri in r],
        pre=" {:2d}  ".format(1+n),
        post="     {}         {:2d}".format(
            "✔" if event_A else "✘", event_A_counter
        )
    ))
    time.sleep(0.1)  # Animation


## Ergebnisse
# #############################################################################

print("─────────────────────────────────────────────────────────")
print("Wahrscheinlichkeit nach LaPlace:")
print("  P(A) = |A|/|Ω| = {}/{} = {:.5}".format(
    event_A_counter, 2 ** 5, event_A_counter / 2 ** 5))
