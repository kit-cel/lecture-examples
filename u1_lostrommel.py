#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
WT-Übung 1: Aufgabe 5

100 Lose, davon 20 Gewinne, davon ein Hauptgewinn.
Es werden vier Lose gezogen.

A = {Hauptgewinn wird gezogen}
B = {mindestens ein Gewinn}
C = {höchstens ein Gewinn}
D = {genau ein Nebengewinn}
"""

from __future__ import division, print_function
import time
from collections import defaultdict
from itertools import product
from random import sample


## Aufgabenstellung
# #############################################################################

print(__doc__)

# Lose und Lostrommel
NIETE, NEBENGEWINN, HAUPTGEWINN = range(3)  # Zahlenwerte egal
LOSTROMMEL = [HAUPTGEWINN] + 19 * [NEBENGEWINN] + 80 * [NIETE]
assert len(LOSTROMMEL) == 100

# Indikatorfunktionen für Ereignisse
event_checks = dict(
    A = lambda: HAUPTGEWINN in draw,
    B = lambda: draw != 4 * [NIETE],
    C = lambda: sum(los != NIETE for los in draw) <= 1,
    D = lambda: sorted(draw) == [NIETE, NIETE, NIETE, NEBENGEWINN]
)

# Anzahl der Versuche
N = 100000
# Zähler für absolute Häufigkeiten der Ereignisse bei N Versuchen
event_counter = defaultdict(int)


## Simulation
# #############################################################################

print(" n     Lose     Ereignisse")
print("───────────────────────────")

format_draw = lambda: ' '.join(
    {HAUPTGEWINN:'H', NEBENGEWINN:'G', NIETE:'N'}[los] for los in draw
)
format_events = lambda: ' '.join(
    event if event in events else '_' for event in sorted(event_checks)
)

for n in range(N):  # Versuch n = 0 ... (N-1)
    draw = sample(LOSTROMMEL, 4)  # zufällige Auswahl von 4 Losen

    # Prüfe ob Ereignis A, B, C, D für Ergebnis 'draw' eintritt
    events = []
    for name, event_check in event_checks.items():
        if event_check():  # Ereignis eingetreten?
            event_counter[name] += 1  # Zähler inkrementieren
            events.append(name)  # für Ausgabe merken

    # Ausgabe
    if n < 10:
        print(" {:d}    {}    {}".format(n, format_draw(), format_events()))
        time.sleep(0.1)
    elif n == 10:
        print(" ...")

print("───────────────────────────")


## Ergebnisse
# #############################################################################

print("Relative Häufigkeiten:")
for name, count in sorted(event_counter.items()):
    print("h({0:}) = {1:5d}/{2:} = {3:.03f}".format(
        name, count, N, count / N
    ))
