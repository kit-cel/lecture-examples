#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
WT-Übung 2: Aufgabe 8

Ein Schimpanse hat zwei Urnen vor sich: Urne 1 enthält drei weiße und zwei
schwarze, Urne 2 eine weiße, zwei grüne und zwei rote Kugeln. Über das
Verhalten des Schimpansen ist bekannt, dass er mit der Wahrscheinlichkeit 0,7
in die erste und mit der Wahrscheinlichkeit 0,3 in die zweite Urne greift.

Der Schimpanse darf nun solange Kugeln (ohne Zurucklegen) ziehen bis er eine
rote Kugel wählt. Wie groß ist die Wahrscheinlichkeit, dass er maximal drei
Kugeln zieht?
"""

from __future__ import division, print_function
import platform
import time
from random import random, choice

## Aufgabenstellung
# ##############################################################################

N = 100000  # Anzahl der Durchläufe

# Kugeln
if platform.system() == "Linux":
    WEISS   = "\033[1;37m●\033[0m"
    SCHWARZ = "\033[1;30m●\033[0m"
    GRUEN   = "\033[1;32m●\033[0m"
    ROT     = "\033[1;31m●\033[0m"
else:
    WEISS, SCHWARZ, GRUEN, ROT = "WBGR"

# Urnen
k_urne1 = [WEISS, WEISS, WEISS, SCHWARZ, SCHWARZ]
k_urne2 = [WEISS, GRUEN, GRUEN, ROT, ROT]

# Wkeiten
p_urne1 = 0.7
p_urne2 = 0.3  # 1 - p_urne1

# Ausgabe der Parameter
print(__doc__.strip())
print("─" * 80)
print("Urne 1:", " ".join(k_urne1), "\t", "P(Urne 1) =", p_urne1)
print("Urne 2:", " ".join(k_urne2), "\t", "P(Urne 2) =", p_urne2)
print("─" * 80)
# format strings
str_len_N = len(str(N))
output_fmt = ("{runde:" + str(str_len_N) + "d}  {verlauf:}",  # normal
              " " * str_len_N + "  ...",                      # spacer
              "{:" + str(str_len_N - 1) + "}")                # results


## Spielverlauf
# ##############################################################################

print(" " * (str_len_N - 6), "Spiel  gezogene Kugeln")
# absolute Häufigkeit der unterschiedlich langen Spiele
haeufigkeit = [0] * (1 + len(k_urne1) + len(k_urne2))
for n in range(N):
    # Spiel initiallisieren
    urne1, urne2 = list(k_urne1), list(k_urne2)  # copy
    runden, kugel, verlauf = 0, None, []
    # Spiel durchführen:
    # - Ziehen ohne Zurücklegen aus einer der Urnen
    # - Abbruch wenn eine rote Kugel gezogen wird
    while kugel != ROT:
        runden += 1
        # Urne auswählen (Urne 1 kann leer sein)
        urne = urne1 if random() < p_urne1 and len(urne1) else urne2
        # Kugel aus der ausgewählten Urne ziehen
        kugel = choice(urne)
        # gezogene Kugel aus der Urne entfernen
        del urne[urne.index(kugel)]
        # gezogene Kugel merken
        verlauf.append(kugel)

    # Anzahl der Runden merken
    haeufigkeit[runden] += 1

    # Ausgabe des Durchlaufs
    if n < 10 or n >= N - 5:
        print(output_fmt[0].format(runde=n+1, verlauf=" ".join(verlauf)))
        time.sleep(0.5)
    elif n == 10:
        print(output_fmt[1])

## Ergebnisse
# ##############################################################################

print("─" * 80)
print("Absolute Häufigkeit von 'genau k Runden' bei {} Durchläufen:".format(N))
print("  Runden:", "  ".join(output_fmt[2].format(runden)
    for runden in range(1, len(haeufigkeit))))
print("  Anzahl:", "  ".join(output_fmt[2].format(h)
    for h in haeufigkeit[1:]))
print()

print("Relative Häufigkeit für 'maximal 3 Runden':")
print("  P(Runden < 4) = {}/{} = {:.2}".format(
    sum(haeufigkeit[:4]), N, sum(haeufigkeit[:4]) / N))
