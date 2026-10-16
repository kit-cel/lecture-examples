Beispielskripte zur Vorlesung "Nachrichtensysteme"
=================================================================================

**Note: An English version of this document can be found below.**

Dieses Repository enthält Beispiele und Simulationen, die in der Vorlesung "Nachrichtensysteme" der Fakultät für Elektrotechnik und Informationstechnik (ETIT) am Karlsruher Institut für Technologie (KIT) verwendet werden. Die zugehörigen Folien und Skripte sind in ILIAS verfügbar (nur für eingeschriebene Studierende).

Die Notebooks liegen im Ordner `vorlesung`:

- `1_quellencodierung`: Huffman-Codierung, gleichmäßige Quantisierung (Sinus und Audio)
- `2_kanalcodierung`: Blockcodierung, Faltungscodierung
- `4_kanaele`: Intersymbolinterferenz, frequenzselektive Kanäle
- `5_ofdm`: OFDM
- `6_mimo`: MIMO, Water-Filling


Python-Abhängigkeiten
---------------------
Es wird Python >= 3.12 benötigt. Die Abhängigkeiten sind in `pyproject.toml` deklariert und in `uv.lock` festgelegt:

- jupyterlab (Jupyter-Notebook-Umgebung)
- numpy
- scipy
- matplotlib
- ipympl
- ipywidgets
- scikit-learn

Verwendung der Python-Notebooks
-------------------------------
Wir empfehlen [uv](https://docs.astral.sh/uv/) zum Einrichten der Python-Umgebung. uv installiert automatisch eine passende Python-Version und alle festgelegten Abhängigkeiten.

1. uv installieren (siehe [Installationsanleitung](https://docs.astral.sh/uv/getting-started/installation/)), z.B. unter Linux/macOS mit `curl -LsSf https://astral.sh/uv/install.sh | sh`, unter Windows mit `powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex""`.
2. Dieses Repository klonen und in den Ordner `nsys/` wechseln.
3. Die Umgebung erstellen und alle Abhängigkeiten mit `uv sync` installieren. Dadurch wird gemäß `uv.lock` eine virtuelle Umgebung in `.venv` angelegt.
4. JupyterLab mit `uv run jupyter lab` starten und die Notebooks im Dateibrowser öffnen. Alternativ kann die `.venv`-Umgebung in einem Editor (z.B. VS Code) als Kernel ausgewählt werden oder die Umgebung mit `source .venv/bin/activate` aktiviert werden.

Einzelne Befehle können mit `uv run <Befehl>` in der Umgebung ausgeführt werden. Ein zusätzliches Paket kann mit `uv add <Paket>` hinzugefügt werden.


Lizenz und Wiederverwendung
---------------------------

Dieser Code steht unter der GPLv2-Lizenz.

Wenn Sie Teile dieses Codes in Ihrer eigenen Forschung oder Lehre verwenden und daraus Publikationen entstehen, zitieren Sie ihn bitte wie folgt:<br>
* L. Schmalen, "Nachrichtensysteme: Vorlesungsbeispiele," online verfügbar unter http://www.github.org/KIT-CEL/lecture-examples/, 2025

Bei Fragen, Kommentaren oder Verbesserungsvorschlägen zu diesem Code wenden Sie sich bitte an Laurent Schmalen (Laurent.Schmalen@kit.edu) oder erstellen Sie einen Pull Request.


---

English version
===============

Example scripts for the "Nachrichtensysteme" (Communication Systems) lecture
=================================================================================

This repository contains examples and simulations used in the lecture "Nachrichtensysteme" (Communication Systems, taught in German) of the faculty for electrical engineering and information technology (ETIT) at Karlsruhe Institute of Technology (KIT). The corresponding slides and lecture notes can be found in ILIAS (only available for enrolled students).

The notebooks are located in the `vorlesung` folder:

- `1_quellencodierung`: Huffman coding, uniform quantization (sine and audio)
- `2_kanalcodierung`: block coding, convolutional coding
- `4_kanaele`: intersymbol interference, frequency-selective channels
- `5_ofdm`: OFDM
- `6_mimo`: MIMO, water-filling


Python Dependencies
---------------------
Python >= 3.12 is required. The dependencies are declared in `pyproject.toml` and pinned in `uv.lock`:

- jupyterlab (Jupyter notebook environment)
- numpy
- scipy
- matplotlib
- ipympl
- ipywidgets
- scikit-learn

Usage of Python Notebooks
-------------------------
We recommend [uv](https://docs.astral.sh/uv/) to set up the Python environment. uv installs a suitable Python version and all pinned dependencies automatically.

1. Install uv (see the [installation guide](https://docs.astral.sh/uv/getting-started/installation/)), e.g., on Linux/macOS with `curl -LsSf https://astral.sh/uv/install.sh | sh`, on Windows with `powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"`.
2. Clone this repository and change into the `nsys/` folder.
3. Create the environment and install all dependencies with `uv sync`. This creates a virtual environment in `.venv` according to `uv.lock`.
4. Start JupyterLab with `uv run jupyter lab` and open the notebooks from the file browser. Alternatively, select the `.venv` environment as kernel in your editor (e.g., VS Code).

Single commands can be run inside the environment with `uv run <command>`. An additional package can be added with `uv add <package>`.


Licensing and Reuse
-------------------

This code is licensed under the GPLv2 license.

If you use, in any way, parts of this code in your own research or teaching that results in publications, please cite it as follows:<br>
* L. Schmalen, "Nachrichtensysteme: Lecture Examples," available online at http://www.github.org/KIT-CEL/lecture-examples/, 2025

In case you have questions, comments, or suggestions regarding this code, please contact Laurent Schmalen (Laurent.Schmalen@kit.edu) or propose a pull request.
