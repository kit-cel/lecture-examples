Example scripts for "Channel Coding - Graph-based Codes" lecture
=================================================================================

This repository contains examples and simulations used in the lecture [Channel Coding - Graph-based Codes](https://www.youtube.com/watch?v=BAwHn95atvo&list=PLLDF9ieaJiSX-qDeYdrZNfpMo9UrUcu46) of the faculty for electrical engineering and information technology (ETIT) at Karlsruhe Institute of Technology (KIT). The corresponding slides and lecture notes can be found in ILIAS (only available for enrolled students). An old version of the lecture is freely available on Youtube (click on link above)


Python Dependencies
---------------------
Python >= 3.12 is required. The dependencies are declared in `pyproject.toml` and pinned in `uv.lock`:

- jupyterlab (Jupyter notebook environment)
- ipywidgets (interactive widgets)
- ipynb (importing functions from other notebooks)
- numpy <= 2.3.5
- scipy
- matplotlib
- sympy
- tqdm
- galois
- scikit-commpy
- cvxpy
- pyscipopt (SCIP solver, used via cvxpy)

Usage of Python Notebooks
-------------------------
We recommend [uv](https://docs.astral.sh/uv/) to set up the Python environment. uv installs a suitable Python version and all pinned dependencies automatically.

1. Install uv (see the [installation guide](https://docs.astral.sh/uv/getting-started/installation/)), e.g., on Linux/macOS with `curl -LsSf https://astral.sh/uv/install.sh | sh`, on Windows with `powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"`.
2. Clone this repository and change into the `ccgbc` folder.
3. Create the environment and install all dependencies with `uv sync`. This creates a virtual environment in `.venv` according to `uv.lock`. 
4. Start JupyterLab with `uv run jupyter lab` and open the notebooks from the file browser. Alternatively, select the `.venv` environment as kernel in your editor (e.g., VS Code).

Single commands can be run inside the environment with `uv run <command>`. An additional package can be added with `uv add <package>`.


Licensing and Reuse
-------------------

This code is licensed under the GPLv2 license. 

If you use, in any way, parts of this code in your own research or teaching that results in publications, please cite it as follows:<br>
* L. Schmalen, "Channel Coding - Graph-based Codes: Lecture Examples," available online at http://www.github.org/KIT-CEL/lecture-examples/, 2025

In case you have questions, comments, or suggestions regarding this code, please contact Laurent Schmalen (Laurent.Schmalen@kit.edu) or propose a pull request.




