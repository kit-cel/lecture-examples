Example scripts for "Channel Coding - Graph-based Codes" lecture
=================================================================================

This repository contains examples and simulations used in the lecture [Channel Coding - Graph-based Codes]([http://www.cel.kit.edu/lehre_2701.php](https://www.youtube.com/watch?v=BAwHn95atvo&list=PLLDF9ieaJiSX-qDeYdrZNfpMo9UrUcu46)) of the faculty for electrical engineering and information technology (ETIT) at Karlsruhe Institute of Technology (KIT). The corresponding slides and lecture notes can be found in ILIAS (only available for enrolled students). The lecture is freely available on Youtube (click on link above)


Python Dependencies
---------------------
- python >= 3.5, tested with python 3.6
- jupyter notebook
- numpy
- scipy version <= 1.2.1
- matplotlib
- lapack
- cvxopt
- cvxpy

cvxpy has sometimes issues of working with newer versions of scipy and will produce cryptic error messages. It is therefore best to install scipy version 1.2.1. Also, it is necessary to install lapack, which doesn't come with the main anaconda channel. If possible, use the environment.yml (see below).

Usage of Python Notebooks
-------------------------
The programming language [Python](http://www.python.org) is usually pre-installed in current Linux distributions and OSX. Additionally required modules need to be installed by hand from the packet sources. Alternatively, we highly recommend to use readily available Python distributions that are tuned for data science. One such distribution is [Anaconda](https://www.anaconda.com/). Anaconda is also the preferred method to install a complete Python environment on a Windows machine. If you are using Anaconda, we advise you to create an environment within you run the notebooks. You can directly create the environment for running the notebooks using the provided environment.yml file using `conda env create -f environment.yml`. You can then activate the envinronment using `conda activate lecture_CC2`. As the provided environment.yml did not work under the Windows version of Anaconda, we have provided additionally a Windows version called environment_win10.yml which you may try on a Windows10 machine.


Passive Viewing of the Notebooks using Pre-Generated, Static Content
--------------------------------------------------------------------
Please try [Jupyter notebook Viewer](https://nbviewer.jupyter.org/github/KIT-CEL/lecture-examples/tree/master/mloc/). 


Licensing and Reuse
-------------------

This code is licensed under the GPLv2 license. 

If you use, in any way, parts of this code in your own research or teaching that results in publications, please cite it as follows:<br>
* L. Schmalen, "Channel Coding - Graph-based Codes: Lecture Examples," available online at http://www.github.org/KIT-CEL/lecture-examples/, 2025

In case you have questions, comments, or suggestions regarding this code, please contact Laurent Schmalen (Laurent.Schmalen@kit.edu) or propose a pull request.



