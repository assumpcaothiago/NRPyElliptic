# NRPyElliptic
``NRPyElliptic`` [[Assumpcao et. al., PhysRevD.105.104037](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.105.104037)] (also [[arXiv:gr-qc/2111.02424](https://arxiv.org/abs/2111.02424)]) is a new, extensible elliptic solver that sets up initial data (ID) for numerical relativity (NR) using the same numerical methods employed for solving hyperbolic evolution equations. Specifically, ``NRPyElliptic`` implements the hyperbolic relaxation method of [[Rüter et. al., arXiv:gr-qc/1708.07358](https://arxiv.org/abs/1708.07358)] to solve complex nonlinear elliptic PDEs for NR ID. The hyperbolic PDEs are evolved forward in (pseudo)time, resulting in an exponential relaxation of the arbitrary initial guess to a steady state that coincides with the solution of the elliptic system. ``NRPyElliptic`` solves these equations on highly efficient numerical grids exploiting underlying symmetries in the physical scenario. To this end, ``NRPyElliptic`` is built within the [SymPy](https://peerj.com/articles/cs-103/)-based [``NRPy+``](https://nrpyplus.net/) framework, which facilitates the solution of hyperbolic PDEs on Cartesian-like, spherical-like, cylindrical-like, or bispherical-like numerical grids. For the purposes of setting up BBH puncture ID, ``NRPyElliptic`` makes use of the latter.  

This repository hosts the following:  

`Tutorial-NRPyElliptic_BasicEquations.ipynb`: documents the hyperbolized Hamiltonian constraint equation under the CTT formalism 

`Tutorial-Start_to_Finish-NRPyElliptic.ipynb`: documents the C code generation for the standalone version of `NRPyElliptic` 

`NRPyElliptic_codegen`: `NRPy+`-based code

Several `NRPy+`files (with minor modifications that will be incorporated into the main NRPy+ repository in the near future)

`NRPyElliptic_Ccodes`: trusted C code directory for the standalone version

`NRPyEllipticET`: Einstein Toolkit thorn based on the standalone version


Citation:

Thiago Assumpção, Leonardo R. Werneck, Terrence Pierre Jacques, and Zachariah B. Etienne. Fast hyperbolic relaxation elliptic solver for numerical relativity: Conformally flat, binary puncture initial data. Phys. Rev. D 105, 104037, 2022 (doi:10.1103/PhysRevD.105.104037)
