=====================================================================================================================
Solving a 1D Reaction Diffusion Equation with operator splitting (strang) and where the reaction is defined in CellML
=====================================================================================================================

This example solves the weak form of the one-dimensional reaction-diffusion equation using the Galerkin finite element method. The classical form of the governing partial differential equation (PDE) in three-dimension is given by,  

|3d_reaction_diffusion_equation|

where |u| is a scalar dependent variable (e.g. concentration of a chemical species), |conductivity_tensor| is a rank-two diffusion tensor, |R| is the reaction term (function of the dependent variable |u|) and |x| and |t| are spatial coordinates and time respectively.

The reaction term is obtained by solving the following ordinary differential equation (ODE).

|rate_equation|

where |a| is a constant. Here two different time steps are used for the time integration of PDE and ODE and this approach is known as the 'splitting' method. 


.. |3d_reaction_diffusion_equation| image:: ./docs/images/3d_reaction_diffusion_equation.svg
   :align: middle

.. |u| image:: ./docs/images/u.svg
   :align: bottom

.. |conductivity_tensor| image:: ./docs/images/conductivity_tensor.svg
   :align: bottom

.. |R| image:: ./docs/images/r.svg
   :align: bottom

.. |x| image:: ./docs/images/x.svg
   :align: bottom
   
.. |t| image:: ./docs/images/t.svg
   :align: bottom   
   
.. |rate_equation| image:: ./docs/images/rate_equation.svg
   :align: middle   
   
.. |a| image:: ./docs/images/a.svg
   :align: bottom
   

Downloading the example 
====================

This example can be downloaded using git clone::

  git clone https://github.com/OpenCMISS-Examples/cellml_split_reaction_diffusion_equation

The python is immediately executable, but the Fortran executable needs to be built first

Building the example in Fortran
===============================

The fortran version of the example can be configured and built with CMake::

  mkdir cellml_split_reaction_diffusion_equation-build
  cd cellml_split_reaction_diffusion_equation-build
  cmake -DOpenCMISSLibs_DIR=/path/to/opencmisslib/install ../cellml_split_reaction_diffusion_equation
  make

This will create the example executable "cellml_split_reaction_diffusion_equation" in ./src/fortran/ directory.


Running the example
===================

Python version::

  cd ./src/python/
  python cellml_split_reaction_diffusion_equation.py

Fortran version::

  cd ./src/fortran/
  ./cellml_split_reaction_diffusion_equation


Verifying the example
=====================

Results can be visualised by running `visualise.cmgui <./src/fortran/visualise.cmgui>`_ with the `Cmgui visualiser <http://physiomeproject.org/software/opencmiss/cmgui/download>`_.

The following figure shows the one-dimensional computational domain and solution of the primary variable, |u|.

.. |figure1a| image:: ./docs/images/mesh.svg
   :width: 400
   :scale: 125

.. |figure1b| image:: ./docs/images/solution_u.svg
   :width: 400
   :scale: 125
   
|figure1a|  |figure1b|   

Figure 1. (a) One-dimensional finite element mesh (b) Primary variable solution

The expected results from this example are available in `expected_results <./src/fortran/expected_results>`_ folder.  

Prerequisites
=============

The ODE that determines the source/reaction term is solved using CellML and the ODE model is input via `constant_rate.xml <./src/fortran/constant_rate.xml>`_.

License
=======

License applicable to this example is described in `LICENSE <./LICENSE>`_.


