Installation
=====================================================
The program is written with Python. Thus, the program is executable in Windows platform or Linux platform.  Please make sure that you install Python 3 including **numpy, matplotlib, scipy, csv, copy, datetime, sys and os package** on your computer.

To download the project, please install git in your computer, and type the following command in git bash environment.

.. code-block:: console

	git clone https://github.com/mingchunNTU/Batch_Crystallizer_Simulator.git
	
Important Script
=======================================================
There're four important scripts that provide different functions, which can be executed by a python core directly. There's a section that defines and explain adjustable parameters in the beginning of each script. 
To avoid editing the source code frequently, the physical properties and operating condition are defined in separated files from the python script. You can define the directory that contains these setting files using a parameters setting_dir in each script.  

* src/run.py: execute the simulation, and output the trajectory and product CSD
* src/plot_CSD.py: plot the product CSD
* src/plot_trajectory.py: plot the trajectory during crystallization
* src/plot_solubility.py: plot the relationship between temperature and solute solubility


File Structure of the example/KNO3
==========================================================
You can modify the examples in examples/KNO3 to fit your need. The kinetic parameters used in the example are from the work of Chung et al. (1). The description of the files in examples/KNO3 are shown below.

* operating_condition.csv: the operating condition of the crystallizer
* physical_property.csv: the crystal property and crystallization kinetic parameters
* Result/CSD(number).csv: the simulated crystal size distriubtion (number basis)
* Result/CSD(volume).csv: the simulated crystal size distriubtion (volume basis)
* Result/trajectory.csv: the trajectory during batch crystallization


Ref: (1) Chung, S. H.; Ma, D. L.; Braatz, R. D. Optimal Seeding in Batch Crystallization. The Canadian journal of chemical engineering 1999, 77 (3), 590–596.

Simulation of Crystallization
====================================================
Please edit the operating_condition.csv and physical_property.csv according to the system you studied. The crystallization kinetic parameter, solubility and crystal property are stored in physical_property.csv file. The operation parameters such as seed condition and cooling trajectory are stored in operation.csv file.

The form of solubility in crystallization_simulator is shown as

.. math::

	Csat=A+B*T+C*T^2 

T is in unit of :math:`^oC` and Csat is in unit of g solute/L solvent

The unit of solubility, seed loading and crystal density should be the same because of the mass balance equation. The time and length unit of kg and kb should be consistent because of moment equation. The birth rate used in the simulation is normallized by the solvent volume. Thus, the third moment is dimensionless. 

The seed CSD is assumed to be parabolic form.

.. math:: 

	f(r)=A*r_0(1+w)(1-w)

for :math:`r_0(1+w)>r>r_0(1-w)`

The temperature trajectory is defined in a polynomial form.

.. math::

	T(t)=T_i-(T_i-T_e)*(\frac{t}{batch time})^j
	
T is temperature. Ti is initial temperature. Te is end temperature. t is time. j is the order of the trajectory 

After specifying all the parameters, you can use run.py to execute the simulation. After simulation finished, you can use plot_CSD.py or plot_trajectory.py scripts to visuallize the results. Please edit setting_dir in each script to fit your environment.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`