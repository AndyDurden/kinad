# kinad
Chemical Kinetics Simulator: solves for concentrations over time in a reaction system



Consider a chemical system with a set of potential energy surface minima (or chemical species) and a set of transformations (reactions) between the minima.

Given a sufficient description of the transformations and an initial state, kinad solves for the concentration of each species at given points in time.

Uncertainty may be calculated via a monte-carlo method.


Input information needed:
-------------------------
List of species/minima <br>
Stoiciometry of transformation <br>
Gibbs Free Activation barrier of transformation (enthalpy barrier may be sufficient but not rigorous [1]) <br>
[Optional Uncertainty of Activation Barriers]
Rate law exponents (order) of transformation <br>

**See the input template when setting up an input file**


Parameters:
-----------
Temperature <br>
Initial State <br>
Points in time solved for <br>


Output:
-------
Concentration at each point in time <br>
Conversion through each tranformation at each point in time <br>


Usage:
-------
python kinad [OPTIONS]

-c, --config=CONFIG_FILE    Specify configuration file, default is config.ini


Configuration:
-------
All parameters are defined in config.ini and described below:

[parameters]
**infile**: Name of Input file located in kinad/. For formatting help see input_template.txt or tosh.in (an example)
**runname**: Name of output folder
**tstart, tend, tstep**: these three parameters define the points in time solved for. at t=0 the system state is equal to your initial state input.
**temp**: Temperature of the simulated reaction, a parameter in the rate law equation.

[montecarlo]
**montecarlo**: 1 to calculate uncertainty via monte-carlo trials, or 0 to run without calculating uncertainties. If 1, your input file must include uncertainty of every activation barrier on the same line separated by a space. [Activation barrier] [Uncertainty]
**monteN**: The number of monte-carlo trials to run. For rigorous statistics typically 1e3-1e6 trials are needed. Note: calculation time should scale linearly with this parameter
**keepmonte**: If 1, the results of each individual monte-carlo trial will be kept. If 0, these trial results are deleted.


To Do:
------
Generate net flow diagram


Foot Notes:
------
[1]: For some systems, frequency calculations to obtain the Gibbs Free Energy barriers rather than Enthalpy barriers are rather expensive, and depending on the system the entropic factor may or may not be significant. For more information on deciding between enthalpy and gibbs barriers see: http://pubs.acs.org/doi/abs/10.1021/acs.jpca.5b09447
