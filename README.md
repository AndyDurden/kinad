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



To Do:
------
Generate net flow diagram


Foot Notes:
------
[1]: For some systems, frequency calculations to obtain the Gibbs Free Energy barriers rather than Enthalpy barriers are rather expensive, and depending on the system the entropic factor may or may not be significant. For more information on deciding between enthalpy and gibbs barriers see: http://pubs.acs.org/doi/abs/10.1021/acs.jpca.5b09447
