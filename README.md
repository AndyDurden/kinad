# kinad
Chemical Kinetics Simulator: solves for concentrations over time in a reaction system


**Finally working!!!**

Consider a chemical system with a set of potential energy surface minima (or chemical species) and a set of transformations (reactions)
Given a sufficient description of the transformations and an initial state, kinad solves for the concentration of each species at given points in time.


Input information needed:
-------------------------
List of species/minima <br>
Stoiciometry of transformation <br>
Gibbs Free Activation barrier of transformation (enthalpy barrier is sufficient but not rigorous) <br>
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
python kinad [input file] [output name] [time start (s)] [time end] [time step] [temperature (K)]




To Do:
------
Create a configuration file which may be used instead of parameters.
