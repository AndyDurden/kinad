# kinad
Chemical Kinetics Simulator: solves for concentrations over time in a reaction system

Consider a chemical system with a set of potential energy surface minima (or chemical species) and a set of transformations (reactions)
Given a sufficient description of the transformations and an initial state, kinad solves for the concentration of each species at given points in time.

Input information needed:
-------------------------
List of species/minima
Stoiciometry of transformation
Gibbs Free Activation barrier of transformation (enthalpy barrier is sufficient but not rigorous)
Rate law exponents (order) of transformation

Parameters:
-----------
Temperature
Initial State
Points in time solved for


Output:
-------
Concentration at each point in time
Conversion through each tranformation at each point in time


Usage:
-------
Under Construction. Program was originally written for only intramolecular reactions (1:1 stoiciometry and 1st order), and the current usage involves an NxN activation barrier matrix in an input file and is executed with the command:
python kinetic.py [input file] [output name] [time start (s)] [time end] [time step] [temperature (K)]
