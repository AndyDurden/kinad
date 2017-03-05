#!/usr/bin/python

import os
import sys
import string
import sys
import time
import scipy.constants
from scipy.integrate import odeint, trapz
import math
import numpy as np
#import matplotlib.pyplot as plt
import ast
import copy

#Vector of (value, unit, uncertainty)
h_unit = scipy.constants.physical_constants["Planck constant"] # 6.62e-34  J*s
h = h_unit[0]
kb_unit = scipy.constants.physical_constants["Boltzmann constant"] # 1.38e-23  J/K
kb = kb_unit[0]
R_unit = scipy.constants.physical_constants["molar gas constant"] # 8.314 J/(mol*K)
R = R_unit[0]
# NOTE: Input energies from gaussian 09 will be in HARTREES, so we must convert to joules.
hartree_unit = scipy.constants.physical_constants["Hartree energy"] # 4.35e-18 J/Hartree
hartree = hartree_unit[0]
joules_per_kcal = 4184



def rate_constant(Ea, T, uncertainty=0):
    if Ea == 0 : return 0  # Setting the rate_constant to 0 means no change to the derivative sum; no reaction
    exponent =  math.exp( -(Ea*joules_per_kcal)/(R * T) ) # Using kcal/mol input
    #exponent =  math.exp( -(Ea*hartree)/(R * T) ) # Using Au (Hartree) input. I think this is missing avogadro's number somewhere.
    #print( str(gibbs_activation) + " -> " + str(exponent))
    if uncertainty == 0:
      rate_constant = ( kb * T  / h ) * exponent
    else:  # Monte Carlo
      rate_constant = np.random.normal( (( kb * T  / h ) * exponent), uncertainty, 1)[0]
    return rate_constant


# Generates T from an input file. T describes all transformations in the system and is an input for dydt
# This function could be cleaned up a lot.
def getT(infile, temp, montecarlo=0):

    # T (Output) is be a list of every reaction/transformation X
    # X should be in the following format:
    # X = [ A, B, C, D, E, k] Where A-E are lists and k is the rate constant.
    # A = y indexs of Reactants, B = Reactant Stoic. Coefficients
    # C = y indexs of Products, D = Product Stoic. Coefficients
    # E = Rate law reactant exponents ( E[0] is the exponent of A[0] )

    f = open(str(infile), 'r')
    l = []
    for line in f: l.append(line)
    f.close()
    names = {}
    print(len(l))
    i = 0
    while (l[i][0] == '#') or (l[i][0] == "\n"): i = i+1; print(i)
    namelist = l[i].split()
    i = i + 1
    n = len(namelist)
    print(namelist)
    Ti = 0
    T = [[]]
    while i < len(l):
        while Ti < 6:
            if (i<len(l)):
                while  (l[i][0] == '#') or (l[i][0] == "\n"): # skip 
                    i = i+1
                    if i >= len(l): break
                    print(i)
                if (i>=len(l)): break
                
                # Labels
                if (Ti == 0 or Ti == 2):
                    print(T)
                    print(len(T))
                    T[len(T)-1].append( map(lambda x: namelist.index(x), l[i].split() ) )
                    Ti = Ti + 1
                    print(T[len(T)-1])

                # Stoiciometries or Rate Orders
                elif (Ti == 1 or Ti == 3 or Ti == 4):
                    T[len(T)-1].append( map(lambda x: float(x), l[i].split()) )
                    Ti = Ti + 1
                    print(T[len(T)-1])
                
                # k, rate constant
                elif Ti == 5:
                  if montecarlo == 0:
                    if len(l[i].split())>1: # In case uncertainty included in input file without monte carlo option
                      T[len(T)-1].append( rate_constant(float(l[i].split()[0]), temp) ) 
                    else: T[len(T)-1].append( rate_constant(float(l[i]), temp) )
                  elif montecarlo == 1:
                    if len(l[i].split()) < 2: raise Exception("Fatal Error: Monte Carlo uncertainty selected, but no uncertainty is included in the input file.")
                    T[len(T)-1].append( rate_constant( float(l[i].split()[0]), temp, float(l[i].split()[1])) ) # Monte Carlo
                  Ti = Ti + 1
                  print(T[len(T)-1])
                i = i+1; print(i)

        if i >= len(l): break
        while (l[i][0] == '#'):
            i = i+1
            if i >= len(l): break
        if l[i][0] != "\n": raise SyntaxWarning("More than 6 lines in a reaction description??")
        if len(T[len(T)-1]) != 6: print(T[len(T)-1]); raise Exception("Sanity Check Failed: Reaction description is not 6 lines long")
        Ti = 0
        T.append([])
    while T[len(T)-1] == []: T.remove([])
    return [T, namelist]







# t is time (does not appear in our equations, but needed for scipy ode solver)
# y is the state of the system, a vector with all concentrations
# output is the vector of derivatives of y


# This function can calculate the dydt of a system with complicated stoiciometry (not 1->1)
def dydt(y, t, T, output_components = 0):
    # T should be a list of all reactions/transformations X
    # X should be in the following format:
    # X = [ A, B, C, D, E, k] Where A-E are lists and k is the rate constant.
    # A = y indexs of Reactants, B = Reactant Stoic. Coefficients
    # C = y indexs of Products, D = Product Stoic. Coefficients
    # E = Rate law reactant exponents ( E[0] is the exponent of A[0] )

    dTdt = np.zeros(len(T)) # List of rate of each transformation
    dydt = np.zeros(len(y))
    dyTdt = np.zeros( (len(y),len(T)) ) # dyTdt[i][j] = contribution of T[j] to dydt[i]
    for i in range(0,len(T)):
        X = T[i]
        if len(X[0]) != len(X[4]): raise Exception("Error: Reactant list is not the same length as rate order list.")
        k = X[5]
        # Reaction Rate = k * A[0]**E[0] * A[1] * E[1] *...
        dTdt[i] = k
        for j in range(0,len(X[0])):
            dTdt[i] = dTdt[i] * y[(X[0][j])]**(X[4][j]) # A[i2]**E[i2]
    # Now for each component we sort through T and calculate dy[i]/dt
    for i in range(0,len(y)):
        for j in range(0,len(T)):
            X = T[j]
            coeff = 0
            if i in X[0]: # Reactants are CONSUMED
                coeff = X[1][ X[0].index(i) ]
                if coeff > 0: coeff=coeff*(-1) # Make sure this is negative
                dyTdt[i][j] = dyTdt[i][j] + (coeff * dTdt[j] )
            dydt[i] = dydt[i] + (coeff * dTdt[j] )
            if i in X[2]: # Products are GENERATED
                coeff = X[1][ X[2].index(i) ]
                if coeff < 0: coeff=coeff*(-1) # Make sure this is positive
                dyTdt[i][j] = dyTdt[i][j] + (coeff * dTdt[j] )
                dydt[i] = dydt[i] + (coeff * dTdt[j] )
    if output_components == 1: return dyTdt
    return dydt
        





# Calculate integral of conversion for a transformation
def sumconversion(y, t, listrep=1):
    return np.trapz(y, t)

