#!/usr/bin/python


# input a file which gives a Ea Matrix for all structures being considered on the PES
# Format of input: First line gives names to structures separated by spaces,
#   the order of names on line #1 will determine the order in which the Ea's are
# Following lines will begin with a name, and then a list of numbers, which are


import string
import sys
import time
import scipy.constants
from scipy.integrate import odeint
from scipy.integrate import trapz
import math
import numpy as np
#import matplotlib.pyplot as plt
import ast
import copy

print(sys.argv)

# Default Args
infile = "in.txt"
lt = time.localtime(time.time())
runname = str(lt[0])+"-"+str(lt[1])+"-"+str(lt[2])+"-"+str(lt[3])+"-"+str(lt[4])+"-"+str(lt[5])

#tstart = 4*10**5
#tend = 4.1*10**5
#tstep = 1*10**3
tstart = 0
tend = 2
tstep = 600
temp = 400


#Arguments:
if len(sys.argv) > 1:
    infile = str(sys.argv[1])
    if len(sys.argv) > 2:
        runname = str(sys.argv[2])
        if len(sys.argv) > 6:
             tstart = float(sys.argv[3])
	     tend = float(sys.argv[4])
	     tstep = float(sys.argv[5])
	     temp = float(sys.argv[6])


#Vector of (value, unit, uncertainty)
h_unit = scipy.constants.physical_constants["Planck constant"] # 6.62e-34  J*s
h = h_unit[0]
kb_unit = scipy.constants.physical_constants["Boltzmann constant"] # 1.38e-23  J/K
kb = kb_unit[0]
R_unit = scipy.constants.physical_constants["molar gas constant"] # 8.314 J/(mol*K)
R = R_unit[0]
# NOTE: Input energies will be in HARTREES, so we must convert to joules.
hartree_unit = scipy.constants.physical_constants["Hartree energy"] # 4.35e-18 J/Hartree
hartree = hartree_unit[0]
joules_per_kcal = 4184


# Put an input file into an array
# Input should be formatted an array with columns separated by spaces and rows separated by new lines
# (Row, Column) is the activation barrier for path (From, To)
def getmatrix(infile):
    f = open(str(infile), 'r')
    matrix = []
    for line in f:
        matrix.append(line.split())
    f.close()
    return matrix


# Formats T from an input file. T describes all transformations in the system and is an input for dydt
# This function is bad don't judge me
def getT(infile):

    # T should be a list of all reactions/transformations X
    # X should be in the following format:
    # X = [ A, B, C, D, E, k] Where A-E are lists and k is the rate constant.
    # A = y indexs of Reactants, B = Reactant Stoic. Coefficients
    # C = y indexs of Products, D = Product Stoic. Coefficients
    # E = Rate law reactant exponents ( E[0] is the exponent of A[0] )

    f = open(str(infile), 'r')
    l = []
    for line in f: l.append(line)
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
                while  (l[i][0] == '#') or (l[i][0] == "\n"):
                    i = i+1
                    if i >= len(l): break
                    print(i)
                if (i>=len(l)): break

                if (Ti == 0 or Ti == 2):
                    print(T)
                    print(len(T))
                    T[len(T)-1].append( map(lambda x: namelist.index(x), l[i].split() ) )
                    Ti = Ti + 1
                    print(T[len(T)-1])

                elif (Ti == 1 or Ti == 3 or Ti == 4):
                    T[len(T)-1].append( map(lambda x: float(x), l[i].split()) )
                    Ti = Ti + 1
                    print(T[len(T)-1])

                elif Ti == 5:
                    T[len(T)-1].append( float(l[i]) )
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
    return T

# Input: Ea matrix, output: reaction rate constant matrix
# 0 if there is no reaction between two structures
def calc_rate_constants(matrix,T):
    outmatrix = matrix[:]
    i1 = 0
    i2 = 0
    while i1 < len(matrix):
        while i2 < len(matrix[i1]):
            if float(matrix[i1][i2]) != 0:
                gibbs_activation = float(matrix[i1][i2])
                exponent =  math.exp( -(gibbs_activation*joules_per_kcal)/(R * T) ) # Using kcal/mol input
		#exponent =  math.exp( -(gibbs_activation*hartree)/(R * T) ) # Using Au (Hartree) input. I think this is missing avogadro's number somewhere.
                #print( str(gibbs_activation) + " -> " + str(exponent))
                rate_constant = ( kb * T  / h ) * exponent
                print( str(gibbs_activation) + " -> " + str(rate_constant)+"\n")
	    # Setting the rate_constant to 0 means no change to the sum calculating the derivative, i.e. there is no path between those two structures
            else: rate_constant = 0
            outmatrix[i1][i2] = rate_constant
            i2 = i2+1
        i2 = 0
        i1 = i1+1
    return outmatrix


# Input, 1-d list of rate constants
def complex_rate_constants(klist):
    outmatrix = klist[:]
    i = 0
    while i < len(klist):
        if klist[i] != 0:
            gibbs_activation = float(klist[i])
            exponent =  math.exp( -(gibbs_activation*joules_per_kcal)/(R * T) ) # Using kcal/mol input
		    #exponent =  math.exp( -(gibbs_activation*hartree)/(R * T) ) # Using Au (Hartree) input. I think this is missing avogadro's number somewhere.
            #print( str(gibbs_activation) + " -> " + str(exponent))
            rate_constant = ( kb * T  / h ) * exponent
            print( str(gibbs_activation) + " -> " + str(rate_constant)+"\n")
	    # Setting the rate_constant to 0 means no change to the sum calculating the derivative, i.e. there is no path between those two structures
        else: rate_constant = 0
        outmatrix[i] = rate_constant
        i = i +1
    return outmatrix



# t is time (does not appear in our equations, but needed for scipy ode solver)
# y is the state of the system, a vector with all concentrations
# output is the vector of derivatives of y

# The extra parameter is because odeint expects to pass two extra arguments
def dydt(y, t, matrix):
    # Array of zeros for derivatives
    output = [0] * len(y)
    # Compute each term of the derivatives
    for i in range(0,len(y)):
        for i2 in range(0,len(y)):
	    # Sum rate of removal of reactant in step
            output[i] = output[i] - y[i]*float(matrix[i][i2])
	    # Sum rate of addition of product in step
            output[i2] = output[i2] + y[i]*float(matrix[i][i2])
    return output

def dydt_detail(y, t, matrix):
    # Should help tell how much of each IC is contributing to each PYR production.
    output = np.zeros( (len(y), len(y) ) )
    for i in range(0,len(y)):
        for i2 in range(0,len(y)):
	    # rate of removal of reactant in step
            output[i][i2] = y[i]*float(matrix[i][i2])
    return output



# This function can calculate the dydt of a system with complicated stoiciometry (not 1->1)
def dydt_complex(y, t, T, output_components = 0):
    # T should be a list of all reactions/transformations X
    # X should be in the following format:
    # X = [ A, B, C, D, E, k] Where A-E are lists and k is the rate constant.
    # A = y indexs of Reactants, B = Reactant Stoic. Coefficients
    # C = y indexs of Products, D = Product Stoic. Coefficients
    # E = Rate law reactant exponents ( E[0] is the exponent of A[0] )
    dTdt = np.zeros(len(T)) # List of rate of each transformation
    dydt = np.zeros(len(y0))
    dyTdt = np.zeros( (len(y0),len(T)) ) # dyTdt[i][j] = contribution of T[j] to dydt[i]
    for i in range(0,len(T)):
        X = T[i]
        k = X[5]
        # Reaction Rate = k * A[0]**E[0] * A[1] * E[1] *...
        dTdt[i] = k
        for j in range(0,len(y)):
            dTdt[i] = dTdt[i] * y[(X[0][j])]**(X[4][j]) # A[i2]**E[i2]
    # Now for each component we sort through T and calculate dy[i]/dt
    for i in range(0,len(y)):
        for j in range(0,len(T)):
            if i in T[0]: # Reactants are CONSUMED
                coeff = T[1][ T[0].index(i) ]
                if coeff > 0: coeff=coeff*(-1) # Make sure this is negative
                dyTdt[i][j] = dyTdt[i][j] + (coeff * dTdt[j] )
                dydt[i] = dydt[i] + (coeff * dTdt[j] )
            if i in T[2]: # Products are GENERATED
                coeff = T[1][ T[0].index(i) ]
                if coeff < 0: coeff=coeff*(-1) # Make sure this is positive
                dyTdt[i][j] = dyTdt[i][j] + (coeff * dTdt[j] )
                dydt[i] = dydt[i] + (coeff * dTdt[j] )
    if output_components == 1: return dyTdt
    return dydt
        


# This function creates an array of one component's concentration over time
# (rather than the concentration vector of the whole system)
def component(soln,i,j=-1):
    out = []
    if j == -1:
        for n in range(0,len(soln)): out.append(soln[n][i])
    else: 
        for n in range(0,len(soln)): out.append(soln[n][i][j])
    return out

#Get the matrix and calculate rate constants from it
#print("Input: \n"+str(getmatrix(str(infile)))+"\n\n")
#rate_constants = calc_rate_constants(getmatrix(str(infile)),273.15 + 150)
rate_constants = calc_rate_constants(getmatrix(str(infile)), temp)
full_output = 1
#print("Rate Constants: \n"+ str(rate_constants)+"\n\n")

#rate_constants = [[0,.01],[0,0]]

# Initial Conditions
#For H: OX, AB_ICNR, AB_IC, AC_IC, BB_IC, BC_IC, AB2BC_INT, A_PYR, B_PYR,
y0 = np.zeros(len(rate_constants)) # Could be any unit of concentration, since the rate constant is 1/s and dydt = y0 * rate_constant
y0[0] = 10

#t = np.logspace(0,10**+0.865,300)
t = np.linspace(tstart, tend, tstep)
soln = odeint(dydt, y0, t, (rate_constants,))
soln_deriv = [[0]*len(soln[0])] *len(soln)
soln_deriv_detail = np.zeros( (len(soln), len(y0), len(y0)) )
for i in range(0,len(soln)):
    soln_deriv[i] = dydt(soln[i], t[i], rate_constants)
    soln_deriv_detail[i] = dydt_detail(soln[i], t[i], rate_constants)


# Print solution
#for point in soln:
#	print(str(point) + "\n")
print("Last Step: ")
print(soln[len(soln)-1])


# Excel can copy/paste stuff separated by newlines, so lets format our data that way
def excelwrite(list, filename, listrep=0):
    outfile = open(str(filename), 'w')
    if listrep == 1:
        outfile.write(str(list))
    else: 
        for item in list: outfile.write(str(item)+'\n')
    return outfile.close()



#Excel write information:

excelwrite(t, str(runname)+"_time.txt")

for i in range(0, len(soln[0])): excelwrite(component(soln,i), str(runname)+"_comp"+str(i)+".txt" )

for i in range(0, len(soln[0])):
    for j in range(0, len(soln[0])):
        if i == j: continue
        excelwrite(component(soln_deriv_detail,i,j), str(runname)+"_conv_"+str(i)+"_to_"+str(j)+".txt", 1)

pcf = open("product_conversion.txt",'w')
for i in range(2,6):
    pcf.write("From "+str(i)+" to "+str(len(y0)-1)+"\n")
    pcf.write(str(np.trapz(component(soln_deriv_detail,i,len(y0)-1), t))+"\n\n")    

    pcf.write("From "+str(i)+" to "+str(len(y0)-2)+"\n")
    pcf.write(str(np.trapz(component(soln_deriv_detail,i,len(y0)-2), t))+"\n\n\n\n")  
pcf.close()    


# Calculate integral of conversion for a transformation

def sumconversion(y, t, listrep=1):
    return np.trapz(y, t)








