# This is an older version of the program which 
# only supports intramolecular transformations,
# which have 1 -> 1 stoiciometry.
# These type of systems have the benefit of more simple input files,
# which can be made as an NxN activation barrier matrix
# Rather than a file with 5 lines per transformation.

from kinetic import *

# Legacy Functions
# -----------------

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

# The extra parameter is because odeint expects to pass two extra arguments
def dydt_i(y, t, matrix):
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
# -----------------



infile = sys.argv[1]
runname = sys.argv[2]
tstart = sys.argv[3]
tend = sys.argv[4]
tstep = sys.argv[5]
temp = sys.argv[6]



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
soln = odeint(dydt_i, y0, t, (rate_constants,))
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


# This function creates an array of one component's concentration over time
# (rather than the concentration vector of the whole system)
# Next time just do map( lambda x: x[i], soln)
def component(soln,i,j=-1):
    out = []
    if j == -1:
        for n in range(0,len(soln)): out.append(soln[n][i])
    else: 
        for n in range(0,len(soln)): out.append(soln[n][i][j])
    return out

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


