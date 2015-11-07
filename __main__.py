#!/usr/bin/python

from kinetic import *
import os

def main(argv):
    if len(argv) < 6:
        print("\nError: Not enough Arguments.\nUsage: python diffeq [inputfile] [runname] [tstart] [tend] [tstep] [temp]\n")
        return 1
    
    infile = str(sys.argv[1])
    runname = str(sys.argv[2])
    tstart = float(sys.argv[3])
    tend = float(sys.argv[4])
    tstep = float(sys.argv[5])
    temp = float(sys.argv[6])
    
    t = np.linspace(tstart, tend, tstep)

    # Read transformation data from input file
    T = getT(infile, temp)[0] # Transformation matrix
    namelist = getT(infile,temp)[1] # Ordered list of structure names
    
    
    y0 = np.zeros(len(namelist)) # Could be any unit of concentration, since the rate constant is 1/s and dydt = y0 * rate_constant
    # Start with 10 units of first structure
    y0[0] = 100

    # Generate solution
    soln = odeint(dydt, y0, t, (T,))

    # Format lists by component for easier output
    soln_deriv = []
    T_conversion = []
    soln_by_component = []
    soln_deriv_by_component = []
    T_deriv_by_trans = []
    for pt in soln:
        soln_deriv.append(dydt(pt, 0, T) )
        T_conversion.append(dydt(pt, 0, T, 1) )
    for i in range(0, len(soln[0])):
        soln_by_component.append(map(lambda x: x[i], soln))
        soln_deriv_by_component.append(map(lambda x: x[i], soln_deriv))
    for i in range(0, len(T)):
        component = []
        for j in range(0, len(namelist)):
            component.append(map( lambda x: x[j][i], T_conversion))
        T_deriv_by_trans.append(component)

    # Write output to files
    if not os.path.exists(runname):
        os.makedirs(runname)
    excelwrite(soln, runname+'/soln.txt', 0)
    excelwrite(soln_deriv, runname+'/soln_deriv.txt', 0)
    excelwrite(T_conversion, runname+'/T_conversion.txt', 0)
    if not os.path.exists(runname+'/soln'):
        os.makedirs(runname+'/soln')
    for i in range(0,len(soln_by_component)):
        excelwrite(soln_by_component[i], runname+'/soln/'+str(i)+'.txt')
    if not os.path.exists(runname+'/soln_deriv'):
        os.makedirs(runname+'/soln_deriv')
    for i in range(0,len(soln_deriv_by_component)):
        excelwrite(soln_by_component[i], runname+'/soln_deriv/'+str(i)+'.txt')
    if not os.path.exists(runname+'/T_deriv'):
        os.makedirs(runname+'/T_deriv')
    for i in range(0,len(T_deriv_by_trans)):
        excelwrite(T_deriv_by_trans[i], runname+'/T_deriv/'+str(i)+'.txt')
    
    excelwrite(t, runname+'/time.txt')
    
    
    miscfile = open(runname+"/run_info.txt",'w')
    miscfile.write("Parameters:\n--------------\ninput file: "+str(infile)+"\ntemp: "+str(temp)+"\ntstart: "+str(tstart)+"\ntend: "+str(tend)+"\ntstep: "+str(tstep))
    miscfile.close()
    

    
    # Done!
    print("\nRun ended successfully!\n")
    return 0

if __name__ == "__main__": main(sys.argv)


