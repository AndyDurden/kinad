import os
from kinetic import dydt
import pickle
import numpy as np

# Excel can copy/paste rows separated by newlines, so lets format our data that way
def excelwrite(list, filename, pickleout=0):
    outfile = open(str(filename), 'w')
    if pickleout == 1:
        pickle.dump(list, outfile)
    else:
        for item in list: outfile.write(str(item)+'\n')
    return outfile.close()
    

def write_soln(soln, p,pickleout=0, i=''):
    runname = p['runname'] + str(i)
    soln_deriv = []
    T_deriv = []
    soln_by_component = []
    soln_deriv_by_component = []
    T_deriv_by_trans = []
    for pt in soln: # back-calculate derivatives from soln, a little redundant but easier than prying open odeint
        soln_deriv.append(dydt(pt, 0, p['T']) )
        T_deriv.append(dydt(pt, 0, p['T'], 1) )
    for i in range(0, len(soln[0])): # Separate data by component
        soln_by_component.append(map(lambda x: x[i], soln))
        soln_deriv_by_component.append(map(lambda x: x[i], soln_deriv))
    for i in range(0, len(p['T'])): # Separate by component and transformation
        component = []
        for j in range(0, len(p['namelist'])):
            component.append(map( lambda x: x[j][i], T_deriv))
        T_deriv_by_trans.append(component)

    # Write output to files
    if not os.path.exists(runname):
        os.makedirs(runname)
    excelwrite(soln, runname+'/soln.txt', pickleout)
    excelwrite(soln_deriv, runname+'/soln_deriv.txt', pickleout)
    excelwrite(T_deriv, runname+'/T_deriv.txt', pickleout)
    if pickleout == 0:
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
            if not os.path.exists(runname+'/T_deriv/'+str(i)+'/'):
              os.makedirs(runname+'/T_deriv/'+str(i)+'/')
            T_total_conversion = []
            for j in range(0,len(p['namelist'])):
              print("\n\ni: "+str(i))
              print("\nj: "+str(j))
              print("\nlen T_deriv_by_trans[i]: " + str(len(T_deriv_by_trans[i])))
              
              print("\n\nlen p['t']: "+str(len(p['t'])))
              print("\nlen T_deriv_by_trans[i][j]: "+str(len(T_deriv_by_trans[i][j])))
              excelwrite(T_deriv_by_trans[i][j], runname+'/T_deriv/'+str(i)+'/'+str(j)+'.txt')
              T_total_conversion.append(np.trapz(T_deriv_by_trans[i][j],p['t']))
            excelwrite(T_total_conversion, runname+'/T_deriv/'+str(i)+'/total_conversion.txt')
    
    excelwrite(p['t'], runname+'/time.txt')
    
    miscfile = open(runname+"/run_info.txt",'w')
    import pprint
    miscfile.write("Parameters:\n--------------\n\n"+pprint.pformat(p,1,1))
    miscfile.close()

from uncertainties import unumpy
# For monte carlo uncertainty calculation
# Finds mean and stdev of every output    
def combine_solns(p): # Should only need 2*tstep*N(components) float64s of memory
    for filename in ['soln.txt', 'soln_deriv.txt', 'T_deriv.txt']:
        i=0
        soln = np.asarray(pickle.load(open(p['runname']+str(i)+'/'+filename))) # input file must be pickled
        soln_sq = soln**2
        i=i+1
        while i<p['monteN']:
          tempsoln = np.asarray(pickle.load(open(p['runname']+str(i)+'/'+filename)))
          soln = soln + tempsoln
          soln_sq = soln_sq + tempsoln**2
          i=i+1
        soln = soln/p['monteN']
        stdev = np.sqrt(abs( (soln_sq/p['monteN']) - soln**2 ))
        if not os.path.exists(p['runname']):
          os.makedirs(p['runname'])
        excelwrite(soln, p['runname']+'/'+filename, 0)
        excelwrite(stdev, p['runname']+'/'+filename[:-4]+'_MCstdev.txt', 0)
        if not os.path.exists(p['runname']+'/'+filename[:-4]):
          os.makedirs(p['runname']+'/'+filename[:-4])
          
        if filename=='soln.txt':
          soln_by_component = []
          stdev_by_component = []
          for i in range(0, len(soln[0])): # Separate data by component
            soln_by_component.append(map(lambda x: x[i], soln))
            stdev_by_component.append(map(lambda x: x[i], stdev))
          for i in range(0,len(soln_by_component)):
            excelwrite(soln_by_component[i], p['runname']+'/soln/'+str(i)+'.txt')
            excelwrite(stdev_by_component[i], p['runname']+'/soln/'+str(i)+'_stdev.txt')
            
        if filename=='soln_deriv.txt':
          soln_deriv_by_component = []
          stdev_deriv_by_component = []
          for i in range(0, len(soln[0])): # Separate data by component
            soln_deriv_by_component.append(map(lambda x: x[i], soln))
            stdev_deriv_by_component.append(map(lambda x: x[i], stdev))
          for i in range(0,len(soln_deriv_by_component)):
            excelwrite(soln_deriv_by_component[i], p['runname']+'/soln_deriv/'+str(i)+'.txt')
            excelwrite(stdev_deriv_by_component[i], p['runname']+'/soln_deriv/'+str(i)+'_stdev.txt')
            
        if filename=='T_deriv.txt':
          T_deriv_by_trans = []
          stdev_deriv_by_trans = []
          for i in range(0, len(p['T'])): # Separate by component and transformation
            component = []
            stdev_component = []
            for j in range(0, len(p['namelist'])):
              component.append(map( lambda x: x[j][i], soln))
              stdev_component.append(map( lambda x: x[j][i], stdev))
            T_deriv_by_trans.append(component)
            stdev_deriv_by_trans.append(stdev_component)
          for i in range(0,len(T_deriv_by_trans)): # Iterate over TRANSFORMATIONS
            if not os.path.exists(p['runname']+'/T_deriv/'+str(i)+'/'): # Make folder for each transformation
                os.makedirs(p['runname']+'/T_deriv/'+str(i)+'/')
            
            u_T_deriv = unumpy.uarray(T_deriv_by_trans[i],stdev_deriv_by_trans[i])
            #excelwrite(u_T_deriv, p['runname']+'/T_deriv/'+str(i)+'.txt')
            u_T_total_conversion = []
            for j in range(0,len(u_T_deriv)): # iterate over COMPONENTS in transformation [i]
              excelwrite(T_deriv_by_trans[i][j], p['runname']+'/T_deriv/'+str(i)+'/'+str(j)+'.txt')
              excelwrite(stdev_deriv_by_trans[i][j], p['runname']+'/T_deriv/'+str(i)+'/'+str(j)+'_stdev.txt')
              u_T_total_conversion.append(np.trapz(u_T_deriv[j],p['t']))
            excelwrite(u_T_total_conversion, p['runname']+'/T_deriv/'+str(i)+'/'+'total_conversion.txt')
    
    excelwrite(p['t'], p['runname']+'/time.txt')

    miscfile = open(p['runname']+"/run_info.txt",'w')
    import pprint
    miscfile.write("Parameters:\n--------------\n\n"+pprint.pformat(p,1,1))
    miscfile.close()
    
    # Now we should clean up all the folders made from each monte-carlo simulation
    if not p['keepmonte']:
      import shutil
      for i in range(0, p['monteN']+1):
        shutil.rmtree(p['runname']+str(i), ignore_errors=False, onerror=None)
    
    return 0

    
    
    
    
    
    