#!/usr/bin/python

from kinetic import *
from formatting import *
import os
import getopt
import ConfigParser

def main(argv):
    
    
    help=str("\nUsage: python kinad [OPTION]\n" +
            "-c, --config=CONFIG_FILE     Specify confile file, default is config.ini\n\n" +
            "Please check config.ini to ensure that all parameters are set properly")
            
    cfg="config.ini"
    try:
      opts, args = getopt.getopt(argv, "hc:", ["config="])
    except getopt.GetoptError:
      print(help)
      sys.exit(2)
    for opt, arg in opts:
      if opt == '-h':
        print(help)
        sys.exit(2)
      elif opt in ("-c", "--config"):
        cfg = arg

    c = ConfigParser.ConfigParser()
    c.readfp(open('kinad/'+cfg))
    
    p = {}
    
    p['infile'] = c.get('parameters','infile')
    p['runname'] = c.get('parameters','runname')
    p['tstart'] = float(c.get('parameters','tstart'))
    p['tend'] = float(c.get('parameters','tend'))
    p['tstep'] = float(c.get('parameters','tstep'))
    p['temp'] = float(c.get('parameters','temp'))
    
    p['montecarlo'] = int(c.get('montecarlo','montecarlo'))
    p['monteN'] = int(c.get('montecarlo','monteN'))
    p['keepmonte'] = int(c.get('montecarlo','keepmonte'))
    
    p['t'] = np.linspace(p['tstart'], p['tend'], p['tstep'])

    # Read transformation data from input file
    p['T'], p['namelist'] = getT(p['infile'], p['temp'], montecarlo=p['montecarlo']) # Transformation matrix, Ordered list of structure names
    
    
    y0 = np.zeros(len(p['namelist'])) # no units on concentration, all relative since the rate constant is 1/s and dydt = y0 * rate_constant
    # Start with 100 units of first structure
    y0[0] = 100

    # Generate solution
    soln = odeint(dydt, y0, p['t'], (p['T'],))

    # Format lists by component for easier output
    if p['montecarlo'] == 0:
      write_soln(soln, p, pickleout=0)
    if p['montecarlo'] == 1:
        i=0
        while i<p['monteN']:
          write_soln(soln, p, pickleout=1, i=i)
          p['T'], p['namelist'] = getT(p['infile'], p['temp'], montecarlo=p['montecarlo'])
          y0 = np.zeros(len(p['namelist']))
          y0[0] = 100
          soln = odeint(dydt, y0, p['t'], (p['T'],))
          i = i + 1
        write_soln(soln, p, pickleout=1, i=i)
        combine_solns(p)
    

    
    # Done!
    print("\nRun ended successfully!\n")
    return 0

if __name__ == "__main__": main(sys.argv)


