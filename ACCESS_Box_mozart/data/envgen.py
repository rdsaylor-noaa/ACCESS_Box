#!/usr/bin/env python
#
# envgen
#
# Creates a test environmental data file for an ACCESS_Box simulation
#
# Rick D. Saylor, NOAA/ARL/ATDD, July 2015
# Updated March 2019
#
import os
import sys
import numpy as np
import math
from datetime import datetime, timedelta

# usage
def usage():
    print("usage: %s DATA_FILENAME" % os.path.basename(sys.argv[0]))

# write the header for the file
def wrtheader(envf, hstr, strdesc, strdate):
    for i in range(8):
       envf.write(hstr)
    envf.write("\n#\n")
    envf.write("#   Environmental data for the\n")
    envf.write("#   Atmospheric Chemistry and Canopy Exchange Simulation System (ACCESS)\n")
    envf.write("#   Box model - version 3\n")
    envf.write("#\n")
    envf.write("#   tac  - air temperature (C)\n")
    envf.write("#   pmb  - air pressure (mb)\n")
    envf.write("#   rh   - relative humidity (%)\n")
    envf.write("#   fj   - actinic flux attenuation ()\n")
    envf.write("#\n")
    for i in range(8):
       envf.write(hstr)
    envf.write("\n")
    envf.write("File: "+strdesc+", "+strdate+"\n") 
    envf.write("#\n")
    envf.write("#  date     time   ")
    hvar = ["tac", "pmb", "rh", "fj"]
    nvars = len(hvar)
    for i in range(nvars):
       envf.write(" {0: ^12}".format(hvar[i]))
    envf.write("\n")
    return

# main
def main(argv=None):
    if argv is None:
        argv = sys.argv

    # enforce proper number of arguments passed
    if len(argv) != 2:
        usage()
        return 2

    # get output filename
    fname = argv[1]
    outName = "./"+fname

    # STATIC or DYNAMIC
    ans = input("\n Simulation type\n 1. STATIC, or\n 2. DYNAMIC ?\n")
    if (ans == "1"):
       simtype = "STATIC"
    else:
       simtype = "DYNAMIC"

    print("Creating %s env file\n" % simtype)

    # get start date/time, timestep size (if DYNAMIC), and number of timesteps
    sdtstart = input("Enter start date/time as 'YYYY-mm-dd HH:MM:SS' : ")
    if (simtype == "DYNAMIC"):
       dthr = float(input("Enter timestep in decimal hours: "))
    else:
       dthr = 0.0
    nts  = int(input("Enter number of time steps: ")) + 1
    dtstart = datetime.strptime(sdtstart, "%Y-%m-%d %H:%M:%S")

    # get data
    tacm = float(input("Enter air temperature (deg C) :\n"))
    rhm  = float(input("Enter relative humidity (%) :\n"))
    pmbm = float(input("Enter air pressure (mb) :\n"))
    fjm  = float(input("Enter actinic flux strength (clear sky = 1): \n"))

    # initialize arrays
    dt     = [dtstart]
    tac    = np.zeros(nts)     # air temperature at zmeas (deg C)
    rh     = np.zeros(nts)     # relative humidity (%)
    pmb    = np.zeros(nts)     # air pressure (mb)
    fj     = np.zeros(nts)     # actinic flux strength (clear sky = 1)
    
    if (simtype == "DYNAMIC"):
       tac[0] = tacm+tacm*0.1*math.cos(0.0)
       rh[0]  = rhm*(1.0-0.1*math.cos(0.0))
       pmb[0] = pmbm
       fj[0]  = fjm 
    else:
       tac[0] = tacm
       rh[0]  = rhm
       pmb[0] = pmbm
       fj[0]  = fjm

    for i in range(1,nts):
       hrs=dthr*i
       dt.append(dtstart + timedelta(hours=hrs))
       if (simtype == "DYNAMIC"):
          radhrs=0.2618*hrs
          diurnfac = math.cos(radhrs) 
          tac[i] = tacm+tacm*0.1*diurnfac
          rh[i]  = rhm*(1.0-0.1*diurnfac)
          pmb[i] = pmbm
          fj[i]  = fjm*(0.5+diurnfac)/1.5
          if (fj[i] < 0.0):
              fj[i] = 0.0
       else:
          tac[i] = tacm
          rh[i]  = rhm
          pmb[i] = pmbm
          fj[i]  = fjm
       
    # open output file and write data
    fout = open(outName, "w")
    hstr = "#########!"
    strtrange = dt[0].strftime("%Y-%m-%d %H:%M:%S")+" to "+dt[nts-1].strftime("%Y-%m-%d %H:%M:%S")
    print("Data for %s" % strtrange)
    print("Writing data to %s" % outName)
    wrtheader(fout, hstr, fname, strtrange)
    
    for i in range(nts):
       fout.write("{0} {1: .5e} {2: .5e} {3: .5e} {4: .5e}\n".format(dt[i], tac[i], pmb[i], rh[i], fj[i])) 

    for l in range(8):
       fout.write(hstr)
    fout.close()

    return 0

if __name__ == "__main__":
    sys.exit(main())
