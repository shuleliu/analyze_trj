#!/usr/bin/python

# print some configs in a lammpstrj to another file

import sys

inputname = sys.argv[1]
start = int(sys.argv[2])
end = int(sys.argv[3])

outname = str(start)+"_"+str(end)+".lammpstrj"

g = open(outname,'w')

start_write=False

with open(inputname,'r') as f:
     for line in f:
         if "TIMESTEP" in line:
            nxtline = f.next()
            tstep = int(nxtline)
            if tstep>=start and tstep<=end:
               start_write=True
               write_line = line.rstrip('\n')
               print>>g, "%s" %write_line
               write_line = nxtline.rstrip('\n')
               print>>g, "%s" %write_line
            elif tstep<start:
               start_write=False
            elif tstep>end:
               start_write=False
               break

         if start_write and "TIMESTEP" not in line:
            write_line = line.rstrip('\n')
            print>>g, "%s" %write_line

      

g.close()
            
