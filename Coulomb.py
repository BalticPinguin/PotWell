#!/usr/bin/python3.4
import numpy as np
import sys

def main(argv):
    #read the file and store in mmap-object
    infile=open(argv[0], 'r')
    for line in infile:
      data=line.split()
      x=float(data[0])
      y=float(data[1])
      z=float(data[2])
      absphi=float(data[5])
      r=np.sqrt(x*x+y*y+z*z)
      phi=np.cos(x/y)
      print(r,phi, absphi*r)


if __name__ == '__main__':
   main(sys.argv[1:])

