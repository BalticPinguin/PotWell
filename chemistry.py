#!/usr/bin/python
# Author: Gael Varoquaux <gael.varoquaux@normalesup.org>
# Copyright (c) 2008-2015, Enthought, Inc.
# License: BSD Style.

# modified and generalised by Hubert.

# Plot the atoms and the bonds ################################################
import sys
import numpy as np
from mayavi import mlab
mlab.figure(1, bgcolor=(0, 0, 0), size=(350, 350))
mlab.clf()

# The position of the atoms
atoms_x=[]
atoms_y=[]
atoms_z=[]

stri = ''.join(file(sys.argv[1]).readlines()[2])
data = np.fromstring(stri, sep=' ')
num_atoms= int(data[0])  # first number in 3-rd row
start_point=[data[1], data[2], data[3]]

stri = ''.join(file(sys.argv[1]).readlines()[3])
data = np.fromstring(stri, sep=' ')
NX = int(data[0])
dx = [data[1], data[2], data[3]]
stri = ''.join(file(sys.argv[1]).readlines()[4])
data = np.fromstring(stri, sep=' ')
NY =  int(data[0])
dy =  [data[1], data[2], data[3]]
stri = ''.join(file(sys.argv[1]).readlines()[5])
data = np.fromstring(stri, sep=' ')
NZ = int(data[0])
dz = [data[1], data[2], data[3]]
for i in range(num_atoms):
   stri = ''.join(file(sys.argv[1]).readlines()[6+i])
   data = np.fromstring(stri, sep=' ')
   atoms_x.append( data[1] ) #line 7+i, 2. row
   atoms_y.append( data[2] ) #line 7+i, 3. row
   atoms_z.append( data[3] ) #line 7+i, 4. row
   # Display the orbital ##################################
   atom_i = mlab.points3d(atoms_x[i], atoms_y[i], atoms_z[i],
                  scale_factor=2,
                  resolution=20,
                  color=(1./float(data[0]), float(data[0])/92, 1.-1./float(data[0])), ## make color dependent on atom-type.
                  scale_mode='none')

# The bounds between the atoms, we use the scalar information to give
# color  --> this is too easy at this point; understand and change later.
mlab.plot3d(atoms_x, atoms_y, atoms_z, atoms_z,
            tube_radius=0.4, colormap='Reds')

stri = ''.join(file(sys.argv[1]).readlines()[6+num_atoms:])
data = np.fromstring(stri, sep=' ')
data.shape = (NX, NY, NZ)
x,y,z=np.mgrid[start_point[0]:start_point[0]+dx[0]*NX+dy[0]*NY+dz[0]*NZ:NX*1j,
               start_point[1]:start_point[1]+dx[1]*NX+dy[1]*NY+dz[1]*NZ:NY*1j,
               start_point[2]:start_point[2]+dx[2]*NX+dy[2]*NY+dz[2]*NZ:NZ*1j]
#source = mlab.pipeline.scalar_field(data)
source = mlab.pipeline.scalar_field(x,y,z, data)
min = data.min()
max = data.max()
vol = mlab.pipeline.volume(source, vmin=min + 0.65 * (max - min),
                                   vmax=min + 0.9 * (max - min))

#mlab.view(132, 54, 45, [21, 20, 21.5])

mlab.show()
