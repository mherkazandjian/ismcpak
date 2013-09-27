#<keywords>
#post processing, visualise, viz, bird eye, perspective, sph, fi
#</keywords>
import numpy

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot

from amuse.units import generic_unit_converter
from amuse.units import units, constants

from amuse.community.fi.interface import Fi

from cooling import global_mu

from amuse.io import read_set_from_file

from amuse.community.tdc.interface import TDC

from time import time

def make_cube(sph,L,N=100):

     x,y,z=numpy.indices( ( N+1,N+1,N+1 ))

     x=L*(x.flatten()-N/2.)/N
     y=L*(y.flatten()-N/2.)/N
     z=L*(z.flatten()-N/2.)/N

     vx=units.kms(numpy.zeros_like(x.number))
     vy=units.kms(numpy.zeros_like(x.number))
     vz=units.kms(numpy.zeros_like(x.number))

     t1=time()

rho,rhovx,rhovy,rhovz,rhoe=sph.get_hydro_state_at_point(x,y,z,vx,vy,vz)
     t2=time()
     print t2-t1,N
#    kin_energy=0.5*(rhovx**2+rhovy**2+rhovz**2)/rho
#    internal_energy=abs(rhoe-kin_energy)/rho
#    temp=global_mu / constants.kB * internal_energy

     rho=rho.reshape((N+1,N+1,N+1))
#    temp=temp.reshape((N+1,N+1))
     dl=L/N
     return dl,rho#,numpy.transpose(temp)

def column_density_cube(sph,L,N=100,axis=2):
   dL,rho_cube=make_cube(sph,L,N=N)
   rho=rho_cube.number
   unit=rho_cube.unit
   column_density=dL*unit(rho.cumsum(axis=axis))
   return column_density


for i in range(1,121,20):
   N=100

   part=read_set_from_file("particles-%6.6i"%i,"amuse")

   species="H2"

   allspecies=TDC.species
   part.mass=part.mass*part.abundances[:,eval("allspecies."+species)]

   print 'i:', i
   print 'H2 max:',part.abundances[:,allspecies.H2].max()
   print 'rho max:',part.rho.amax().in_(units.amu/units.cm**3)

   label="mc"

   dt= 100| units.yr
   L=10.| units.parsec
   UnitLength=L
   UnitMass=part.mass.sum()
   UnitVelocity=units.kms

convert=generic_unit_converter.ConvertBetweenGenericAndSiUnits(UnitLength, 
UnitMass, constants.G)


sph=Fi(convert,mode='periodic',redirection="none")#number_of_workers=1,use_gl=False,debugger='gdb')

   sph.parameters.use_hydro_flag=True
   sph.parameters.radiation_flag=False
   sph.parameters.self_gravity_flag=False
   sph.parameters.gamma=1
   sph.parameters.isothermal_flag=True
   sph.parameters.integrate_entropy_flag=False
   sph.parameters.timestep=dt/2
   sph.parameters.verbosity=0

   sph.gas_particles.add_particles(part)

   mass_column_density=column_density_cube(sph,L=2*L,N=N)

number_column_density=(mass_column_density/(1.|units.amu)).in_(units.cm**-2)

   print number_column_density[:,:,N].amax()

   projected=number_column_density[:,:,N]

   LL=L.number
   f=pyplot.figure(figsize=(8,8))

   vmax=5.e20 #projected.value_in(1/units.cm**2).max()
   vmin=0
   print vmax
   print vmin

   pyplot.imshow(projected.value_in(1/units.cm**2),
          extent=[-LL,LL,-LL,LL],vmin=vmin,vmax=vmax,origin='lower')
   pyplot.xlabel("parsec")
   pyplot.savefig('projected_'+species+'-'+label+'-%6.6i.png'%i)
   f.clear()
   pyplot.close(f)
