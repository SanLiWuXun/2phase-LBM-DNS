#******lattice unit setup******#
LatticeLength=6.23e-6		#lattice length in SI unit
LatticeTime=0.126e-6		#lattice time step in SI unit
LatticeRho=1.0				#lattice density in SI unit
#******lattice unit setup******#

def ConvertToLU(value, type):
	CFL=1.0/LatticeLength		#conversion factor for length
	CFT=1.0/LatticeTime			#conversion factor for time
	CFRho=1.0/LatticeRho	 	#conversion factor for density

	#Conversion of length#
	if type=='length':
		return value*CFL

	#Conversion of velocity#
	if type=='velocity':
		return value*CFL/CFT

	#Conversion of acceleration#
	if type=='acceleration':
		return value*CFL/CFT/CFT

	#Conversion of viscosity#
	if type=='viscosity':
		return value*CFL*CFL/CFT

#convert some value needed
Radius=62.3e-6
temp=ConvertToLU(Radius,'length')
print('Radius in lattice unit is: '+str(temp))

g=9.8
temp=ConvertToLU(g,'acceleration')
print('Gravity in lattice unit is: '+str(temp))

vis=1.54e-5
temp=ConvertToLU(vis,'viscosity')
print('viscosity in lattice unit is: '+str(temp))
tao=3.0*temp+0.5
print('Relaxation time tao is: '+str(tao))
