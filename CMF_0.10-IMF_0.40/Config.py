# -*- coding: utf-8 -*-
#!/usr/bin/env python

#*****************************************************************************#
# PRE-RUN SETTINGS
#*****************************************************************************#
  #Multiple planets
Multiple = True
  #Compute thermal evolution
th_ev = True
  #Ice Planet
ice_planet = True

  #Conditional to use a stop time after solid core arise
stop_after = False
  #Computation time after solid core arises			[yr]
t_after = 3e9

  #Conditional to use multiple time steps
mult_step = True
  #Time step after solid core arise				[yr]
t_step_1 =  1.0E7
  #Inner core radius after solid core arise, in order 
  # to apply the second time step 				[Rc]
Ri_stop = 0.1
  #Time step after inner solid core reach the Ri_stop value	[yr]
t_step_2 =  1.0E8

  #Turn off evolution when Qconv becomes negative
Qconv_condition	= True
  #Plot data in the end
plot_data = False



#*****************************************************************************#
# PLYNET CONFIGURATION (planet creation)
#*****************************************************************************#

#-----------------------------------------------------------------------------#
#PLANET PROPERTIES
#-----------------------------------------------------------------------------#
  #Planet mass								[ME]  
if Multiple == True:
    #Multiple runs: in this case the user must run <Multiple_runs.py>
    Mp = np.loadtxt('M')
else:
    #Individual run: in this case the user must run <Thermal_Model.py>
    Mp = 1.0
    #Rp = 1.61

#Update radius or mass after integration				[-]
update = 'radius'
  #Core Mass Fraction CMF						[%]
CMF = 0.10
  #Ice Mass Fraction IMF						[%]
IMF = 0.40
  #Density in surface							[kg/m^3]
rho_srf = 1000.0
  #Pressure in surface							[Pa]
P_srf = 0.0
  #Core radius guess							[Rp]
Rc = 0.55
  #Mass-Radius scaled index (guess)                                     [-]
index_scaled = 0.31			#0.27
  #Radius of a planet with M = 1 Mp (with the same properties)		[Rearth]
R1Me = 1.30


  #Core Properties ( Several references )
Core={
    "comp" :	'Fe08FeS02',
    "rho1" :	7171.0,
    "K1"   : 	150.2E9,
    "K1p"  : 	5.675,
    "cp"   :	850.0, 
    "k"	   :	40.,
    "q"	   :	0.91,
    "alpha":	1.11E-5,
    "alphc" :	0.64,
    "grun" : 	1.36,
    "dS"   :	118.0,
    "LH"   :	750E3,
    "lambd":	2.0,
    "kappa":	6.5E-6,
    "sigma":	5E5}
  #Mantle Properties ( Several references )  
Mantle={
    "comp" :	'pv_fmw',
    "rho1" :	4152.,
    "K1"   : 	223.6e9,
    "K1p"  : 	4.274,
    "cp"   :	1250.0, 
    "k"	   :	6.0,
    "q"	   :	1.4,
    "alpha":	3.0E-5,
    "grun" : 	1.48,
    "dS"   :	130.0,
    "lambd":	0.0,
    "kappa":	7.5E-7,
    "sigma":	0.0}
    
  #Radiative elements in mantle
    #"El" : Element name		"Ci" : Initial cocentration
    #"shp": Specefic heat production	"hlt": Half life time
vec_rad = [ {"El":"K",    "Ci":30.7E-9, "shp":2.92E-5, "hlt":1.26E9},
	    {"El":"Th",   "Ci":84.1E-9, "shp":2.64E-5, "hlt":14.0E9},
	    {"El":"U235", "Ci":0.15E-9, "shp":56.9E-5, "hlt":0.704E9},
	    {"El":"U238", "Ci":21.0E-9, "shp":9.46E-5, "hlt":4.4704E9}]

#-----------------------------------------------------------------------------#
#LOAD-ITERATE PROPERTIES
#-----------------------------------------------------------------------------#
  #Load of Iterate planet properties ('load', 'iterate')
mode = 'load'
  #File name to load (*.txt) ['planet-state000']
filename = 'planets/M%1.2f-STRUC'%(Mp)
  #File format ( None, plynet default )
fmt = None

#-----------------------------------------------------------------------------#
#NUMERIC PROPERTIES
#-----------------------------------------------------------------------------#
  #Resolution of urvec
N_ur = 200.
  #Integration Scheme
ply.numeric.confnum.scheme = 'rk4'
  #Integratio step
ply.numeric.confnum.h_step = 1./110
  #Error in residual mass condition
ply.numeric.confnum.accuracy_mr = 1e-11
  #Number maxim of iterations in residual mass condition
ply.numeric.confnum.n_max_mr = 60
  #Number of bisections in planet radius optimization
ply.numeric.confnum.n_section = 3
  #Adimensional minim radius of integration.
  # where criterion mass convergence will be performance
ply.numeric.confnum.r_min_int = 0.0



#*****************************************************************************#
# THERMAL CONFIGURATION
#*****************************************************************************#

#-----------------------------------------------------------------------------#
#INITIAL CONDITIONS
#-----------------------------------------------------------------------------#
  #Superficial Temperature Function [K]
Ts = lambda t: 290.0
  #Initial CMB delta temperature
eps_dTc = None

#-----------------------------------------------------------------------------#
# NUMERIC PARAMETERS
#-----------------------------------------------------------------------------#
  #Time step								[yr]
t_step = 1.0E7
  #Initial time								[yr]
t_min = 0
  #Time of evolution							[yr]
t_max = 1.0e10
  #Thermal evolution integration scheme
num.confnum.scheme = 'rk4'

#-----------------------------------------------------------------------------#
#FREE THERMAL PARAMETERS
#-----------------------------------------------------------------------------#
  #Fe Light element density deficit ( Gaidos, 2010 )			[%]
d_rho = .06
  #Reference Rayleigh ( Gaidos, 2010 )					[-]
Ra_ref = 1000.
  #Coefficient of Radiactive Heat production in mantle			[-]
Qr_fac = 1.2531
  #Type of lid ('mobile', 'stagnant')					[-]
type_lid = 'mobile'
  #Viscosity model in upper mantle ('Arrenhius', 'Nabarro-Herring')	[-]
vismodel_up = 'Arrenhius'
  #Viscosity model in upper mantle ('Arrenhius', 'Nabarro-Herring')	[-]
vismodel_lw = 'Nabarro-Herring'
  #Fe melting point (Solidus temperature reference)			[K]
Tao_c0 = 1808.
  #Correction to melting point of Fe					[-]
Tao_cr = 1.0
  #Potential temperature in zero pressure				[K]
T_adb0 = 1700.0
  #Fraction of contribution of core properties in film CMB zone		[-]
frac_l = 0.4
  #Power law exponent of Nusselt-Rayleigh relation			[-]
delta_Nu = 1/3.