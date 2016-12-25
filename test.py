from planetarygrid import *
import ternary

#========================================
#LOAD PLANETARY GRID
#========================================
loadPlanetaryGrid()

#========================================
#GET GRID POINT
#========================================
Mp=2.9
CMF=0.28
IMF=0.15

print "Interpolating planet: Mp = ",Mp,", CMF = ",CMF,", IMF = ",IMF

#Cell object
cell=loadPlanetCell(Mp,CMF,IMF,verbose=False)

#Neighbors 
print "Grid neighbors: ",cell.sig

#Structure data of fist neighbor
print "Structure of neighbor (CMF,IMF,Mp)=",cell.sig[0],":\n",cell.struct[0]

#Planetary property
Rp=planetProperty(cell,"Radius",data="struct")
print "Planetary radius (R_Earth) = ",Rp

#========================================
#INTERPOLATE PROPERTY
#========================================

#STRUCTURE
P=planetProperties(Mp,CMF,IMF,
                   properties=[
                       'Radius',
                       'CentralPressure',
                       'CentralDensity',
                       'SurfaceGravitationalField',
                       'CoreDensity',
                       'CoreRadius',
                       'CoreGravitationalField',
                       'MantleAverageDensity',
                       'MantleThick',
                       'MantleGravitationalField',
                       'IceThick',
                       'IcePressure'
                   ],
                   data='struct',
                   verbose=False,test=False)
print "Structure Properties:\n",P

#THERMAL EVOLUTION
P=planetProperties(Mp,CMF,IMF,
                   properties=[
                       'DynamoLifetime',
                       'InitialTCMB',
                       'InitialDeltaTCMB',
                       'TimeInnerCore',
                       'AverageQconv',
                       'AverageQsurface',
                       'MaximumQconv',
                       'MaximumRic'
                   ],
                   data='tevol',
                   verbose=False,test=False)
print "Thermal evolution:\n",P

#MAGNETIC PROPERTIES
P=planetProperties(Mp,CMF,IMF,
                   properties=[
                       'MinimumMagneticField',
                       'AverageMagneticField',
                       'MaximumMagneticField',
                       'MinimumDipoleMoment',
                       'AverageDipoleMoment',
                       'MaximumDipoleMoment'
                   ],
                   data='full',
                   verbose=False,test=False)
print "Magnetic properties:\n",P

