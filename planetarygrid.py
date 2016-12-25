#######################################################################
#EXTERNAL LIBRARIES
#######################################################################
import commands
import numpy as np

#######################################################################
#PARAMETERS
#######################################################################

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#BEHAVIOR
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DIRPLGRID="./"
SIMULATE=False
VERBOSE=False
GPARAMS=dict(
    Prot=1.0, #Days
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#PLANET PROPERTIES
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""
Structure ('struct'):

    Radius (REarth)
    CentralPressure (Pa)
    CentralDensity
    SurfaceGravitationalField
    CoreDensity
    CoreRadius
    CoreGravitationalField
    MantleAverageDensity
    MantleThick
    MantleGravitationalField
    IceThick
    IcePressure

Thermal evolution ('tevol'):

    DynamoLifetime
    InitialTCMB
    InitialDeltaTCMB
    TimeInnerCore
    AverageQconv
    AverageQsurface
    MaximumQconv
    MaximumRic

Both ('full'):

    MinimumMagneticField
    AverageMagneticField
    MaximumMagneticField
    MinimumDipoleMoment
    AverageDipoleMoment
    MaximumDipoleMoment
"""

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#DATA FIELDS
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
"""
STRUC FIELDS:
0:ur, 1:r, 2:mr, 3:rho, 4:P, 5:g, 6:phi, 7:T, 8:composition
"""
IUR=0
IR=1
IMR=2
IRHO=3
IP=4
IG=5
IPHI=6
ITEMP=7
ICOMP=8

"""
TEVOL FIELDS:
0:t, 1:Qconv[W], 2:Ri[Rp], 3:R*[Rp], 4:Qc[W], 5:Qm[W], 6:Qr[W], 7:Tcmb[K], 
8:Tl[K], 9:Tup[K], 10:RiFlag, 11:Bs[T]
"""
IT=0
IQCONV=1
IRIC=2
IRSTAR=3
IQC=4
IQM=5
IQR=6
ITCMB=7
ITL=8
ITUP=9
IRIFLAG=10

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#PHYSICAL
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
REARTH=6.371E6 #m

#######################################################################
#AUXILIARY VARIABLES
#######################################################################
MPS=[]
MODELS=[]
CMFS=[]
IMFS=[]

#######################################################################
#ROUTINES
#######################################################################
class dict2obj(object):
    """Object like dictionary
    
    Parameters:
    ----------
    dict:
       Dictionary with the attributes of object
    
    Examples:
    --------
    >>> c=dictobj({'a1':0,'a2':1})
    
    Addition:

    >>> c+=dictobj({'a3':2})

    """
    def __init__(self,dic={}):self.__dict__.update(dic)
    def __add__(self,other):
        for attr in other.__dict__.keys():
            exec("self.%s=other.%s"%(attr,attr))
        return self

def System(cmd,out=False,sim=SIMULATE):
    """
    Execute a command and return standard output
    """
    if VERBOSE or sim:print "CMD:\n\t%s"%cmd
    if sim:return ""
    if not out:
        system(cmd)
        output=""
    else:
        output=commands.getoutput(cmd)
    return output

def loadPlanetaryGrid(dirplgrid=DIRPLGRID):
    """
    LOAD PLANETARY GRID PROPERTIES:
    Masses (MPS), Core mass fractions (CMFS), Ice mass fractions (IMFS)
    """
    global MODELS,CMFS,IMFS,MPS
    list=System("cd %s;ls -d CMF*"%dirplgrid,out=True)
    MODELS+=list.split()
    cmfs=[]
    imfs=[]
    for model in MODELS:
        cmf,imf=model.split('-')
        cmfs+=[float(cmf[4:])]
        imfs+=[float(imf[4:])]
    CMFS+=np.unique(np.array(cmfs)).tolist()
    IMFS+=np.unique(np.array(imfs)).tolist()
    list=System("cd %s/%s;ls -d *STRUC*"%(dirplgrid,MODELS[0]),out=True)
    files=list.split()
    Mps=[]
    for file in files:
        parts=file.split('-')
        Mps+=[float(parts[0][1:])]
    MPS+=np.unique(np.array(Mps)).tolist()

    CMFS=np.array(CMFS)
    IMFS=np.array(IMFS)
    MPS=np.array(MPS)

    del((list,cmfs,imfs,Mps))

def loadPlanetCell(Mp=1.0,CMF=0.325,IMF=0.00,dirplgrid=DIRPLGRID,verbose=True):
    """
    LOAD PLANETARY CELL AROUND A MASS AND COMPOSITION.

    A planetary grid sample a 3D space: MP, CMF, IMF.  Grid is sampled in
    discrete values of this space.  If you provide a point in this space
    you will have 8 neighbors:

        ML < MP < MU

        for ML and MU:

                  IMFL

           CMFL   POINT   CMFU

                  IMFU


    Get a cell object:

       model : dictionary having values of the corners (Mp,CMF,IMF)

       sig : signature of the cell.  Array with the values of CMF:IMF:Mp
             at the corners of the cell

       cmfs : range of cmf enclosing data point

       imfs : range of imf enclosing data point 

       mps : range of masses enclosing data point 

       struct : Matrix of structure data - model x points x fields
                0:ur, 1:r, 2:mr, 3:rho, 4:P, 5:g, 6:phi, 7:T, 8:composition

       tevol : Matrix of thermal evolution data - model x points x fields
               0:t, 1:Qconv[W], 2:Ri[Rp], 3:R*[Rp], 4:Qc[W], 5:Qm[W], 6:Qr[W], 7:Tcmb[K], 
               8:Tl[K], 9:Tup[K], 10:RiFlag, 11:Bs[T]

    Usage:

       cell=loadPlanetCell(Mp=2.9,CMF=0.28,IMF=0.15,verbose=False)
       print cell.sig
       print cell.model.Mp
       print len(cell.struct)
       print cell.struct[0] 

    """
    #GLOBALS
    global CMFS,IMFS,MPS

    #ARRAYS
    PGSTRUC=[]
    PGTEVOL=[]
    PGFULL=[]
    PGVAL=[]

    #CMF
    cond=abs(CMFS-CMF)<1E-5
    if np.size(CMFS[cond])>0:
        lcmf=ucmf=CMF
        if verbose:print "Point in the Grid."
    else:
        lcmfs=CMFS[CMFS<=CMF]
        ucmfs=CMFS[CMFS>=CMF]
        lcmf=lcmfs[-1]
        ucmf=ucmfs[0]
    if verbose:print "LCMF,UCMF:",lcmf,ucmf

    #IMF
    cond=abs(IMFS-IMF)<1E-5
    if np.size(IMFS[cond])>0:
        limf=uimf=IMF
        if verbose:print "Point in the Grid."
    else:
        limfs=IMFS[IMFS<=IMF]
        uimfs=IMFS[IMFS>=IMF]
        limf=limfs[-1]
        uimf=uimfs[0]
    if verbose:print "LIMF,UIMF:",limf,uimf

    #LCMF,LIMF
    mps=[]
    for c in ['l','u']:
        for i in ['l','u']:
            cmd="'CMF_%.2f-IMF_%.2f'"
            cmd+="%%(%scmf,%simf)"%(c,i)
            model=eval(cmd)
            cmd="'%.2f:%.2f'"
            cmd+="%%(%scmf,%simf)"%(c,i)
            modelsig=eval(cmd)
            if verbose:print "Reading model:%s"%model
            list=System("cd %s/%s;ls *STRUC*"%(dirplgrid,model),out=True)
            if 'No' in list:continue
            files=list.split()
            Mps=[]
            for file in files:
                parts=file.split('-')
                Mps+=[float(parts[0][1:])]
            Mps=np.array(Mps)
            Mps.sort()
            lcmps=Mps[Mps<=Mp]
            ucmps=Mps[Mps>=Mp]
            lcmp=lcmps[-1];ucmp=ucmps[0]
            mps+=[lcmp,ucmp]
            if verbose:print "UCMP,LCMP:",ucmp,lcmp
            file="%s/%s/M%.2f-STRUC.dat"%(dirplgrid,model,lcmp)
            strl=np.loadtxt(file)
            file="%s/%s/M%.2f-STRUC.dat"%(dirplgrid,model,ucmp)
            stru=np.loadtxt(file)
            PGSTRUC+=[strl,stru]
            file="%s/%s/M%.2f-TEVOL.dat"%(dirplgrid,model,lcmp)
            tevl=np.loadtxt(file)
            file="%s/%s/M%.2f-TEVOL.dat"%(dirplgrid,model,ucmp)
            tevu=np.loadtxt(file)
            PGTEVOL+=[tevl,tevu]
            PGFULL+=[(strl,tevl),(stru,tevu)]
            PGVAL+=["%s:%s"%(modelsig,lcmp),"%s:%s"%(modelsig,ucmp)]

    if np.size(np.unique(PGVAL))==1:
        PGVAL=np.unique(PGVAL)
        nsig=1
    else:
        nsig=np.size(PGVAL)/2
        if nsig<4:
            #IF IT'S A 3 POINT INTERPOLATION CHECK THAT IT IS IN HALF-CELL
            dC=CMF-lcmf
            rdC=ucmf-CMF
            dI=IMF-limf
            ratIC=dI/rdC
            if dI/rdC-1>=1E-5:1/0

    cell=dict2obj(dict(
            model=dict2obj(dict(Mp=Mp,CMF=CMF,IMF=IMF)),
            nsig=nsig,
            sig=np.array(PGVAL),
            cmfs=[lcmf,ucmf],
            imfs=[limf,uimf],
            mps=mps,
            struct=PGSTRUC,
            tevol=PGTEVOL,
            full=PGFULL,
            ))

    return cell

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#STRUCTURE PROPERTIES
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
def getRadius(struc,**args):
    R=struc[-1,1]
    return R

def getCentralPressure(struc,**args):
    P=struc[0,4]
    return P

def getCentralDensity(struc,**args):
    rho=struc[0,3]
    return rho

def getSurfaceGravitationalField(struc,**args):
    g=struc[-1,5]
    return g

def getCoreDensity(struc,**args):
    core=struc[:,8]==0
    rhocs=struc[core,3]
    return rhocs[-3]

def getCoreRadius(struc,**args):
    core=struc[:,8]==0
    Rs=struc[core,1]
    return Rs[-1]

def getCoreGravitationalField(struc,**args):
    core=struc[:,8]==0
    gs=struc[core,5]
    g=gs.mean()
    return g

def getMantleAverageDensity(struc,**args):
    mantle=struc[:,8]==1
    rhocs=struc[mantle[3:-3],3]
    return rhocs.mean()

def getMantleThick(struc,**args):
    mantle=struc[:,8]==1
    Rs=struc[mantle,0]
    return Rs[-1]-Rs[0]

def getMantleGravitationalField(struc,**args):
    mantle=struc[:,8]==1
    gs=struc[mantle,5]
    g=gs.mean()
    return g

def getIceThick(struc,**args):
    ice=struc[:,8]==2
    Rs=struc[ice,0]
    if np.size(Rs)==0:return 0
    else:return Rs[-1]-Rs[0]

def getIcePressure(struc,**args):
    mantle=struc[:,8]==1
    Ps=struc[mantle,4]
    return Ps[-1]

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#THERMAL EVOLUTION PROPERTIES
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
def getDynamoLifetime(tevol,**args):
    try:tau=tevol[-1,0]
    except IndexError:tau=1E8
    return tau

def getInitialTCMB(tevol,**args):
    try:TCMB=tevol[1,7]
    except:TCMB=tevol[7]
    return TCMB

def getInitialDeltaTCMB(tevol,**args):
    try:DeltaT=tevol[1,7]-tevol[1,8]
    except:DeltaT=tevol[7]-tevol[8]
    return DeltaT

def getTimeInnerCore(tevol,**args):
    try:tdyn=tevol[-1,0]
    except:tdyn=1E8
    try:
        Riflag=tevol[:,10]==0
        ts=tevol[Riflag,0]
        tic=ts[-1]
    except:
        tic=tevol[0]
    return tic/tdyn

TPROTECTION=1.0E9
def getAverageQconv(tevol,**args):
    try:
        itav=tevol[:,0]<TPROTECTION
        Qconv=tevol[itav,1]
        Qconv_mean=Qconv.mean()
    except:
        Qconv_mean=tevol[1]
    return Qconv_mean 

def getAverageQsurface(tevol,**args):
    try:
        itav=tevol[:,0]<TPROTECTION
        Qm=tevol[itav,5]
        Qm_mean=Qm.mean()
    except:
        Qm_mean=tevol[5]
    return Qm_mean

def getMaximumQconv(tevol,**args):
    try:Qmax=tevol[:,1].max()
    except:Qmax=tevol[1]
    return Qmax

def getMaximumRic(tevol,**args):
    try:Ricmax=tevol[:,2].max()
    except:Ricmax=tevol[2]
    return Ricmax

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#MAGNETIC PROPERTIES
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
"""
Buffet (Nature, 2010): 
'Numerical models and theoretical consideration
suggest an internal magnetic field of 1-4 mT'
"""
def scaleBrms(Qconv,Ric,Rc,rhoc,verbose=False,**args):
    Chi=Ric/Rc
    D=Rc*(1-Chi)
    V=Rc**3*(1-Chi**3)
    Brms=0.24*np.sqrt(MU0)*rhoc**(1./6)*(D/V)**(1./3)*Qconv**(1./3)
    return Brms

def scaleRolm(Qconv,Ric,Rc,rhoc,P,sigma=6E5,kappa=8E-6,verbose=False,**args):
    CROLM=(1-Chi_E)**(1./3)*(1-Chi_E**3)**(1./2)
    QconvE=3E12
    sigma_E=6E5
    kappa_E=8E-6

    Chi=Ric/Rc
    D=Rc*(1-Chi)
    V=Rc**3*(1-Chi**3)

    if verbose:
        print "Rolm scaling Factors:"
        print "\tQconv = ",(Qconv/QconvE)
        print "\tP = ",(P/DAY)
        print "\tRc = ",(Rc/Rc_E)
        print "\trhoc = ",(rhoc/rhoc_E)
        print "\tsigma = ",(sigma/sigma_E)
        print "\tkappa = ",(kappa/kappa_E)
        print "\tCROLM = ",CROLM/((1-Chi)**(1./3)*(1-Chi**3)**(1./2))

    Rolm=Rosl_E*CROLM*(Qconv/QconvE)**(1./2)*(P/DAY)**(7./6)/((rhoc/rhoc_E)**(1./6)*(Rc/Rc_E)**(11./6)*(1-Chi)**(1./3)*(1-Chi**3)**(1./2)*((sigma/sigma_E)/(kappa/kappa_E))**(1./5))

    return Rolm

def scaleDipfacs(Rolm,**args):
    cmul=funcScalecmul(Rolm)
    fdip,bdip=funcScalefdip(Rolm)
    return cmul,bdip

def surfaceField(Brms,Rolm,Rc,Rp,**args):
    cmul,bdip=scaleDipfacs(Rolm)
    Bs=cmul*Brms/bdip*(Rc/Rp)**3.0
    Ms=DipolarMoment(Bs,Rp)
    return Bs,Ms

def gridSurfaceMagneticField(full,**args):

    defargs=dict(P=GPARAMS['Prot'],prop='Bsurf')
    args=dict2obj(argapp(defargs,args))

    struc=full[0]
    tevol=full[1]

    #Radius
    Rp=struc[-1,1]*Rp_E
    
    #Core density
    core=struc[:,8]==0
    rhocs=struc[core,3]
    rhoc=rhocs[-3]

    #Core radius
    Rs=struc[core,1]
    Rc=Rs[-1]*Rp_E

    #Inner core radius
    try:
        Ric=tevol[:,IRIC]*Rp_E
        Qconv=tevol[:,IQCONV]
        iQ=(Qconv>0)*(-np.isnan(Qconv))
        Ric=Ric[iQ]
        Qconv=Qconv[iQ]
        if np.size(Qconv)<=1:return 0
    except:
        Ric=tevol[IRIC]*Rp_E
        Qconv=tevol[IQCONV]
        if Qconv<0:return 0

    #Brms
    Brms=scaleBrms(Qconv,Ric,Rc,rhoc)
    
    #Rolm
    Rolm=scaleRolm(Qconv,Ric,Rc,rhoc,args.P*DAY)

    #cmul,bdip
    cmul,bdip=scaleDipfacs(Rolm)

    #Bsurf,Mdip
    Bsurf,Mdip=surfaceField(Brms,Rolm,Rc,Rp)
    
    if args.prop=='Mdip':
        P=Mdip/MDIPE
    else:
        P=Bsurf
    
    return P

def getMinimumMagneticField(full,**args):
    struct=full[0]
    tevol=full[1]
    Bs=gridSurfaceMagneticField(full,**args)
    try:
        itav=tevol[:,0]<TPROTECTION
        Bsurf=Bs[itav]
        Bsurf_min=Bsurf.min()
    except:
        try:
            Bsurf_min=Bs[0]
        except:
            Bsurf_min=Bs
    return Bsurf_min 

def getAverageMagneticField(full,**args):
    struct=full[0]
    tevol=full[1]
    Bs=gridSurfaceMagneticField(full,**args)
    try:
        itav=tevol[:,0]<TPROTECTION
        Bsurf=Bs[itav]
        Bsurf_m=Bsurf.mean()
    except:
        try:
            Bsurf_m=Bs[0]
        except:
            Bsurf_m=Bs
    return Bsurf_m 

def getMaximumMagneticField(full,**args):
    struct=full[0]
    tevol=full[1]
    Bs=gridSurfaceMagneticField(full,**args)
    try:
        itav=tevol[:,0]<TPROTECTION
        Bsurf=Bs[itav]
        Bsurf_max=Bsurf.max()
    except:
        try:
            Bsurf_max=Bs[0]
        except:
            Bsurf_max=Bs

    return Bsurf_max

def getMinimumDipoleMoment(full,**args):
    struct=full[0]
    tevol=full[1]
    args['prop']='Mdip'
    Bs=gridSurfaceMagneticField(full,**args)
    try:
        itav=tevol[:,0]<TPROTECTION
        Bsurf=Bs[itav]
        Bsurf_min=Bsurf.min()
    except:
        try:
            Bsurf_min=Bs[0]
        except:
            Bsurf_min=Bs
    return Bsurf_min 

def getAverageDipoleMoment(full,**args):
    struct=full[0]
    tevol=full[1]
    args['prop']='Mdip'
    Bs=gridSurfaceMagneticField(full,**args)
    try:
        itav=tevol[:,0]<TPROTECTION
        Bsurf=Bs[itav]
        Bsurf_m=Bsurf.mean()
    except:
        try:
            Bsurf_m=Bs[0]
        except:
            Bsurf_m=Bs
    return Bsurf_m 

def getMaximumDipoleMoment(full,**args):
    struct=full[0]
    tevol=full[1]
    args['prop']='Mdip'
    Bs=gridSurfaceMagneticField(full,**args)
    try:
        itav=tevol[:,0]<TPROTECTION
        Bsurf=Bs[itav]
        Bsurf_max=Bsurf.max()
    except:
        try:
            Bsurf_max=Bs[0]
        except:
            Bsurf_max=Bs

    return Bsurf_max

#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#INTERPOLATION ROUTINE
#&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
def planetProperty(cell,prop,data='struct',verbose=False,**args):
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #NO INTERPOLATION
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if cell.nsig==1:
        #VALUE
        R=eval("get%s(cell.%s[0],**args)"%(prop,data))
        return R

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #INTERPOLATION WITH 3 POINTS
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if cell.nsig<4:
        R=0
        dcmf=cell.model.CMF-cell.cmfs[0]
        dimf=cell.model.IMF-cell.imfs[0]
        cmfr=cell.model.CMF+dimf
        imfr=cell.model.IMF+dcmf
        
        #INFERIOR MASS
        R1=eval("get%s(cell.%s[0])"%(prop,data))
        if verbose:print "R1=",R1
        R2=eval("get%s(cell.%s[4])"%(prop,data))
        if verbose:print "R2=",R2
        RC=(R2-R1)/(cell.cmfs[1]-cell.cmfs[0])*(cmfr-cell.cmfs[0])+R1
        if verbose:print "RC=",RC
        R2=eval("get%s(cell.%s[2])"%(prop,data))
        if verbose:print "R2=",R2
        RI=(R2-R1)/(cell.imfs[1]-cell.imfs[0])*(imfr-cell.imfs[0])+R1
        if verbose:print "RI=",RI

        dRIC=RI-RC
        if verbose:print "dRIC = ",dRIC

        dIC=np.sqrt((cmfr-cell.cmfs[0])**2+(imfr-cell.imfs[0])**2)
        dmIC=np.sqrt((cell.model.CMF-cell.cmfs[0])**2+(imfr-cell.model.IMF)**2)
        DIC=dIC-dmIC

        if verbose:print "dIC,dmIC,DIC=",dIC,dmIC,DIC

        dRf=DIC*(dRIC/dIC)
        Rml=RC+dRf
        if verbose:print "dRf,Rml=",dRf,Rml

        if cell.mps[0]==cell.mps[1]:return Rml

        #SUPERIOR MASS
        R1=eval("get%s(cell.%s[1])"%(prop,data))
        if verbose:print "R1=",R1
        R2=eval("get%s(cell.%s[5])"%(prop,data))
        if verbose:print "R2=",R2
        RC=(R2-R1)/(cell.cmfs[1]-cell.cmfs[0])*(cmfr-cell.cmfs[0])+R1
        if verbose:print "RC=",RC
        R2=eval("get%s(cell.%s[3])"%(prop,data))
        if verbose:print "R2=",R2
        RI=(R2-R1)/(cell.imfs[1]-cell.imfs[0])*(imfr-cell.imfs[0])+R1
        if verbose:print "RI=",RI

        dRIC=RI-RC
        if verbose:print "dRIC = ",dRIC

        dIC=np.sqrt((cmfr-cell.cmfs[0])**2+(imfr-cell.imfs[0])**2)
        dmIC=np.sqrt((cell.model.CMF-cell.cmfs[0])**2+(imfr-cell.model.IMF)**2)
        DIC=dIC-dmIC

        if verbose:print "dIC,dmIC,DIC=",dIC,dmIC,DIC

        dRf=DIC*(dRIC/dIC)
        Rmu=RC+dRf
        if verbose:print "dRf,Rmu=",dRf,Rmu

        R=(Rmu-Rml)/(cell.mps[1]-cell.mps[0])*(cell.model.Mp-cell.mps[0])+Rml

        if verbose:print "Final R:",R
        
        return R

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    #INTERPOLATION WITH 4 POINTS
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Rs=[]
    ir=0

    #GET CORNER RADIUS
    for i in xrange(0,8,2):
        try:
            R1=eval("get%s(cell.%s[i])"%(prop,data))
        except NameError:
            print "Property '%s' not recognized."%prop
            exit(1)
        R2=eval("get%s(cell.%s[i+1])"%(prop,data))
        if cell.mps[i+1]==cell.mps[i]:
            Rs+=[R2]
        else:
            Rs+=[(R2-R1)/(cell.mps[i+1]-cell.mps[i])*(cell.model.Mp-cell.mps[i])+R1]
        if verbose:print "Property point %d:"%ir,R1,R2,Rs[ir]
        ir+=1

    #INTERPOLATE BETWEEN IMF
    if cell.imfs[1]==cell.imfs[0]:
        Rlcmf=Rs[0]
        Rucmf=Rs[2]
    else:
        Rlcmf=(Rs[1]-Rs[0])/(cell.imfs[1]-cell.imfs[0])*(cell.model.IMF-cell.imfs[0])+Rs[0]
        Rucmf=(Rs[3]-Rs[2])/(cell.imfs[1]-cell.imfs[0])*(cell.model.IMF-cell.imfs[0])+Rs[2]

    #INTERPOLATE BETWEEN CMF
    if cell.cmfs[1]==cell.cmfs[0]:
        R=Rlcmf
    else:
        R=(Rucmf-Rlcmf)/(cell.cmfs[1]-cell.cmfs[0])*(cell.model.CMF-cell.cmfs[0])+Rlcmf

    return R

PG_Property=planetProperty

def planetProperties(Mp=1.0,CMF=0.3,IMF=0.0,
                     properties=['Radius'],data='struct',
                     verbose=False,test=False):
    MMF=1-(CMF+IMF)
    if verbose:print "Retrieving information for: (IMF = %.2f, CMF = %.2f, MMF = %.2f)"%(IMF,CMF,MMF)
    Ps=np.ones(np.size(properties))
    if IMF<0 or CMF<0.1 or MMF<0.1:
        sig=[]
        Ps*=0
    else:
        if test:
            if verbose:print "Mp=%.2f, CMF=%.2f, MMF=%.2f, IMF=%.2f"%(Mp,CMF,MMF,IMF)
            planetcell=loadPlanetCell(Mp=Mp,CMF=CMF,IMF=IMF)
            Rp=planetProperty(planetcell,'Radius',data='struct')
            if verbose:print "Cell signature:",planetcell.sig
        try:
            planetcell=loadPlanetCell(Mp=Mp,CMF=CMF,IMF=IMF,verbose=False)
            sig=planetcell.sig
            Ps=[planetProperty(planetcell,property,data=data) for property in properties]
        except:
            sig=[]
            Ps*=0
    return Ps
