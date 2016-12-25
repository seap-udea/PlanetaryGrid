import matplotlib
from matplotlib import pyplot as plt,cm,path,patches,_cntr as cntr
from numpy import *
from sys import exit

#################################################################################
# CONSTANTS 
#################################################################################
SQRT3OVER2 = sqrt(3) / 2.
COS60 = cos(60*pi/180)
SIN60 = sin(60*pi/180)

#Default colormap, other options here: http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps
DEFAULT_COLOR_MAP = plt.get_cmap('jet')

#################################################################################
# UTIL
#################################################################################
class dict2obj(object):
    def __init__(self,dic={}):self.__dict__.update(dic)
    def __add__(self,other):
        for attr in other.__dict__.keys():
            exec("self.%s=other.%s"%(attr,attr))
        return self
glob=dict2obj(dict())

def unzip(l):
    #return [x for (x,y) in l], [y for (x,y) in l]
    return zip(*l)

def normalize(xs):
    s = float(sum(xs))
    return [x / s for x in xs]

def argapp(args,app):
    for key in app:args[key]=app[key]
    return args

def colch(val,min,max):
    icol=(val-min)/(max-min)
    return icol

#################################################################################
#PROJECTION ROUTINES
#################################################################################
## Curve Plotting ##
def project_point(p):
    """Maps (x,y,z) coordinates to planar-simplex."""
    x = 0.5 * (2 * p[0] + p[1])
    y = SQRT3OVER2 * p[1]
    return (x, y)

def project(s):
    """Maps (x,y,z) coordinates to planar-simplex."""
    # Is s an appropriate sequence or just a single point?
    try:
        return unzip(map(project_point, s))
    except TypeError:
        return project_point(s)
    except IndexError: # for numpy arrays
        return project_point(s)

def simplex_points(steps=100, boundary=True):
    """Systematically iterate through a lattice of points on the simplex."""
    steps = steps - 1
    start = 0
    if not boundary:
        start = 1
    for x1 in range(start, steps + (1-start)):
        for x2 in range(start, steps + (1-start) - x1):
            x3 = steps - x1 - x2
            yield (x1, x2, x3)

def colormapper(x, a=0, b=1, cmap=None):
    """Maps color values to [0,1] and obtains rgba from the given color map for triangle coloring."""
    if b - a == 0:
        rgba = cmap(0)
    else:
        rgba = cmap((x - a) / float(b - a))
    hex_ = matplotlib.colors.rgb2hex(rgba)
    return hex_

def globalCmap(cmap_name,range=(0.0,1.0),glob=(0.0,1.0)):
    cmap=plt.get_cmap(cmap_name)
    Ncmap=cmap.N
    
    dglob=glob[1]-glob[0]
    
    ini=(range[0]-glob[0])/dglob*cmap.N
    end=(range[1]-glob[1])/dglob*cmap.N+Ncmap
    
    cmaplist=[cmap(i) for i in xrange(int(ini),int(end))]
    Nmap=len(cmaplist)
    
    cmap=cmap.from_list('Custom',cmaplist,Nmap)
    return cmap

def triangle_coordinates(i, j, alt=False):
    """Returns the ordered coordinates of the triangle vertices for i + j + k = N. Alt refers to the averaged triangles; the ordinary triangles are those with base parallel to the axis on the lower end (rather than the upper end)"""
    # N = i + j + k
    if not alt:
        return [(i/2. + j, i * SQRT3OVER2), 
                (i/2. + j + 1, i * SQRT3OVER2), 
                (i/2. + j + 0.5, (i + 1) * SQRT3OVER2)]
    else:
        # Alt refers to the inner triangles not covered by the default case
        return [(i/2. + j + 1, i * SQRT3OVER2), 
                (i/2. + j + 1.5, (i + 1) * SQRT3OVER2), 
                (i/2. + j + 0.5, (i + 1) * SQRT3OVER2)]

def mesh(xvec,yvec,func,excluded=None,verbose=False):
    X,Y=meshgrid(xvec,yvec)
    nvecx=size(xvec)
    nvecy=size(yvec)
    if verbose:print "Generating Mesh Grid: nx = %d, ny = %d"%(nvecx,nvecy)
    Z=zeros((nvecy,nvecx))
    S=[]
    for i in xrange(0,nvecx):
        for j in xrange(0,nvecy):
            x=xvec[i]
            y=yvec[j]
            if verbose:print "\ti = %d/%d, j = %d/%d: x = %e, y = %e"%(i,nvecx,j,nvecy,
                                                                       x,y)
            Z[j,i]=func(xvec[i],yvec[j])
            
            if isnan(Z[j,i]):
                print "Last values:",xvec[i],yvec[j]
                print func(xvec[i],yvec[j])
                exit(0)
            
            if verbose:print "\t\tValue = %e"%(Z[j,i])
            try:
                if Z[j,i]!=excluded:
                    S+=[Z[j,i]]
                else:
                    if verbose:print "\t\t\tValue excluded"
            except:
                S+=[Z[j,i]]
    if verbose:print "%d values of %d not excluded"%(size(S),nvecx*nvecy)

    L=[min(S),max(S),mean(S),std(S)]
    if verbose:print "Limits: ",L
    return X,Y,Z,L

#################################################################################
#MAIN CLASS
#################################################################################
class init(object):

    #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    #INITIALIZATION
    #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    def __init__(self,fig,ax,**args):
        ax.set_aspect('equal')
        self.fig=fig
        self.ax=ax
        self.__dict__.update(args)

    #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    #BOUNDARY
    #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    def boundary(self,**args):
        defargs=dict(linewidth=3.0, color='black')
        scale=1.0
        # Note that the sqrt term is such to prevent noticable roundoff on the top corner point.
        self.ax.plot([0, scale, scale / 2, 0], [0, 0, sqrt(scale * scale * 3.) / 2, 0], 
                     zorder=1000,
                     **argapp(defargs,args))
        self.ax.set_ylim([-0.05 * scale, .90 * scale])
        self.ax.set_xlim([-0.05 * scale, 1.05 * scale])

        self.ax.set_frame_on(False)
        self.ax.axes.get_xaxis().set_visible(False)
        self.ax.set_xticks([]);self.ax.set_yticks([])

    #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    #PLOT TRAJECTORY
    #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    def trajectory(self,t,**args):
            scale=1.0
            defargs=dict(linestyle='-',linewidth=1,color='k')
            xs,ys=project(t)
            self.ax.plot(array(xs)*scale,array(ys)*scale,**argapp(defargs,args))

    def trajectories(self,trajs,**args):
        for t in trajs:
            self.trajectory(t,**args)

    #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    #CONTOURS
    #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    def contour(self,X,Y,Z,levels,**args):
        C=cntr.Cntr(X,Y,Z)
        paths=[]
        values=[]
        i=0
        for level in levels:
            nlist=C.trace(level,level,0)
            nseg = len(nlist)//2
            segs = nlist[:nseg]
            kinds = nlist[nseg:]
            for seg in segs:
                paths+=[seg]
                values+=[level]

        return paths,values

    def contourpath(self,paths,values,cbar=None,cmap=cm.jet,color=-1):
        #READ CBAR
        if cbar is not None:
            cbar.ax.text(-0.5,1.02,"Contours",fontsize=8,
                         ha='right',va='bottom')
            yl=cbar.ax.get_yticklabels()
            yt=cbar.ax.get_yticks()
            valini=float(yl[0].get_text())
            valend=float(yl[-1].get_text())
            fac=1/(valend-valini)
            
        i=0
        for path in paths:
            self.trajectory(path,color=cmap(color))
            if cbar is not None:
                y=fac*(values[i]-valini)
                if y<=1:
                    cbar.ax.text(-0.5,y+0.01,"%.2f"%values[i],
                                  fontsize=8,
                                  ha='right',va='center')
                    cbar.ax.plot([0.0,1.0],[y,y],linewidth=1,color='k')
            i+=1

    def text(self,x=0.0,y=0.0,trajectory=None,text='Text',t=0,center=(0.5,0.5),toff=0.01,**args):
        defargs=dict(fontsize=10,rotation=0,rotation_mode='anchor',
                     horizontalalignment='center',verticalalignment='bottom',
                     bbox=dict(boxstyle='square,pad=0.1',fc="none",ec="none"))
        farg=argapp(defargs,args)
        if trajectory is None:
            xs,ys=project_point([x,y])
        else:
            ntraj=trajectory.shape[0]
            if t>0:
                ipos=int(t*ntraj)
            else:
                dmin=1
                for i in xrange(0,ntraj,ntraj/5):
                    xe,ye=trajectory[i]
                    d=sqrt((xe-center[0])**2+(ye-center[1])**2)
                    if d<=dmin:
                        ipos=i
                        dmin=d
            x,y=trajectory[ipos]
            try:
                xn,yn=trajectory[ipos+1]
            except:
                xn,yn=trajectory[ipos-1]
            xs,ys=project_point([x,y])
            xsn,ysn=project_point([xn,yn])
            try:
                angle=arctan((ysn-ys)/(xsn-xs))*180/pi
            except:
                angle=90
            farg['rotation']=angle
            xs=xs+toff*cos((angle+90)*pi/180)
            ys=ys+toff*sin((angle+90)*pi/180)

        self.ax.text(xs,ys,text,**farg)
            
    def contourf(self,X,Y,Z,L,cmap=cm.jet,colorbar=True,verbose=False,**args):
        defargs=dict(cbar_fontsize=8,cbar_fmt="%.2f",
                     cbar_label='Contour',cmin=L[0],cmax=L[1],
                     cbar_labelsize=12,cbar_labeloff=3,
                     )
        farg=argapp(defargs,args)

        #EXTREMES
        min=farg['cmin']
        max=farg['cmax']
        
        #COORDINATE
        xvec=X[0,:]
        dx=xvec[1]-xvec[0]
        yvec=Y[:,0]
        dy=yvec[1]-yvec[0]
        N=size(xvec)
        N=50
        #COLORBAR
        levels=linspace(min,max,N)
        cset=self.ax.contourf(X,Y,Z,levels=levels,cmap=cmap);self.ax.cla()
        if colorbar:
            cbar=plt.colorbar(cset)
            yt=cbar.ax.get_yticks()
            yl=cbar.ax.get_yticklabels()
            nticks=size(yt)
            ynl=[]
            for i in xrange(0,nticks):
                yval=float(yl[i].get_text())
                ynl+=[farg['cbar_fmt']%yval]
            cbar.ax.set_yticklabels(ynl,fontsize=farg['cbar_fontsize'])
            cbar.ax.set_position([0.78+0.04,0.17,0.9,0.68])
            cbar.ax.text(3.0,0.5,farg['cbar_label'],
                         horizontalalignment='center',
                         verticalalignment='bottom',
                         fontsize=farg['cbar_labelsize'],
                         rotation_mode='anchor',
                         rotation=-90)
                         

        else:cbar=None
    
        #CREATE TRIANGLES
        axpath=self.ax
        i=0
        for x in xvec:
            j=0
            for y in yvec:
                if y+x>1:break
                #FIRST TRIANGLE
                value=Z[j,i]
                color=cmap(colch(value,min,max))
                if (y+x+dx<=1) and (y+x+dy<=1):
                    x1,y1=project_point([x,y])
                    x2,y2=project_point([x+dx,y])
                    x3,y3=project_point([x,y+dy])
                    axpath.add_patch(patches.Polygon([[x1,y1],[x2,y2],[x3,y3]],
                                                     closed=True,
                                                     fill=True,
                                                     color=color,
                                                     edgecolor=None))
                    
                #SECOND TRIANGLE
                if i+1<N and j+1<N and Z[j,i]>0:
                    Zhalf=(Z[j,i+1]+Z[j+1,i]+Z[j,i])/3
                    if verbose:print "\t\tNeighbors: j,i+1=%d,%d,R=%e and j+1,i=%d,%d,R=%e"%(j,i+1,Z[j,i+1],
                                                                                             j+1,i,Z[j+1,i])
                else:Zhalf=Z[j,i]
                value=Zhalf
                color=cmap(colch(value,min,max))
                if y+x+dx+dy<=1:
                    x1,y1=project_point([x+dx,y])
                    x2,y2=project_point([x+dx,y+dy])
                    x3,y3=project_point([x,y+dy])
                    axpath.add_patch(patches.Polygon([[x1,y1],[x2,y2],[x3,y3]],
                                                     closed=True,
                                                     fill=True,
                                                     color=color,
                                                     edgecolor=None))
                    if verbose:print "\tvalue second half=%e"%(value)
                    if verbose:raw_input()
                j+=1
            i+=1
        return cbar

    #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    #DECORATION
    #&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    def grid(self,ut=None,**args):
        defargs=dict(linestyle='--',linewidth=0.1,color=cm.gray(0.5))

        if ut is None:
            try:
                ut=self.ut
            except:
                ut=arange(0.0,1.1,0.1)
        
        for x in ut:
            t=[]
            for y in ut:
                if x+y>1:break
                t+=[[x,y]]
            self.trajectory(t,**argapp(defargs,args))
        for y in ut:
            t=[]
            for x in ut:
                if x+y>1:break
                t+=[[x,y]]
            self.trajectory(t,**argapp(defargs,args))
        for z in ut:
            t=[]
            for y in ut:
                if y+z>1:break
                x=1-y-z
                t+=[[x,y]]
            self.trajectory(t,**argapp(defargs,args))

    def set_title(self,text,**args):
        defargs=dict(position=(0.5,1.05))
        """
        self.ax.text(defargs['position'][0],defargs['position'][1],text,
                     transform=self.ax.transAxes,
                     **argapp(defargs,args))
        """
        self.ax.set_title(text,**argapp(defargs,args))

    def set_xlabel(self,text,loff=None,**args):
        scale=1.0
        defargs=dict(color='k',
                     horizontalalignment='center',
                     verticalalignment='top',
                     rotation_mode='anchor',
                     fontsize=10,
                     rotation=0)
        farg=argapp(defargs,args)
        if loff is None:
            try:
                loff=self.loff
            except:
                loff=0.06

        x=0.5
        y=0.0
        xs,ys=project_point([x,y])
        ys-=loff
        self.ax.text(xs*scale,ys*scale,text,**farg)

    def set_ylabel(self,text,loff=None,**args):
        scale=1.0
        defargs=dict(color='k',
                     horizontalalignment='center',
                     verticalalignment='bottom',
                     rotation_mode='anchor',
                     fontsize=10,
                     rotation=-60)
        farg=argapp(defargs,args)
        if loff is None:
            try:
                loff=self.loff
            except:
                loff=0.06
                
        x=0.5
        y=0.5
        xs,ys=project_point([x,y])
        xs+=loff
        self.ax.text(xs*scale,ys*scale,text,**farg)
        
    def set_zlabel(self,text,loff=None,**args):
        scale=1.0
        defargs=dict(color='k',
                     horizontalalignment='center',
                     verticalalignment='bottom',
                     fontsize=10,
                     rotation_mode='anchor',
                     rotation=60)
        farg=argapp(defargs,args)
        if loff is None:
            try:
                loff=self.loff
            except:
                loff=0.06
                
        x=0.0
        y=0.5
        xs,ys=project_point([x,y])
        xs-=loff
        self.ax.text(xs*scale,ys*scale,text,**farg)
        
    def set_majorticks(self,ut=arange(0.0,1.1,0.1),tick_size=0.02,tick_args=dict(),aoff=0.01,ulabels=[],**args):
        scale=1.0

        defargs=dict(color='k',fontsize=10)
        fargl=argapp(defargs,dict())

        deftickargs=dict(color='k',linewidth=1,linestyle='-')
        fargt=argapp(deftickargs,args)
            
        fargl['rotation_mode']='anchor'
        for x in ut:
            y=0
            xs,ys=project_point([x,y])
            fargl['horizontalalignment']='center'
            fargl['verticalalignment']='top'
            fargl['rotation']=0
            xp=xs
            yp=ys-tick_size*COS60-aoff
            xt=xs-tick_size*COS60
            yt=ys-tick_size*SIN60
            self.ax.text((xt-aoff*COS60)*scale,(yt-aoff*SIN60)*scale,"%.1f"%x,**fargl)
            self.ax.plot([xs*scale,xt*scale],[ys*scale,yt*scale],**fargt)
        for y in ut:
            x=1-y
            xs,ys=project_point([x,y])
            fargl['horizontalalignment']='center'
            fargl['verticalalignment']='bottom'
            fargl['rotation']=-60
            xp=xs+tick_size+aoff*COS60
            yp=ys
            xt=xs+tick_size
            yt=ys
            self.ax.text((xt+aoff)*scale,yt*scale,"%.1f"%y,**fargl)
            self.ax.plot([xs*scale,xt*scale],[ys*scale,yt*scale],**fargt)
        for z in ut:
            y=1-z
            x=1-y-z
            xs,ys=project_point([x,y])
            fargl['horizontalalignment']='center'
            fargl['verticalalignment']='bottom'
            fargl['rotation']=60
            xp=xs-tick_size
            yp=ys
            xt=xs-tick_size*COS60
            yt=ys+tick_size*SIN60
            self.ax.text((xt-aoff*COS60)*scale,(yt+aoff*SIN60)*scale,"%.1f"%z,**fargl)
            self.ax.plot([xs*scale,xt*scale],[ys*scale,yt*scale],**fargt)
        
        self.ut=ut
        self.loff=tick_size+fargl['fontsize']/200.0
