from numpy import *
import pyfits
import healpy
import plot
import sys
import scipy
import scipy.signal

f=pyfits.open(sys.argv[1])
ijd=f[1].data['IJD']

try:
    hp=f[2].data[0,:]
except:
    hp=f[2].data

fl=f[3].data
vr=f[4].data
sg=fl/vr**0.5

mp=zeros((healpy.nside2npix(1024)));mp[hp.astype(int64)]=sg[0,:]
healpy.mollview(mp)
plot.plot()

th,ph=healpy.pix2ang(1024,hp.astype(int64))
dec=th/pi*180-90
ra=ph/pi*180

m=(vr>0) & (~isinf(vr)) & (~isnan(vr))
vr[~m]=NaN
fl[~m]=NaN

w_vr=nansum(vr**-1,0)**-1
w_fl=nansum(fl/vr,0)*w_vr
w_sg=w_fl/w_vr**0.5

print w_fl,w_vr

plot.p.clf()

totsig=nansum(sg**2,0)**0.5
maxsig=nanmax(sg,0)

avsig=nansum(sg,0)
variance=(totsig**2-avsig**2)**0.5

plot.p.clf()
from matplotlib.pyplot import cm 

lombs=[]

for ipx in range(64**2):
    freq=logspace(-2,2,10)
    i,a,b,c=ijd.astype('float64'),fl[:,ipx].astype('float64'),vr[:,ipx].astype('float64'),freq.astype('float64')
    m=~isnan(a) & (b>0)
    pgrm=scipy.signal.lombscargle(i[m],a[m],c)

 #   print i,nansum(a)

    plot.p.plot(c,pgrm,alpha=abs(nansum(a)/300),lw=1)

    lombs.append(pgrm)

plot.p.loglog()
plot.plot()
    
def plotit(it,vmin,vmax,name):
    plot.p.clf()
    plot.gshow=False

 #   print it

    plot.p.scatter(ra,dec,c=it,vmin=vmin,vmax=vmax,lw=0)
   # plot.p.xlim([81,88])
   # plot.p.ylim([19,25])

    b1=linspace(ra.min(),ra.max(),64)
    b2=linspace(dec.min(),dec.max(),64)
    h=histogram2d(ra,dec,weights=it,bins=(b1,b2))
    pyfits.PrimaryHDU(h[0]).writeto(name.replace(".png",".fits"),clobber=True)

    plot.p.colorbar()
    plot.plot(name)

ev=nansum(fl**2/vr,0)-nansum(vr*0+1,0)

lombs=array(lombs)

print lombs.shape

for l in range(10):
    plotit(lombs[:,l],nanmin(lombs[:,l]),nanmax(lombs[:,l]),"lomb_sky_%i.png"%l)

plotit(w_sg,-2,6,"sig.png")
plotit(maxsig,0,15,"maxsig.png")
plotit(ev,1,600,"ev.png")


m=~isnan(variance)

print w_sg.min(),w_sg.max()

nscw=1
for i in range(ijd.shape[0]/nscw):
    cfl=fl[i*nscw:i*nscw+nscw,:]
    cvr=vr[i*nscw:i*nscw+nscw,:]

    w_vr=nansum(cvr**-1,0)**-1
    w_fl=nansum(cfl/cvr,0)*w_vr
    w_sg=w_fl/w_vr**0.5

    plot.p.clf()
    plot.gshow=False
    plot.p.title("IJD %.15lg"%float(ijd[i*nscw]))
    #plot.p.scatter(ra,dec,c=log10(w_sg),vmin=-1,vmax=2,lw=0)
    plot.p.scatter(ra,dec,c=w_sg,vmin=-2,vmax=6,lw=0)
 #   plot.p.xlim([81,88])
 #   plot.p.ylim([19,25])
    plot.p.colorbar()
    plot.plot("sig_%i.png"%i)

