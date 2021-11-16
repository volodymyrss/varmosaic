from numpy import *
import pyfits
import healpy
import plot
import os

class Merged:
    obj=None

    def add(self,obj):
        if self.obj is None:
            self.obj=obj
        else:
            self.obj[1].data=concatenate((self.obj[1].data,obj[1].data))
            self.obj[2].data=vstack((self.obj[2].data,obj[2].data))
            self.obj[3].data=vstack((self.obj[3].data,obj[3].data))
            self.obj[4].data=vstack((self.obj[4].data,obj[4].data))


if True:
    #ra=83
    #dec=22

    ra,dec=266.416817,-26.007825

    theta=(dec+90)/180.*pi
    phi=ra/180.*pi
    hp1024=healpy.ang2pix(1024,theta,phi)
    hp16=healpy.ang2pix(16,theta,phi)
    hp4=healpy.ang2pix(4,theta,phi)
    hp1=healpy.ang2pix(1,theta,phi)

    merged=Merged()

    for r in range(1500):
        try:
            f=pyfits.open("hpdata/%.4i/hp_%i/hp_%i/hp_%i.fits"%(r,hp1,hp4,hp16))
            hp=f[2].data
            fl=f[3].data
            vr=f[4].data
            sg=fl/vr**0.5

            mp=zeros((healpy.nside2npix(1024)));mp[hp.astype(int64)]=sg[0,:]
           # healpy.mollview(mp)
           # plot.plot()

            th,ph=healpy.pix2ang(1024,hp.astype(int64))

            m=(vr>0) & (~isinf(vr)) & (~isnan(vr))
            vr[~m]=NaN
            fl[~m]=NaN

            w_vr=nansum(vr**-1,0)**-1


            w_fl=nansum(fl/vr,0)*w_vr
            w_sg=w_fl/w_vr**0.5

            print w_fl,w_vr

            #plot.p.clf()

            totsig=nansum(sg**2,0)**0.5

            avsig=nansum(sg,0)
            variance=(totsig**2-avsig**2)**0.5

            m=~isnan(variance)

            print w_sg.min(),w_sg.max()

        #    plot.p.scatter(th,ph,c=w_sg,vmin=-3,vmax=5)

            print variance[m].min()
            print variance[m].max()

            merged.add(f)

        #    plot.p.figure()
        #    plot.p.scatter(th,ph,c=variance,vmin=-1,vmax=1)

            if r%2==0:
                path="hpdata/merged/hp_%i/hp_%i/"%(hp1,hp4)
                try:
                    os.makedirs(path)
                except Exception as e:
                    print e
                fn="hpdata/merged/hp_%i/hp_%i/hp_%i.fits"%(hp1,hp4,hp16)
                print "fn:",fn
                merged.obj.writeto(fn,clobber=True)

        except:
            pass
