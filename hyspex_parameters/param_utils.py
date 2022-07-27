#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  9 15:34:25 2022
utility functions for parameter calculations

@author: phillms1
"""
import numpy as np

        
def getBand(cube,wvt,wl,kwidth = 5):
    delta = [q-wl for q in wvt]
    bindex = delta.index(min(delta,key=abs))
    if kwidth == 1:
        r = cube[:,:,bindex]
    else:
        w = (kwidth-1)/2
        r = np.median(cube[:,:,int(bindex-w):int(bindex+w)],axis=2)
    return r

def getClosestWavelength(wl,band_list):
    delta = [q-wl for q in band_list]
    return band_list[delta.index(min(delta,key=abs))]

def buildSummary(p1,p2,p3):
    shp = np.shape(p1)
    shp = np.append(shp,3)
    a = np.empty(shp)
    a[:,:,0]=p1
    a[:,:,1]=p2
    a[:,:,2]=p3
    return a

def getBandDepth(cube, wvt,low,mid,hi,lw=5,mw=5,hw=5):
    # retrieve bands from cube
    Rlow = getBand(cube,wvt,low,kwidth=lw)
    Rmid = getBand(cube, wvt,mid,kwidth=mw)
    Rhi = getBand(cube,wvt,hi,kwidth=hw)
    
    # determine wavelengths for low, mid, hi
    WL = getClosestWavelength(low,wvt)
    WM = getClosestWavelength(mid,wvt)
    WH = getClosestWavelength(hi,wvt)
    
    a = (WM-WL)/(WH-WL)     #a gets multipled by the longer band
    b = 1.0-a               #b gets multiplied by the shorter band
    
    # compute the band depth using precomputed a and b
    img = 1.0 - (Rmid/(b*Rlow + a*Rhi))
    nmin = np.nanmin(np.where(img>-np.inf,img,np.nan))
    img = np.where(img>-np.inf,img,nmin)
    return img

def getBandDepthInvert(cube,wvt,low,mid,hi,lw=5,mw=5,hw=5):
    # retrieve bands from cube
    Rlow = getBand(cube,wvt,low,kwidth=lw)
    Rmid = getBand(cube,wvt,mid,kwidth=mw)
    Rhi = getBand(cube,wvt,hi,kwidth=hw)

    # determine wavelength values for closest channels
    WL = getClosestWavelength(low,wvt)
    WM = getClosestWavelength(mid,wvt)
    WH = getClosestWavelength(hi,wvt)
    a = (WM-WL)/(WH-WL)     # a gets multipled by the longer band
    b = 1.0-a               # b gets multiplied by the shorter band

    # compute the band depth using precomputed a and b
    img = 1.0 - ((b*Rlow + a*Rhi)/Rmid)
    nmin = np.nanmin(np.where(img>-np.inf,img,np.nan))
    img = np.where(img>-np.inf,img,nmin)
    return img

def getBandRatio(cube,wvt,num_l,denom_l,num_w=5,denom_w=5):
    num = getBand(cube, wvt, num_l,kwidth=num_w)
    denom = getBand(cube, wvt, denom_l,kwidth=denom_w)
    img = num/denom
    nmin = np.nanmin(np.where(img>-np.inf,img,np.nan))
    img = np.where(img>-np.inf,img,nmin)
    return img

def normalizeParameter(p):
    return (p-np.nanmean(p))/np.nanstd(p)

def cropZeros(p):
    return np.where(p<0,0.0,p)

def browse2bit(B):
    A=B
    for i in range(np.shape(B)[2]):
        b = B[:,:,i]
        A[:,:,i] = np.array((255*((b-np.nanmin(b))/(np.nanmax(b)-np.nanmin(b)))),dtype='int')
    return A

def stretchBand(p,stype='linear',perc = 0.02):
    if stype == 'linear':
        pn = np.where(p==0.0,np.nan,p)
        pn = np.where(p==2.0,np.nan,p)
        sp = (pn-(1-perc)*np.nanmin(pn))/((1-perc)*(np.nanmax(pn)-np.nanmin(pn)))
        sp = np.where(np.isnan(sp) is True,0.0,sp)
        sp = np.where(sp<0,0.0,sp)
        sp = np.where(sp>1,1.0,sp)
    # if flag is '1std':
    #     sp = 
    return sp

def cropNZeros(img):
    img2=img
    n=np.shape(img)[2]
    for b in range(n):
        img2[:,:,b] = cropZeros(img[:,:,b])
    return img2

def stretchNBands(img,perc=0.02):
    img2=img
    n = np.shape(img)[2]
    for b in range(n):
        img2[:,:,b]=stretchBand(img[:,:,b],perc=perc)
    return img2
    


    