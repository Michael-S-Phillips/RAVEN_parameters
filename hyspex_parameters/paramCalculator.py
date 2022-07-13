#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 10 13:33:42 2022

@author: phillms1
"""
import numpy as np
import sys
import re
import os
import timeit
import IPython.display as Disp
import scipy.integrate as it
from param_utils import *
from matplotlib.widgets import LassoSelector
from matplotlib import pyplot as plt 
from ipywidgets import widgets
from spectral import *
import spectral.io.envi as envi



class paramCalculator:
    '''
    this class handles an input hyperspectral image and returns browse product 
    summary images
    
    I/O
    vfile: /path/to/vis.hdr
        this is a .hdr file, not the .img
    sfile: /path/to/swir.hdr
        this is a .hdr file, not the .img
        
    '''
    
    def __init__(self, vfile, sfile,outdir):
        self.vfile = vfile
        self.sfile = sfile
        self.outdir = outdir
        print('loading objects using spectral')
        self.v_ = envi.open(vfile)
        self.s_ = envi.open(sfile)
        print('objects loaded\nloading data')
        tic = timeit.default_timer()
        self.v = np.flip(np.transpose(self.v_.load(),(1,0,2)),axis=0)
        self.s = np.flip(np.transpose(self.s_.load(),(1,0,2)),axis=0)
        toc = timeit.default_timer()-tic
        print(f'{np.round(toc/60,2)} minutes to load data')
        # load wave tables
        self.v_bands = [float(b) for b in self.v_.metadata['wavelength']]
        self.s_bands = [float(b) for b in self.s_.metadata['wavelength']]
        self.v_preview_bands = [int(float(b)) for b in self.v_.metadata['default bands']]
        self.s_preview_bands = [int(float(b)) for b in self.s_.metadata['default bands']]
        # self.HCPINDEX2()
        # self.LCPINDEX2()
        # self.OLINDEX3()
    # -------------------------------------------------------------------------
    # utilities
    # -------------------------------------------------------------------------
    def previewData(self):
        # preview data
        plt.subplot(2,1,1)
        plt.title('VIS Preview')
        v_preview = self.v[:,:,self.v_preview_bands]
        v_preview = [(r-np.min(r))/(np.max(r)-np.min(r)) for r in v_preview]
        plt.imshow(v_preview)
        plt.subplot(2,1,2)
        plt.title('SWIR Preview')
        s_preview = self.s[:,:,self.s_preview_bands]
        s_preview = [(r-np.min(r))/(np.max(r)-np.min(r)) for r in s_preview]
        plt.imshow(s_preview)
        plt.show()

    # -------------------------------------------------------------------------
            
    # -------------------------------------------------------------------------
    # parameter library
    # -------------------------------------------------------------------------
    def HCPINDEX2(self):
        cube = self.s
        wvt = self.s_bands
        # extract data from image cube
        R2120 = getBand(cube, wvt,2120)
        R2140 = getBand(cube, wvt,2140)
        R2230 = getBand(cube, wvt,2230)
        R2250 = getBand(cube, wvt,2250)
        R2430 = getBand(cube, wvt,2430)
        R2460 = getBand(cube, wvt,2460)
        R2530 = getBand(cube, wvt,2530)
        R1810 = getBand(cube, wvt,1810)
        R1690 = getBand(cube, wvt,1690)
    
        W2120 = getClosestWavelength(2120,wvt)
        W2140 = getClosestWavelength(2140,wvt)
        W2230 = getClosestWavelength(2230,wvt)
        W2250 = getClosestWavelength(2250,wvt)
        W2430 = getClosestWavelength(2430,wvt)
        W2460 = getClosestWavelength(2460,wvt)
        W2530 = getClosestWavelength(2530,wvt)
        W1690 = getClosestWavelength(1690,wvt)
        W1810 = getClosestWavelength(1810,wvt)
    
    
        # compute the corrected reflectance interpolating 
        slope = (R2530 - R1690) / (W2530 - W1690)      
        intercept = R2530 - slope * W2530
    
        # weighted sum of relative differences
        Rc2120 = slope*W2120 + intercept
        Rc2140 = slope*W2140 + intercept
        Rc2230 = slope*W2230 + intercept
        Rc2250 = slope*W2250 + intercept
        Rc2430 = slope*W2430 + intercept
        Rc2460 = slope*W2460 + intercept
    
        img=((1-(R2120/Rc2120))*0.1) + ((1-(R2140/Rc2140))*0.1) + ((1-(R2230/Rc2230))*0.15) + ((1-(R2250/Rc2250))*0.3) + ((1-(R2430/Rc2430))*0.2) + ((1-(R2460/Rc2460))*0.15)
        nmin = np.nanmin(np.where(img>-np.inf,img,np.nan))
        img = np.where(img>-np.inf,img,nmin)
        return img
    def LCPINDEX2(self):
        cube=self.s
        wvt=self.s_bands
        # extract data from image cube
        R1690 = getBand(cube, wvt, 1690)
        R1750 = getBand(cube, wvt, 1750)
        R1810 = getBand(cube, wvt, 1810)
        R1870 = getBand(cube, wvt, 1870)
        R1560 = getBand(cube, wvt, 1560)
        R2450 = getBand(cube, wvt, 2450)
    
        W1690 = getClosestWavelength(1690,wvt)
        W1750 = getClosestWavelength(1750,wvt)
        W1810 = getClosestWavelength(1810,wvt)
        W1870 = getClosestWavelength(1870,wvt)
        W1560 = getClosestWavelength(1560,wvt)
        W2450 = getClosestWavelength(2450,wvt)
    
        # compute the corrected reflectance interpolating 
        slope = (R2450 - R1560)/(W2450 - W1560)
        intercept = R2450 - slope * W2450
    
        # weighted sum of relative differences
        Rc1690 = slope*W1690 + intercept
        Rc1750 = slope*W1750 + intercept
        Rc1810 = slope*W1810 + intercept
        Rc1870 = slope*W1870 + intercept
    
        img=((1-(R1690/Rc1690))*0.2) + ((1-(R1750/Rc1750))*0.2) + ((1-(R1810/Rc1810))*0.3) + ((1-(R1870/Rc1870))*0.3)
        nmin = np.nanmin(np.where(img>-np.inf,img,np.nan))
        img = np.where(img>-np.inf,img,nmin)
        return img
    def OLINDEX3(self):
        cube=self.s
        wvt=self.s_bands
        # extract data from image cube
        R1210 = getBand(cube, wvt,1210)
        R1250 = getBand(cube, wvt,1250)
        R1263 = getBand(cube, wvt,1263)
        R1276 = getBand(cube, wvt,1276)
        R1330 = getBand(cube, wvt,1330)
        R1750 = getBand(cube, wvt,1750)
        R1862 = getBand(cube, wvt,1862)
    
        # find closest Hyspex wavelength
        W1210 = getClosestWavelength(1210,wvt)
        W1250 = getClosestWavelength(1250,wvt)
        W1263 = getClosestWavelength(1263,wvt)
        W1276 = getClosestWavelength(1276,wvt)
        W1330 = getClosestWavelength(1330,wvt)
        W1750 = getClosestWavelength(1750,wvt)
        W1862 = getClosestWavelength(1862,wvt)
    
        # ; compute the corrected reflectance interpolating 
        slope = (R1862 - R1750)/(W1862 - W1750)   #;slope = ( R2120 - R1690 ) / ( W2120 - W1690 )
        intercept = R1862 - slope*W1862               #;intercept = R2120 - slope * W2120
    
        Rc1210 = slope * W1210 + intercept
        Rc1250 = slope * W1250 + intercept
        Rc1263 = slope * W1263 + intercept
        Rc1276 = slope * W1276 + intercept
        Rc1330 = slope * W1330 + intercept
    
        img = (((Rc1210-R1210)/(abs(Rc1210)))*0.1) + (((Rc1250-R1250)/(abs(Rc1250)))*0.1) + (((Rc1263-R1263)/(abs(Rc1263)))*0.2) + (((Rc1276-R1276)/(abs(Rc1276)))*0.2) + (((Rc1330-R1330)/(abs(Rc1330)))*0.4)
        nmin = np.nanmin(np.where(img>-np.inf,img,np.nan))
        img = np.where(img>-np.inf,img,nmin)
        return img
    def RPEAK1(self):
        cube = self.v
        wvt = self.v_bands
        rp_wv = [442,533,600,710,740,775,800,833,860,892]
        rp_i = [wvt.index(getClosestWavelength(i,wvt)) for i in rp_wv]
        rp_w = [getClosestWavelength(i,wvt) for i in rp_wv]
        rp_ = cube[:,:,rp_i]
        x_ = np.linspace(rp_w[0],rp_w[-1],num=5000)
        rp_l = np.empty(np.shape(cube)[0:2])
        rp_r = np.empty(np.shape(cube)[0:2])
        for i in range(np.shape(rp_)[0]):
            for j in range(np.shape(rp_)[1]):
                coefs = np.polyfit(rp_w,rp_[i,j,:],4)
                poly = np.poly1d(coefs)
                y_ = list(poly(x_))
                rp_l[i,j] = x_[y_.index(np.max(y_))]/1000
                rp_r[i,j] = np.max(y_)
        return rp_l,rp_r
    def BDI1000VIS(self,rp_r=False):
        cube = self.v
        wvt = self.v_bands
        if rp_r is False:
            rp_l, rp_r = self.RPEAK1()
        bdi_wv = [833,860,892,925,951,984,1023] 
        vi = [wvt.index(getClosestWavelength(i,wvt)) for i in bdi_wv]
        wv_um = [getClosestWavelength(i,wvt)/1000 for i in bdi_wv]
        bdi1000_cube = cube[:,:,vi]
        bdi_norm = np.empty(np.shape(bdi1000_cube))
        for b in range(len(vi)):
            bdi_norm[:,:,b] = bdi1000_cube[:,:,b]/rp_r
        
        bdi1000vis_value = np.empty((np.shape(cube)[0],np.shape(cube)[1]))
        for i in range(np.shape(cube)[0]):
            for j in range(np.shape(cube)[1]):
                spec_vec = bdi_norm[i,j,:]
                coefs = np.polyfit(wv_um,spec_vec,3)
                poly = np.poly1d(coefs)
                pint = np.poly1d(poly.integ())
                bdi1000vis_value[i,j] = pint(wv_um[-1])-pint(wv_um[0])#it.(wv_um,1.0-spec_vec)
        return bdi1000vis_value
    def BD530_2(self):
        cube = self.v
        wvt = self.v_bands
        img = getBandDepth(cube,wvt,440,530,614)
        return img
    def BD920_2(self):
        cube = self.v
        wvt = self.v_bands
        img = getBandDepth(cube,wvt,807,920,984)
        return img
    def BD2210_2(self):
        cube = self.s
        wvt = self.s_bands
        img = getBandDepth(cube,wvt,2165,2210,2290) #(kaolinite group)
        return img
    def BD2190(self):
        cube = self.s
        wvt = self.s_bands
        img = getBandDepth(cube,wvt,2120,2185,2250,mw=3,hw=3) #(Beidellite, Allophane)
        return img
    def BD2250(self):
        cube = self.s
        wvt = self.s_bands
        img = getBandDepth(cube,wvt,2120, 2245, 2340,mw=7,hw=3) 
        return img
    def BD2165(self):
        cube = self.s
        wvt = self.s_bands
        img = getBandDepth(cube,wvt,2120,2165,2230,mw=3,hw=3) #(kaolinite group)
        return img
    def BD2355(self):
        cube = self.s
        wvt = self.s_bands
        img = getBandDepth(cube,wvt,2300, 2355, 2450) #(fe/mg phyllo group)
        return img
    def BD2290(self):
        cube = self.s
        wvt = self.s_bands
        img = getBandDepth(cube,wvt,2250, 2290, 2350) #(fe/mg phyllo group)
        return img
    def D2200(self):
        cube = self.s
        wvt = self.s_bands
        # extract individual channels
        R1815 = getBand(cube, wvt,1815, kwidth=7)
        R2165 = getBand(cube, wvt,2165)
        R2210 = getBand(cube, wvt,2210, kwidth=7)
        R2230 = getBand(cube, wvt,2230, kwidth=7)
        R2430 = getBand(cube, wvt,2430, kwidth=7)
    
    
        # retrieve the CRISM wavelengths nearest the requested values
        W1815 = getClosestWavelength(1815, wvt)
        W2165 = getClosestWavelength(2165, wvt) 
        W2210 = getClosestWavelength(2210, wvt)
        W2230 = getClosestWavelength(2230, wvt)
        W2430 = getClosestWavelength(2430, wvt)
        
        slope = (R2430 - R1815)/(W2430 - W1815)
    
        CR2165 = R1815 + slope*(W2165 - W1815)    
        CR2210 = R1815 + slope*(W2210 - W1815)
        CR2230 = R1815 + slope*(W2230 - W1815)
    
        # compute d2300 with IEEE NaN values in place of CRISM NaN
        img = 1 - (((R2210/CR2210) + (R2230/CR2230))/(2*(R2165/CR2165)))
        nmin = np.nanmin(np.where(img>-np.inf,img,np.nan))
        img = np.where(img>-np.inf,img,nmin)
        return img
    def D2300(self):
        cube = self.s
        wvt = self.s_bands
        # extract individual channels
        R1815 = getBand(cube,wvt,1815)
        R2120 = getBand(cube,wvt,2120)
        R2170 = getBand(cube,wvt,2170)
        R2210 = getBand(cube,wvt,2210)
        R2290 = getBand(cube,wvt,2290,kwidth=3)
        R2320 = getBand(cube,wvt,2320,kwidth=3)
        R2330 = getBand(cube,wvt,2330,kwidth=3)
        R2530 = getBand(cube,wvt,2530)
        
        # get closestgetClosestWavelengthngth
        W1815 = getClosestWavelength(1815,wvt)
        W2120 = getClosestWavelength(2120,wvt)
        W2170 = getClosestWavelength(2170,wvt)
        W2210 = getClosestWavelength(2210,wvt)
        W2290 = getClosestWavelength(2290,wvt)
        W2320 = getClosestWavelength(2320,wvt)
        W2330 = getClosestWavelength(2330,wvt)
        W2530 = getClosestWavelength(2530,wvt)
        
        # compute the interpolated continuum values at selected wavelengths between 1815 and 2530
        slope = (R2530 - R1815)/(W2530 - W1815)
        CR2120 = R1815 + slope*(W2120 - W1815)
        CR2170 = R1815 + slope*(W2170 - W1815)
        CR2210 = R1815 + slope*(W2210 - W1815)
        CR2290 = R1815 + slope*(W2290 - W1815)
        CR2320 = R1815 + slope*(W2320 - W1815)
        CR2330 = R1815 + slope*(W2330 - W1815)
    
        # compute d2300 with IEEE NaN values in place of CRISM NaN
        img = 1 - (((R2290/CR2290) + (R2320/CR2320) + (R2330/CR2330))/((R2120/CR2120) + (R2170/CR2170) + (R2210/CR2210)))
        nmin = np.nanmin(np.where(img>-np.inf,img,np.nan))
        img = np.where(img>-np.inf,img,nmin) 
        return img
    def BD1900r2(self):
        cube = self.s
        wvt = self.s_bands
        # extract individual channels, replacing CRISM_NANs with IEEE_NaNs
        R1908 = getBand(cube,wvt,1908, kwidth = 1) 
        R1914 = getBand(cube,wvt,1914, kwidth = 1) 
        R1921 = getBand(cube,wvt,1921, kwidth = 1) 
        R1928 = getBand(cube,wvt,1928, kwidth = 1) 
        R1934 = getBand(cube,wvt,1934, kwidth = 1) 
        R1941 = getBand(cube,wvt,1941, kwidth = 1) 
        R1862 = getBand(cube,wvt,1862, kwidth = 1) 
        R1869 = getBand(cube,wvt,1869, kwidth = 1) 
        R1875 = getBand(cube,wvt,1875, kwidth = 1) 
        R2112 = getBand(cube,wvt,2112, kwidth = 1) 
        R2120 = getBand(cube,wvt,2120, kwidth = 1) 
        R2126 = getBand(cube,wvt,2126, kwidth = 1) 
        
        R1815 = getBand(cube, wvt, 1815);
        R2132 = getBand(cube, wvt, 2132); 
        
        # retrieve the CRISM wavelengths nearest the requested values
        W1908 = getClosestWavelength(1908, wvt)
        W1914 = getClosestWavelength(1914, wvt)
        W1921 = getClosestWavelength(1921, wvt)
        W1928 = getClosestWavelength(1928, wvt)
        W1934 = getClosestWavelength(1934, wvt)
        W1941 = getClosestWavelength(1941, wvt)
        W1862 = getClosestWavelength(1862, wvt)#
        W1869 = getClosestWavelength(1869, wvt)
        W1875 = getClosestWavelength(1875, wvt)
        W2112 = getClosestWavelength(2112, wvt)
        W2120 = getClosestWavelength(2120, wvt)
        W2126 = getClosestWavelength(2126, wvt)
        W1815 = getClosestWavelength(1815, wvt)#  
        W2132 = getClosestWavelength(2132, wvt) 
        
        # compute the interpolated continuum values at selected wavelengths between 1815 and 2530
        slope = (R2132 - R1815)/(W2132 - W1815)
        CR1908 = R1815 + slope *(W1908 - W1815)
        CR1914 = R1815 + slope *(W1914 - W1815)
        CR1921 = R1815 + slope *(W1921 - W1815) 
        CR1928 = R1815 + slope *(W1928 - W1815)
        CR1934 = R1815 + slope *(W1934 - W1815)
        CR1941 = R1815 + slope *(W1941 - W1815) 
           
        CR1862 = R1815 + slope*(W1862 - W1815)
        CR1869 = R1815 + slope*(W1869 - W1815)
        CR1875 = R1815 + slope*(W1875 - W1815)    
        CR2112 = R1815 + slope*(W2112 - W1815)
        CR2120 = R1815 + slope*(W2120 - W1815)
        CR2126 = R1815 + slope*(W2126 - W1815)
        img= 1.0-((R1908/CR1908+R1914/CR1914+R1921/CR1921+R1928/CR1928+R1934/CR1934+R1941/CR1941)/(R1862/CR1862+R1869/CR1869+R1875/CR1875+R2112/CR2112+R2120/CR2120+R2126/CR2126))
        nmin = np.nanmin(np.where(img>-np.inf,img,np.nan))
        img = np.where(img>-np.inf,img,nmin)
        return img
    def MIN2295_2480(self):
        cube = self.s
        wvt = self.s_bands
        img1 = getBandDepth(cube, wvt,2165,2295,2364)
        img2 = getBandDepth(cube,wvt,2364,2480,2570)
        img3 = np.empty((np.shape(img1)[0],np.shape(img1)[1],2))
        img3[:,:,0] = img1
        img3[:,:,1] = img2
        img = np.min(img3,axis=2)
        nmin = np.nanmin(np.where(img>-np.inf,img,np.nan))
        img = np.where(img>-np.inf,img,nmin)
        return img
    def MIN2250(self):
        cube = self.s
        wvt = self.s_bands
        img1 = getBandDepth(cube, wvt,2165, 2210, 2350)
        img2 = getBandDepth(cube,wvt,2165, 2265, 2350)
        img3 = np.empty((np.shape(img1)[0],np.shape(img1)[1],2))
        img3[:,:,0] = img1
        img3[:,:,1] = img2
        img = np.min(img3,axis=2)
        nmin = np.nanmin(np.where(img>-np.inf,img,np.nan))
        img = np.where(img>-np.inf,img,nmin)
        return img
    def MIN2345_2537(self):
        cube = self.s
        wvt = self.s_bands
        img1 = getBandDepth(cube, wvt,2250, 2345, 2430)
        print('check')
        img2 = getBandDepth(cube,wvt,2430, 2537, 2602)
        img3 = np.empty((np.shape(img1)[0],np.shape(img1)[1],2))
        img3[:,:,0] = img1
        img3[:,:,1] = img2
        img = np.min(img3,axis=2)
        nmin = np.nanmin(np.where(img>-np.inf,img,np.nan))
        img = np.where(img>-np.inf,img,nmin)
        return img    
    def BDCARB(self):
        cube = self.s
        wvt = self.s_bands
        # extract channels, replacing CRISM_NAN with IEEE NAN
        R2230 = getBand(cube, wvt, 2230)
        R2320 = getBand(cube, wvt, 2320)
        R2330 = getBand(cube, wvt, 2330)
        R2390 = getBand(cube, wvt, 2390)
        R2520 = getBand(cube, wvt, 2520)
        R2530 = getBand(cube, wvt, 2530)
        R2600 = getBand(cube, wvt, 2600)
    
        # identify nearest CRISM wavelengths
        WL1 = getClosestWavelength(2230,wvt)
        WC1 = (getClosestWavelength(2330,wvt)+getClosestWavelength(2320,wvt))*0.5
        WH1 = getClosestWavelength(2390,wvt)
        a =  (WC1 - WL1)/(WH1 - WL1)  # a gets multipled by the longer (higher wvln)  band
        b = 1.0-a                     # b gets multiplied by the shorter (lower wvln) band
    
        WL2 =  getClosestWavelength(2390,wvt)
        WC2 = (getClosestWavelength(2530,wvt) + getClosestWavelength(2520,wvt))*0.5
        WH2 =  getClosestWavelength(2600,wvt)
        c = (WC2 - WL2)/(WH2 - WL2)   # c gets multipled by the longer (higher wvln)  band
        d = 1.0-c                           # d gets multiplied by the shorter (lower wvln) band
    
        # compute bdcarb
        img = 1.0 - (np.sqrt((((R2320 + R2330)*0.5)/(b*R2230 + a*R2390))*(((R2520 + R2530)*0.5)/(d*R2390 + c*R2600))))
        nmin = np.nanmin(np.where(img>-np.inf,img,np.nan))
        img = np.where(img>-np.inf,img,nmin)
        return img
    def SINDEX2(self):
        cube = self.s
        wvt = self.s_bands
        img = getBandDepthInvert(cube,wvt,2120, 2290, 2400,mw=7,hw=3)
        return img
    def BD2100_2(self):
        cube = self.s
        wvt = self.s_bands
        img = getBandDepth(cube,wvt,1930, 2132, 2250,lw=3,hw=3)
        return img
    def BD1900_2(self):
        cube = self.s
        wvt = self.s_bands
        img = getBandDepth(cube,wvt,1850, 1930, 2067)
        return img
    def ISLOPE(self):
        cube = self.s
        wvt = self.s_bands
        # extract individual bands
        R1815 = getBand(cube, wvt, 1815)
        R2530 = getBand(cube, wvt, 2530)
    
        W1815 = getClosestWavelength(1815,wvt)
        W2530 = getClosestWavelength(2530,wvt)
    
        # want in units of reflectance / um
        img = 1000.0 * ( R1815 - R2530 )/ ( W2530 - W1815 )
        nmin = np.nanmin(np.where(img>-np.inf,img,np.nan))
        img = np.where(img>-np.inf,img,nmin)
        return img
    def BD1400(self):
        cube = self.s
        wvt = self.s_bands
        img = getBandDepth(cube,wvt,1330, 1395, 1467,mw=3)
        return img
    def IRR2(self):
        cube = self.s
        wvt = self.s_bands
        img = getBandRatio(cube,wvt,2530,2210)
        return img
    # -------------------------------------------------------------------------
    # Browse Products
    # -------------------------------------------------------------------------
    def MAF(self,norm=False):
        print('calculating MAF browse product')
        tic = timeit.default_timer()
        p1 = self.OLINDEX3()
        p2 = self.LCPINDEX2()
        p3 = self.HCPINDEX2()
        if norm is True:
            p1 = normalizeParameter(p1)
            p2 = normalizeParameter(p2)
            p3 = normalizeParameter(p3)
        img = buildSummary(p1, p2, p3)
        toc = timeit.default_timer()-tic
        print(f'MAF finished in {round(toc,2)} seconds')
        return img
    def FM2(self,norm=False):
        print('calculating FM2 browse product')
        tic = timeit.default_timer()
        p1 = self.BD530_2()
        p2 = self.BD920_2()
        p3 = self.BDI1000VIS()
        if norm is True:
            p1 = normalizeParameter(p1)
            p2 = normalizeParameter(p2)
            p3 = normalizeParameter(p3)
        img = buildSummary(p1, p2, p3)
        toc = timeit.default_timer()-tic
        print(f'FM2 finished in {round(toc/60,2)} minutes')
        return img
    def FAL(self,norm=True):
        print('calculating FAL browse product')
        tic = timeit.default_timer()
        p1 = getBand(self.s,self.s_bands,2529)
        p2 = getBand(self.s,self.s_bands,1506)
        p3 = getBand(self.s,self.s_bands,1080)
        if norm is True:
            p1 = normalizeParameter(p1)
            p2 = normalizeParameter(p2)
            p3 = normalizeParameter(p3)
        img = buildSummary(p1, p2, p3)
        toc = timeit.default_timer()-tic
        print(f'FAL finished in {round(toc,2)} seconds')
        return img
    def PAL(self, norm=False): 
        print('calculating PAL browse product')
        tic = timeit.default_timer()
        p1 = self.BD2210_2()
        p2 = self.BD2190()
        p3 = self.BD2165()
        if norm is True:
            p1 = normalizeParameter(p1)
            p2 = normalizeParameter(p2)
            p3 = normalizeParameter(p3)
        img = buildSummary(p1, p2, p3)
        toc = timeit.default_timer()-tic
        print(f'PAL finished in {round(toc,2)} seconds')
        return img   
    def PHY(self, norm=False): 
        print('calculating PHY browse product')
        tic = timeit.default_timer()
        p1 = self.D2200()
        p2 = self.D2300()
        p3 = self.BD1900r2()
        if norm is True:
            p1 = normalizeParameter(p1)
            p2 = normalizeParameter(p2)
            p3 = normalizeParameter(p3)
        img = buildSummary(p1, p2, p3)
        toc = timeit.default_timer()-tic
        print(f'PHY finished in {round(toc,2)} seconds')
        return img
    def PFM(self, norm=False): 
        print('calculating PFM browse product')
        tic = timeit.default_timer()
        p1 = self.BD2355()
        p2 = self.D2300()
        p3 = self.BD2290()
        if norm is True:
            p1 = normalizeParameter(p1)
            p2 = normalizeParameter(p2)
            p3 = normalizeParameter(p3)
        img = buildSummary(p1, p2, p3)
        toc = timeit.default_timer()-tic
        print(f'PFM finished in {round(toc,2)} seconds')
        return img
    def CR2(self, norm=False): 
        print('calculating CR2 browse product')
        tic = timeit.default_timer()
        p1 = self.MIN2295_2480()
        p2 = self.MIN2345_2537()
        p3 = self.BDCARB()
        if norm is True:
            p1 = normalizeParameter(p1)
            p2 = normalizeParameter(p2)
            p3 = normalizeParameter(p3)
        img = buildSummary(p1, p2, p3)
        toc = timeit.default_timer()-tic
        print(f'CR2 finished in {round(toc,2)} seconds')
        return img
    def HYD(self, norm=False): 
        print('calculating HYD browse product')
        tic = timeit.default_timer()
        p1 = self.SINDEX2()
        p2 = self.BD2100_2()
        p3 = self.BD1900_2()
        if norm is True:
            p1 = normalizeParameter(p1)
            p2 = normalizeParameter(p2)
            p3 = normalizeParameter(p3)
        img = buildSummary(p1, p2, p3)
        toc = timeit.default_timer()-tic
        print(f'HYD finished in {round(toc,2)} seconds')
        return img
    def CHL(self, norm=False): 
        print('calculating CHL browse product')
        tic = timeit.default_timer()
        p1 = self.ISLOPE()
        p2 = self.BD1400()
        p3 = self.IRR2()
        if norm is True:
            p1 = normalizeParameter(p1)
            p2 = normalizeParameter(p2)
            p3 = normalizeParameter(p3)
        img = buildSummary(p1, p2, p3)
        toc = timeit.default_timer()-tic
        print(f'CHL finished in {round(toc,2)} seconds')
        return img
    def HYS(self, norm=False): 
        print('calculating HYS browse product')
        tic = timeit.default_timer()
        p1 = self.MIN2250()
        p2 = self.BD2250()
        p3 = self.BD1900r2()
        if norm is True:
            p1 = normalizeParameter(p1)
            p2 = normalizeParameter(p2)
            p3 = normalizeParameter(p3)
        img = buildSummary(p1, p2, p3)
        toc = timeit.default_timer()-tic
        print(f'HYS finished in {round(toc,2)} seconds')
        return img
    
    def MNF_(self):
        s = self.s
        v = self.v
        tic = timeit.default_timer()
        s_signal = calc_stats(s)
        s_noise = noise_from_diffs(s)
        s_mnfr = mnf(s_signal, s_noise)
        s_mnf10 = s_mnfr.reduce(s, num=10)
        
        v_signal = calc_stats(v)
        v_noise = noise_from_diffs(v)
        v_mnfr = mnf(v_signal, v_noise)
        v_mnf10 = v_mnfr.reduce(v, num=10)
        toc = timeit.default_timer()-tic
        print(f'{round(toc/60,2)} minutes to calculate MNF')
        return s_mnf10, v_mnf10
    
    
    
    
    
    