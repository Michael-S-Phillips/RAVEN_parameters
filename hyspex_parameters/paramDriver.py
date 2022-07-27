#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 10 14:40:28 2022
script to drive the paramCalculator, generates spectral parameter prodcuts from
HySpex data

@author: phillms1
"""
import os
# import sys
import cv2
from param_utils import *
from paramCalculator import paramCalculator
from spectral import *
import spectral.io.envi as envi
# sys.path.append('/Users/phillms1/Documents/Work/RAVEN/RAVEN_parameters/hyspex_parameters/')
data_path = os.path.abspath(os.path.join('/Volumes/HySpex_Back/RAVEN/FieldSeason_2022/TriPod/DryRun220723'))
shdr = data_path + '/Crop/sTargetRocks.hdr'
vhdr = data_path + '/Crop/vTargetRocks.hdr'
# shdr = data_path + '/TrainingData/SWIR_OUTPUT_training/boresight_26_april_1_Mjolnir_S620_SN7064_raw_rad_keystone_smile_bsq_float32_geo.hdr'
# vhdr = data_path + '/TrainingData/VNIR_OUTPUT_training/boresight_26_april_1_Mjolnir_V1240_SN5037_raw_rad_keystone_smile_spectralbinningx2_bsq_float32_geo.hdr'

p = paramCalculator(vhdr,shdr,data_path,crop=True)

date_path='220723' #YYMMDD
savepath = '/Volumes/HySpex_Back/RAVEN/FieldSeason_2022/TriPod/DryRun220723/ParameterOutput/Crop/'+date_path
if os.path.isdir(savepath) is False:
    os.mkdir(savepath)

#%% SWIR parameter block
'''
calculate spectral parameters from SWIR image data, save as .img and as .png
'''
#------------------------------------------------------------------------------
#------------- SWIR Parameters ------------------------------------------------
#------------------------------------------------------------------------------
paramList = ['OLINDEX3','LCPINDEX2','HCPINDEX2','BD1400','BD1900_2','BD1900r2','BD2100_2','BD2165','BD2190','BD2210_2','BD2250','BD2290','BD2355','BDCARB','D2200','D2300','IRR2','ISLOPE','MIN2250','MIN2295_2480','MIN2345_2537','R2529', 'R1506', 'R1080','SINDEX2']
sMeta = p.s_.metadata.copy()
sMeta['wavelength'] = paramList
sMeta['wavelength units'] = 'parameters'
sMeta['default bands'] = ['R2529', 'R1506', 'R1080']
print('calculating SWIR parameters')
sParams=p.calculateSwirParams()
envi.save_image(savepath+'/swir_parameters.hdr',sParams,metadata=sMeta,dtype=np.float32)

#--------- SWIR Browse PNG ----------------------------------------------------
i0 = paramList.index('OLINDEX3')
i1 = paramList.index('LCPINDEX2')
i2 = paramList.index('HCPINDEX2')
MAF = buildSummary(sParams[:,:,i0], sParams[:,:,i1], sParams[:,:,i2])
MAF = np.flip(browse2bit(stretchNBands(cropNZeros(MAF))),axis=2)
n = '/MAF.png'
cv2.imwrite(savepath+n, MAF)

i0 = paramList.index('R2529')
i1 = paramList.index('R1506')
i2 = paramList.index('R1080')
FAL = buildSummary(sParams[:,:,i0], sParams[:,:,i1], sParams[:,:,i2])
FAL = np.flip(browse2bit(stretchNBands(cropNZeros(FAL))),axis=2)
n = '/FAL.png'
cv2.imwrite(savepath+n, FAL)

i0 = paramList.index('BD2210_2')
i1 = paramList.index('BD2190')
i2 = paramList.index('BD2165')
PAL = buildSummary(sParams[:,:,i0], sParams[:,:,i1], sParams[:,:,i2])
PAL = np.flip(browse2bit(stretchNBands(cropNZeros(PAL))),axis=2)
n = '/PAL.png'
cv2.imwrite(savepath+n, PAL)

i0 = paramList.index('D2200')
i1 = paramList.index('D2300')
i2 = paramList.index('BD1900r2')
PHY = buildSummary(sParams[:,:,i0], sParams[:,:,i1], sParams[:,:,i2])
PHY = np.flip(browse2bit(stretchNBands(cropNZeros(PHY))),axis=2)
n = '/PHY.png'
cv2.imwrite(savepath+n, PHY)

i0 = paramList.index('BD2355')
i1 = paramList.index('D2300')
i2 = paramList.index('BD2290')
PFM = buildSummary(sParams[:,:,i0], sParams[:,:,i1], sParams[:,:,i2])
PFM = np.flip(browse2bit(stretchNBands(cropNZeros(PFM))),axis=2)
n = '/PFM.png'
cv2.imwrite(savepath+n, PFM)

i0 = paramList.index('MIN2295_2480')
i1 = paramList.index('MIN2345_2537')
i2 = paramList.index('BDCARB')
CR2 = buildSummary(sParams[:,:,i0], sParams[:,:,i1], sParams[:,:,i2])
CR2 = np.flip(browse2bit(stretchNBands(cropNZeros(CR2))),axis=2)
n = '/CR2.png'
cv2.imwrite(savepath+n, CR2)

i0 = paramList.index('SINDEX2')
i1 = paramList.index('BD2100_2')
i2 = paramList.index('BD1900_2')
HYD = buildSummary(sParams[:,:,i0], sParams[:,:,i1], sParams[:,:,i2])
HYD = np.flip(browse2bit(stretchNBands(cropNZeros(HYD))),axis=2)
n = '/HYD.png'
cv2.imwrite(savepath+n, HYD)

i0 = paramList.index('ISLOPE')
i1 = paramList.index('BD1400')
i2 = paramList.index('IRR2')
CHL = buildSummary(sParams[:,:,i0], sParams[:,:,i1], sParams[:,:,i2])
CHL = np.flip(browse2bit(stretchNBands(cropNZeros(CHL))),axis=2)
n = '/CHL.png'
cv2.imwrite(savepath+n, CHL)

i0 = paramList.index('MIN2250')
i1 = paramList.index('BD2250')
i2 = paramList.index('BD1900r2')
HYS = buildSummary(sParams[:,:,i0], sParams[:,:,i1], sParams[:,:,i2])
HYS = np.flip(browse2bit(stretchNBands(cropNZeros(HYS))),axis=2)
n = '/HYS.png'
cv2.imwrite(savepath+n, HYS)

#%% VIS parameter block
'''
calculate spectral parameters from VIS image data, save as .img and as .png
'''
#------------------------------------------------------------------------------
#------------- VIS Parameters -------------------------------------------------
#------------------------------------------------------------------------------
paramList = ['R637','R550','R463','BD530_2','BD920_2','RPEAK1','BDI1000VIS']
vMeta = p.v_.metadata.copy()
vMeta['wavelength'] = paramList
vMeta['wavelength units'] = 'parameters'
vMeta['default bands'] = ['R637', 'R550', 'R463']
print('calculating VIS parameters')
vParams=p.calculateVisParams()
envi.save_image(savepath+'/vis_parameters.hdr',vParams,metadata=vMeta,dtype=np.float32)

#-- VIS Browse PNG ------------------------------------------------------------
i0 = paramList.index('BD530_2')
i1 = paramList.index('BD920_2')
i2 = paramList.index('RPEAK1')
FM2 = buildSummary(vParams[:,:,i0], vParams[:,:,i1], vParams[:,:,i2])
FM2 = np.flip(browse2bit(stretchNBands(cropNZeros(FM2))),axis=2)
n = '/FM2.png'
cv2.imwrite(savepath+n, FM2)

i0 = paramList.index('R637')
i1 = paramList.index('R550')
i2 = paramList.index('R463')
TRU = buildSummary(vParams[:,:,i0], vParams[:,:,i1], vParams[:,:,i2])
TRU = np.flip(browse2bit(stretchNBands(cropNZeros(TRU))),axis=2)
# TRU = np.flip(browse2bit(TRU),axis=2)
n = '/TRU.png'
cv2.imwrite(savepath+n, TRU)
#%% MNF block
'''
calculate minimum noise fraction transform images and save as PNG and ENVI .img
'''
#SWIR
# s_mnf10 = p.SWIR_MNF()
bs = [4,3,2] #because cv2 writes BGR for some reason
sName = '/SWIR_MNF_234.png'
cv2.imwrite(savepath+sName,s_mnf10[:,:,bs])
sName = '/SWIR_MNF_234_8bit.png'
cv2.imwrite(savepath+sName,browse2bit(s_mnf10[:,:,bs]))

# bandList = ['1','2','3','4','5','6','7','8','9','10']
# sMeta = p.s_.metadata.copy()
# sMeta['wavelength'] = bandList
# sMeta['wavelength units'] = 'MNF Band'
# sMeta['default bands'] = ['1', '2', '3']
# envi.save_image(savepath+'/SWIR_MNF.hdr', s_mnf10, metadata=sMeta,dtype=np.float32)


#Vis
# v_mnf10 = p.VIS_MNF()
bs = [4,3,2]
vName = '/VNIR_MNF_234.png'
cv2.imwrite(savepath+vName,v_mnf10[:,:,bs])
vName = '/VNIR_MNF_234_8bit.png'
cv2.imwrite(savepath+vName,browse2bit(v_mnf10[:,:,bs]))

# bandList = ['1','2','3','4','5','6','7','8','9','10']
# vMeta = p.v_.metadata.copy()
# vMeta['wavelength'] = bandList
# vMeta['wavelength units'] = 'MNF Band'
# vMeta['default bands'] = ['1', '2', '3']
# envi.save_image(savepath+'/VIS_MNF.hdr', v_mnf10, metadata=vMeta,dtype=np.float32)

#%% browse block
'''
generate one-off browse product images (as PNG files)
'''
savepath = '/Users/phillms1/Documents/Work/RAVEN/RAVEN_parameters/hyspex_parameters/BrowseProducts/'+date_path
if os.path.isdir(savepath) is False:
    os.mkdir(savepath)

# calculate and save browse products as .png
bp = p.MAF()
n = 'MAF'+'.png'
img = np.flip(browse2bit(stretchNBands(cropNZeros(bp))),axis=2)
cv2.imwrite(savepath+n,img)
del(img)

bp = p.FM2()
n = 'FM2'+'.png'
img = np.flip(browse2bit(stretchNBands(cropNZeros(bp))),axis=2)
cv2.imwrite(savepath+n,img)
del(img)

bp = p.FAL()
n = 'FAL'+'.png'
img = np.flip(browse2bit(stretchNBands(cropNZeros(bp))),axis=2)
cv2.imwrite(savepath+n,img)
del(img)

bp = p.PAL()
n = 'PAL'+'.png'
img = np.flip(browse2bit(stretchNBands(cropNZeros(bp))),axis=2)
cv2.imwrite(savepath+n,img)
del(img)

bp = p.PFM()
n = 'PFM'+'.png'
img = np.flip(browse2bit(stretchNBands(cropNZeros(bp))),axis=2)
cv2.imwrite(savepath+n,img)
del(img)

bp = p.PHY()
n = 'PHY'+'.png'
img = np.flip(browse2bit(stretchNBands(cropNZeros(bp))),axis=2)
cv2.imwrite(savepath+n,img)
del(img)

bp = p.CR2()
n = 'CR2'+'.png'
img = np.flip(browse2bit(stretchNBands(cropNZeros(bp))),axis=2)
cv2.imwrite(savepath+n,img)
del(img)

bp = p.HYD()
n = 'HYD'+'.png'
img = np.flip(browse2bit(stretchNBands(cropNZeros(bp))),axis=2)
cv2.imwrite(savepath+n,img)
del(img)

bp = p.CHL()
n = 'CHL'+'.png'
img = np.flip(browse2bit(stretchNBands(cropNZeros(bp))),axis=2)
cv2.imwrite(savepath+n,img)
del(img)

bp = p.HYS()
n = 'HYS'+'.png'
img = np.flip(browse2bit(stretchNBands(cropNZeros(bp))),axis=2)
cv2.imwrite(savepath+n,img)
del(img)


#%% utilities

# crop metadata
vCropMeta = p.v_.metadata.copy()
vCropMeta['lines']=np.shape(vCrop)[1]
vCropMeta['samples']=np.shape(vCrop)[0]

sCropMeta = p.s_.metadata.copy()
sCropMeta['lines']=np.shape(sCrop)[1]
sCropMeta['samples']=np.shape(sCrop)[0]
