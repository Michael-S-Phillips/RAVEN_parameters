#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 10 14:40:28 2022

@author: phillms1
"""
import os
# import sys
import cv2
from param_utils import browse2bit
from paramCalculator import paramCalculator
# sys.path.append('/Users/phillms1/Documents/Work/RAVEN/RAVEN_parameters/hyspex_parameters/')
data_path = os.path.abspath(os.path.join('/Volumes/Utumno/HySpex_training_data/'))
shdr = data_path + '/TrainingData/SWIR_OUTPUT_training/boresight_26_april_1_Mjolnir_S620_SN7064_raw_rad_keystone_smile_bsq_float32_geo.hdr'
vhdr = data_path + '/TrainingData/VNIR_OUTPUT_training/boresight_26_april_1_Mjolnir_V1240_SN5037_raw_rad_keystone_smile_spectralbinningx2_bsq_float32_geo.hdr'

p = paramCalculator(vhdr,shdr,data_path)

#%% browse block
savepath = '/Users/phillms1/Documents/Work/RAVEN/RAVEN_parameters/hyspex_parameters/BrowseProducts/Aerial/'
if os.path.isdir(savepath) is False:
    os.mkdir(savepath)

# calculate and save browse products
bp = p.MAF()
n = 'MAF'+'.png'
cv2.imwrite(savepath+n,browse2bit(bp))
del(bp)

bp = p.FM2()
n = 'FM2'+'.png'
cv2.imwrite(savepath+n,browse2bit(bp))

bp = p.FAL()
n = 'FAL'+'.png'
cv2.imwrite(savepath+n,browse2bit(bp))
del(bp)

bp = p.PAL()
n = 'PAL'+'.png'
cv2.imwrite(savepath+n,browse2bit(bp))
del(bp)

bp = p.PFM()
n = 'PFM'+'.png'
cv2.imwrite(savepath+n,browse2bit(bp))
del(bp)

bp = p.PHY()
n = 'PHY'+'.png'
cv2.imwrite(savepath+n,browse2bit(bp))

bp = p.CR2()
n = 'CR2'+'.png'
cv2.imwrite(savepath+n,browse2bit(bp))

bp = p.HYD()
n = 'HYD'+'.png'
cv2.imwrite(savepath+n,browse2bit(bp))

bp = p.CHL()
n = 'CHL'+'.png'
cv2.imwrite(savepath+n,browse2bit(bp))

bp = p.HYS()
n = 'HYS'+'.png'
cv2.imwrite(savepath+n,browse2bit(bp))

#%% MNF block

s_mnf10, v_mnf10 = p.MNF_()
bs = [2,1,0]
sName = 'SWIR_MNF.png'
cv2.imwrite(savepath+sName,s_mnf10[:,:,bs])
sName = 'SWIR_MNF_8bit.png'
cv2.imwrite(savepath+sName,browse2bit(s_mnf10[:,:,bs],cropZero=True,norm=True))
vName = 'VNIR_MNF.png'
cv2.imwrite(savepath+vName,v_mnf10[:,:,bs])
vName = 'VNIR_MNF_8bit.png'
cv2.imwrite(savepath+vName,browse2bit(v_mnf10[:,:,bs],cropZero=True,norm=True))
