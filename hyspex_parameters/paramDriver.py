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
data_path = os.path.abspath(os.path.join('/Volumes/Utumno/HySpex_training_data/RockGarden2022/RAD'))
shdr = data_path + '/RAD_SWIR/rockGardenTest_Mjolnir_S620_SN7062_39781us_2022-06-01T195010_raw_rad_keystone_smile_float32.hdr'
vhdr = data_path + '/RAD_VNIR/rockGardenTest_Mjolnir_V1240_SN5011_19784us_2022-06-01T195010_raw_rad_keystone_smile_float32.hdr'

p = paramCalculator(vhdr,shdr,data_path)

#%% browse block
savepath = '/Users/phillms1/Documents/Work/RAVEN/RAVEN_parameters/hyspex_parameters/BrowseProducts/'
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
