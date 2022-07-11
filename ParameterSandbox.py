#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  9 13:11:50 2022

@author: phillms1
"""

import rasterio as rio
import numpy as np
import sys
import re
import os
import timeit
import IPython.display as Disp
import cv2
import scipy.integrate as it
from pylab import *
from matplotlib.widgets import LassoSelector
from matplotlib import pyplot as plt 
from ipywidgets import widgets
from spectral import *
import spectral.io.envi as envi

data_path = os.path.abspath(os.path.join('/Volumes/Utumno/HySpex_training_data/RockGarden2022/RAD'))
# plt.rcParams["figure.figsize"] = (20,10)


tic = timeit.default_timer()
# get file paths
sfile = data_path + '/RAD_SWIR/rockGardenTest_Mjolnir_S620_SN7062_39781us_2022-06-01T195010_raw_rad_keystone_smile_float32.img'
vfile = data_path + '/RAD_VNIR/rockGardenTest_Mjolnir_V1240_SN5011_19784us_2022-06-01T195010_raw_rad_keystone_smile_float32.img'
shdr = data_path + '/RAD_SWIR/rockGardenTest_Mjolnir_S620_SN7062_39781us_2022-06-01T195010_raw_rad_keystone_smile_float32.hdr'
vhdr = data_path + '/RAD_VNIR/rockGardenTest_Mjolnir_V1240_SN5011_19784us_2022-06-01T195010_raw_rad_keystone_smile_float32.hdr'

# load data
s_ = envi.open(shdr)
v_ = envi.open(vhdr)
s = s_.load()
s = np.flip(np.transpose(s,(1,0,2)),axis=0)
v = v_.load()
v = np.flip(np.transpose(v,(1,0,2)),axis=0)
toc = timeit.default_timer()-tic
print(f'{np.round(toc/60,2)} minutes')


# retrieve bands
s_bands = s_.metadata['wavelength']
s_bands = [float(b) for b in s_bands]
v_bands = v_.metadata['wavelength']
v_bands = [float(b) for b in v_bands]

# default preview bands
s_preview_bands = s_.metadata['default bands']
s_preview_bands = [int(b) for b in s_preview_bands]
v_preview_bands = v_.metadata['default bands']
v_preview_bands = [int(b) for b in v_preview_bands]

# preview data
plt.subplot(2,1,1)
plt.title('VIS Preview')
v_preview = v[:,:,v_preview_bands]
v_preview = [(r-np.min(r))/(np.max(r)-np.min(r)) for r in v_preview]
plt.imshow(v_preview)
plt.subplot(2,1,2)
plt.title('SWIR Preview')
s_preview = s[:,:,s_preview_bands]
s_preview = [(r-np.min(r))/(np.max(r)-np.min(r)) for r in s_preview]
plt.imshow(s_preview)
plt.show()

#%% bounding box select

# cv bounding box select tool
class bbox_select():
    # %matplotlib notebook 


    def __init__(self,im):
        self.im = im
        self.selected_points = []
        self.fig,ax = plt.subplots()
        self.img = ax.imshow(self.im.copy())
        self.ka = self.fig.canvas.mpl_connect('button_press_event', self.onclick)
        disconnect_button = widgets.Button(description="Disconnect mpl")
        Disp.display(disconnect_button)
        disconnect_button.on_click(self.disconnect_mpl)


        
    def poly_img(self,img,pts):
        pts = np.array(pts, np.int32)
        pts = pts.reshape((-1,1,2))
        cv2.polylines(img,[pts],True,(np.random.randint(0,255),np.random.randint(0,255),np.random.randint(0,255)),7)
        return img

    def onclick(self, event):
    #display(str(event))
        self.selected_points.append([event.xdata,event.ydata])
        if len(self.selected_points)>1:
            self.fig
            self.img.set_data(self.poly_img(self.im.copy(),self.selected_points))
    def disconnect_mpl(self,_):
        self.fig.canvas.mpl_disconnect(self.ka)
        
bs = bbox_select(v_preview)
arr = np.array([bs.selected_points],'int')
