#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 12 12:45:04 2022

@author: phillms1
"""
import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from oreXpressParamCalculator import oreXpressParamCalculator

rootPath = '/Users/phillms1/Documents/Work/RAVEN/RAVEN_parameters/OreXpressParameters/'
dataFile = 'CompiledFieldSpecta_111219.csv'
df = pd.read_csv(rootPath + dataFile)
wvt = df.iloc[:,0]
spectra = df.iloc[:,1:]/100
specNames = list(spectra.columns)
paramNames = ['OLINDEX3','LCPINDEX2','HCPINDEX2','BDI1000VIS','BD530_2','BD920_2','BD2210_2','BD2190','D2165','BD2250','BD2355','BD2290','D2300','D2200','BD1900r2','MIN2295_248','MIN2345_253','BDCARB','SINDEX2','BD2100_2','BD1900_2','ISLOPE','BD1400','IRR2','MIN2250','BD2250','BD1900r2']

pdf_ = np.empty((27,len(specNames)))
# pd.DataFrame(paramNames)
# paramDF.columns = list(specNames)
i = 0
for spectrum in specNames:
    
    opc = oreXpressParamCalculator(wvt,spectra.loc[:,spectrum])
    pdf_[0,i] = opc.OLINDEX3()
    pdf_[1,i] = opc.LCPINDEX2()
    pdf_[2,i] = opc.HCPINDEX2()
    # FM2
    pdf_[3,i] = opc.BDI1000VIS()
    pdf_[4,i] = opc.BD530_2()
    pdf_[5,i] = opc.BD920_2()
    # PAL
    pdf_[6,i] = opc.BD2210_2()
    pdf_[7,i] = opc.BD2190()
    pdf_[8,i] = opc.BD2165()
    # PFM
    pdf_[9,i] = opc.BD2250()
    pdf_[10,i] = opc.BD2355()
    pdf_[11,i] = opc.BD2290()
    # PHY
    pdf_[12,i] = opc.D2300()
    pdf_[13,i] = opc.D2200()
    pdf_[14,i] = opc.BD1900r2()
    # CR2
    pdf_[15,i] = opc.MIN2295_2480()
    pdf_[16,i] = opc.MIN2345_2537()
    pdf_[17,i] = opc.BDCARB()
    # HYD
    pdf_[18,i] = opc.SINDEX2()
    pdf_[19,i] = opc.BD2100_2()
    pdf_[20,i] = opc.BD1900_2()
    # CHL
    pdf_[21,i] = opc.ISLOPE()
    pdf_[22,i] = opc.BD1400()
    pdf_[23,i] = opc.IRR2()
    # HYS
    pdf_[24,i] = opc.MIN2250()
    pdf_[25,i] = opc.BD2250()
    pdf_[26,i] = opc.BD1900r2()
    
    i=i+1

paramDF = pd.DataFrame(pdf_,columns=specNames)


#%% plot cell
import os
dateStr = '220712' #YYMMDD
figSavePath = rootPath+'Output/'+dateStr+'/'

if os.path.isdir(figSavePath) is False:
    os.mkdir(figSavePath)
    
i=101
plt.figure(num=1,figsize=(6.5,4),dpi=300)
plt.bar(np.linspace(1,len(paramNames),len(paramNames)),paramDF.iloc[:,i])
plt.title(specNames[i]+' Parameter Values')
plt.xticks(np.linspace(1,len(paramNames),len(paramNames)),labels=paramNames,rotation=90)
plt.subplots_adjust(bottom=0.3)
plt.savefig(figSavePath+specNames[i]+'_Bar.png', dpi=300)

plt.figure(num=2,figsize=(6.5,4),dpi=300)
plt.plot(wvt,spectra.iloc[:,i])
plt.title('Spectrum '+specNames[i])
plt.xlabel('Wavelength (nm)')
plt.savefig(figSavePath+specNames[i]+'_Spectrum.png', dpi=300)







