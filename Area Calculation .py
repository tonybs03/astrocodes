# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 19:13:10 2021

@author: tonyb
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
import photutils as pt
from scipy.optimize import curve_fit
from numpy.linalg import inv
from matplotlib.colors import LogNorm
from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from astropy import units as u
from astropy.coordinates import SkyCoord
from reproject import reproject_interp, reproject_exact
from reproject.mosaicking import reproject_and_coadd
from reproject.mosaicking import find_optimal_celestial_wcs
import FITS_tools
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import matplotlib.patches as patches
from astropy.coordinates import Angle
import xlsxwriter

RA0, DEC0, C, Quadrant  = np.genfromtxt ('ALL COORDS IN DEG.cvs', unpack = True)
Q0,Q1,Q2,Q3,Q4,Q5,Q6,Q7 = np.genfromtxt ('Chip Edges.cvs', unpack = True)

Xo= 146.41250
Yo = -31.19110
COSINE = np.cos(np.deg2rad(Yo))

edge = SkyCoord('9h45m17.5s -31d23m59.55s', frame='icrs')





lonR0, lonL0, latd0, latu0 = Q0
lonR0, lonL0 = (lonR0-Xo)*COSINE, (lonL0-Xo)*COSINE
latd0, latu0 = latd0-Yo, latu0-Yo

lonR1, lonL1, latd1, latu1 = Q1
lonR1, lonL1 = (lonR1-Xo)*COSINE, (lonL1-Xo)*COSINE
latd1, latu1 = latd1-Yo, latu1-Yo

lonR2, lonL2, latd2, latu2 = Q2
lonR2, lonL2 = (lonR2-Xo)*COSINE, (lonL2-Xo)*COSINE
latd2, latu2 = latd2-Yo, latu2-Yo

lonR3, lonL3, latd3, latu3 = Q3
lonR3, lonL3 = (lonR3-Xo)*COSINE, (lonL3-Xo)*COSINE
latd3, latu3 = latd3-Yo, latu3-Yo

lonR4, lonL4, latd4, latu4 = Q4
lonR4, lonL4 = (lonR4-Xo)*COSINE, (lonL4-Xo)*COSINE
latd4, latu4 = latd4-Yo, latu4-Yo

lonR5, lonL5, latd5, latu5 = Q5
lonR5, lonL5 = (lonR5-Xo)*COSINE, (lonL5-Xo)*COSINE
latd5, latu5 = latd5-Yo, latu5-Yo

lonR6, lonL6, latd6, latu6 = Q6
lonR6, lonL6 = (lonR6-Xo)*COSINE, (lonL6-Xo)*COSINE
latd6, latu6 = latd6-Yo, latu6-Yo

lonR7, lonL7, latd7, latu7 = Q7
lonR7, lonL7 = (lonR7-Xo)*COSINE, (lonL7-Xo)*COSINE
latd7, latu7 = latd7-Yo, latu7-Yo



fig, ax = plt.subplots(figsize=(15, 15))
#CREATE RECTANGLE PATCHES
x0 = lonR0
y0 = latd0
w0 = np.abs(lonL0-lonR0)
h0 = np.abs(latd0-latu0)
rect = patches.Rectangle((x0, y0), w0, h0, linewidth=2, edgecolor='r', facecolor='none')
ax.add_patch(rect)

x1 = lonR1
y1 = latd1
w1 = np.abs(lonL1-lonR1)
h1 = np.abs(latd1-latu1)
rect = patches.Rectangle((x1, y1), w1, h1, linewidth=2, edgecolor='r', facecolor='none')
ax.add_patch(rect)

x2 = lonR2
y2 = latd2
w2 = np.abs(lonL2-lonR2)
h2 = np.abs(latd2-latu2)
rect = patches.Rectangle((x2, y2), w2, h2, linewidth=2, edgecolor='r', facecolor='none')
ax.add_patch(rect)

x3 = lonR3
y3 = latd3
w3 = np.abs(lonL3-lonR3)
h3 = np.abs(latd3-latu3)
rect = patches.Rectangle((x3, y3), w3, h3, linewidth=2, edgecolor='r', facecolor='none')
ax.add_patch(rect)

x4 = lonR4
y4 = latd4
w4 = np.abs(lonL4-lonR4)
h4 = np.abs(latd4-latu4)
rect = patches.Rectangle((x4, y4), w4, h4, linewidth=2, edgecolor='r', facecolor='none')
ax.add_patch(rect)

x5 = lonR5
y5 = latd5
w5 = np.abs(lonL5-lonR5)
h5 = np.abs(latd5-latu5)
rect = patches.Rectangle((x5, y5), w5, h5, linewidth=2, edgecolor='r', facecolor='none')
ax.add_patch(rect)

x6 = lonR6
y6 = latd6
w6 = np.abs(lonL6-lonR6)
h6 = np.abs(latd6-latu6)
rect = patches.Rectangle((x6, y6), w6, h6, linewidth=2, edgecolor='r', facecolor='none')
ax.add_patch(rect)

x7 = lonR7
y7 = latd7
w7 = np.abs(lonL7-lonR7)
h7 = np.abs(latd7-latu7)
rect = patches.Rectangle((x7, y7), w7, h7, linewidth=2, edgecolor='r', facecolor='none')
ax.add_patch(rect)

plt.rcParams.update({'font.size':13})
ax.set_xlim(lonL2+0.1, lonR7-0.1)
ax.set_ylim(latd3-0.1, latu4+0.1)

ax.plot(0, 0, 'D', mfc='none', markeredgewidth=2.0, 
         color='purple', markersize=10)

circle0 = plt.Circle((0, 0), 0.5, color='lime', fill=False)
circle1 = plt.Circle((0, 0), 1.0, color='lime', fill=False)
circle2 = plt.Circle((0, 0), 1.5, color='lime', fill=False)
circle3 = plt.Circle((0, 0), 2.0, color='lime', fill=False)
circle4 = plt.Circle((0, 0), 2.5, color='lime', fill=False)

ax.add_patch(circle0)
ax.add_patch(circle1)
ax.add_patch(circle2)
ax.add_patch(circle3)
ax.add_patch(circle4)
plt.xlabel('RA (DEG)')
plt.ylabel('DEC (DEG)')
plt.rcParams.update({'font.size':15})



#WORKING WITH UNITS
NofP = 241
NofD = 1.75
XLEFT = NofD
XRIGHT = -NofD
YDOWN = -NofD
YUP = NofD

PSize = (np.abs(YUP-YDOWN)/NofP)


XXX = np.linspace(XLEFT, XRIGHT, NofP)
XXX = (XXX - XXX[int(NofP/2)])
fctrX = np.abs(((NofP-1)/2)/XXX[0])
XXX = XXX*(((NofP-1)/2)/XXX[0])

YYY = np.linspace(YDOWN, YUP, NofP)
YYY = (YYY - YYY[int(NofP/2)])
fctrY = np.abs(((NofP-1)/2)/YYY[0])
YYY = YYY*(((NofP-1)/2)/YYY[0])


fig2, ax2 = plt.subplots(figsize=(15, 15))    



#UPDATE PATCHES
Q0 = np.array([lonR0*fctrY, lonL0*fctrY, latd0*fctrY, latu0*fctrY])
lonR0, lonL0, latd0, latu0 = Q0

Q1 = np.array([lonR1*fctrY, lonL1*fctrY, latd1*fctrY, latu1*fctrY])
lonR1, lonL1, latd1, latu1 = Q1

Q2 = np.array([lonR2*fctrY, lonL2*fctrY, latd2*fctrY, latu2*fctrY])
lonR2, lonL2, latd2, latu2 = Q2

Q3 = np.array([lonR3*fctrY, lonL3*fctrY, latd3*fctrY, latu3*fctrY])
lonR3, lonL3, latd3, latu3 = Q3

Q4 = np.array([lonR4*fctrY, lonL4*fctrY, latd4*fctrY, latu4*fctrY])
lonR4, lonL4, latd4, latu4 = Q4

Q5 = np.array([lonR5*fctrY, lonL5*fctrY, latd5*fctrY, latu5*fctrY])
lonR5, lonL5, latd5, latu5 = Q5

Q6 = np.array([lonR6*fctrY, lonL6*fctrY, latd6*fctrY, latu6*fctrY])
lonR6, lonL6, latd6, latu6 = Q6

Q7 = np.array([lonR7*fctrY, lonL7*fctrY, latd7*fctrY, latu7*fctrY])
lonR7, lonL7, latd7, latu7 = Q7



#CREATE RECTANGLE PATCHES
x0 = lonR0
y0 = latd0
w0 = np.abs(lonL0-lonR0)
h0 = np.abs(latd0-latu0)
rect = patches.Rectangle((x0, y0), w0, h0, linewidth=2, edgecolor='r', facecolor='none')
ax2.add_patch(rect)

x1 = lonR1
y1 = latd1
w1 = np.abs(lonL1-lonR1)
h1 = np.abs(latd1-latu1)
rect = patches.Rectangle((x1, y1), w1, h1, linewidth=2, edgecolor='r', facecolor='none')
ax2.add_patch(rect)

x2 = lonR2
y2 = latd2
w2 = np.abs(lonL2-lonR2)
h2 = np.abs(latd2-latu2)
rect = patches.Rectangle((x2, y2), w2, h2, linewidth=2, edgecolor='r', facecolor='none')
ax2.add_patch(rect)

x3 = lonR3
y3 = latd3
w3 = np.abs(lonL3-lonR3)
h3 = np.abs(latd3-latu3)
rect = patches.Rectangle((x3, y3), w3, h3, linewidth=2, edgecolor='r', facecolor='none')
ax2.add_patch(rect)

x4 = lonR4
y4 = latd4
w4 = np.abs(lonL4-lonR4)
h4 = np.abs(latd4-latu4)
rect = patches.Rectangle((x4, y4), w4, h4, linewidth=2, edgecolor='r', facecolor='none')
ax2.add_patch(rect)

x5 = lonR5
y5 = latd5
w5 = np.abs(lonL5-lonR5)
h5 = np.abs(latd5-latu5)
rect = patches.Rectangle((x5, y5), w5, h5, linewidth=2, edgecolor='r', facecolor='none')
ax2.add_patch(rect)

x6 = lonR6
y6 = latd6
w6 = np.abs(lonL6-lonR6)
h6 = np.abs(latd6-latu6)
rect = patches.Rectangle((x6, y6), w6, h6, linewidth=2, edgecolor='r', facecolor='none')
ax2.add_patch(rect)

x7 = lonR7
y7 = latd7
w7 = np.abs(lonL7-lonR7)
h7 = np.abs(latd7-latu7)
rect = patches.Rectangle((x7, y7), w7, h7, linewidth=2, edgecolor='r', facecolor='none')
ax2.add_patch(rect)

GoodX = []
GoodY = []
for i in XXX:
    for j in YYY:
        if (Q0[0] < i < Q0[1] and Q0[2] < j < Q0[3]) or (
                Q1[0] < i < Q1[1] and Q1[2] < j < Q1[3]) or (
                    Q2[0] < i < Q2[1] and Q2[2] < j < Q2[3]) or (
                        Q3[0] < i < Q3[1] and Q3[2] < j < Q3[3]) or (
                            Q4[0] < i < Q4[1] and Q4[2] < j < Q4[3]) or (
                                Q5[0] < i < Q5[1] and Q5[2] < j < Q5[3]) or (
                                    Q6[0] < i < Q6[1] and Q6[2] < j < Q6[3]) or (
                                        Q7[0] < i < Q7[1] and Q7[2] < j < Q7[3]):
            GoodX.append(i)
            GoodY.append(j)

GoodX = np.array(GoodX)
GoodY = np.array(GoodY)                            
         


BIN0 = 0
BIN1 = 0
BIN2 = 0
BIN3 = 0   
BIN4 = 0     
BIN5 = 0           
for i in np.arange(0, len(GoodX)):
    R = np.sqrt((GoodX[i]**2)+(GoodY[i])**2)*PSize
    if 0.1 <= R < 0.4:
        plt.scatter(GoodX[i], GoodY[i], s=1, marker='s', c='black')
        BIN0=BIN0+1  
    if 0.4 <= R < 0.7:
        plt.scatter(GoodX[i], GoodY[i], s=1, marker='s', c='red')
        BIN1=BIN1+1   
    if 0.7 <= R < 1.0:
        plt.scatter(GoodX[i], GoodY[i], s=1, marker='s', c='green')
        BIN2=BIN2+1
    if 1.0 <= R < 1.3:
        plt.scatter(GoodX[i], GoodY[i], s=1, marker='s', c='cyan')
        BIN3=BIN3+1 
    if 1.3 <= R < 1.6:
        plt.scatter(GoodX[i], GoodY[i], s=1, marker='s', c='purple')
        BIN4=BIN4+1
    if 1.6 <= R < 2.2:
        plt.scatter(GoodX[i], GoodY[i], s=1, marker='s', c='yellow')
        BIN5=BIN5+1

circleo = plt.Circle((0, 0), 0.1/PSize, color='lime', fill=False)
circle0 = plt.Circle((0, 0), 0.4/PSize, color='lime', fill=False)
circle1 = plt.Circle((0, 0), 0.7/PSize, color='lime', fill=False)
circle2 = plt.Circle((0, 0), 1.0/PSize, color='lime', fill=False)
circle3 = plt.Circle((0, 0), 1.3/PSize, color='lime', fill=False)
circle4 = plt.Circle((0, 0), 1.6/PSize, color='lime', fill=False)
circle5 = plt.Circle((0, 0), 2.2/PSize, color='lime', fill=False)

ax2.add_patch(circleo)
ax2.add_patch(circle0)
ax2.add_patch(circle1)
ax2.add_patch(circle2)
ax2.add_patch(circle3)
ax2.add_patch(circle4)
ax2.add_patch(circle5)

ax2.set_xlim(XXX[0], XXX[int(len(XXX)-1)])
ax2.set_ylim(YYY[int(len(YYY)-1)], YYY[0])
    
ax2.plot(0, 0, 'x', mfc='none', markeredgewidth=2.0, 
          color='cyan', markersize=10)

A0 = BIN0*PSize**2
A1 = BIN1*PSize**2
A2 = BIN2*PSize**2
A3 = BIN3*PSize**2
A4 = BIN4*PSize**2
A5 = BIN5*PSize**2

AREA = [A0, A1, A2, A3, A4, A5]
txtdata = np.array([AREA])
txtdata = txtdata.T
with open('AREA.cvs', 'w+') as datafile_id0:
    np.savetxt(datafile_id0, txtdata, fmt='%1.5f')
    
    
RA = (RA0-Xo)*COSINE
DEC = DEC0-Yo
RA, DEC = RA*fctrY, DEC*fctrY

for i in range(0,len(RA)):
    if C[i] == 1:
        plt.plot(RA[i],DEC[i],'o', mfc='none', markeredgewidth=2.0, 
              color='red', markersize=6)
    if C[i] == 2:
        plt.plot(RA[i],DEC[i],'o', mfc='none', markeredgewidth=2.0, 
              color='red', markersize=6)     
    if C[i] == 3:
        plt.plot(RA[i],DEC[i],'o', mfc='none', markeredgewidth=2.0, 
              color='black', markersize=6)

ax2.tick_params(labelbottom=False, labelleft=False, bottom=False, left=False)

# Convert celsius to Fahrenheit
T_f = lambda T_c: (T_c*PSize)
# Convert Fahrenheit to Celsius
T_c = lambda T_f: (T_f*PSize)

axb = ax2.secondary_xaxis("bottom", functions=(T_f, T_c))
axb.tick_params(which='major', direction='in', length=16)
axb.tick_params(which='minor', direction='in', length=8)
axb.xaxis.set_minor_locator(MultipleLocator(0.1))
axb.set_xlabel('R.A. [deg]', fontsize=24, family='Times New Roman', labelpad=12)
axb.tick_params(labelsize=18, pad=3)

axl = ax2.secondary_yaxis("left", functions=(T_f, T_c))
axl.tick_params(which='major', direction='in', length=16)
axl.tick_params(which='minor', direction='in', length=8)
axl.yaxis.set_minor_locator(MultipleLocator(0.1))
axl.set_ylabel('Decl. [deg]', fontsize=24, family='Times New Roman', labelpad=12)
axl.tick_params(labelsize=18, pad=3)





# Convert celsius to Fahrenheit
T_f = lambda T_c: (T_c*PSize)*0.21235
# Convert Fahrenheit to Celsius
T_c = lambda T_f: (T_f)/(0.21235*PSize)

axr = ax2.secondary_yaxis("right", functions=(T_f, T_c))
axr.tick_params(which='major', direction='in', length=16)
axr.tick_params(which='minor', direction='in', length=8)
axr.yaxis.set_minor_locator(MultipleLocator(0.05))
axr.set_ylabel('Decl. [Mpc]', fontsize=24, family='Times New Roman', labelpad=12)
axr.tick_params(labelsize=18, pad=3)

axt = ax2.secondary_xaxis("top", functions=(T_f, T_c))
axt.tick_params(which='major', direction='in', length=16)
axt.tick_params(which='minor', direction='in', length=8)
axt.set_xlabel('R.A. [Mpc]', fontsize=24, family='Times New Roman', labelpad=12)
axt.tick_params(labelsize=18, pad=3)
axt.minorticks_on()
axt.yaxis.set_minor_locator(MultipleLocator(0.05))

plt.savefig('Area Plot.pdf', bbox_inches='tight')























