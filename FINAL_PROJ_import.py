#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 14:31:51 2019
    EART/ASTR 119

@author: Andrew Quartuccio
"""

# =============================================================================
#                           Module Imports
# =============================================================================
from __future__ import division
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import Modules.opt_utils as opt

# =============================================================================
#                           File/Directory Parmeters for Fig 1,2
# =============================================================================


file_in  = './Data/hygdata_v3.csv'
file_out = './Data/star_location.gif'
dir_in   = './Data'
pi       = np.pi

hyg_data   = np.genfromtxt( file_in, dtype = float, delimiter = ',', skip_header = 1, usecols = (16, 33)).T
hyg_ptdata = np.genfromtxt( file_in, dtype = float, delimiter = ',', skip_header = 1, usecols = (17,18,19,13)).T
hyg_vptdata= np.genfromtxt( file_in, dtype = float, delimiter = ',', skip_header = 1, usecols = (20,21,22)).T

m_X = hyg_ptdata[0]
m_Y = hyg_ptdata[1]
m_Z = hyg_ptdata[2]


m_vX = hyg_vptdata[0]
m_vY = hyg_vptdata[1]
m_vZ = hyg_vptdata[2]

rho= ((m_X**2+m_Y**2+m_Z**2)**.5)
mag= hyg_ptdata[3]


#=======================================================================================
def lin_LS( aX, aY):
    """
    - linear least squares assuming normal distributed errors in Y, no errors in X

    :param aX: - independent variable
    :param aY: - measured dependent variable

    :return: float(<slope>)
    """
    meanX = aX.mean()
    meanY = aY.mean()
    # variance and co-variance - 1./N term omitted because it cancels in the following ratio
    VarX  = ( (aX - meanX)**2).sum()
    #VarY  = ( (aY - meanY)**2).sum()
    CovXY = ( (aY-meanY)*(aX-meanX)).sum()
    slope = CovXY/VarX
    a     = meanY - meanX*slope
    return slope, a
#=======================================================================================



ci= hyg_data[0]
Lum= hyg_data[1]


sigma= 5.670*1e-8
T= 4600*(1/((.92*ci)+1.7)+(1/((.92*ci)+.62)))
R= (np.sqrt(Lum/(4*pi*sigma*(T**4))))/1000

# Power Law fit
tmin, tmax = 1, 1e5
sel   = np.logical_and( T[0::] >= tmin, T[0::] <= tmax)
slope, f_a = lin_LS( np.log10( T[sel]), np.log10( Lum[sel]))

aX_fit = np.linspace( tmin*.5, tmax*5, 100)
aPLfit = 10**(f_a)*aX_fit**slope

# Power Law Plot
plt.figure(1)
ax = plt.subplot(111)
ax.set_title( 'Power-Law Fit')
ax.loglog( T, Lum, 'ko', label = 'data')
ax.loglog( aX_fit, aPLfit, 'r--', label = 'L ~ T^(%.2f)'%( round( slope, 2)))
ax.legend( loc = 'lower left')
plt.ax.set_xlim(ax, 0, 100000)
plt.ax.set_ylim(ax, 0, 100000)
ax.set_xlabel( 'Temperature [degree C]')
ax.set_ylabel( 'Luminosity [solar units]')
plt.show()

#plt.figure(1)
#ax = plt.subplot()
#plt.plot(rho, mag, 'r.')
#plt.Axes.set_xlim(ax, 0, 1100)
#plt.Axes.set_ylim(ax, 0, 15)
#plt.show()



#for t in range(0, 10000001, 250000):
#    m_X= m_vX*t+m_X
#    m_Y= m_vY*t+m_Y
#    m_X= m_vZ*t+m_Z
#
#    fig2 = plt.figure(2, figsize=(20,20))
#    print 'Time Step: ', int( t), 'years'  #,(n+1)/plot_step,(n+1)%plot_step
#    ax = axes3d.Axes3D(fig2)
#
#    for i in range(0,len(m_X)):
#        if (abs(m_X[i]) <= 1000):
#            m_X[i] = 0
#            m_Y[i] = 0
#            m_Z[i] = 0
#
#
#    ax.scatter( m_X, m_Y, m_Z, c='r', marker='o')
#    ax.scatter(0,0,0, color='blue', marker='o')
#    plt.pause(0.08)
#    ax.set_title( 'Time Step: %.1f [yrs]'%( t))
##    plt.savefig( file_out)

#
## =============================================================================
##           File/Directory Parmeters for Hertz-Sprung
## =============================================================================
#
#file_in = './Data/hygdata_v3.csv'
#dir_in  = './data'
#
#df = pd.read_table(file_in, delimiter=',', header=0, index_col = 0, usecols = ( 0, 14, 9, 16, 15),
#                   names = ['ID', 'dist', 'M_V', 'SpType', 'B-V'])
#df_clean = df.applymap(lambda x: np.nan if isinstance(x, basestring)
#                       and x.isspace() else x)
#df_clean= df_clean.dropna()
#
#
##create new row with first two characters of spectral class
#f = lambda s: (len(s) >= 2)  and (s[0].isalpha()) and (s[1].isdigit())
#i  = df_clean['SpType'].apply(f)
#df_clean = df_clean[i]
#f = lambda s: s[0:2]
#df_clean['SpType2'] = df_clean['SpType'].apply(f)
#
## Check spectral classes f = lambda s: s[0] #clases = df_clean['SpType'].map(f) #clases.value_counts()
#
##remove special classes C,N,R,S
#f = lambda s: s[0] in 'OBAFGKM'
#df_clean = df_clean[df_clean['SpType'].map(f)]
#
## order presicely
#orden = {'O':'0', 'B':'1', 'A':'2', 'F':'3', 'G':'4', 'K':'5', 'M':'6'}
#f = lambda s: orden[s[0]]+s[1]
#df_clean['SpType2'] = df_clean['SpType2'].apply(f)
#df_clean.head()
#
#def plot_lum_class(b,c, label):
#    ''' b: boolean Series to make the selection
#        c: Color
#        label: for the legend
#    '''
#    x = df_clean['B-V'][b]
#    y = df_clean['M_V'][b]
#    ax.scatter(x, y, c = c, s=5, edgecolors='none', label = label)
#
#fig = plt.figure(figsize=(8,10))
#ax = fig.add_subplot(111, facecolor='1.00')
#
#ax.set_xlim(-0.5, 2.75)
#ax.set_ylim(15, -15)
#ax.grid()
#ax.set_title('H-R Diagram /n HYD Star Database')
#ax.title.set_fontsize(15)
#ax.set_xlabel('Color index B-V')
#ax.xaxis.label.set_fontsize(15)
#ax.set_ylabel('Absolute magnitude')
#ax.yaxis.label.set_fontsize(15)
#
##code in luminosity class
#f = lambda s: 'VII' in s
#b = df_clean['SpType'].map(f)
#plot_lum_class(b,'black', 'VII: white dwarfs')
#
#f = lambda s: ('VI' in s) and ('VII' not in s)
#b = df_clean['SpType'].map(f)
#plot_lum_class(b,'darkblue', 'VI: subdwarfs')
#
#f = lambda s: ('V' in s) and ('VI' not in s) and ('IV' not in s)
#b = df_clean['SpType'].map(f)
#plot_lum_class(b,'teal', 'V: main-sequence')
#
#f = lambda s: 'IV' in s
#b = df_clean['SpType'].map(f)
#plot_lum_class(b,'hotpink', 'IV: subgiants')
#
#f = lambda s: 'III' in s
#b = df_clean['SpType'].map(f)
#plot_lum_class(b,'green', 'III: giants')
#
#f = lambda s: ('II' in s) and ('III' not in s) and ('VII' not in s)
#b = df_clean['SpType'].map(f)
#plot_lum_class(b,'orange', 'II: bright giants')
#
#f = lambda s: ('I' in s) and ('II' not in s) and ('V' not in s)
#b = df_clean['SpType'].map(f)
#plot_lum_class(b,'yellow', 'I: supergiants')
#
#ax.tick_params(axis='both', labelsize=10)
#ax.legend(loc = 'best', scatterpoints=1,markerscale = 4, shadow=True)
