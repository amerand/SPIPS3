import spips
import numpy as np
import matplotlib.pyplot as plt
import os

def computeUD(filename, wavelength=0.8, B=100.0, plot=False):
    """
    compute UD for different phases (0.00-0.99 in 0.01 increment) at a given
    wavelength, baseline for a SPIPS model from a FITS file

    filename is the exported FITS file
    wavelength: in microns
    B: baseline in meters

    return  phases, UD

    """
    global _p
    # -- load parameters and data from FITS file
    _p, _obs = spips.importFits(filename)
    # -- create fake obs:
    phases = np.linspace(0,1,101)[:-1]
    obs = [(p, 'UDdiam', (wavelength, B), 0., 0.) for p in phases]
    diams = np.array(spips.model(obs, _p, plot=False, verbose=0))
    _p = spips.dephaseParam(_p, spips.phaseOffset)
    diams = np.array(spips.model(obs, _p, plot=False, verbose=0))
    if plot:
        plt.clf()
        plt.plot(phases, diams, 'ok')
        plt.title(filename+', B=%.0fm, wl=%.2fum'%(B, wavelength))
    return phases, diams

stars = ['T Vul', 'zeta Gem', 'X Cyg', 'T Mon', 'SV Vul']
data = {}
wavelength, B = 0.8, 100.
for s in stars:
    d = s.replace(' ', '_').upper()
    f = s.replace(' ', '_').lower()
    filename = '../SPIPS_STARS/FINAL/%s/%s.fits'%(d,f)
    data[s] = computeUD('../SPIPS_STARS/FINAL/%s/%s.fits'%(d,f),
                        wavelength=wavelength, B=B)

print('# WL=%.2fum B=%.0fm'%(wavelength, B))
print('# Phase ,', end=' ')
for s in stars:
    print(s ,',', end=' ')
print('')
for i in range(len(data[stars[0]][0])):
    for j,s in enumerate(stars):
        if j==0:
            print('%.3f'%data[s][0][i], ',', end=' ')
        print('%.3f'%data[s][1][i], end=' ')
        if j!=len(stars)-1:
            print(',', end=' ')
    print('')
