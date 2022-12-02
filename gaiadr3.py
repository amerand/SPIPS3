from astroquery.gaia import Gaia
# -- update astroquery to query DR2! e.g.:
# /Applications/anaconda2/lib/python2.7/site-packages/astroquery/gaia/core.py

from astroquery.simbad import Simbad
import astropy.units as u
from astropy.coordinates import SkyCoord
import urllib.request, urllib.error, urllib.parse
import numpy as np
import matplotlib.pyplot as plt
import pickle # Python 2!
import os

datafile = 'gaiadr2.cpickle'

def getDr2Phot(starName):
    """
    - query simbad to get sky coordinates for "starName"
    - query Gaia DR2 to get ID
    - query ESA to get epoch photometry

    returns dictionnary with all measurements, key 'SPIPS' contains the SPIPS data

    """
    if os.path.exists(datafile):
        f = open(datafile, 'rb')
        database = pickle.load(f)
        f.close()
        if starName in list(database.keys()):
            return database[starName]
    else:
        database = {}


    s = Simbad.query_object(starName)
    s.pprint()
    coord = SkyCoord(ra=s['RA'][0], dec=s['DEC'][0], unit=(u.hour, u.degree), frame='icrs')
    width = u.Quantity(1, u.arcsec)
    height = u.Quantity(1, u.arcsec)
    r = Gaia.query_object(coordinate=coord, width=width, height=height)
    data = {'SPIPS':[], 'MJD':[], 'RA':s['RA'][0], 'DEC':s['DEC'][0]}
    try:
        data['parallax'] = float(r['parallax'])
        data['parallax_error'] = float(r['parallax_error'])
    except:
        data['parallax'] = np.nan
        data['parallax_error'] = np.nan
        # -- save result
        database[starName] = data
        f = open(datafile, 'wb')
        pickle.dump(database, f)
        f.close()
        return data

    r['source_id', 'ra', 'dec', 'phot_g_mean_mag', 'parallax', 'parallax_error'].pprint()
    id = r['source_id']
    url = 'http://geadata.esac.esa.int/data-server/data?RETRIEVAL_TYPE=epoch_photometry&ID=%d&VALID_DATA=false&FORMAT=CSV'
    f = urllib.request.urlopen(url%(int(id)))
    for l in f.readlines():
        if not 'mag' in list(data.keys()):
            cols =  l.split(',')
            data.update({c:[] for c in cols})
        elif len(l)>10 and l.split(',')[3]!='':
            for i,v in enumerate(l.split(',')):
                try:
                    v = float(v)
                except:
                    pass
                data[cols[i]].append(v)
            data['MJD'].append(data['time'][-1]+55197.0)

            # -- make SPIPS data points:
            f = {'G':'G_GAIA_GAIA2',
                'BP':'Gbp_GAIA_GAIA2',
                'RP':'Grp_GAIA_GAIA2'}

            data['SPIPS'].append([data['MJD'][-1], 'mag;Gaia DR2',
                                  f[data['band'][-1]], data['mag'][-1],
                                  2.5/np.log(10)*1/data['flux_over_error'][-1],])
        else:
            pass
    for c in cols+['MJD']:
        data[c] = np.array(data[c])
    # -- save result
    database[starName] = data
    f = open(datafile, 'wb')
    pickle.dump(database, f)
    f.close()
    return data

def checkColorExcess(s, oplot=False):
    """
    s = getDr2Phot('...')
    """
    if not oplot:
        plt.figure(0, figsize=(12,5))
        plt.clf()
    if isinstance(s, str):
        if not oplot:
            plt.suptitle(s)
        s = getDr2Phot(s)
    if not 'band' in list(s.keys()):
        return
    F, M = {}, {}
    if not 'G' in s['band'] or not 'RP' in s['band'] or not 'BP' in s['band']:
        return
    for b in set(s['band']):
        w = np.where(s['band']==b)
        F[b] = (s['MJD'][w], s['flux'][w])
        M[b] = (s['MJD'][w], s['mag'][w])

    #     plt.plot(s['MJD'][w], s['mag'][w], 'o-', label=b)
    # plt.legend()
    RP = np.interp(M['G'][0], M['RP'][0], M['RP'][1])
    IRP = np.interp(F['G'][0], F['RP'][0], F['RP'][1])
    BP = np.interp(M['G'][0], M['BP'][0], M['BP'][1])
    IBP = np.interp(F['G'][0], F['BP'][0], F['BP'][1])
    IG = F['G'][1]
    color = BP-RP
    excess = (IRP+IBP)/IG
    plt.plot(color, excess, '.k', alpha=0.2)
    if not oplot:
        plt.xlabel('Bp-Rp')
        plt.ylabel('$(I_{BP}+I_{Rp})/I_{G}$')
        c = np.linspace(-1, 2, 100)
        plt.plot(c, 1.3+0.06*c**2, '--r', label='upper limit well behaved source')
        plt.legend()
    return

stars = ['alf UMi', 'v* L Car', 'eta Aql', 'bet Dor', 'Y Oph', 'zet Gem',
        'del Cep', 'X Sgr', 'W Sgr', 'RS Pup', 'T Mon', 'X Cyg', 'FF Aql',
        'Y Sgr', 'U Aql', 'SV Vul', 'V0636 Cas', 'U Car', 'AX Cir', 'U Vul',
         'S Vul', 'U Sgr', 'S Sge', 'RZ Vel', 'S Mus', 'BG Cru', 'TT Aql',
         'V0440 Per', 'SU Cas', 'MY Pup', 'AW Per', 'RX Cam', 'V0636 Sco',
         'T Vul', 'AH Vel', 'V1334 Cyg', 'EW Sct', 'GY Sge', 'TX Cyg',
         'V0470 Sco', 'SU Cru', 'RY Sco', 'V1496 Aql', 'U Nor', 'RT Aur',
         'KQ Sco', 'S Nor', 'SZ Tau', 'RW Cam', 'RU Sct', 'RY Vel', 'KN Cen',
         'R Mus', 'BB Sgr', 'WZ Sgr', 'V Cen', 'AQ Pup', 'S Cru', 'V0496 Aql',
         'DT Cyg', 'SW Vel', 'V0473 Lyr', 'FM Aql', 'AV Sgr', 'CK Cam', 'BM Per',
         'V0659 Cen', 'YZ Sgr', 'X Pup', 'SZ Aql', 'VY Car', 'VX Cyg', 'AV Cir',
         'V0350 Sgr', 'R TrA', 'FN Aql', 'RX Aur', 'V Car', 'BF Oph', 'V0340 Nor',
         'CD Cyg', 'IR Cep', 'BP Cir', 'XX Cen', 'XX Sgr', 'DL Cas', 'YZ Car',
         'T Vel', 'WZ Car', 'AS Per', 'UZ Sct', 'SU Cyg', 'CO Aur', 'Z Sct',
         'V0340 Ara', 'VZ Pup', 'Y Car', 'VZ Cyg', 'Y Lac', 'GH Lup', 'X Vul',
         'W Gem', 'VY Sgr', 'V0386 Cyg', 'Y Sct', 'BZ Cyg', 'SZ Cyg', 'CR Cep',
         'TW Nor', 'Z Lac', 'MW Cyg', 'V0600 Aql', 'V0495 Cyg', 'RS Cas', 'TY Sct',
         'V1162 Aql', 'RS Ori', 'CP Cep', 'SS Sct', 'CR Ser', 'RY CMa', 'CK Sct',
         'X Lac', 'SV Mon', 'SV Per', 'SX Vel', 'VY Cyg', 'ST Tau', 'VX Per',
         'RY Cas', 'CV Mon', 'RW Cas', 'VW Cen', 'YZ Aur', 'V0459 Cyg', 'BG Lac',
         'SW Cas', 'RR Lac', 'QZ Nor', 'FM Cas', 'V0538 Cyg', 'DD Cas', 'RZ Gem',
         'RZ CMa', 'UU Mus', 'BN Pup', 'V0402 Cyg', 'SY Cas', 'LS Pup', 'BE Mon',
         'TZ Mon', 'CS Vel', 'CF Cas', 'V1344 Aql', 'V0737 Cen', 'T Cru', 'QY Cen',
         'BG Vel', 'V0609 Cyg', 'RV Sco', 'S TrA', 'ER Car', 'AP Sgr', 'V0500 Sco',
         'DG Vul', 'XZ Car', 'R Cru', 'SY Nor', 'DR Vel', 'V0950 Sco', 'V0482 Sco',
         'AP Pup', 'BQ Ser', 'OR Cam', 'AY Sgr', 'SZ Cas', 'V0339 Cen', 'TX Cen',
         'VW Cru', 'IQ Nor', 'V0378 Cen', 'V5567 Sgr', 'GS Lup', 'SV Vel', 'CH Cas',
         'V0367 Sct', 'V0381 Cen', 'CY Cas', 'OX Cam', 'XY Car', 'V Vel', 'IT Car',
         'V0383 Cyg', 'AT Pup', 'V0379 Cas', 'X Cru', 'LR TrA', 'V0532 Cyg',
         'V0336 Aql', 'V0411 Lac', 'ST Vel', 'RS Nor', 'HO Vul', 'GH Cyg', 'AY Cen',
         'V0458 Sct', 'VY Per', 'V0391 Nor', 'SY Aur', 'CD Cas', 'GX Car', 'SS CMa',
         'V0419 Cen', 'AC Mon', 'AE Vel', 'EU Tau', 'WX Pup', 'V0397 Nor', 'V0496 Cen',
         'UW Car', 'HW Car', 'BP Cas', 'AG Cru', 'V0389 Sct', 'V Lac', 'V0901 Cep',
         'XX Car', 'VZ CMa', 'CF Cam', 'X Sct', 'U TrA', 'V0397 Car', 'GI Car',
         'FR Car', 'V1210 Cen', 'V0898 Cen', 'AX Vel', 'UY Per', 'GQ Ori', 'BY Cas',
         'TU Cas', 'GH Car', 'AQ Car', 'BR Vul', 'SX Car', 'GU Nor', 'EV Sct',
         'V2475 Cyg', 'RT Mus', 'V1726 Cyg', 'BK Aur', 'V0335 Pup', 'BD Cas',
         'AD Pup', 'FN Vel', 'UX Car', 'DW Cas', 'UY Car', 'EX Vel', 'VX Pup',
         'XY Cas', 'FO Car', 'AP Vel', 'V1154 Cyg', 'FI Car', 'V0351 Cep',
         'AA Gem', 'TYC 4034-222-1', 'V0824 Cas', 'V1019 Cas', 'UZ Cen',
         'CM Sct', 'MZ Cen', 'TW CMa', 'UZ Car', 'AZ Cen', 'V0701 Car', 'MN Cam',
         'V0520 Cyg', 'TV CMa', 'V0526 Mon', 'VX Cru', 'CE Pup', 'BB Cen', 'Y Aur',
         'WW Car', 'XX Vel', 'CY Car', 'GZ Car', 'VV Cas', 'BB Her', 'BM Pup',
         'VW Cas', 'EK Mon', 'GI Cyg', 'V0493 Aql', 'CZ Cas', 'CN Car', 'CS Mon',
         'OO Pup', 'MZ Cam', 'CR Car', 'EY Car', 'HK Car', 'BK Cen', 'V0395 Cas',
         'DF Cas', 'T Ant', 'V1397 Cyg', 'TX Mon', 'CG Cas', 'KK Cen', 'V0733 Aql',
         'V0637 Aur', 'TYC 8308-2055-1', 'NO Cas', 'IO Car', 'HL Pup', 'SX Per',
         'WW Pup', 'AY Cas', 'KL Aql', 'TZ Mus', 'V0508 Mon', 'DX Gem', 'LL Pup',
         'AD Gem', 'FN Car', 'VW Pup', 'V1100 Cas', 'UY Mon', 'EK Pup', 'DW Per',
         'WY Pup', 'MS Mus', 'V0363 Cas', 'WZ Pup', 'V0914 Mon', 'BV Mon', 'V1154 Cas',
         'DK Vel', 'XX Mon', 'V0465 Mon', 'UZ Cas', 'DY Car', 'V0371 Gem', 'V0924 Cyg',
         'TV Cam', 'OP Pup', 'UX Per', 'V1048 Cen', 'CS Ori', 'V1345 Cen', 'LR Pup',
         'BB Gem', 'V0572 Aql', 'V0371 Per', 'V0901 Mon', 'BQ Pup', 'IX Cas', 'DQ And',
         'V0720 Car', 'CE Cas A', 'CE Cas B',]
def getAll():
    global stars
    for i,s in enumerate(stars):
        if i%10==0:
            print('='*100)
            print(i, '/', len(stars))
            print('='*100)
        tmp = getDr2Phot(s)
    return

def plotAll():
    if os.path.exists(datafile):
        f = open(datafile, 'rb')
        database = pickle.load(f)
        f.close()
    else:
        return
    star = list(database.keys())
    p = np.array([database[s]['parallax'] for s in star])
    ep = np.array([database[s]['parallax_error'] for s in star])
    plt.figure(1)
    ax = plt.subplot(111)
    plt.clf()
    #plt.plot(p, 100*ep/p, '.k', alpha=0.5)
    emax = 10
    for i in range(len(star)):
        if not np.isnan(p[i]) and p[i]>0.1 and 100*ep[i]/p[i]<emax:
            plt.text(p[i], 100*ep[i]/p[i], star[i], size=5,
                    ha='center', va='center')

    plt.xlabel('parallax (mas)')
    plt.ylabel('parallax error (%)')
    #plt.xlim(0.1, 5); plt.xscale('log')
    #plt.ylim(1, 20); plt.yscale('log')
    plt.xlim(0,5); plt.ylim(0,emax); plt.grid()
    #Y = [1,2,5,10,20,50]
    #ax.set_yticks(Y)
    #ax.set_yticklabels([str(y)+'%' for y in Y])

    return
