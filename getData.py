import numpy as np
import urllib.request, urllib.error, urllib.parse
import pickle
import os

"""

"""

# -- database is a pickeld dict indexed by HIP number (integer) pointing at data in SPIPS format
_dBfile_hip = './DATA/PHOTOM/hip.dpy'

def readDb_hip(HIP):
    global _dBfile_hip
    if not os.path.exists(_dBfile_hip):
        return None
    f = open(_dBfile_hip, 'rb')
    db = pickle.load(f)
    f.close()
    if HIP in list(db.keys()):
        return db[HIP]
    else:
        print(list(db.keys()))
        return None

def writeDb_hip(HIP, data):
    global _dBfile_hip
    if not os.path.exists(_dBfile_hip):
        print(' > getData: creating', _dBfile_hip)
        db = {HIP:data}
    else:
        print(' > getData: updating', _dBfile_hip)
        f = open(_dBfile_hip, 'rb')
        db = pickle.load(f)
        f.close()
    db[HIP] = data
    f = open(_dBfile_hip, 'wb')
    pickle.dump(db, f)
    f.close()
    return

def getHipPhotom(HIP, errSys=0.0):
    """
    HIP == HIPPARCOS NUMBER
    """
    global _dBfile_hip
    # -- get data from saved file
    tmp = readDb_hip(HIP)
    if not tmp is None:
        print(' > getData: HIP%d already in database %s'%(HIP, _dBfile_hip))
        for t in tmp:
            t[-1] = np.sqrt(t[-1]**2+errSys**2)
        return tmp
    print(' > getData: fetching data for HIP%d'%HIP)
    # -- query ESA page
    url = 'http://www.rssd.esa.int/hipparcos_scripts/HIPcatalogueSearch.pl?hipepId=%d'%HIP
    lines = urllib.request.urlopen(url).read().decode('utf-8').split('\n')
    # -- keep only data
    lines = [l for l in lines if '|' in l and not 'HT' in l]
    # -- keep only good points:
    lines = [l for l in lines if len(l.split('|'))>=4]
    lines = [l for l in lines if '0' in l.split('|')[3]]
    # -- see http://www.rssd.esa.int/SA-general/Projects/Hipparcos/CATALOGUE_VOL1/sect2_05.pdf
    data = [[float(l.split('|')[0])+2440000-2400000.5,
            'mag; Hipparcos', 'HP_B2000',
            float(l.split('|')[1]),
            float(l.split('|')[2])] for l in lines ]
    writeDb_hip(HIP, data)
    return data
