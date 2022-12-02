import numpy as np
import urllib.request, urllib.error, urllib.parse
import pickle
import os
import simbad


_TycHip = './siphon_tychip.dpy'

def clean():
    try:
        os.remove(_TycHip)
    except:
        pass
    return

def TycHipData(target, errSys=None):
    """
    http://cdsarc.u-strasbg.fr/viz-bin/nph-Plot/Vgraph/htm?I/239/110991&3995-1479-1
    """
    global _TycHip
    if os.path.exists(_TycHip):
        f = open(_TycHip)
        db = pickle.load(f)
        f.close()
    else:
        db = {}

    if target in list(db.keys()):
        res = db[target]
        if not errSys is None:
            for r in res:
                r[-1] = np.sqrt(r[-1]**2+errSys**2)
        return res

    s = simbad.query(target)[0]
    url = 'http://cdsarc.u-strasbg.fr/viz-bin/nph-Plot/Vgraph/htm?I/239/'
    url += str(s['HIP'])+'&'+s['TYC']

    lines = urllib.request.urlopen(url, timeout=20).readlines()
    tyc, hip = [], []
    store = 0
    for l in lines:
        if not '|' in l:
            store = 0
        if store == 1:
            tyc.append(l)
        if store == 2:
            hip.append(l)
        if 'JD-2440000| BTmag' in l:
            store = 1
        if 'JD-2440000| Hpmag' in l:
            store = 2

    res = []
    for d in tyc:
        try:
            if len(d.split('|')[-1].strip())<2:
                res.append([float(d.split('|')[0])+40000-0.5,
                       'mag; Tycho', 'B_MVB_TYCHO',
                       float(d.split('|')[1]),
                       float(d.split('|')[2])])

                res.append([float(d.split('|')[0])+40000-0.5,
                       'mag; Tycho', 'V_MVB_TYCHO',
                       float(d.split('|')[4]),
                       float(d.split('|')[5])])
        except:
            pass
    for d in hip:
        try:
            if len(d.split('|')[-1].strip())<2:
                res.append([float(d.split('|')[0])+40000-0.5,
                       'mag; Hipparcos', 'HP_MVB_HIPPARCOS',
                       float(d.split('|')[1]),
                       float(d.split('|')[2])])
        except:
            pass
    db[target] = res
    f = open(_TycHip, 'wb')
    pickle.dump(db, f)
    f.close()

    if not errSys is None:
        for r in res:
            r[-1] = np.sqrt(r[-1]**2+errSys**2)

    return res

def McMaster():
    """
    list of targets:
    http://crocus.physics.mcmaster.ca/Cepheid/Classical.html
    """
    # -- list stars:
    global _mcMaster
    _mcMaster={}
    url = 'http://crocus.physics.mcmaster.ca/Cepheid/Classical.html'
    lines = urllib.request.urlopen(url, timeout=20).readlines()
    lines = [l for l in lines if '<a href=' in l]
    for l in lines:
        print(l)
        name = l.split('<tt>')[1].split('</tt>')[0].strip()
        name = name.split('.')
        test = True
        while test:
            try:
                name.remove('')
            except:
                test = False
        name = ' '.join(name)
        _mcMaster[name] = {'url':l.split('"')[1].strip()}
    return _mcMaster

def _mcmPage(url):
    lines = urllib.request.urlopen('http://crocus.physics.mcmaster.ca/'+url, timeout=20).readlines()
    source = ''
    for i in range(len(lines)):
        if '[data]' in lines[i]:
            if '[abstract]' in lines[i]:
                data = lines[i].split('[abstract]')[1]
            else:
                data = lines[i]
            data = data.split('[data]')[0].split('<a href="')[1].split('"')[0]
            print('    ', data)

        if '[notes]' in lines[i]:
            if '[data]' in lines[i]:
                notes = lines[i].split('[data]')[1]
            else:
                notes = lines[i]
            notes = notes.split('<a href="')[1].split('"')[0]
            print('    ', notes)
            # -- once whe have the notes, we have everything to star reading!
            _mcmData(notes, data)

        if '<li>' in lines[i]:
            source = lines[i].split('<li>')[1]
            while not 'br>' in lines[i]:
                i += 1
                source += lines[i]
            if '<br>' in source:
                source = source.split('<br>')[0]
            if '</br>' in source:
                source = source.split('</br>')[0]
            source = source.replace('\n', ' ')
            source.replace(';', '')
            print('> ', source)
    return

def _mcmData(notes, data):
    n = urllib.request.urlopen('http://crocus.physics.mcmaster.ca/'+notes, timeout=20).readlines()
    # -- figure out the columns:
    cols = [l for l in n if ')' in l and l.split(')')[0].split()[-1].isdigit()]
    #print cols
    ncol = [int(c.split(')')[0].split()[-1]) for c in cols]
    cols = [')'.join(c.split(')')[1:]) for c in cols]
    cols = [c.replace('\n', '').strip() for c in cols]
    #print zip(ncol, cols)
    #print ncol
    # -- detect if 2 sets of data are present in the file
    if len([n for n in ncol if n==1])>1:
        print('More than one data set!')

    n = urllib.request.urlopen('http://crocus.physics.mcmaster.ca/'+notes, timeout=20).readlines()
    # -- case with only one data set
    # -- find columns
    tmp = [l.replace(' ','').replace('.','') for l in n]

    return
