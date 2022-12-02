#!/usr/bin/env python2.7

"""
simbad.py: query simbad astronomical database by name

example:

import simbad

# get the info from simbad:
s = simbad.query(['Polaris', 'delta Cep', 'eta Aql'])

# pretty print for org mode
simbad.prettyPrint(s)
"""

__author__ = 'Antoine Merand'
__email__ = "amerand@eso.org"
__version__ = "1.0"
__date__ = 'Mon  6 May 2013 05:09:07 UTC'

import urllib.request, urllib.error, urllib.parse
import numpy as np
import pickle
import os

#_default_site = 'http://simbad.u-strasbg.fr/' # default site is Strasbourg
_default_site = 'http://130.79.128.4/'
#_alternate_site = 'http://simbad.cfa.harvard.edu/' # alternate is Harvard
_alternate_site = 'http://131.142.185.22/'
_simbad_site = _default_site
_dbfile = 'simbad.dpy'

def query(identifiers, debug=False, closest=False, around=0):
    """
    identifiers is list of SIMBAD identifier

    if around is set, search around target for the specified value (in minutes)

    returns a list of dictionnaries with parameters (self
    explanatory). DIAM is the angular diameter in mas, estimated from
    photometry. TRANSIT MONTH is the moonth for which the objects
    transit at midnight.

    Note:
    - if magnitudes are not present, -99. is returned
    - if VSINI (km/s) not found, set to -1
    - if PM(A,D) or PLX not found, set to 0
    - IRAS FLUXES are in Jy, for 12, 25, 60, 100
    """
    global _simbad_site, _dbfile

    if not isinstance(identifiers, list): # make a scalar into a list
        identifiers = [identifiers]
    if closest and len(identifiers)>1:
        print('ID:', identifiers, len(identifiers))
        print('closest=True only for single query...')
        return

    ngroup = 40 # max number of objects to query at the same time
    if len(identifiers)>ngroup: # slice list
        # group by ngroup
        res = []
        while len(identifiers)>0:
            print(len(identifiers), len(res))
            tmp = query(identifiers[:ngroup], debug=debug, closest=closest, around=around)
            res.extend(tmp)
            identifiers = identifiers[ngroup:]
        return res

    ###########################################################
    ### from here, assumes identifiers is a list of strings ###
    ###########################################################

    # check if it is in the DataBase
    if os.path.isfile(_dbfile) and around==0:
        dbf = open(_dbfile)
        db = pickle.load(dbf)
        dbf.close()
        res = []
        for i in identifiers:
            if i in db:
                #print i, 'found in local DB'
                res.append(db[i])
            else:
                #print i, 'NOT found in local DB'
                res.append({})
    else:
        res = [{} for i in identifiers]

    # -- all target found in database
    if all(['IDENTIFIER' in r for r in res]):
        return res

    rt_ = '%0D%0A' # cariage return
    plus_ = '%2B' # + in the URL
    separator = ';'
    format_ = "format+object+form1+\""+separator+"+%25IDLIST(1,HD,HIP,TYC)+"+separator+\
            "+%25COO(A+D)+"+separator+"+%25OTYPE+"+separator+"+%25SP+"+\
            separator+"+%25PM(A+D)+"+separator+"+%25PLX(V+E)+"+separator+\
            "+%25FLUXLIST(B)+"+separator+"+%25FLUXLIST(V)+"+separator+\
            "+%25FLUXLIST(R)+"+separator+"+%25FLUXLIST(J)+"+separator+\
            "+%25FLUXLIST(H)+"+separator+"+%25FLUXLIST(K)+"+separator+\
            "+%25MEASLIST(rot;|F)+"+\
            separator+"%25MEASLIST(JP11;|F)"+"\""
            #separator+"%25MEASLIST(iras)"+\


    url = 'simbad/sim-script?submit=submit+script&script='+format_

    Nquery = 0
    IDquery = []
    for k,i in enumerate(identifiers):
        if 'IDENTIFIER' not in res[k]:
            Nquery+=1
            IDquery.append(i)
            obj = i.replace('+', plus_)
            obj = obj.replace('_', ' ')
            obj = obj.replace(' ', '+')
            if ':' in i: # these must be coordinates!
                url = url+rt_+'query+coo+'+obj+'+radius%3D5s'
            elif around>0:
                url = url+rt_+'query+around+'+obj+'+radius%3D'+str(around)+'m'
            else:
                url = url+rt_+'query+id+'+obj

    if debug:
        print(_simbad_site+url)
    try:
        lines = urllib.request.urlopen(_simbad_site+url, timeout=20).read()
    except:
        _simbad_site = _alternate_site
        print('switching to alternate server...')
        try:
            lines = urllib.request.urlopen(_simbad_site+url, timeout=20).read()
        except:
            raise NameError('servers do not respond OR no internet connection')

    if debug:
        print(lines)
    lines = lines.split('\n')

    # go to data
    for k, l in enumerate(lines):
        if ':error:' in l:
            #print '  ERROR:', lines[k+2]
            #print '------------------------------'
            #print lines
            # -- save result to avoid querying it again
            if len(IDquery)==1:
                try:
                    if not isinstance(db, dict):
                        db = {}
                except:
                    db = {}

                for k,i in enumerate(IDquery):
                    db[i]= {}
                dbf = open(_dbfile, 'w')
                pickle.dump(db, dbf)
                dbf.close()
            return None
        if ':data:' in l:
            lines = lines[k+1:]
            break

    lines = [x for x in lines if len(x)>0]

    if len(lines)!=Nquery and not closest and around==0:
        print('  ERROR: too many/few results!')
        return None
    if debug:
        print(lines)

    # read every line which is a different object
    for k, l in enumerate(lines):
        obj = {}
        if around>0:
            obj['IDENTIFIER'] = 'around: '+identifiers[0]
        else:
            obj['IDENTIFIER'] = IDquery[k]

        obj['NAME'] = l.split(separator)[1].split(',')[0].strip()
        if 'HD' in l.split(separator)[1]:
            tmp = l.split(separator)[1].split('HD')[1].split(',')[0]
            if not tmp.isdigit():
                tmp = tmp[:-1]
            if tmp.isdigit():
                obj['HD'] = int(tmp)
            else:
                obj['HD'] = -1
        else:
            obj['HD'] = -1
        if 'HIP' in l.split(separator)[1]:
            obj['HIP'] = int(l.split(separator)[1].split('HIP')[1].split(',')[0])
        else:
            obj['HIP'] = -1
        if 'TYC' in l.split(separator)[1]:
            obj['TYC'] = l.split(separator)[1].split('TYC')[1].split(',')[0].strip()
        else:
            obj['TYC'] = ''

        if '-' in l.split(separator)[2]:
            l_ra = l.split(separator)[2].split('-')[0]
            l_dec = '-'+l.split(separator)[2].split('-')[1]
        else:
            l_ra = l.split(separator)[2].split('+')[0]
            l_dec = '+'+l.split(separator)[2].split('+')[1]

        obj['RA'] = l_ra.strip()
        obj['DEC'] = l_dec.strip()

        if len(l_ra.split())==3:
            obj['RA.h'] = (float(l_ra.split()[0])+
                           float(l_ra.split()[1])/60.+
                           float(l_ra.split()[2])/3600.)
        elif len(l_ra.split())==2:
            obj['RA.h'] = (float(l_ra.split()[0])+
                           float(l_ra.split()[1])/60.)
        else:
            obj['RA.h'] = float(l_ra.split()[0])

        obj['RA.d'] = obj['RA.h']*15

        if len(l_dec.split())==3:
            obj['DEC.d'] = abs(float(l_dec.split()[0]))+\
                           float(l_dec.split()[1])/60.+\
                           float(l_dec.split()[2])/3600.
        elif len(l_dec.split())==2:
            obj['DEC.d'] = abs(float(l_dec.split()[0]))+\
                           float(l_dec.split()[1])/60.
        else:
            obj['DEC.d'] = abs(float(l_dec.split()[0]))

        obj['DEC.d'] = np.copysign(obj['DEC.d'],
                                     float(l_dec.split()[0]))

        # 15th Jan at midnight is ~ LST 6:00
        obj['TRANSIT MONTH'] = int(round((obj['RA.h']-6.00)/2.-1, 0))%12+1
        obj['TYPE'] = l.split(separator)[3].split('~')[0].strip()
        obj['SPTYPE'] = l.split(separator)[4].strip().split()[0]

        try:
            obj['PMA'] = float(l.split(separator)[5].split()[0])/1000.
            obj['PMD'] = float(l.split(separator)[5].split()[1])/1000.
        except:
            obj['PMA'] = 0.0
            obj['PMD'] = 0.0

        try:
            obj['PLX'] = float(l.split(separator)[6].split()[0])/1000.
            obj['EPLX'] = float(l.split(separator)[6].split()[1])/1000.
        except:
            obj['PLX'] = 0.0
            obj['EPLX'] = 0.0

        mags = ['B','V','R','J','H','K']
        for j, m in enumerate(mags):
            try:
                obj[m+'MAG'] = float(l.split(separator)[7+j].split()[1])
            except:
                try:
                    # take first number
                    tmp = l.split(separator)[7+j]
                    for i in range(len(tmp)):
                        if tmp[i].isdigit() or tmp[i]=='-':
                            break
                    obj[m+'MAG'] = float(tmp[i:].split()[0])
                except:
                    obj[m+'MAG'] = np.nan
        try:
            obj['VSINI'] = float(l.split(separator)[13].split('|')[0].split()[0])
        except:
            obj['VSINI'] = -1 # failed
        iras_wl = ['12um', '24um', '60um', '100um']

        obj['IRAS'] = dict(list(zip(iras_wl, np.zeros(len(iras_wl)))))
        for i,j in enumerate(iras_wl):
            try:
                _tmp = l.split(separator)[14].split('|')[2].split()[i]
            except:
                _tmp = ''
            try:
                obj['IRAS'][j] = float(_tmp)
            except:
                obj['IRAS'][j] = _tmp

        JP11_wl = ['U', 'B', 'V', 'R', 'I', 'J', 'K', 'L', 'M', 'N', 'H']
        obj['JP11'] = dict(list(zip(JP11_wl, np.zeros(len(JP11_wl)))))
        for i,j in enumerate(JP11_wl):
            try:
                obj['JP11'][j] = float(l.split(separator)[15].split('|')[i].split()[0])
            except:
                obj['JP11'][j] = np.nan
        if np.isnan(obj['KMAG']) and not np.isnan(obj['JP11']['K']):
            obj['KMAG']= obj['JP11']['K']
        if np.isnan(obj['RMAG']) and not np.isnan(obj['JP11']['R']):
            obj['RMAG']= obj['JP11']['R']

        res[identifiers.index(IDquery[k])] = obj
        if closest:
            break

    if around>0:
        for k in range(len(res)):
            res[k]['DIST D'] = np.sqrt( (res[0]['DEC.d']-res[k]['DEC.d'])**2+
                                          np.cos(res[0]['DEC.d']*3.1415/180)**2*
                                          (res[0]['RA.d']-res[k]['RA.d'])**2)
            res[k]['DIST S'] =  res[k]['DIST D']*3600
    res = addApproxDiam(res, verbose=False)

    if around==0:
        try:
            if not isinstance(db, dict):
                db = {}
        except:
            db = {}

        for k,i in enumerate(IDquery):
            if not res[k] is None:
                #print 'adding ->', k, i
                db[i] = res[k]
        if db!={}:
            #print 'saving', _dbfile
            dbf = open(_dbfile, 'w')
            pickle.dump(db, dbf)
            dbf.close()
    return res

def prettyPrint(dics, keys=['IDENTIFIER', 'RA', 'DEC', 'SPTYPE', 'VMAG', 'HMAG', 'KMAG', 'DIAM']):
    """
    dics is the result of 'query' (a list of dict). keys is the list of keys
    of the dict you want to display in the table
    """
    # -- make each line with list of strings
    res = [keys]
    for d in dics: # for each object
        line = []
        for k in keys:
            line.append(str(d[k]))
        res.append(line)
    # -- get the longest string for each column
    n = [max([len(r[k]) for r in res]) for k in range(len(res[0]))]
    # -- make each line
    form = '| '+' | '.join([('%' if 'MAG' in keys[i] else '%-')+str(k)+'s' for i,k in enumerate(n)])+' |\n'
    # -- make the intermediate line
    _line = (form.replace('|', '+').replace(' ','-'))%tuple(['-'*k for k in n])
    # -- print all lines
    _res = _line
    _res += form%tuple(res[0])
    _res += _line
    for r in res[1:]:
        _res += form%tuple(r)
    _res += _line
    return _res

def addApproxDiam(dics, verbose=True):
    """
    add the approximated diameter estimated using V-K

    uses 'BMAG', 'VMAG', 'JMAG', 'HMAG', 'KMAG' keyword of the
    dictionnary
    """
    # surface brightness relations for dwarf stars
    # from Kervella et al. 2004
    k04 = {}
    #           coef0  coef1  error
    k04['BV']=[.9095, .4889, .0918]
    #k04['BR']=[.4771, .5116, .0475]
    k04['BJ']=[.3029, .5216, .0307]
    k04['BH']=[.2630, .5134, .0189]
    k04['BK']=[.2538, .5158, .0100]
    k04['VJ']=[.3547, .5310, .0475]
    #k04['VR']=[.7900, .5217, .0853]
    k04['VH']=[.2893, .5148, .0185]
    k04['VK']=[.2753, .5175, .0101]
    #k04['RJ']=[.4647, .5392, .0843]
    #k04['RH']=[.3405, .5119, .0288]
    #k04['RK']=[.3833, .5140, .0384]
    k04['JH']=[.6280, .4990, .1044]
    k04['JK']=[.5256, .5097, .0575]

    for k, d in enumerate(dics): # for each star
        diams = []
        errs = []
        couls = []
        for coul in list(k04.keys()): # for each color
            # -- check magnitudes are valid, compute diameter and error
            if coul[0]+'MAG' in d and d[coul[0]+'MAG']>-90 and\
            coul[1]+'MAG' in d and d[coul[1]+'MAG']>-90:
                diams.append(diamSurfBri(d[coul[0]+'MAG'], d[coul[1]+'MAG'],
                                         k04[coul]))
                errs.append(k04[coul][2]*diams[-1])
                couls.append(coul)
        if len(diams)>1:
            diams = np.array(diams)
            errs  = np.array(errs)
            _tmp = np.sum(diams/errs)/np.sum(1./errs)
            dics[k]['DIAM chi2'] = np.mean((diams-_tmp)**2/errs**2)
            dics[k]['DIAM all'] = [(couls[i], diams[i], errs[i]) for i in range(len(diams))]
            dics[k]['DIAM'] = round(_tmp, int(-np.log10(_tmp) + 3))
        elif len(diams)==1:
            dics[k]['DIAM'] = round(diams[0], int(-np.log10(diams[0])+3))
        else:
            dics[k]['DIAM'] = 0
        if verbose:
            print(dics[k]['NAME'], '|', dics[k]['DIAM'])
    return dics

def diamSurfBri(c0, c1, coef):
    """
    surface brightness formula for magnitudes c0 and c1 (color c1-c0)

    see Kervella & Fouque 2008
    """
    return 10**(coef[1] + coef[0]*(c0-c1) - 0.2*c0)

def createEdb(s):
    """
    convert simbad dic into a xEphem format (list of stringes)

    exple: Polaris,f|M|F7,2:31:48.704,89:15:50.72,2.02,2000
    """
    if isinstance(s, list):
        return [createEdb(x) for x in s]

    res = "%s,f|S|%s,%s,%s,%3.1f,2000" % (s['IDENTIFIER'],
                                          s['SPTYPE'][:4],
                                          s['RA'].strip().replace(" ", ":"),
                                          s['DEC'].strip().replace(" ", ":"),
                                          s['VMAG'])
    return res

def createRdb(s):
    """
    create RDB line for APES input catalog. Fields are tab separated, missing
    values are replaced by ~.

    stype:  'T' or 'R'
    system_id: same for 2 stars
    star_id
    alpha, delta in dd/hh:mm:ss.ss
    dalpha, ddelta: ?
    epoch: J2000
    equinox: 2000
    coord_syst: ircs
    mualpha, mudelta: propermotion, in mas/yr
    Tint_max: ?
    stdev_Dphi: ?
    parallax: in mas
    SP_type: spectral type
    Teff: in K
    lambda_eff: in microns?
    magV, magK, magH:
    MS_tar, r, MP_1, T0_1, period_1, ecc_1, a_1, inc_1, omega_1, OMEGA_1,
        MP_2, T0_2, period_2, ecc_2, a_2, inc_2, omega_2, OMEGA_2, MP_3, T0_3,
        period_3, ecc_3, a_3, inc_3, omega_3, OMEGA_3: parameters of astrometric
        signal

    """
    pass
    return
