import spips
import delta_cep_data # contains all the observations
import os
import numpy as np
import matplotlib.pyplot as plt
import time

def fitFromScratch(obs, firstGuess={'PERIOD':5.3663, 'd_kpc':.25, 'DIAMAVG':1.5, 'TEFFAVG':6000}, 
                    vradSplines=False, teffSplines=False,
                    changeOfPeriodOrder=1, withIrExcess=True, oT=5, oV=8):
    """
    find SPIPS Fourier model from dataset "obs"
    
    firstGuess: provide initial guess as dict. Mandatory is "PERIOD" (in days). Recommended are
        'd_kpc': distance in kpc
        'DIAMAVG': average angular diameter in mas
        'TEFFAVG': average effective temperature in K
    
    changeOfPeriodOrder: 0 for no change of period, 1 for linear etc (default 1)
    
    withIrExcess: True or False (default True)

    oT, oV: number of Fourier componenents  for Teff and Vrad (default 5, 8)
    """

    # -- data without Teff and color, just for test
    #obs = [o for o in obs if not any([o[1].startswith(x) for x in ['teff', 'color']])]

    p = {'P-FACTOR':1.27,
         'd_kpc': 1.0,
         'METAL': 0.0,
         'E(B-V)':0.0, 
         'DIAMAVG':1.0,
      }
    if 'TEFFAVG' in firstGuess:
        teff = firstGuess['TEFFAVG']
        firstGuess.pop('TEFFAVG')
    else:
        teff = 6000.

    p.update(firstGuess)
    execution = {}
    # == fit teff and colors
    print('fitting Teff profile to teff data...', end=' ')

    if vradSplines:
        p.update({'VRAD POW': 1.5,
                 'VRAD PHI0': 0.8,
                 'VRAD VAL0': 0.,})
    else:
        p['VRAD A0'] = 0

    if teffSplines:
        p.update({'TEFF VAL0': teff,
                  'TEFF POW': 1.2,         
                  'TEFF PHI0': .8})

        for k in range(oT-1):
            p['TEFF DVAL%d'%(k+1)] = np.interp(k, [0, oT/3, oT/2, 2*oT/3, oT], [0.25, 1, 0, -1, 0])
            p['TEFF DVAL%d'%(k+1)] *= 500
        follow = ['TEFF VAL0']
    else:    
        p['TEFF A0'] = teff
        for k in range(oT):
            if k==0:
                p['TEFF A%d'%(k+1)] = 500/(k+1)**2 
                p['TEFF PHI%d'%(k+1)] = np.pi/2*(k+1)
            else:
                p['TEFF R%d'%(k+1)] = 1/(k+1)**2 
                p['TEFF PHI%d'%(k+1)] = 0.01
        follow = ['TEFF A0']
    
    tech = ['teff']
    _obs = [o for o in obs if any([o[1].startswith(t) for t in tech])]
    maxCores = 1
    if len(_obs)<=(2*oT+1):
        print('not enough data')
        print('fitting Teff profile to teff and color data...', end=' ')
        # -- try including colors
        tech = ['teff', 'color']
        _obs = [o for o in obs if any([o[1].startswith(t) for t in tech])]
        maxCores = 4

    if len(_obs)>(2*oT+1):
        mjd0 = np.mean([o[0] for o in _obs if not o[0] is None ])
        p['MJD0'] = mjd0
        doNotFit = ['MJD0', 'P-FACTOR', 'METAL', 'd_kpc', 'DIAMAVG', 'PERIOD','PERIOD']
        doNotFit += list(filter(lambda x: x.startswith('VRAD '), p.keys()))
        if not (any([o[1].startswith('teff') for o in _obs]) and
                any([o[1].startswith('color') for o in _obs])):
            doNotFit.append('E(B-V)')
        t0 = time.time()

        fit = spips.fit(_obs, p, doNotFit=doNotFit, verbose=maxCores>1, plot=False, maxCores=maxCores,
                        normalizeErrors='techniques', ftol=1e-3,
                        follow=['E(B-V)']+follow)
        p = fit['best']
        print('chi2=%.3f [in %.1fs]'%(fit['chi2'], time.time()-t0))
        execution['Teff'] = time.time()-t0
        # -- 

    # == only Vrad =====================================
    print('fitting VRAD only...', end=' ')
    tech = ['vrad']
    _obs = [o for o in obs if any([o[1].startswith(t) for t in tech])]

    p['MJD0'] += np.round((np.mean([o[0] for o in _obs])-p['MJD0'])/p['PERIOD'], 0)*p['PERIOD']
    if vradSplines:
        p['VRAD VAL0'] = np.mean([o[-2] for o in _obs]) 
        p['VRAD PHI0'] = 0.8
        p['VRAD POW'] = 1.8
        for k in range(oV-1):
            p['VRAD DVAL%d'%(k+1)] = np.interp(k, [0, oV/3, oV/2, 2*oV/3, oV], [0.25, 1, 0.01, -1, -0.01])
            p['VRAD DVAL%d'%(k+1)] *= np.ptp([o[-2] for o in _obs])/2
        follow = []
    else:
        p['VRAD A0'] = np.mean([o[-2] for o in _obs])
        p['VRAD A1'] = np.ptp([o[-2] for o in _obs])/oV
        p['VRAD PHI1'] = 0
        for k in range(oV-1):
            p['VRAD R%d'%(k+2)] = 1.0
            p['VRAD PHI%d'%(k+2)] = 0 
        follow = []

    doNotFit = ['MJD0', 'P-FACTOR', 'METAL', 'd_kpc', 'DIAMAVG','E(B-V)']
    for k in p:
        if k.startswith('TEFF '):
            doNotFit.append(k)
    t0 = time.time()
    fit = spips.fit(_obs, p, doNotFit=doNotFit, verbose=False, plot=False, maxCores=1, ftol=1e-3, 
                    epsfcn=1e-8)
    p = fit['best']
    print('chi2=%.3f [in %.1fs]'%(fit['chi2'], time.time()-t0))
    execution['Vrad'] = time.time()-t0

    if True and 'Teff' in execution:
        # -- not really needed?

        print('fitting all data with fixed Vrad and Teff profiles', end=' ')
        doNotFit = ['MJD0', 'P-FACTOR', 'METAL']
        # -- keep Teff and Vrad profile fixed
        for k in p:
            if k.startswith('VRAD') and not k in ['VRAD PHI0', 'VRAD PHI1']:
                doNotFit.append(k)
            if k.startswith('TEFF') and not k in ['TEFF PHI0', 'TEFF PHI1']:
                doNotFit.append(k) 
        mjd0 = np.mean([o[0] for o in obs if not o[0] is None ])
        p['MJD0'] += np.round((mjd0-p['MJD0'])/p['PERIOD'], 0)*p['PERIOD']
        t0 = time.time()
        fit = spips.fit(obs, p, doNotFit=doNotFit, verbose=True, plot=False, maxCores=4,
                        normalizeErrors='techniques', ftol=1e-2,
                        follow=['E(B-V)', 'd_kpc'])
        p = fit['best']
        print('chi2=%.3f [in %.1fs]'%(fit['chi2'], time.time()-t0))
        execution['All, fixed Vrad & Teff'] = time.time()-t0

    # -- fitting everything
    print('fitting all parameters to all data ', end=' ')
    doNotFit = ['MJD0', 'P-FACTOR', 'METAL']
    if teffSplines:
        follow = ['TEFF VAL0']
    else:
        follow = ['TEFF A0']
    if withIrExcess:
        p.update({
                  'EXCESS EXP':   0.4 , # exponential law for IR Excess
                  'EXCESS SLOPE': 0.0 , # slope for IR excess
                  'EXCESS WL0':   1.2 , # starting WL, in um, for IR Excess
                  })
        doNotFit += ['EXCESS EXP', 'EXCESS WL0']
    if changeOfPeriodOrder>0:
        for k in range(changeOfPeriodOrder):
            if not 'PERIOD%d'%(k+1) in p:
                p['PERIOD%d'%(k+1)] = 1e-3    
    mjd0 = np.mean([o[0] for o in obs if not o[0] is None ])
    p['MJD0'] += np.round((mjd0-p['MJD0'])/p['PERIOD'], 0)*p['PERIOD']
    t0 = time.time()
    fit = spips.fit(obs, p, doNotFit=doNotFit, verbose=True, plot=False, maxCores=4,
                    normalizeErrors='techniques', ftol=1e-5, epsfcn=1e-8,
                    follow=['d_kpc', 'E(B-V)', 'DIAMAVG']+follow)
    p = fit['best']
    print('chi2=%.3f [in %.0fs]'%(fit['chi2'], time.time()-t0))
    execution['All'] = time.time()-t0
    print(execution)
    print('total execution %.1fmin'%(np.sum([execution[k] for k in execution])/60))
    return fit


def fitFromScratchFourier(obs, firstGuess={'PERIOD':5.3663, 'd_kpc':.25, 'DIAMAVG':1.5, 'TEFF A0':6000}, 
                            changeOfPeriodOrder=1, withIrExcess=True, oT=5, oV=8):
    """
    find SPIPS Fourier model from dataset "obs"
    
    firstGuess: provide initial guess as dict. Mandatory is "PERIOD" (in days). Recommended are
        'd_kpc': distance in kpc
        'DIAMAVG': average angular diameter in mas
        'TEFF A0': average effective temperature in K
    
    changeOfPeriodOrder: 0 for no change of period, 1 for linear etc (default 1)
    
    withIrExcess: True or False (default True)

    oT, oV: number of Fourier componenents  for Teff and Vrad (default 5, 8)
    """

    # -- data without Teff and color, just for test
    #obs = [o for o in obs if not any([o[1].startswith(x) for x in ['teff', 'color']])]

    p = {'P-FACTOR':1.27,
         'd_kpc': 1.0,
         'TEFF A0': 6000,
         'METAL': 0.0,
         'E(B-V)':0.0, 
         'DIAMAVG':1.0,
         'VRAD A0': 0.0,
      }
    p.update(firstGuess)
    execution = {}
    # == fit teff and colors
    print('fitting Teff profile to teff data...', end=' ')

    for k in range(oT):
        if k==0:
            p['TEFF A%d'%(k+1)] = 500/(k+1)**2 
            p['TEFF PHI%d'%(k+1)] = np.pi/2*(k+1)
        else:
            p['TEFF R%d'%(k+1)] = 1/(k+1)**2 
            p['TEFF PHI%d'%(k+1)] = 0.01
    
    tech = ['teff']
    _obs = [o for o in obs if any([o[1].startswith(t) for t in tech])]
    maxCores = 1
    if len(_obs)<=(2*oT+1):
        print('not enough data')
        print('fitting Teff profile to teff and color data...', end=' ')
        # -- try including colors
        tech = ['teff', 'color']
        _obs = [o for o in obs if any([o[1].startswith(t) for t in tech])]
        maxCores = 4

    if len(_obs)>(2*oT+1):
        mjd0 = np.mean([o[0] for o in _obs if not o[0] is None ])
        p['MJD0'] = mjd0
        doNotFit = ['MJD0', 'P-FACTOR', 'METAL', 'd_kpc', 'DIAMAVG', 'PERIOD', 'VRAD A0', 'PERIOD']
        doNotFit += list(filter(lambda x: x.startswith('VRAD '), p.keys()))
        if not (any([o[1].startswith('teff') for o in _obs]) and
                any([o[1].startswith('color') for o in _obs])):
            doNotFit.append('E(B-V)')
        t0 = time.time()
        fit = spips.fit(_obs, p, doNotFit=doNotFit, verbose=maxCores>1, plot=False, maxCores=maxCores,
                        normalizeErrors='techniques', ftol=1e-3,
                        follow=['E(B-V)', 'TEFF A0'])
        p = fit['best']
        print('chi2=%.3f [in %.1fs]'%(fit['chi2'], time.time()-t0))
        execution['Teff'] = time.time()-t0
        # -- 

    # == only Vrad =====================================
    print('fitting VRAD only...', end=' ')
    tech = ['vrad']
    _obs = [o for o in obs if any([o[1].startswith(t) for t in tech])]

    p['MJD0'] += np.round((np.mean([o[0] for o in _obs])-p['MJD0'])/p['PERIOD'], 0)*p['PERIOD']
    p['VRAD A0'] = np.mean([o[-2] for o in _obs])
    p['VRAD A1'] = np.ptp([o[-2] for o in _obs])/oV
    p['VRAD PHI1'] = 0
    for k in range(oV-1):
        p['VRAD R%d'%(k+2)] = 1.0
        p['VRAD PHI%d'%(k+2)] = 0 

    doNotFit = ['MJD0', 'P-FACTOR', 'METAL', 'd_kpc', 'DIAMAVG', 'TEFF A0', 'E(B-V)']
    for k in p:
        if k.startswith('TEFF '):
            doNotFit.append(k)
    t0 = time.time()
    fit = spips.fit(_obs, p, doNotFit=doNotFit, verbose=False, plot=False, maxCores=1, ftol=1e-3, 
                    epsfcn=1e-8)
    p = fit['best']
    print('chi2=%.3f [in %.1fs]'%(fit['chi2'], time.time()-t0))
    execution['Vrad'] = time.time()-t0

    if True and 'Teff' in execution:
        # -- not really needed?

        print('fitting all data with fixed Vrad and Teff profiles', end=' ')
        doNotFit = ['MJD0', 'P-FACTOR', 'METAL', 'EXCESS EXP', 'EXCESS WL0']
        # -- keep Teff and Vrad profile fixed
        for k in p:
            if k.startswith('VRAD R') or k.startswith('VRAD A'):
                doNotFit.append(k)
            if k.startswith('TEFF R') or k.startswith('TEFF A'):
                doNotFit.append(k)
            if ' PHI' in k and not k.endswith('PHI1'):
                doNotFit.append(k)
        mjd0 = np.mean([o[0] for o in obs if not o[0] is None ])
        p['MJD0'] += np.round((mjd0-p['MJD0'])/p['PERIOD'], 0)*p['PERIOD']
        t0 = time.time()
        fit = spips.fit(obs, p, doNotFit=doNotFit, verbose=True, plot=False, maxCores=4,
                        normalizeErrors='techniques', ftol=1e-2,
                        follow=['E(B-V)', 'd_kpc'])
        p = fit['best']
        print('chi2=%.3f [in %.1fs]'%(fit['chi2'], time.time()-t0))
        execution['All, fixed Vrad & Teff'] = time.time()-t0

    # -- fitting everything
    print('fitting all parameters to all data ', end=' ')
    doNotFit = ['MJD0', 'P-FACTOR', 'METAL']

    if withIrExcess:
        p.update({
                  'EXCESS EXP':   0.4 , # exponential law for IR Excess
                  'EXCESS SLOPE': 0.0 , # slope for IR excess
                  'EXCESS WL0':   1.2 , # starting WL, in um, for IR Excess
                  })
        doNotFit += ['EXCESS EXP', 'EXCESS WL0']
    if changeOfPeriodOrder>0:
        for k in range(changeOfPeriodOrder):
            if not 'PERIOD%d'%(k+1) in p:
                p['PERIOD%d'%(k+1)] = 1e-3    
    mjd0 = np.mean([o[0] for o in obs if not o[0] is None ])
    p['MJD0'] += np.round((mjd0-p['MJD0'])/p['PERIOD'], 0)*p['PERIOD']
    t0 = time.time()
    fit = spips.fit(obs, p, doNotFit=doNotFit, verbose=True, plot=False, maxCores=4,
                    normalizeErrors='techniques', ftol=1e-5, epsfcn=1e-8,
                    follow=['d_kpc', 'E(B-V)', 'DIAMAVG', 'TEFF A0'])
    p = fit['best']
    print('chi2=%.3f [in %.0fs]'%(fit['chi2'], time.time()-t0))
    execution['All'] = time.time()-t0
    print(execution)
    print('total execution %.1fmin'%(np.sum([execution[k] for k in execution])/60))
    return fit

def fitFromScratchSplines(obs, firstGuess={'PERIOD':5.3663, 'd_kpc':.25, 'DIAMAVG':1.5, 'TEFF VAL0':6000}, 
                          changeOfPeriodOrder=1, withIrExcess=True, oT=4, oV=6):
    """
    find SPIPS Spline model from dataset "obs"
    
    firstGuess: provide initial guess as dict. Mandatory is "PERIOD" (in days). Recommended are
        'd_kpc': distance in kpc
        'DIAMAVG': average angular diameter in mas
        'TEFF VAL0': average effective temperature in K
    
    changeOfPeriodOrder: 0 for no change of period, 1 for linear etc (default 1)
    
    withIrExcess: True or False (default True)

    oT, oV: number of spline Nodes for Teff and Vrad (default 4, 6)

    """

    # -- data without Teff and color, just for test
    #obs = [o for o in obs if not any([o[1].startswith(x) for x in ['teff', 'color']])]

    p = {'P-FACTOR':1.27,
         'd_kpc': 1.0,
         'TEFF VAL0': 6000,
         'TEFF POW': 1.2,         
         'TEFF PHI0': .8,
         'VRAD POW': 1.5,
         'VRAD PHI0': 0.8,
         'VRAD VAL0': 0.,
         'METAL': 0.0,
         'E(B-V)':0.0, 
         'DIAMAVG':1.0,
      }
    p.update(firstGuess)
    execution = {}

    # == fit teff and colors
    print('fitting Teff profile to teff data...', end=' ')

    for k in range(oT-1):
        p['TEFF DVAL%d'%(k+1)] = np.interp(k, [0, oT/3, oT/2, 2*oT/3, oT], [0.25, 1, 0, -1, 0])
        p['TEFF DVAL%d'%(k+1)] *= 500
    
    tech = ['teff']
    _obs = [o for o in obs if any([o[1].startswith(t) for t in tech])]
    maxCores=1
    if len(_obs)<=(2*oT+1):
        print('not enough data')
        print('fitting Teff profile to teff and color data...', end=' ')
        # -- try including colors
        tech = ['teff', 'color']
        _obs = [o for o in obs if any([o[1].startswith(t) for t in tech])]
        maxCores = 4

    if len(_obs)>(2*oT+1):
        mjd = [o[0] for o in _obs if not o[0] is None ]
        print('[MJD ptp/period= %.1f]'%(np.mean(mjd)/p['PERIOD']), end=' ')
        #p['MJD0'] += np.round((np.mean(mjd)-p['MJD0'])/p['PERIOD'], 0)*p['PERIOD']
        p['MJD0'] = np.mean(mjd)
        
        doNotFit = ['MJD0', 'P-FACTOR', 'METAL', 'd_kpc', 'DIAMAVG', 'PERIOD']

        doNotFit += list(filter(lambda x: x.startswith('VRAD '), p.keys()))
        if not (any([o[1].startswith('teff') for o in _obs]) and
                any([o[1].startswith('color') for o in _obs])):
            doNotFit.append('E(B-V)')
        t0 = time.time()
        fit = spips.fit(_obs, p, doNotFit=doNotFit, verbose=maxCores>1, plot=False, maxCores=maxCores,
                        normalizeErrors='techniques', ftol=1e-3, 
                        follow=['E(B-V)', 'TEFF VAL0'])
        p = fit['best']
        print('chi2=%.3f [in %.1fs]'%(fit['chi2'], time.time()-t0))
        execution['Teff'] = time.time()-t0

    # == only Vrad =====================================
    print('fitting VRAD only...', end=' ')
    tech = ['vrad']
    _obs = [o for o in obs if any([o[1].startswith(t) for t in tech])]
    mjd = [o[0] for o in _obs if not o[0] is None ]
    print('[MJD ptp/period=%.0f]'%(np.mean(mjd)/p['PERIOD']), end=' ')
    p['MJD0'] += np.round((np.mean(mjd)-p['MJD0'])/p['PERIOD'], 0)*p['PERIOD']
    p['VRAD VAL0'] = np.mean([o[-2] for o in _obs]) 
    p['VRAD PHI0'] = 0.8
    p['VRAD POW'] = 1.8

    for k in range(oV-1):
        p['VRAD DVAL%d'%(k+1)] = np.interp(k, [0, oV/3, oV/2, 2*oV/3, oV], [0.25, 1, 0.01, -1, -0.01])
        p['VRAD DVAL%d'%(k+1)] *= np.ptp([o[-2] for o in _obs])/2
    doNotFit = ['MJD0', 'P-FACTOR', 'METAL', 'd_kpc', 'DIAMAVG',  'E(B-V)','PERIOD']
    for k in p:
        if k.startswith('TEFF '):
            doNotFit.append(k)

    t0 = time.time()
    fit = spips.fit(_obs, p, doNotFit=doNotFit, verbose=False, plot=False, maxCores=1)
    p = fit['best']
    print('chi2=%.3f [in %.1fs]'%(fit['chi2'], time.time()-t0))
    execution['Vrad'] = time.time()-t0

    if True and 'Teff' in execution:
        # -- not really needed?
        print('fitting all data with fixed Vrad and Teff profiles', end=' ')
        doNotFit = ['MJD0', 'P-FACTOR', 'METAL', 'EXCESS EXP', 'EXCESS WL0']
        # -- keep Teff and Vrad profile fixed
        for k in p:
            if k.startswith('VRAD') and k!='VRAD PHI0':
                doNotFit.append(k)
            if k.startswith('TEFF') and k!='TEFF PHI0':
                doNotFit.append(k)

        mjd0 = np.mean([o[0] for o in obs if not o[0] is None ])
        p['MJD0'] += np.round((mjd0-p['MJD0'])/p['PERIOD'], 0)*p['PERIOD']
        t0 = time.time()
        fit = spips.fit(obs, p, doNotFit=doNotFit, verbose=True, plot=False, maxCores=4,
                        normalizeErrors='techniques', ftol=1e-2, 
                        follow=['E(B-V)', 'd_kpc', 'DIAMAVG'])
        p = fit['best']
        print('chi2=%.3f [in %.1fs]'%(fit['chi2'], time.time()-t0))
        execution['All, fixed Vrad & Teff'] = time.time()-t0

    # -- fitting everything
    print('fitting all parameters to all data ', end=' ')
    doNotFit = ['MJD0', 'P-FACTOR', 'METAL']

    if withIrExcess:
        p.update({
                  'EXCESS EXP':   0.4 , # exponential law for IR Excess
                  'EXCESS SLOPE': 0.0 , # slope for IR excess
                  'EXCESS WL0':   1.2 , # starting WL, in um, for IR Excess
                  })
        doNotFit += ['EXCESS EXP', 'EXCESS WL0']
    if changeOfPeriodOrder>0:
        for k in range(changeOfPeriodOrder):
            if not 'PERIOD%d'%(k+1) in p:
                p['PERIOD%d'%(k+1)] = 1e-3    

    mjd0 = np.mean([o[0] for o in obs if not o[0] is None ])
    p['MJD0'] += np.round((mjd0-p['MJD0'])/p['PERIOD'], 0)*p['PERIOD']
    t0 = time.time()
    fit = spips.fit(obs, p, doNotFit=doNotFit, verbose=True, plot=False, maxCores=4,
                    normalizeErrors='techniques', ftol=1e-4, 
                    follow=['d_kpc', 'E(B-V)', 'DIAMAVG', 'TEFF VAL0', 'VRAD POW'])
    p = fit['best']
    print('chi2=%.3f [in %.0fs]'%(fit['chi2'], time.time()-t0))
    execution['All'] = time.time()-t0
    print(execution)
    print('total execution %.1fmin'%(np.sum([execution[k] for k in execution])/60))
    return fit

def findPeriod(obs, period=5, relativeRange=0.1, Np=1000, plot=None):
    """

    """
    typs = [o[1]+str(o[2]) if o[1].startswith('mag') or o[1].startswith('color') else o[1] for o in obs ]

    res = 0
    p = np.linspace(1-relativeRange, 1+relativeRange, Np)*period
    # -- get periodogram for each observable separatly
    for typ in set(typs):
        t = np.array([o[0] for i,o in enumerate(obs) if typs[i]==typ and not o[0] is None and o[0]>1]) 
        v = np.array([o[-2] for i,o in enumerate(obs) if typs[i]==typ and not o[0] is None and o[0]>1]) 
        res += periodogram(t, v, p)
    # -- find maximum
    i0 = np.argmax(res)
    ## d/dx ax**2+bx+c == 0 -> 2ax + b = 0 -> x = -b/2a
    if i0>0 and i0<len(res)-1:
        c = np.polyfit(p[i0-1:i0+1]-p[i0], res[i0-1:i0+1], 2)
        pmax = p[i0]-c[1]/2/c[0]
        #print('period:', pmax)
    else:
        #print("can't find maximum")
        pmax = None
    if plot is True:
        plot = 1
    if type(plot)==int:
        plt.close(plot)
        plt.figure(plot, figsize=(9,7))
        plt.subplot(121)
        plt.plot(p, res, label='power')
        plt.vlines(pmax, 0, max(res), linestyle='dotted', color='g', 
                label='best period')
        plt.xlabel('period (days)')
        plt.legend()
        mjd0 = np.mean([o[0] for o in obs if not o[0] is None and o[0]>1])
        for j,typ in enumerate(set(typs)):
            t = np.array([o[0] for i,o in enumerate(obs) if typs[i]==typ and not o[0] is None and o[0]>1]) 
            v = np.array([o[-2] for i,o in enumerate(obs) if typs[i]==typ and not o[0] is None and o[0]>1]) 
            ax = plt.subplot(len(set(typs)),2,2*j+2)
            ax.yaxis.set_visible(False)
            plt.plot(((t-mjd0)/pmax)%1, v, '.k')
            plt.title(typ, x=0.5, y=0.5, fontsize=6)
            plt.xlim(-0.05,1.05)
        plt.tight_layout()
        plt.subplots_adjust(hspace=0)
    return pmax

def periodogram(t, v, p):
    """
    t: times (1D ndarray)
    v: values (1D ndarray, same length as t)
    p: periods (1D array, same units as t)
    
    return power spectrum: 1D ndarray, same length as p
    """
    _v = (v-np.mean(v))/np.std(v)
    return np.abs(np.mean(np.cos(2*np.pi*t[:,None]/p[None,:])*_v[:,None] + 
                       1j*np.sin(2*np.pi*t[:,None]/p[None,:])*_v[:,None], axis=0))**2


    



