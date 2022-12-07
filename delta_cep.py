import spips
import delta_cep_data # contains all the observations
import os
import numpy as np
import time
# -- load data from aux file
obs = delta_cep_data.data

#import gaiadr3
#obs += gaiadr3.getDr3Phot('delta Cep')['SPIPS']

# -- this is the model using splines
# Parameters description:
# 'DIAMAVG':      1.45573 , # Average angular diameter (mas)
# 'E(B-V)':       0.09109 , # reddenning
# 'EXCESS EXP':   0.4 , # exponential law for IR Excess
# 'EXCESS SLOPE': 0.06278 , # slope for IR excess
# 'EXCESS WL0':   1.2 , # starting WL, in um, for IR Excess
# 'METAL':        0.06 , # metalicity / solar
# 'MJD0':         48304.7362421 , # 0-phase
# 'P-FACTOR':     1.2687 , # projection factor
# 'PERIOD':       5.36626201 , # period in days
# 'PERIOD1':      -0.0086 , # period change, in s/yrs
# 'TEFF DVAL1':   -219.91 , # offset of split point to VAL0
# 'TEFF DVAL2':   163.7 , # offset of split point to VAL0
# 'TEFF DVAL3':   709.69 , # offset of split point to VAL0
# 'TEFF DVAL4':   571.61 , # offset of split point to VAL0
# 'TEFF PHI0':    0.8844 , # ref phase for spline comb
# 'TEFF POW':     1.5856 re, # spline comb spread (1 is regular)
# 'TEFF VAL0':    5704.094 , # first spline node value
# 'VRAD DVAL1':   14.002 , #
# 'VRAD DVAL2':   22.084 , #
# 'VRAD DVAL3':   18.935 , #
# 'VRAD DVAL4':   -1.221 , #
# 'VRAD DVAL5':   -13.712 , #
# 'VRAD PHI0':    0.84471 , #
# 'VRAD POW':     1.9098 , #
# 'VRAD VAL0':    -21.371 , #
# 'd_kpc':        0.274 , # distance in kilo-pc

p_splines = {'DIAMAVG':      1.45616 , # +/- 0.00105
            'E(B-V)':       0.0908 , # +/- 0.00192
            'EXCESS EXP':   0.4 ,
            'EXCESS SLOPE': 0.06187 , # +/- 0.0015
            'EXCESS WL0':   1.2 ,
            'METAL':        0.06 ,
            'MJD0':         48304.7362421 ,
            'P-FACTOR':     1.2712 , # +/- 0.0177
            'PERIOD':       5.36627863 , # +/- 5.5e-06
            'PERIOD1':      -0.0851 , # +/- 0.0293
            'TEFF DVAL1':   -221.253 , # +/- 4.327
            'TEFF DVAL2':   167.46 , # +/- 18.02
            'TEFF DVAL3':   711.8 , # +/- 13.76
            'TEFF DVAL4':   577.53 , # +/- 12.55
            'TEFF PHI0':    0.88491 , # +/- 0.00161
            'TEFF POW':     1.5952 , # +/- 0.0372
            'TEFF VAL0':    5702.675 , # +/- 6.041
            'VRAD DVAL1':   13.882 , # +/- 0.275
            'VRAD DVAL2':   22.02 , # +/- 0.162
            'VRAD DVAL3':   19.021 , # +/- 0.406
            'VRAD DVAL4':   -1.405 , # +/- 0.524
            'VRAD DVAL5':   -13.67 , # +/- 0.321
            'VRAD PHI0':    0.84753 , # +/- 0.00197
            'VRAD POW':     1.8918 , # +/- 0.0478
            'VRAD VAL0':    -21.379 , # +/- 0.17
            'd_kpc':        0.274 ,
            }


# Alternatively, the TEFF and VRAD profiles can be described using FOURIER parameters:
# 'TEFF A0':      5887.886 , # average Teff
# 'TEFF A1':      469.915 , # amplitude of first harmonic
# 'TEFF PHI1':    -0.3581 , # phase of first harmonic
# 'TEFF PHI2':    -0.2403 , # etc.
# 'TEFF PHI3':    0.3564 ,
# 'TEFF PHI4':    0.7853 ,
# 'TEFF PHI5':    1.71 ,
# 'TEFF R2':      0.39556 , # amp1/amp2
# 'TEFF R3':      0.15563 , # etc.
# 'TEFF R4':      0.06821 ,
# 'TEFF R5':      0.02028 ,

p_fourier = {    'DIAMAVG'     : 1.45108, # +/- 0.00144
        'E(B-V)'      : 0.09549, # +/- 0.00212
        'EXCESS EXP'  : 0.4 ,
        'EXCESS SLOPE': 0.05929, # +/- 0.00184
        'EXCESS WL0'  : 1.2 ,
        'METAL'       : 0.06 ,
        'MJD0'        : 48304.732663390016 ,
        'P-FACTOR'    : 1.2375, # +/- 0.0211
        'PERIOD'      : 5.36627561, # +/- 0.00000588
        'PERIOD1'     : -0.0600, # +/- 0.0271
        'TEFF A0'     : 5909.59, # +/- 6.17
        'TEFF A1'     : 472.07, # +/- 3.66
        'TEFF PHI1'   : 5.91673, # +/- 0.00580
        'TEFF PHI2'   : 5.4450, # +/- 0.0154
        'TEFF PHI3'   : 1.3499, # +/- 0.0415
        'TEFF PHI4'   : 0.5966, # +/- 0.0866
        'TEFF PHI5'   : 2.338, # +/- 0.237
        'TEFF R2'     : 0.39779, # +/- 0.00588
        'TEFF R3'     : -0.15503, # +/- 0.00544
        'TEFF R4'     : -0.06561, # +/- 0.00523
        'TEFF R5'     : 0.02219, # +/- 0.00500
        'VRAD A0'     : -18.4306, # +/- 0.0798
        'VRAD A1'     : 15.483, # +/- 0.130
        'VRAD PHI1'   : 2.17859, # +/- 0.00796
        'VRAD PHI2'   : -1.5853, # +/- 0.0212
        'VRAD PHI3'   : 6.1507, # +/- 0.0449
        'VRAD PHI4'   : 4.5615, # +/- 0.0590
        'VRAD PHI5'   : 2.7519, # +/- 0.0905
        'VRAD PHI6'   : 10.574, # +/- 0.163
        'VRAD PHI7'   : 2.280, # +/- 0.302
        'VRAD PHI8'   : 0.719, # +/- 0.602
        'VRAD R2'     : -0.41554, # +/- 0.00678
        'VRAD R3'     : -0.21487, # +/- 0.00535
        'VRAD R4'     : 0.12098, # +/- 0.00736
        'VRAD R5'     : -0.05963, # +/- 0.00760
        'VRAD R6'     : -0.03824, # +/- 0.00568
        'VRAD R7'     : 0.01914, # +/- 0.00526
        'VRAD R8'     : -0.00986, # +/- 0.00602
        'd_kpc'       : 0.274 ,
        }

def fit(p=None):
    """
    p: dictionnary containgin the model
    """
    if p is None:
        p = p_splines
    # - list parameters which we do not wish to fit
    doNotFit= ['MJD0','METAL', 'd_kpc', 'EXCESS WL0', 'EXCESS EXP']
    fitOnly = None
    # - alternatively, we can list the only parameters we wich to fit
    # fitOnly = filter(lambda x: x.startswith('TEFF ') or x.startswith('VRAD '), p.keys())
    # obs = filter(lambda o: 'vrad' in o[1] or 'teff' in o[1], obs)
    #fitOnly=['P-FACTOR', 'DIAMAVG']
    f = spips.fit(obs, p, doNotFit=doNotFit, fitOnly=fitOnly,
            normalizeErrors='techniques', # 'observables' is the alternative
            ftol=5e-4, # stopping tolerance on chi2
            epsfcn=1e-8, # by how much parameters will vary
            maxfev=500, # maximum number of iterations
            maxCores=4, # max number of CPU cores, None will use all available
            starName='delta Cep',
            follow=['P-FACTOR'], # list here parameters you want to see during fit
            exportFits=True,
            )
    spips.dispCor(f) # show the correlation matrix between parameters
    return

def show(p=None):
    """
    p: dictionnary containgin the model
    """
    if p is None:
        p = p_splines
    Y = spips.model(obs, p, starName='delta Cep', verbose=True, plot=True)

def fitsDemo(mode='export', p=None):
    if p is None:
        p = p_splines

    if mode=='export':
        Y = spips.model(obs, p, starName='delta Cep',
                        exportFits=True, verbose=True)
    elif mode=='import':
        filename = os.path.join('DELTA_CEP', 'delta_cep.fits')
        if os.path.exists(filename):
            tmp = spips.importFits(filename)
            Y = spips.model(tmp[1], tmp[0], starName='delta Cep', verbose=True,
                            plot=True)
        else:
            print('ERROR:', filename, 'does not exist!')
    else:
        print("use: fitsDemo(mode='export')")
        print("  or fitsDemo(mode='import')")


def fitFromScratchFourier(obs, firstGuess={'PERIOD':5.366, 'd_kpc':.25, 'DIAMAVG':1.5, 'TEFF A0':6000}, 
                            changeOfPeriodOrder=1, withIrExcess=True, oT=5, oV=8):
    """
    find SPIPS Fourier model from dataset "obs"
    
    firstGuess: provide initial guess as dict. Mandatory is "PERIOD" (in days). Recommended are
        'd_kpc': distance in kpc
        'DIAMAVG': average angular diameter in mas
        'TEFF A0': average effective temperature in K
    
    changeOfPeriodOrder: 0 for no change of period, 1 for linear etc (default 1)
    
    withIrExcess: True or False (default True)

    oT, oV: number of Fourier componenents (default 5, 8)
    """

    # -- data without Teff and color, just for test
    #obs = [o for o in obs if not any([o[1].startswith(x) for x in ['teff', 'color']])]

    p = {'P-FACTOR':1.27,
         'd_kpc': 1.0,
         'TEFF A0': 6000,
         'METAL': 0.0,
         'E(B-V)':0.0, 
         'DIAMAVG':1.0,
      }
    p.update(firstGuess)
    execution = {}
    oV = 8 # Fourier Components for first
    oT = 6
    # == only Vrad =====================================
    print('fitting VRAD only...', end=' ')
    tech = ['vrad']
    _obs = [o for o in obs if any([o[1].startswith(t) for t in tech])]

    p['MJD0'] = np.mean([o[0] for o in _obs])
    p['VRAD A0'] = np.mean([o[-2] for o in _obs])
    p['VRAD A1'] = np.ptp([o[-2] for o in _obs])/oV
    p['VRAD PHI1'] = 0
    for k in range(oV-1):
        p['VRAD R%d'%(k+2)] = 1.0
        p['VRAD PHI%d'%(k+2)] = 0 

    doNotFit = ['MJD0', 'P-FACTOR', 'METAL', 'd_kpc', 'DIAMAVG', 'TEFF A0', 'E(B-V)']
    t0 = time.time()
    fit = spips.fit(_obs, p, doNotFit=doNotFit, verbose=False, plot=False, maxCores=1)
    p = fit['best']
    print('chi2=%.3f [in %.1fs]'%(fit['chi2'], time.time()-t0))
    execution['Vrad'] = time.time()-t0
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
    if len(_obs)<=(2*oT+1):
        print('not enough data')
        print('fitting Teff profile to teff and color data...', end=' ')
        # -- try including colors
        tech = ['teff', 'color']
        _obs = [o for o in obs if any([o[1].startswith(t) for t in tech])]

    if len(_obs)>(2*oT+1):
        mjd0 = np.mean([o[0] for o in _obs if not o[0] is None ])
        p['MJD0'] += np.round((mjd0-p['MJD0'])/p['PERIOD'], 0)*p['PERIOD']
        doNotFit = ['MJD0', 'P-FACTOR', 'METAL', 'd_kpc', 'DIAMAVG', 'PERIOD']
        doNotFit += list(filter(lambda x: x.startswith('VRAD '), p.keys()))
        if not (any([o[1].startswith('teff') for o in _obs]) and
                any([o[1].startswith('color') for o in _obs])):
            doNotFit.append('E(B-V)')
        t0 = time.time()
        fit = spips.fit(_obs, p, doNotFit=doNotFit, verbose=False, plot=False, maxCores=1,
                        normalizeErrors='techniques', ftol=1e-3)
        p = fit['best']
        print('chi2=%.3f [in %.1fs]'%(fit['chi2'], time.time()-t0))
        execution['Teff'] = time.time()-t0
        # -- 
        if True:
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
                            normalizeErrors='techniques', ftol=1e-2)
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
                    normalizeErrors='techniques', ftol=1e-4)
    p = fit['best']
    print('chi2=%.3f [in %.0fs]'%(fit['chi2'], time.time()-t0))
    execution['All'] = time.time()-t0
    print(execution)
    print('total execution %.1fmin'%(np.sum([execution[k] for k in execution])/60))
    return p

def fitFromScratchSplines(obs, firstGuess={'PERIOD':5.366, 'd_kpc':.25, 'DIAMAVG':1.5, 'TEFF VAL0':6000}, 
                          changeOfPeriodOrder=1, withIrExcess=True, oT=4, oV=6):
    """
    find SPIPS Spline model from dataset "obs"
    
    firstGuess: provide initial guess as dict. Mandatory is "PERIOD" (in days). Recommended are
        'd_kpc': distance in kpc
        'DIAMAVG': average angular diameter in mas
        'TEFF VAL0': average effective temperature in K
    
    changeOfPeriodOrder: 0 for no change of period, 1 for linear etc (default 1)
    
    withIrExcess: True or False (default True)

    oT, oV: number of spline Nodes (default 4, 6)

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
    oV = 6 # Spline nodes
    oT = 4

    # == fit teff and colors
    print('fitting Teff profile to teff data...', end=' ')

    for k in range(oT-1):
        p['TEFF DVAL%d'%(k+1)] = np.interp(k, [0, oT/3, oT/2, 2*oT/3, oT], [0.25, 1, 0, -1, 0])
        p['TEFF DVAL%d'%(k+1)] *= 500
    
    tech = ['teff']
    _obs = [o for o in obs if any([o[1].startswith(t) for t in tech])]
    if len(_obs)<=(2*oT+1):
        print('not enough data')
        print('fitting Teff profile to teff and color data...', end=' ')
        # -- try including colors
        tech = ['teff', 'color']
        _obs = [o for o in obs if any([o[1].startswith(t) for t in tech])]

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
        fit = spips.fit(_obs, p, doNotFit=doNotFit, verbose=False, plot=False, maxCores=1,
                        normalizeErrors='techniques', ftol=1e-3)
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
                        normalizeErrors='techniques', ftol=1e-2)
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
                    normalizeErrors='techniques', ftol=1e-4)
    p = fit['best']
    print('chi2=%.3f [in %.0fs]'%(fit['chi2'], time.time()-t0))
    execution['All'] = time.time()-t0
    print(execution)
    print('total execution %.1fmin'%(np.sum([execution[k] for k in execution])/60))
    return p



    



