import spips
import numpy as np


def derivedUncer(f, o, N=128):
    """
    f: results from spips.fits()
    o: observations used in the fit, to compute uncertainties
    N: number of randomisation of parameters 
    uncertainties for derived parameters: <Teff>, <Lbol> and magnitudes
    """
    _phi = np.linspace(0,1,101)[:-1]
    _o = []
    # -- all photometric bands
    bands = set([x[2] for x in o if 'mag' in x[1]])
    for b in bands:
        _o.extend([(f['best']['MJD0']+x*f['best']['PERIOD'], 
                    'mag', b, 0.0, 0.0) for x in _phi])
    _o.extend([(f['best']['MJD0']+x*f['best']['PERIOD'], 
                'teff', 0.0, 0.0) for x in _phi])
    _o.extend([(f['best']['MJD0']+x*f['best']['PERIOD'], 
                 'diam', 0.0, 0.0) for x in _phi])
    _o.extend([(f['best']['MJD0']+x*f['best']['PERIOD'], 
                 'vrad', 0.0, 0.0) for x in _phi])

    f['func'] = spips.model
    f = spips.dpfit.randomParam(f, x=_o, N=N, multi=True)

    # -- retrieve results
    ref = {'phi':_phi}
    for b in bands:
        ref[b] = np.array([f['r_y'][i] for i,x in enumerate(_o) if 'mag' in x[1] and b in x[2]])
        ref[b+' +1sigma'] = np.array([f['r_yp1s'][i] for i,x in enumerate(_o) if 'mag' in x[1] and b in x[2]])
        ref[b+' -1sigma'] = np.array([f['r_ym1s'][i] for i,x in enumerate(_o) if 'mag' in x[1] and b in x[2]])
        print('<'+b+'> = ', '%.4f'%np.mean(ref[b]), 
            '+- %.4f'%(0.5*np.mean(ref[b+' +1sigma']-ref[b+' -1sigma'])))

    ref['teff'] = np.array([f['r_y'][i] for i,x in enumerate(_o) if 'teff' in x[1]])
    ref['teff +1sigma'] = np.array([f['r_yp1s'][i] for i,x in enumerate(_o) if 'teff' in x[1]])
    ref['teff -1sigma'] = np.array([f['r_ym1s'][i] for i,x in enumerate(_o) if 'teff' in x[1]])
    print('<TEFF> =', '%.1f'%np.mean(ref['teff']),
        '+- %.1f'%(0.5*np.mean(ref['teff +1sigma']-ref['teff -1sigma'])), '(K)')

    ref['diam'] = np.array([f['r_y'][i] for i,x in enumerate(_o) if 'diam' in x[1]])
    ref['diam +1sigma'] = np.array([f['r_yp1s'][i] for i,x in enumerate(_o) if 'diam' in x[1]])
    ref['diam -1sigma'] = np.array([f['r_ym1s'][i] for i,x in enumerate(_o) if 'diam' in x[1]])
    print('<DIAM> =', '%.5f'%np.mean(ref['diam']),
        '+- %.5f'%(0.5*np.mean(ref['diam +1sigma']-ref['diam -1sigma'])), '(mas)')

    ref['vrad'] = np.array([f['r_y'][i] for i,x in enumerate(_o) if 'vrad' in x[1]])
    ref['vrad +1sigma'] = np.array([f['r_yp1s'][i] for i,x in enumerate(_o) if 'vrad' in x[1]])
    ref['vrad -1sigma'] = np.array([f['r_ym1s'][i] for i,x in enumerate(_o) if 'vrad' in x[1]])
    print('<VRAD> =', '%.2f'%np.mean(ref['vrad']),
        '+- %.2f'%(0.5*np.mean(ref['vrad +1sigma']-ref['vrad -1sigma'])), '(km/s)')


    return ref
