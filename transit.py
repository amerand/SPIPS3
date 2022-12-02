import numpy as np
from matplotlib import pyplot as plt
import os
import dpfit
from scipy import special


# ===  points on a disk: =======================================
_disk = {'x':np.array([]), 'y':np.array([])}
Nr = 200
for r in np.linspace(0,1,Nr):
    n = max(int(2*Nr*np.sqrt(r)),1)
    t = np.linspace(0,2*np.pi,n+1)[1:]
    _disk['x'] = np.append(_disk['x'], r*np.cos(t))
    _disk['y'] = np.append(_disk['y'], r*np.sin(t))
# -- surface element
_disk['ds'] = np.pi/float(len(_disk['x']))

# ====  read Neilson linear laws for dwarfs =================================
col2 = ['Teff', 'logg', 'mass']
for b in ['B','V','R','I','H','K','CoR','Kep']:
    for i in ['2']:
        col2.append('f%s(%s)'%(i,b))
nei2 = {c:[] for c in col2}
# -- read the 4 parameters for each band
for l in open('./DATA/LD_NEILSON/DWARFS/table2.dat').readlines():
    for k,w in enumerate(l.split()):
        nei2[col2[k]].append(float(w))
for c in col2:
    nei2[c] = np.array(nei2[c])

# ====  read Neilson 4P laws for dwarfs =================================
col5 = ['Teff', 'logg', 'mass']
for b in ['B','V','R','I','H','K','CoR','Kep']:
    for i in ['1','2','3','4']:
        col5.append('f%s(%s)'%(i,b))
nei5 = {c:[] for c in col5}
# -- read the 4 parameters for each band
for l in open('DATA/LD_NEILSON/DWARFS/table5.dat').readlines():
    for k,w in enumerate(l.split()):
        nei5[col5[k]].append(float(w))
for c in col5:
    nei5[c] = np.array(nei5[c])


def test2(filename=None, noisePpm=None, showV2=False):
    """
    this should go to 0 at r = 0.0:
    """
    global nei5, _tot
    # === setup transit with SATLAS KEPLER profile
    if filename is None:
        #filename = 'DATA/LD_NEILSON/DWARFS/spheric/ld_satlas_surface.2t5000g400m8.dat'
        #filename = 'DATA/LD_NEILSON/DWARFS/spheric/ld_satlas_surface.2t3900g450m5.dat'
        filename = 'DATA/LD_NEILSON/DWARFS/spheric/ld_satlas_surface.2t6000g450m11.dat'
    # -- read SATLAS file:
    f = open(filename)
    cols = ['mu','B','V','R','I','H','K','CoR','Kep']
    data = {c:[] for c in cols}
    for l in f.readlines():
        for k,c in enumerate(cols):
            data[c].append(float(l.split()[k]))
    for c in cols:
        data[c] = np.array(data[c])
    f.close()
    laws = {'sigmo power':{'mu0':0.04, 'c':.8,  'ap':0.7, 'p':.2, 'e':1.6},
            #'sigmo power':{'mu0':0.04, 'c':.8, 'ap':0.7, 'p':.3},
            }

    Teff = int(filename.split('surface.2t')[1].split('g')[0])
    logg = float(filename.split('surface.2t')[1].split('g')[1].split('m')[0])/100.
    mass = float(filename.split('surface.2t')[1].split('g')[1].split('m')[1].split('.')[0])/10
    # -- Find P4 coef:

    i5 = np.where((nei5['Teff']==Teff)*(nei5['logg']==logg)*(nei5['mass']==mass))[0][0]
    # laws['P4'] = {'f1': nei5['f1(Kep)'][i5],
    #       'f2': nei5['f2(Kep)'][i5],
    #       'f3': nei5['f3(Kep)'][i5],
    #       'f4': nei5['f4(Kep)'][i5]}

    i2 = np.where((nei2['Teff']==Teff)*(nei2['logg']==logg)*(nei2['mass']==mass))[0][0]
    laws['Lin'] = {'f2': nei2['f2(Kep)'][i2]}
    colors= {'P4':'r', 'sigmo power':'b', 'Lin':'g', 'sigmo lin':'y'}

    # -- Figures: ----
    plt.close(1)
    plt.figure(1, figsize=(10,8))
    plt.subplots_adjust(wspace=0.3, right=0.95, left=0.1, top=0.93)
    if showV2:
        # - ax1, ax2: CLD profile and resiuals:
        ax1 = plt.subplot(231)
        ax1.set_title(os.path.basename(filename))
        ax2 = plt.subplot(234, sharex=ax1)
        # - ax3, ax4: transit and resiuals:
        ax3 = plt.subplot(232)
        ax4 = plt.subplot(235, sharex=ax3)
        # - ax5, ax6: visibility and resiuals:
        ax5 = plt.subplot(233)
        ax5.set_title('interferometry')
        ax6 = plt.subplot(236, sharex=ax5)
    else:
        # - ax1, ax2: CLD profile and resiuals:
        ax1 = plt.subplot(221)
        ax1.set_title(os.path.basename(filename))
        ax2 = plt.subplot(223, sharex=ax1)
        # - ax3, ax4: transit and resiuals:
        ax3 = plt.subplot(222)
        ax4 = plt.subplot(224, sharex=ax3)

    r = np.sqrt(1-data['mu']**2)[::-1]
    Ir = data['Kep'][::-1]
    r0 = maxVar(Ir,r)
    r /= r0
    func(r, None, Ir=Ir)

    ax1.plot(data['mu'], data['Kep'], '-k', linewidth=3, label='SATLAS', alpha=0.2)
    #ax1.plot(r*r0, Ir, '-k', linewidth=3, label='SATLAS', alpha=0.4)

    # -- compute transit SATLAS
    p = {'b':0.43, 'rp/r*':0.18, 'r*':1} # HATS-6b
    #p = {'b':0.5, 'rp/r*':0.10, 'r*':1} # random

    rt = np.linspace(0,1.1*p['r*']*(1+2*p['rp/r*']),100)
    f = func(rt, p)
    if not noisePpm is None:
        noise = np.random.randn(len(rt))*noisePpm*1e-6
    else:
        # -- no noise
        noise = np.zeros(len(rt))
        noisePpm = 1.0
    ax3.set_title(r'Transit: b=%4.2f, $R_p/R_\star$=%3.1f%%'%(p['b'], 100*p['rp/r*']))
    ax3.plot( rt, -1000*2.5*np.log10(f), '-k', linewidth=3, alpha=0.4, label='SATLAS')
    ax3.plot(-rt, -1000*2.5*np.log10(f), '-k', linewidth=3, alpha=0.4)

    # -- visibility:
    x = np.linspace(0,7,100)
    V2_satlas = np.trapz(special.jv(0,x[None,:]*r[:,None])*
                 (Ir[:,None]*r[:,None]), r[:,None], axis=0)**2
    V2_satlas /= np.trapz(special.jv(0,0*r)*(Ir*r), r)**2
    if showV2:
        ax5.plot(x, V2_satlas, '-k', linewidth=3, alpha=0.4, label='SATLAS')

    for k in list(laws.keys()):
        print('='*12, k, '='*12)
        if True: #'sigmo' in k: # -- fit to I(mu) profile
            fit = dpfit.leastsqFit(Imu_an, data['mu'], laws[k],
                                data['Kep'], verbose=True)
            laws[k] = fit['best']
            dpfit.dispCor(fit)
            print('fit to Imu:', fit['best'])
            ax1.plot(data['mu'], fit['model'], '.', color=colors[k],
                        linestyle='', alpha=0.5)

        r = np.sqrt(1-data['mu']**2)[::-1]
        Ir = Imu_an(data['mu'][::-1], laws[k])
        r0 = maxVar(Ir,r)
        r /= r0
        ax1.plot(data['mu'], Ir[::-1], '-', label=k, color=colors[k])
        ax2.plot(data['mu'], Ir[::-1] - data['Kep'] , '-', label=k, color=colors[k])

        # -- transit curve
        tmp = p.copy()
        tmp.update(laws[k])
        f_ = func_p(rt, tmp)
        ax3.plot( rt, -1000*2.5*np.log10(f_), '-', label=k, color=colors[k])
        ax3.plot(-rt, -1000*2.5*np.log10(f_), '-', color=colors[k])
        res = (f_-f)*1e6
        ax4.plot( rt, res, '-', color=colors[k])
        ax4.plot(-rt, res, '-', color=colors[k])
        rstar = p['r*']
        V2_ = np.trapz(special.jv(0,rstar*x[None,:]*r[:,None])*
                         (Ir[:,None]*r[:,None]), r[:,None], axis=0)**2
        V2_ /= np.trapz(special.jv(0,0*r)*(Ir*r), r)**2
        if showV2:
            ax5.plot(x, V2_, '-', label=k, color=colors[k])
            ax6.plot(x, 1000*(V2_-V2_satlas), '-', label=k, color=colors[k])

        # -- fitting in the transit space:
        # --------------------------------
        doNotFit=['mu0'] # mu0 and r* are perfectly correlated...
        fit = dpfit.leastsqFit(func_p, rt, tmp, f+noise, noisePpm*np.ones(len(rt))*1e-6,
                               verbose=True, doNotFit=doNotFit, ftol=1e-3, maxfev=1000)
        dpfit.dispCor(fit)
        print('-->', k)
        for c in ['r*', 'rp/r*', 'b']:
            if fit['uncer'][c]>0:
                print('  %6s : %3.1f sigmas'%(c, (fit['best'][c] - p[c])/fit['uncer'][c]))
        f_tf = fit['model']
        Ir_tf = Imu_an(data['mu'], fit['best'])[::-1]
        #ax1.plot(data['mu'], Ir_tf[::-1], '--', label=k+' fit transit', color=colors[k])
        #ax2.plot(data['mu'], Ir_tf[::-1]-data['Kep'],'--', label=k, color=colors[k])
        #ax3.plot( rt, -1000*2.5*np.log10(f_tf), '--', label=k, color=colors[k])
        #ax3.plot(-rt, -1000*2.5*np.log10(f_tf), '--', color=colors[k])
        res = (f_tf-f)*1e6
        #ax4.plot( rt, res, '--', label=k+' fitted, r*=%5.4f'%fit['best']['r*'],color=colors[k])
        #ax4.plot(-rt, res, '--', color=colors[k])

        rstar = fit['best']['r*']
        V2_tf = np.trapz(special.jv(0,rstar*x[None,:]*r[:,None])*
                         (Ir_tf[:,None]*r[:,None]), r[:,None], axis=0)**2
        V2_tf /= np.trapz(special.jv(0,0*r)*(Ir_tf*r), r)**2
        #if showV2:
            #ax5.plot(x, V2_tf, '--', label=k, color=colors[k])
            #ax6.plot(x, 1000*(V2_tf-V2_satlas), '--', label=k, color=colors[k])

        # --

    # -- finalize plots:
    ax4.hlines([-5, 5], -rt.max(), rt.max(), linestyle='dotted', label='+- 5 ppm')

    ax1.set_ylabel('Kepler: I / Imax')
    ax2.set_ylabel('O - C')
    ax2.set_xlabel(r'$\mu$')
    ax4.set_ylabel('O - C (ppm)')
    ax3.set_ylabel('Kepler transit (mmag)')
    if showV2:
        ax5.set_ylabel('V2')
        ax6.set_ylabel('O - C (1e-3)')
        ax6.set_xlabel(r'$\pi B \theta/\lambda$' )

    ax1.legend(loc='lower right', prop={'size':10})
    ax1.set_ylim(0.005, 1.05)

    ax1.set_ylim(0.0001, 1.5)
    ax1.semilogx()
    #ax1.semilogy()
    ax2.semilogx()

    ax4.set_xlabel(r'd, with R$_\star$=%5.3f'%p['r*'])
    ax4.legend(loc='lower left', prop={'size':10})
    ax4.set_ylim(-np.abs(ax4.get_ylim()).max(),
                 np.abs(ax4.get_ylim()).max())
    ax2.set_ylim(-np.abs(ax2.get_ylim()).max(),
                 np.abs(ax2.get_ylim()).max())
    ax3.set_ylim(ax3.get_ylim()[1], ax3.get_ylim()[0])
    if showV2:
        ax5.semilogy()
        ax5.set_ylim(0.001,1)
        ax6.set_ylim(-np.abs(ax6.get_ylim()).max(),
                     np.abs(ax6.get_ylim()).max())
        ax6.hlines(0.0, x.min(), x.max(), linestyle='dotted')
    return

def maxVar(I,r):
    """
    maximum variation
    """
    g = np.abs(np.gradient(I)/np.gradient(r))
    i_0 = np.argmax(g)
    i_0 = min(len(I)-2, i_0)
    i_0 = max(1, i_0)

    i = np.array([i_0-1,i_0,i_0+1])
    try:
        c = np.polyfit(r[i]-1, g[i], 2)
    except:
        print(I, r)
        print(r[i], g[i])
    res = -c[1]/(2*c[0])+1
    return res

def Imu_an(mu,p):
    if 'mu0' in list(p.keys()):
        # -- sigmoid using 3 parameters
        res = (1+np.exp(-(1-p['mu0'])/np.abs(p['c'])))/(1+np.exp(-(mu-p['mu0'])/np.abs(p['c'])))
        if 'e' in list(p.keys()):
            res = res**p['e']
        mup = (mu-p['mu0'])/(1-p['mu0'])
    else:
        res = np.ones(len(mu))
        mup = mu
    c = 1.
    for k in list(p.keys()):
        # -- f1...f4 gives 4 parameters law, a la Claret
        # -- f2 only for linear CLD
        if k.startswith('f'):
            c -= p[k]*(1-mu**(float(k[1:])/2.))

    if 'ap' in list(p.keys()) and 'p' in list(p.keys()):
        c -= p['ap']*(1-np.sign(mup)*np.abs(mup)**p['p'])
    res *= c
    return res

def func_p(d,p):
    global _Ir, _r, _r2, _tot
    mu = np.linspace(1,0,1000)
    _r = np.sqrt(1-mu**2)
    _Ir = Imu_an(mu, p)
    _Ir = np.maximum(_Ir, 0.0)

    if 'mu0' in list(p.keys()):
        r0 = maxVar(_Ir,_r)
        _r /= r0
    _r2 = _r**2
    _tot = np.trapz(2*np.pi*_r*_Ir, _r)
    return func(d,p)

def func(d, p, Ir=None):
    """
    impact parameter and planet / star radii ratio
    p = {'b': impact parameter,
         'rp/r*': planet relait radius,
         'r*': stellar radius}

    d is the position in the transit (d==0, closest to star's center),
    in unit of star's radius

    Ir should be initialized first, with r going from 0 to rmax (could be >1)
    """
    global _Ir, _r, _r2, _tot
    if not Ir is None:
        # -- setup
        _Ir = Ir
        _r = d
        _r2 = d**2
        # -- total flux
        _tot = np.trapz(2*np.pi*_r*_Ir, _r)
        return
    _x = _disk['x'][:,None]*p['rp/r*'] + d[None,:]/p['r*']
    _y = _disk['y'][:,None]*p['rp/r*'] + p['b']
    r2 = (_x**2+_y**2)
    # -- result
    I = np.interp(r2, _r2, _Ir, left=_Ir[0], right=0.0)
    return (_tot - np.sum(I, axis=0)*_disk['ds']*p['rp/r*']**2)/_tot
