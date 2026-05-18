#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import sys
import os
import numpy as np
import pandas as pd
import math
import time
from tools.tools import load_config
from qcdlib.aux import AUX
from tools.config import conf
from scipy.integrate import quad, dblquad, fixed_quad
from numba import jit
from qcdlib import core,pdf0,ff0,pdf1,ff1


#AN_theory0_hadinjet.py - program to calculate A_UT hadron-in-jet in pp -> hX

#flavor = ['g','u','ub','d','db','s','sb']
         #  0   1   2    3    4   5    6

# Common color factors and fractions
# c = {'r3': 1. / 3., 'r4': 0.25, 'r6': 1. / 6., 'r8': 0.125,
#         'r9': 1. / 9., 'r18': 1. / 18., 'r24': 1. / 24., 'r27': 1. / 27.}

try: zdepCw=conf['zdep Collins widths']
except: zdepCw=False

if 'basis' not in conf:
    conf['basis'] = 'default'

def get_f(x, Q2): # Collinear unpolarized PDF
    return conf['pdf'].get_C(x, Q2)

def get_ft(x, Q2): # Collinear unpolarized PDF
    return conf['pdf'].get_C(x, Q2)

def get_d(z, jperp2, Q2, had): # D1(z,jperp^2) in Gaussian ansatz
    
    if 'pi' in had:
        w_had = np.abs(conf['ffpi'].get_widths(Q2))
        
        if jperp2==0: return conf['ffpi'].get_C(z, Q2) #This is for jperp-integrated observables
        else: return conf['ffpi'].get_C(z, Q2) * np.exp(-jperp2/w_had)/(np.pi*w_had)
    
    elif 'k' in had:
        w_had = np.abs(conf['ffk'].get_widths(Q2))
        
        if jperp2==0: return conf['ffk'].get_C(z, Q2) #This is for jperp-integrated observables
        else: return conf['ffk'].get_C(z, Q2) * np.exp(-jperp2/w_had)/(np.pi*w_had)

def get_h(x, Q2): # Collinear transversity
    return conf['transversity'].get_C(x, Q2)

def get_H1p(z, jperp2, Q2, had): #H_1^{\perp}(z,jperp^2) in Gaussian ansatz
    jperp = np.sqrt(jperp2)
    
    if 'pi' in had:
        Mh=0.134
        if zdepCw==True: w_had=np.abs(conf['collinspi'].get_widths(Q2,z))
        elif zdepCw==False: w_had=np.abs(conf['collinspi'].get_widths(Q2))
        
        if jperp2==0: return z * Mh * np.sqrt(np.pi)/np.sqrt(w_had) * conf['collinspi'].get_C(z, Q2) #This is for jperp-integrated observables
        else: return (jperp/z/Mh) * 2.* z**2 * Mh**2 * conf['collinspi'].get_C(z, Q2) * np.exp(-jperp2/w_had)/(np.pi*w_had**2)
    
    elif 'k' in had:
        Mh = 0.497
        if zdepCw==True: w_had=np.abs(conf['collinsk'].get_widths(Q2,z))
        elif zdepCw==False: w_had=np.abs(conf['collinsk'].get_widths(Q2))
        
        if jperp2==0: return z * Mh * np.sqrt(np.pi)/np.sqrt(w_had) * conf['collinsk'].get_C(z, Q2) #This is for jperp-integrated observables
        else: return (jperp/z/Mh) * 2.* z**2 * Mh**2 * conf['collinsk'].get_C(z, Q2) * np.exp(-jperp2/w_had)/(np.pi*w_had**2)


@jit(nopython=True)
def get_Hupol(s,t,u):
  # Hard parts for the unpolarized cross section
    Hupol = np.zeros(15)
    Hupol[1] = 4. * (s**2. / t**2. + u**2. / t**2.) * 1/9
    Hupol[2] = 4. * (s**2. / u**2. + t**2. / u**2.) * 1/9
    Hupol[3] = -8. * s/t * s/u * 1/27
    Hupol[4] = 4. * (t**2. / s**2. + u**2. / s**2.) * 1/9
    Hupol[5] = 4. * (s**2. / t**2. + u**2. / t**2.) * 1/9
    Hupol[6] = -8. * u/t * u/s * 1/27
    Hupol[7] = 4. * (t**2. / s**2. + u**2. / s**2.) * 1/9
    Hupol[8] = 4. * (s**2. / u**2. + t**2. / u**2.) * 1/9
    Hupol[9] = -8. * t/u * t/s * 1/27
    Hupol[10] = 32. * 1/27 * (t/u + u/t) * (1. - 9. * t/s * u/s * 1/4)
    Hupol[11] = 4. * 1/9 * (- s/u - u/s) * (1. - 9. * s/t * u/t * 1/4)
    Hupol[12] = 4. * 1/9 * (- s/t - t/s) * (1. - 9. * s/u * t/u * 1/4)
    Hupol[13] = (t/u + u/t) * 1/6 - 3. * (t**2. / s**2. + u**2. / s**2.) * 1/8
    Hupol[14] = 4.5 * (3. - t/s * u/s - s/t * u/t - s/u * t/u)

    return Hupol

@jit(nopython=True)
def get_HTffa(s, t, u):
  # Hard parts for the transversely polarized fragmentation term
    HTffa = np.zeros(13)
    HTffa[0] = 0
    HTffa[1] = 8./9. - 2*s*u/t**2 #qg --> qg
    HTffa[2] = 8./9. * (-s*u/t**2 + 1./3. * s/t) #qq --> qq
    HTffa[3] = 8./9. * (-s*u/t**2 + 1./3. * u/t) #qqbar --> qqbar
    HTffa[4] = -8./27. #qbar q  --> qqbar
    HTffa[5] = -8./9. * s * u / t**2 #qq' --> qq'
    HTffa[6] = -8./9. * s * u / t**2 #qqbar' --> qqbar'
    HTffa[7] = HTffa[1]
    HTffa[8] = HTffa[2]
    HTffa[9] = HTffa[3]
    HTffa[10] = HTffa[4]
    HTffa[11] = HTffa[5]
    HTffa[12] = HTffa[6]
    
    return HTffa


@jit(nopython=True)
def get_Hxxpz(z, Q2, had, s, t, u, HTffa, H1p):
    if had.endswith('-'):
        H1p = core.charge_conj(H1p)

    elif had.endswith('0'):
        H1pp=H1p
        H1pm=core.charge_conj(H1p)

        H1p=0.5*(H1pp+H1pm)

    return HTffa,H1p


#  @profile
# Calculation of the unpolarized cross section
@jit(nopython=True)
def get_var(x, z, jperp, pT, eta, rs):

    if pT > 1.:
        Q = pT
    else:
        Q = 1.

    Q2 = Q * Q

    xT = 2. * pT / rs
    xT2 = xT * xT
    xF = xT*np.sinh(eta)
    xF2 = xF * xF

  # Mandelstam variables at the hadron level
    ss = rs * rs
    tt = -0.5 * ss * (np.sqrt(xF2 + xT2) - xF)
    uu = -0.5 * ss * (np.sqrt(xF2 + xT2) + xF)

    oz = 1. / z

    xp = -x * tt / (x * ss + uu)

  # Mandelstam variables at the parton level
    s = x * xp * ss
    t = x * tt
    u = xp * uu

  # Prefactor
    denfac = 1. / ((x * ss + uu) * x * xp)
    numfac = 1. / ((x * ss + uu) * x * xp)
    
    return denfac,numfac,s,t,u,Q2,xp

@jit(nopython=True)
def _get_dsig(x, z, jperp, pT, eta, rs, tar, had, f, ft, d): 
  # Get arrays of the nonperturbative functions
  #print x,xp
    
    var=get_var(x, z, jperp, pT, eta, rs)
    denfac=var[0]
    s,t,u=var[2],var[3],var[4]
    Q2,xp=var[5],var[6]
    
    Hupol=get_Hupol(s,t,u)

    if had.endswith('-'):
        d = core.charge_conj(d)

    elif had.endswith('0'):
        dp=d
        dm=core.charge_conj(d)
        d=0.5*(dp+dm)

    fg = f[0]
    fu = f[1]
    fub = f[2]
    fd = f[3]
    fdb = f[4]
    fs = f[5]
    fsb = f[6]

    ftg = ft[0]
    ftu = ft[1]
    ftub = ft[2]
    ftd = ft[3]
    ftdb = ft[4]
    fts = ft[5]
    ftsb = ft[6]

    dg = d[0]
    du = d[1]
    dub = d[2]
    dd = d[3]
    ddb = d[4]
    ds = d[5]
    dsb = d[6]

    upol = 0

    upol += (fu * du + fd * dd + fs * ds) * (ftu + ftd + fts) * Hupol[1] \
    + (fu + fd + fs) * (ftu * du + ftd * dd + fts * ds) * Hupol[2] \
    + (fu * ftu * du + fd * ftd * dd + fs * fts * ds) * Hupol[3] 

    upol += (fub * dub + fdb * ddb + fsb * dsb) * (ftub + ftdb + ftsb) * Hupol[1] \
    + (fub + fdb + fsb) * (ftub * dub + ftdb * ddb + ftsb * dsb) * Hupol[2] \
    + (fub * ftub * dub + fdb * ftdb * ddb + fsb * ftsb * dsb) * Hupol[3]  

    upol += (fu * du + fd * dd + fs * ds) * (ftub + ftdb + ftsb) * Hupol[5] \
    + (fu * ftub + fd * ftdb + fs * ftsb) * (du + dd + ds) * Hupol[4] \
    + (fu * ftub * du + fd * ftdb * dd + fs * ftsb * ds) * Hupol[6] 

    upol += (fub * dub + fdb * ddb + fsb * dsb) * (ftu + ftd + fts) * Hupol[5] \
    + (fub * ftu + fdb * ftd + fsb * fts) * (dub + ddb + dsb) * Hupol[4] \
    + (fub * ftu * dub + fdb * ftd * ddb + fsb * fts * dsb) * Hupol[6] 

    upol += (fu + fd + fs) * (ftub * dub + ftdb * ddb + ftsb * dsb) * Hupol[8] \
    + (fu * ftub + fd * ftdb + fs * ftsb) * (dub + ddb + dsb) * Hupol[7] \
    + (fu * ftub * dub + fd * ftdb * ddb + fs * ftsb * dsb) * Hupol[9] 

    upol += (fub + fdb + fsb) * (ftu * du + ftd * dd + fts * ds) * Hupol[8] \
    + (fub * ftu + fdb * ftd + fsb * fts) * (du + dd + ds) * Hupol[7] \
    + (fub * ftu * du + fdb * ftd * dd + fsb * fts * ds) * Hupol[9] 

    upol += (fu * ftub + fd * ftdb + fs * ftsb) * dg * Hupol[10]

    upol += (fub * ftu + fdb * ftd + fsb * fts) * dg * Hupol[10]

    upol += (fu * du + fd * dd + fs * ds) * ftg * Hupol[11]

    upol += (fub * dub + fdb * ddb + fsb * dsb) * ftg * Hupol[11]

    upol += (fu + fd + fs) * ftg * dg * Hupol[12]

    upol += (fub + fdb + fsb) * ftg * dg * Hupol[12]

    upol += fg * ftg * (du + dd + ds) * Hupol[13]

    upol += fg * ftg * (dub + ddb + dsb) * Hupol[13]

    upol += fg * ftg * dg * Hupol[14]

    upol += fg * (ftu * du + ftd * dd + fts * ds) * Hupol[12]

    upol += fg * (ftub * dub + ftdb * ddb + ftsb * dsb) * Hupol[12]

    upol += fg * (ftu + ftd + fts) * dg * Hupol[11]

    upol += fg * (ftub + ftdb + ftsb) * dg * Hupol[11]

    return denfac * upol

#  @profile
# Calculation of the fragmentation term in the transversely polarized cross section
@jit(nopython=True)
def _get_dsigST(x, z, jperp, pT, eta, rs, tar, had, ft, h, Hxxpz):  
    
    var=get_var(x, z, jperp, pT, eta, rs)
    numfac=var[1]
    s,t,u=var[2],var[3],var[4]
    Q2,xp=var[5],var[6]
    
    if had.startswith('pi'): Mh=0.134
    elif had.startswith('k'): Mh = 0.497

    hg = h[0]
    hu = h[1]
    hub = h[2]
    hd = h[3]
    hdb = h[4]
    hs = h[5]
    hsb = h[6]

    ftg = ft[0]
    ftu = ft[1]
    ftub = ft[2]
    ftd = ft[3]
    ftdb = ft[4]
    fts = ft[5]
    ftsb = ft[6]


    Hxxpz1u = Hxxpz[1][1]
    Hxxpz1d = Hxxpz[1][3]
    Hxxpz1s = Hxxpz[1][5]
    Hxxpz2u = Hxxpz[2][1]
    Hxxpz2d = Hxxpz[2][3]
    Hxxpz2s = Hxxpz[2][5]
    Hxxpz3u = Hxxpz[3][1]
    Hxxpz3d = Hxxpz[3][3]
    Hxxpz3s = Hxxpz[3][5]
    Hxxpz4u = Hxxpz[4][1]
    Hxxpz4d = Hxxpz[4][3]
    Hxxpz4s = Hxxpz[4][5]
    Hxxpz5u = Hxxpz[5][1]
    Hxxpz5d = Hxxpz[5][3]
    Hxxpz5s = Hxxpz[5][5]
    Hxxpz6u = Hxxpz[6][1]
    Hxxpz6d = Hxxpz[6][3]
    Hxxpz6s = Hxxpz[6][5]

    Hxxpz7ub = Hxxpz[7][2]
    Hxxpz7db = Hxxpz[7][4]
    Hxxpz7sb = Hxxpz[7][6]
    Hxxpz8ub = Hxxpz[8][2]
    Hxxpz8db = Hxxpz[8][4]
    Hxxpz8sb = Hxxpz[8][6]
    Hxxpz9ub = Hxxpz[9][2]
    Hxxpz9db = Hxxpz[9][4]
    Hxxpz9sb = Hxxpz[9][6]
    Hxxpz10ub = Hxxpz[10][2]
    Hxxpz10db = Hxxpz[10][4]
    Hxxpz10sb = Hxxpz[10][6]
    Hxxpz11ub = Hxxpz[11][2]
    Hxxpz11db = Hxxpz[11][4]
    Hxxpz11sb = Hxxpz[11][6]
    Hxxpz12ub = Hxxpz[12][2]
    Hxxpz12db = Hxxpz[12][4]
    Hxxpz12sb = Hxxpz[12][6]

    ffcs = 0

    ffcs += ftg * (hu * Hxxpz1u + hd * Hxxpz1d + hs * Hxxpz1s)

    ffcs += (hu * ftu * Hxxpz2u) + (hd * ftd * Hxxpz2d) + (hs * fts * Hxxpz2s)

    ffcs += (hu * ftub * Hxxpz3u) + (hd * ftdb * Hxxpz3d) + (hs * ftsb * Hxxpz3s)

    ffcs += (hub * ftu * Hxxpz4u) + (hdb * ftd * Hxxpz4d) + (hsb * fts * Hxxpz4s)

    ffcs += ftu * (hd * Hxxpz5d + hs * Hxxpz5s) + ftd * (hu * Hxxpz5u + hs * Hxxpz5s) + fts * (hu * Hxxpz5u + hd * Hxxpz5d)

    ffcs += ftub * (hd * Hxxpz6d + hs * Hxxpz6s) + ftdb * (hu * Hxxpz6u + hs * Hxxpz6s) + ftsb * (hu * Hxxpz6u + hd * Hxxpz6d)

    ffcs += ftg * (hub * Hxxpz7ub + hdb * Hxxpz7db + hsb * Hxxpz7sb)

    ffcs += (hub * ftub * Hxxpz8ub) + (hdb * ftdb * Hxxpz8db) + (hsb * ftsb * Hxxpz8sb)

    ffcs += (hub * ftu * Hxxpz9ub) + (hdb * ftd * Hxxpz9db) + (hsb * fts * Hxxpz9sb)

    ffcs += (hu * ftub * Hxxpz10ub) + (hd * ftdb * Hxxpz10db) + (hs * ftsb * Hxxpz10sb)

    ffcs += ftub * (hdb * Hxxpz11db + hsb * Hxxpz11sb) + ftdb * (hub * Hxxpz11ub + hsb * Hxxpz11sb) + ftsb * (hub * Hxxpz11ub + hdb * Hxxpz11db)

    ffcs += ftu * (hdb * Hxxpz12db + hsb * Hxxpz12sb) + ftd * (hub * Hxxpz12ub + hsb * Hxxpz12sb) + fts * (hub * Hxxpz12ub + hdb * Hxxpz12db)


    return numfac * ffcs
    
    
def get_dsig(x, z, jperp, pT, eta, rs, tar, had):

    var=get_var(x, z, jperp, pT, eta, rs)
    s,t,u=var[2],var[3],var[4]
    Q2,xp=var[5],var[6]

  # Get arrays of the nonperturbative functions
    f=np.zeros(11)
    ft=np.zeros(11)
    d=np.zeros(11)
    
    jperp2 = jperp**2
    
    f = get_f(x, Q2)
    ft = get_ft(xp, Q2)
    d = get_d(z, jperp2, Q2, 'pi+')
    
    return _get_dsig(x, z, jperp, pT, eta, rs, tar, had, f, ft, d)

def get_dsigST(x, z, jperp, pT, eta, rs, tar, had):

    var=get_var(x, z, jperp, pT, eta, rs)
    s,t,u=var[2],var[3],var[4]
    Q2,xp=var[5],var[6]

  # Get arrays of the nonperturbative functions
    ft=np.zeros(11)
    h=np.zeros(11)
    Hxxpz = np.zeros((13, 7))
    d=np.zeros(11)
    
    HTffa = get_HTffa(s, t, u)
    
    jperp2 = jperp**2
    
    H1p=get_H1p(z, jperp2, Q2, 'pi+')

    ft = get_ft(xp, Q2)
    h = get_h(x, Q2)
    d = get_d(z, jperp2, Q2, 'pi+')
    Hxxpz = get_Hxxpz(z, Q2, had, s, t, u, HTffa, H1p)
    
    Hxxpz = np.einsum('i,j->ij', Hxxpz[0], Hxxpz[1])
    
    return _get_dsigST(x, z, jperp, pT, eta, rs, tar, had, ft, h, Hxxpz)

def get_sig(z, jperp, pT, eta, rs, tar, had, mode='gauss', nx=10):
    xT = 2. * pT / rs
    xT2 = xT * xT
    xF = xT*np.sinh(eta)
    xF2 = xF * xF

    # Mandelstam variables at the hadron level
    ss = rs * rs
    tt = -0.5 * ss * (np.sqrt(xF2 + xT2) - xF)
    uu = -0.5 * ss * (np.sqrt(xF2 + xT2) + xF)

    # Lower limits of the x integration

    xmin = xT * np.exp(eta) / (2. - xT * np.exp(-eta))

    if mode == 'gauss':
        dsigdx = np.vectorize(lambda x: get_dsig(x, z, jperp, pT, eta, rs, tar, had))
        sig = fixed_quad(dsigdx, xmin, 1., n=nx)[0]
    elif mode == 'quad':
        sig = quad(lambda x: get_dsig(x, z, jperp, pT, eta, rs, tar, had), xmin, 1.)[0]
    return sig

def get_sigST(z, jperp, pT, eta, rs, tar, had, mode='gauss', nx=10):
    xT = 2. * pT / rs
    xT2 = xT * xT
    xF = xT*np.sinh(eta)
    xF2 = xF * xF

    # Mandelstam variables at the hadron level
    ss = rs * rs
    tt = -0.5 * ss * (np.sqrt(xF2 + xT2) - xF)
    uu = -0.5 * ss * (np.sqrt(xF2 + xT2) + xF)

    # Lower limits of the x integration
    xmin = xT * np.exp(eta) / (2. - xT * np.exp(-eta))

    if mode == 'gauss':
        dsigdx = np.vectorize(lambda x: get_dsigST(x, z, jperp, pT, eta, rs, tar, had))
        sig = fixed_quad(dsigdx, xmin, 1., n=nx)[0]
    elif mode == 'quad':
        sig = quad(lambda x: get_dsigST(x, z, jperp, pT, eta, rs, tar, had), xmin, 1.)[0]
    return sig

    return sig


if __name__ == '__main__':


    from qcdlib.ff0 import FF as FF0
    from qcdlib.ff1_test import FF as FF1
    from qcdlib.pdf0 import PDF as PDF0
    from qcdlib.pdf1_test import PDF as PDF1
    conf['aux']= AUX()
    conf['pdf']=PDF0()
    conf['collinspi']=FF1('pi')
    conf['collinsk']=FF1('k')
    conf['dcollinspi']=FF1('pi','deriv')
    conf['dcollinsk']=FF1('k','deriv')
    conf['Htildepi']=FF1('pi')
    conf['Htildek']=FF1('k')
    conf['transversity']=PDF1()
    conf['sivers']=PDF1()
    conf['dsivers']=PDF1('deriv')
    conf['ffpi']=FF0('pi')
    conf['ffk']=FF0('k')

    rs = 200.
    tar = 'p'
    #had = 'pi+'
    had='pi+'
    pT = 2.0
    z = 0.6
    jperp = 0.5
    eta = 3

    def test():
    
        den = get_sig(z, jperp, pT, eta, rs, tar, had, mode='gauss', nx=10)
        num = get_sigST(z, jperp, pT, eta,rs, tar, had, mode='gauss', nx=10)

        AN = num / den
        print(AN)

#  test()

# from timeit import Timer
# t = Timer("test()", "from __main__ import test")
# print 't elapsed ',t.timeit(number=1)

# def test2():
#  den = anthy.get_dsig(0.3,0.6,xF,pT,rs,tar,had)
#  num = anthy.get_dsigST(0.3,0.6,xF,pT,rs,tar,had)
#
#  print den,num
#
# from timeit import Timer
# t = Timer("test2()", "from __main__ import test2")
# print 't elapsed ',t.timeit(number=1)

# start = time.time()
# print anthy.get_dsig(0.3,0.6,xF,pT,rs,tar,had)
# print anthy.get_dsigST(0.3,0.6,xF,pT,rs,tar,had)
# end = time.time()
# print 'time=',(end-start)

    for k in range(100):
        start = time.time()
        test()
        end = time.time()
        print('time=', (end - start))

  # Integration of the numerator from xmin to 1 and from zmin to 1 (the values for xmin and zmin are above)

  # Integration of the denominator from xmin to 1 and from zmin to 1 (the values for xmin and zmin are above)
  #den = dblquad(lambda x,z: ANTHEORY().get_dsig(x,z,xF,pT,rs,tar,had),zmin,1.,xmin,lambda x: 1.)

  #AN = num[0]/den[0]
  # print(AN)
