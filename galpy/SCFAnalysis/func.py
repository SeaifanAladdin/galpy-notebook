import numpy as nu

from Constants import *
from astropy import units
from galpy.util.bovy_coords import rect_to_cyl
from galpy import potential

def meanSquares(raw_dec, fit_dec):
    return (nu.sum((raw_dec - fit_dec)**2)/ len(raw_dec))**.5

def ra_decFit(ra, dec, minRa, maxRa, size=100):
    poly = nu.polyfit(ra,dec,3)
    ra_array = nu.linspace(minRa, maxRa, size)
    p = nu.poly1d(poly)
    dec_array = p(ra_array)
    return ra_array, dec_array


def getForces(folder):
    def _getForces(folder):
        loc = "SCFAnalysis/Orbits/{0}".format(folder)
        Acos = nu.load(loc + "/Acos.npy")
        a_nfw = 200./8
        orbits_pos =nu.load(TIMEORBITS.format(loc, TIME,0.004,1*units.Myr))
        scf = potential.SCFPotential(Acos=Acos, a=a_nfw, normalize=.35)
        ps= potential.PowerSphericalPotentialwCutoff(alpha=1.8,rc=1.9/8.,normalize=0.05)
        mn= potential.MiyamotoNagaiPotential(a=3./8.,b=0.28/8.,normalize=.6)
        MWPotential= [ps,mn,scf]
        x,y,z = orbits_pos[-1,1:4]
        R, phi, z = rect_to_cyl(x,y,z)
        print "(R,z,phi) = ({0}, {1}, {2})".format(R,z,phi)
        Rforce = potential.evaluateRforces(MWPotential, R,z,phi)
        zforce = potential.evaluatezforces(MWPotential, R,z,phi)
        phiforce = potential.evaluatephiforces(MWPotential, R,z,phi)
        return Rforce, zforce, phiforce
    if type(folder) is list:
        forces = nu.empty((len(folder), 3))
        for i in range(len(folder)):
            forces[i,:] = _getForces(folder[i])
        return forces
    else:
        return _getForces(folder)
    

    
