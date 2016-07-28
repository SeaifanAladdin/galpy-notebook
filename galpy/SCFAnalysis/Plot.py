import numpy as nu
import matplotlib.pyplot as plt
from seaborn import heatmap

from Constants import *
from func import *

def __heatmap(coeffs, scale=1, scalePower=1.):
    heatmap(scale*nu.fabs(coeffs)**scalePower)

def plotAxiCoeffs(Acos,scalePower=1.):
    plt.figure()
    __heatmap(Acos[:,:,0],1.,scalePower)
    plt.xlabel("l")
    plt.ylabel("n")

def plotAxiCoeffsVarying_n(Acos, n_start=0, n_end=5, l_step=2,log=False):
    r = nu.arange(0, Acos.shape[1], l_step)
    for n in range(n_start, n_end, 1):
        Acos_line = Acos[n,::l_step,0]
        plt.plot(r, Acos_line, "-", label="{0}".format(n))
    plt.xlabel("l")
    plt.ylabel(r"$Acos_{nl}$")
    plt.xticks(r)
    plt.legend()

    if log:
        plt.xscale("symlog")
        plt.yscale("symlog")
    
def plotOrbitPos_RA_DEC(orbits_pos, fit=True, fig=None):
    ra = nu.array(orbits_pos[:,7], dtype=float)
    dec = nu.array(orbits_pos[:,8], dtype=float)
    if fig==None:
        fig = plt.figure()
    if fit:
        ax = fig.add_subplot(211)
    else:
        ax=fig.add_subplot(111)
    ax.set_xlabel("RA")
    ax.set_ylabel("DEC")
    ax.plot(ra[:-1], dec[:-1], '.', label="Scatter")
    ax.plot(ra[-1], dec[-1], 'o', label="Orginal")
    if fit:
        plotRA_DEC_fit(orbits_pos)
    ax.legend(loc=4)
    if fit:
        ax = fig.add_subplot(212)
        #plt.subplot(1,2,2)
        plotRA_DEC_fit_scatter(ra, dec, ax)
    fig.tight_layout()
    

def plotRA_DEC(directory, std, dt,fit=True):
    orbits_pos = nu.load(TIMEORBITS.format(directory, TIME,std,dt))
    fig = plt.figure()
    fig.suptitle(r"$\sigma = {0}$".format(std))
    plotOrbitPos_RA_DEC(orbits_pos, fit,fig)
    

def plotRA_DEC_fit(orbits_pos):
    ra = nu.array(orbits_pos[:,7], dtype=float)
    dec = nu.array(orbits_pos[:,8], dtype=float)
    
    ra_array, dec_array = ra_decFit(ra,dec, min(ra), max(ra))
    plt.plot(ra_array, dec_array, label="Polynomial fit")

def plotRA_DEC_fit_scatter(ra, dec, ax=None):
    poly = nu.polyfit(ra,dec,3)
    p = nu.poly1d(poly)
    new_dec = p(ra)
    if ax==None:
        ax = fig.add_subplot(111)
    ax.set_xlabel("RA")
    ax.set_ylabel(r'$DEC_{raw} - DEC_{fit}$')
    ax.plot(ra, dec-new_dec, ".")
    ax.set_title('Mean Squares Error = {0}'.format(meanSquares(dec, new_dec)), fontsize=9)
    
def animateOrbit(orbits_pos, view_radius=50):
    fig = plt.figure()
    ax = fig.add_subplot("111")
    line, = ax.plot(orbits_pos[:-1,7,0], orbits_pos[:-1,8,0], ".")
    center, = ax.plot(orbits_pos[-1,7,0], orbits_pos[-1,8,0], "+w")
    ax.set_xlabel("RA")
    ax.set_ylabel("DEC")


    plt.ion()
    plt.show()
    for i in range(1, orbits_pos.shape[2], 5):
        line.set_xdata(orbits_pos[:-1,7,i])
        line.set_ydata(orbits_pos[:-1,8,i])

        c_x = orbits_pos[-1,7,i]
        c_y = orbits_pos[-1,8,i]
        center.set_xdata(c_x)
        center.set_ydata(c_y)
        
        ax.set_xbound(c_x - view_radius, c_x + view_radius)
        ax.set_ybound(c_y - view_radius, c_y + view_radius)
        ax.set_title("Time {0} Myr".format(i))
        
        
        plt.pause(.01)
    plt.ioff()
