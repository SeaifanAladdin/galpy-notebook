import numpy as nu

from galpy.util.bovy_coords import  cyl_to_rect
from galpy import potential
from galpy.orbit import Orbit
from Plot import *
from mpl_toolkits.mplot3d import axes3d, art3d
import matplotlib.pyplot as plt
from astropy import units

from galpy.util import bovy_conversion

from time import clock
folder = "SCF_n=1_l=1_Acos=-.66"


##Creating the adjusted nfw density
def newNFWcutoffDensity(R, z, phi=0, a=2.):
    x,y,z = cyl_to_rect(R,phi,z)
    xi = (x**2 + y**2 + .5*z**2)**.5
    return (4*nu.pi*(a**3)*(xi/a)*(1. + (xi/a))**2)**-1 * nu.e**(-(xi/100./a)**2)
def newNFWDensity(R, z, phi=0, a=2.):
    x,y,z = cyl_to_rect(R,phi,z)
    xi = (x**2 + y**2 + .5*z**2)**.5
    return (4*nu.pi*(a**3)*(xi/a)*(1. + (xi/a))**2)**-1 

a_nfw = 200./8
Acos, Asin = potential.scf_compute_coeffs_axi(newNFWcutoffDensity, 30,10, a=a_nfw, costheta_order=40)
Acos[1,1,0] = -.66
scf = potential.SCFPotential(Acos=Acos, Asin=Asin, a=a_nfw)

##plotAxiCoeffs(scf._Acos, .3)
##plt.figure()
##R = nu.linspace(.1,10,1000)
##plt.loglog(R, scf.dens(R,R,0), label="SCF")
##plt.loglog(R, newNFWDensity(R,R,0), label="new NFW")
##plt.loglog(R, newNFWcutoffDensity(R,R,0), label="new NFW with cutoff")
##plt.legend()
##plt.figure()
##r = nu.arange(0, Acos.shape[1], 2)
##for n in range(min(Acos.shape[1], 5)):
##    Acos_line = Acos[n,::2,0]
##    plt.plot(r, Acos_line, "--", label="{0}".format(n))
##plt.legend()    
##plt.show()

nu.save("Orbits/{0}/Acos.npy".format(folder), Acos)
    

TIME = 3000*units.Myr
dt = 1 * units.Myr
div = 1.
step_size = dt/div

##Creating the Milky Way Potential
ps= potential.PowerSphericalPotentialwCutoff(alpha=1.8,rc=1.9/8.,normalize=0.05)
mn= potential.MiyamotoNagaiPotential(a=3./8.,b=0.28/8.,normalize=.6)
MWPotential= [ps,mn,scf]

##Creating the orbit
o = Orbit([229.018,-0.124,23.2,-2.296,-2.257,-58.7],radec=True,ro=8.,vo=220.,solarmotion=[-11.1,24.,7.25])
o.turn_physical_off()

##Integrate backwards in time
o = o.flip()
ts= nu.arange(0,(TIME + step_size).value,step_size.value)*units.Myr

o.integrate(ts, MWPotential, method='dopr54_c')

##Integrating Forward in time
newOrbit = Orbit([o.R(TIME), -o.vR(TIME), -o.vT(TIME), o.z(TIME), -o.vz(TIME), o.phi(TIME)],ro=8.,vo=220.)
newOrbit.turn_physical_off()
newOrbit.integrate(ts, MWPotential, method='dopr54_c')


def randomVelocity(std=.001):
    if type(std).__name__ == "Quantity":
        return nu.random.normal(scale=std.value)*std.unit
    return nu.random.normal(scale=std)

time1 = nu.arange(0, TIME.value, dt.value)*units.Myr
orbits = nu.empty(len(time1), dtype=Orbit)
orbits_pos = nu.empty((len(time1) + 1,9,len(ts)), dtype=units.quantity.Quantity)
orbits_pos[0, :, :] = ts, newOrbit.x(ts), newOrbit.y(ts), newOrbit.z(ts), newOrbit.vx(ts), newOrbit.vy(ts), newOrbit.vz(ts), newOrbit.ra(ts), newOrbit.dec(ts)
orbits_pos[:,:,:] = orbits_pos[0,:,:]
i = 0
std = 0.008
stdR = std
stdT = std
stdz = std
for t in time1:
    print t
    dvR = randomVelocity(stdR)
    dvT = randomVelocity(stdT)
    dvz = randomVelocity(stdz)
    #dvR, dvT, dvz = 0,0,0
    tempOrbit = Orbit([newOrbit.R(t), newOrbit.vR(t) + dvR, newOrbit.vT(t) + dvT, newOrbit.z(t), newOrbit.vz(t) + dvz, newOrbit.phi(t)],ro=8.,vo=220.)
    tempOrbit.turn_physical_off()
    time = nu.arange(0,(TIME + step_size - t).value,step_size.value)*units.Myr
    tempOrbit.integrate(ts, MWPotential, method='dopr54_c')
    orbits[i] = tempOrbit
    #orbits_pos[:,i] = (orbits[i].x(time[-1]), orbits[i].y(time[-1]), orbits[i].z(time[-1]))
    orbits_pos[i, :, i*div:] = ts[ts>=t], tempOrbit.x(time), tempOrbit.y(time), tempOrbit.z(time), tempOrbit.vx(time), tempOrbit.vy(time), tempOrbit.vz(time), tempOrbit.ra(time), tempOrbit.dec(time)
    i +=1

##Saves the orbit at every step
nu.save("Orbits/{0}/Orbits_std{1}_dt{2}.npy".format(folder,std, dt),orbits_pos)
##Saves the Orbit at the very end
nu.save("Orbits/{0}/{1}/Orbits_std{2}_dt{3}.npy".format(folder, TIME, std, dt),orbits_pos[:,:,-1])    
##plt.plot(newOrbit.ra(TIME), newOrbit.dec(TIME), 'o')
##plt.plot(orbits_pos[:,7,-1], orbits_pos[:,8,-1], '.')

##orbits[0].plot(overplot=False)
##for o in orbits[1:]:
##    o.plot(overplot=True)


def plot3D():
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    c_x,c_y,c_z = newOrbit.x(TIME), newOrbit.y(TIME), newOrbit.z(TIME)
    ax.autoscale(enable=False,axis='both')  #you will need this line to change the Z-axis
    ax.set_xbound(c_x-2, c_x+2)
    ax.set_ybound(c_y-2, c_y+2)
    ax.set_zbound(c_z-2, c_z+2)
    
    
    i = 0
    for o in orbits:
        t = time1[i]
        time = ts[ts >= t]
        time2 = time - time[0]
        x,y,z = o.x(time2),o.y(time2), o.z(time2)
        ax.plot(x,y,z, "-")
        ax.plot(x[-1:], y[-1:], z[-1:], "o")
        i+=1
    ax.plot(nu.array([c_x]), nu.array([c_y]), nu.array([c_z]), 'Dk')


#plot3D()
del orbits_pos

