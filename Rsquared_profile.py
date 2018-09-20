
# coding: utf-8

# We use a density distribution that drops off as $r^{-2}$, following a power law profile. This relationship is more common in density structure created by a Type II presupernova mass-loss.

# In[11]:


import numpy as np
import matplotlib.pyplot as plt
import astropy.units as unt
import scipy.integrate as integrate
from matplotlib import rc, rcParams


# In[12]:


# Setting up the radius from the center of the explosion.

pos_x = np.linspace(2.2546e-5, 1.5, 500)*unt.pc
pos_x = pos_x.to(unt.cm).value
pos_y = np.linspace(2.2546e-5, 1.5, 500)*unt.pc
pos_y = pos_y.to(unt.cm).value
pos_z = np.linspace(2.2546e-5, 1.5, 500)*unt.pc
pos_z = pos_z.to(unt.cm).value

def radius(a, b, c):
    r = np.sqrt(a**2 + b**2 + c**2)
    return r

rad = radius(pos_x, pos_y, pos_z) # cm


# density of the circumstellar medium is given by
# \begin{equation}
# \rho_{cs} = \frac{\dot{M}}{4\pi u_w r^2}
# \end{equation}

# In[22]:


dM = 6e-5 #Msol/yr
dM = 5.6658e21 #g/s

u_w = 10*(unt.km/unt.s) #km/s
u_w = u_w.to(unt.cm/unt.s).value #cm/s

den = dM/(4*np.pi*u_w*rad**2) #g/cm^3


fig, ax = plt.subplots(1, figsize=(7, 5))
ax.plot(r, den)
ax.set_xlabel(r'r [cm]')
ax.set_ylabel(r'$\rho_{cs} \, [\frac{g}{cm^3}]$')
ax.ticklabel_format(style='sci', scilimits=(0.0, 0.0), axis='y')


# All the observation dependent parameters

beta = -0.81
alpha = -0.73
delta = -1.88
m = -delta/3

gamma = 2


time = np.linspace(30, 3500, 500)*unt.day #days
time = time.to(unt.s).value

#freq = [0.3, 0.55, 1.4, 4.9, 8.5, 15, 22, 90] #GHz
f = 30.3*unt.GHz
f = f.to(unt.Hz).value

eta = (gamma - 7 + 12*m) / 4
#b = -(gamma - 5 + 6*m) / 2
#a = (1 - gamma) / 2


# In[37]:


# Calculating the optical depth from the emission of the shocked material.
def int_optical_depth(t, mu):
    
    OPTICAL_DEPTH = (mu**(-2.1)) * (t**delta) * (den)**(5+delta)

    return(OPTICAL_DEPTH)


# Calculating the radio flux density from the emission of the shocked material.
def int_flux_density(t, mu, tau):

    FLUX_DENSITY = (mu**alpha) * (t**beta) * (den)**eta * np.exp(-tau) * (2*np.pi*(1-np.cos(0.00203622)))
    
    return(FLUX_DENSITY)



# The synchrotron spectrum of a single electron has a logarithmic slope of 1/3 at low frequencies, a broad peak near the critical frequency $\nu$ max, and falls off sharply at higher frequencies.

# \begin{equation}
# 1 \mbox{Jy} = 10^{-23} \frac{\mbox{erg}}{\mbox{s} \, \mbox{cm}^2 \mbox{Hz}}
# \end{equation}
# 

# \begin{equation}
# 1 \mbox{day} = 86400 \, \mbox{sec}
# \end{equation}
# 

# In[38]:


opt_sum = int_optical_depth(time, f)
flux_sum = int_flux_density(time, f, opt_sum)
# intergrating over the entire solid angle and time
S_v = integrate.cumtrapz(flux_sum, time, initial=0)



fig, ax = plt.subplots(1, figsize=(15, 5))
ax.plot(time/86400, S_v*(1/10**(-23)), 'o')
#ax.legend(freq, loc=1)
ax.set_xlabel("days")
ax.set_ylabel(r"$S(\nu) \, [mJy]$")
plt.savefig("light_curve.png")





fig, ax = plt.subplots(1, figsize=(15, 5))

for f in freq:
    opt_sum = int_optical_depth(time, f)
    tau_v = integrate.cumtrapz(opt_sum, time, initial=0)
    
    flux_sum = int_flux_density(time, f, tau_v)
    S_v = integrate.cumtrapz(flux_sum, time, initial=0)
    
    ax.plot(time, S_v/10000)
    ax.legend(freq, loc=1)
ax.set_xlabel("days")
ax.set_ylabel(r"$S(\nu) [mJy]$")
plt.savefig("light_curve.png")

