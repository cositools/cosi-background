#!/usr/bin/env python

""" Create a .dat file containing the spectrum of the primary
    and secondary protons plus albedo neutrons from the class
    LEOBackgroundGenerator to be used as input for the Step1
    *.source for Activation Simulations with MEGAlib.
"""

import numpy as np
import pandas as pd

from scipy.integrate import quad
import sys
import argparse

from LEOBackgroundGenerator import LEOBackgroundGenerator as LEO

# Instantiate the parser
pars = argparse.ArgumentParser(description='Create a .dat file containing '
                               + 'the spectrum of the primary and secondary '
                               + 'protons plus albedo neutrons from the class '
                               + 'LEOBackgroundGenerator to be used as input '
                               + 'for the Step1 *.source for Activation '
                               + 'Simulations with MEGAlib.')

pars.add_argument('-i', '--inclination', type=float, nargs='?',
                  default=0., help='Inclination of the orbit in degree [0.]')

pars.add_argument('-a', '--altitude', type=float, nargs='?',
                  default=550., help='Altitude of the orbit in km [550.]')

pars.add_argument('-el', '--elow', type=float, nargs='?',
                  default=2, help='Log10 of the lowest energy limit in keV [2]')
                
                                    

pars.add_argument('-eh', '--ehigh', type=float, nargs='?',
                  default=8, help='Log10 of the highest energy limit in keV [8]')


pars.add_argument('-c', '--cutoff', type=float, nargs='?',
                  default=None, help='Value of the geocutoff [Average cutoff]')
                  
pars.add_argument('-s', '--solarmodulation', type=float, nargs='?',
                  default=650., help='solar modulation (550 min and 1100 max) [650]')


pars.add_argument('-o', '--outputpath', type=str, nargs='?',
                  default="OUTPUT", help='output path')


pars.add_argument('-f','--components',type=str,nargs='?',default=None,help=
		'''list of the components you want, separated by a comma : 
		AtmosphericNeutrons\n
		CosmicPhotons\n
		PrimaryProtons\n
		SecondaryProtonsUpward\n
		SecondaryProtonsDownward\n
		PrimaryAlphas\n
		PrimaryElectrons\n
		PrimaryPositrons\n
		SecondaryElectrons\n
		SecondaryPositrons\n
		AlbedoPhotons\n
		
		default : all components''')


args = pars.parse_args()

Inclination = args.inclination
Altitude = args.altitude

Elow = args.elow
Ehigh = args.ehigh

Geocutoff = args.cutoff
outputpath = args.outputpath

solarmod = args.solarmodulation

components = args.components


LEOClass = LEO(1.0*Altitude, 1.0*Inclination,Geocutoff,solarmod)

ViewAtmo = 2*np.pi * (np.cos(np.deg2rad(LEOClass.HorizonAngle)) + 1)
ViewSky = 2*np.pi * (1-np.cos(np.deg2rad(LEOClass.HorizonAngle)))



if components == None :


    Particle = ["AtmosphericNeutrons", 
         "CosmicPhotons", 
         "PrimaryProtons",
         "SecondaryProtonsUpward","SecondaryProtonsDownward", "PrimaryAlphas", "PrimaryElectrons",
         "PrimaryPositrons", "SecondaryElectrons", "SecondaryPositrons",
         "AlbedoPhotons"
         ]

else :


    Particle = components.split(",")
    if Geocutoff is None :
        for f in Particle :
            if f in ["AtmosphericNeutrons", 
                    "SecondaryProtonsUpward",
                    "SecondaryProtonsDownward", 
                    "SecondaryElectrons", 
                    "SecondaryPositrons",
                    "AlbedoPhotons"] :
                        print("Error : If cutoff not set, you need to choose a Primary components (proton,alpha,etc...) or Cosmic photon")
                        sys.exit()
                    

# if cutoff is not set, we can do only the Primary components
if Geocutoff and components is None :
     Particle = [ 
         "CosmicPhotons", 
         "PrimaryProtons",
         "PrimaryAlphas", "PrimaryElectrons",
         "PrimaryPositrons" 
         ]    



fac = [ViewAtmo, ViewSky,ViewSky, 2*np.pi, 2*np.pi, ViewSky, ViewSky, ViewSky,4*np.pi,4*np.pi,ViewAtmo]

for i in range(0, len(Particle)):

    Energies = np.logspace(Elow, Ehigh, num=100, endpoint=True, base=10.0)
    if Geocutoff==None :
        Output = "%s/%s_Spec_%skm_%sdeg_%ssolarmod.dat" % (outputpath,Particle[i], float(Altitude), float(Inclination),float(solarmod))
    else :
        Output = "{0}/{1}_Spec_{2}km_{3}deg_{4:.3f}cutoff_{5}solarmod.dat".format(outputpath,Particle[i], float(Altitude), float(Inclination),float(Geocutoff),float(solarmod))
    #print(Megalibfunc[i](Energies))
    IntSpectrum = np.trapz(getattr(LEOClass, Particle[i])(Energies),Energies)
    print(Particle[i], IntSpectrum*fac[i], " #/cm^2/s")
    with open(Output, 'w') as f:
        print('# %s spectrum ' % Particle[i], file=f)
        print('# Format: DP <energy in keV> <shape of differential spectrum [XX/keV]>', file=f)
        print('# Although cosima doesn\'t use it the spectrum here is given as a flux in #/cm^2/s/keV', file=f)
        print('# Integrated over %s sr' % fac[i], file=f)
        #print('# Flux: %s #/cm^2/s/str' % (IntSpectrum), file=f)
        print('# Integral Flux: %s #/cm^2/s' % (IntSpectrum*fac[i]), file=f)
        print('', file=f)
        print('IP LOGLOG', file=f)
        print('', file=f)
        for j in range(0, len(Energies)):
            print('DP', Energies[j], getattr(LEOClass, Particle[i])(Energies[j]), file=f)
        print('EN', file=f)
