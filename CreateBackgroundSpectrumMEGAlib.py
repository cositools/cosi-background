#!/usr/bin/env python

""" Create a .dat file containing the spectrum of the primary
    and secondary protons plus albedo neutrons from the class
    LEOBackgroundGenerator to be used as input for the Step1
    *.source for Activation Simulations with MEGAlib.
"""

import numpy as np
import pandas as pd

from scipy.integrate import quad

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


args = pars.parse_args()

Inclination = args.inclination
Altitude = args.altitude

Elow = args.elow
Ehigh = args.ehigh

Geocutoff = args.cutoff
outputpath = args.outputpath

solarmod = args.solarmodulation

LEOClass = LEO(1.0*Altitude, 1.0*Inclination,Geocutoff,solarmod)

ViewAtmo = 2*np.pi * (np.cos(np.deg2rad(LEOClass.HorizonAngle)) + 1)
ViewSky = 2*np.pi * (1-np.cos(np.deg2rad(LEOClass.HorizonAngle)))

Particle = ["AtmosphericNeutrons", 
         "CosmicPhotons", 
         "PrimaryProtons",
         "SecondaryProtonsUpward","SecondaryProtonsDownward", "PrimaryAlphas", "PrimaryElectrons",
         "PrimaryPositrons", "SecondaryElectrons", "SecondaryPositrons",
         "AlbedoPhotons"
         ]

Megalibfunc = [LEOClass.AtmosphericNeutrons,
                LEOClass.CosmicPhotons,
               LEOClass.PrimaryProtons, LEOClass.SecondaryProtonsUpward,LEOClass.SecondaryProtonsDownward,
               LEOClass.PrimaryAlphas, LEOClass.PrimaryElectrons,
               LEOClass.PrimaryPositrons, LEOClass.SecondaryElectrons,
               LEOClass.SecondaryPositrons, LEOClass.AlbedoPhotons
              ]

fac = [ViewAtmo, ViewSky,ViewSky, 2*np.pi, 2*np.pi, ViewSky, ViewSky, ViewSky,4*np.pi,4*np.pi,ViewAtmo]

for i in range(0, len(Megalibfunc)):

    Energies = np.logspace(Elow, Ehigh, num=100, endpoint=True, base=10.0)
    if Geocutoff==None :
        Output = "%s/%s_Spec_%skm_%sdeg_%ssolarmod.dat" % (outputpath,Particle[i], float(Altitude), float(Inclination),float(solarmod))
    else :
        Output = "%s/%s_Spec_%skm_%sdeg_%scutoff_%ssolarmod.dat" % (outputpath,Particle[i], float(Altitude), float(Inclination),float(Geocutoff),float(solarmod))
    #print(Megalibfunc[i](Energies))
    IntSpectrum = np.trapz(Megalibfunc[i](Energies),Energies)
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
            print('DP', Energies[j], Megalibfunc[i](Energies[j]), file=f)
        print('EN', file=f)
