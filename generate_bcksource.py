import numpy as np
import os
import shutil
import pandas as pd
import argparse
import aacgmv2
from astropy.constants import R_earth
from astropy.time import Time
import datetime


# Generate the source files for the background, depending on the cutoff value, solar activity
# ,altitude and inclination


# Instantiate the parser
pars = argparse.ArgumentParser(
    description='Generate the source files for the background, depending on the cutoff value, solar activity,altitude, inclination and other parameters')


pars.add_argument(
    'Outputpath', help='full path where the source files will be created')
pars.add_argument('geopath', help='full path of the massmodel')
pars.add_argument('beamfilepath', help='full path of the beam.dat file')
pars.add_argument('Lightcurvepath', help='full path of the light curve file')
pars.add_argument("time", type=float, help='Exposure time (s)')
pars.add_argument('component', type=str,
                   help='''the component you want : 
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
		
		''')

pars.add_argument("-o", '--orientationfile', default=None,
                  help='full path of orientation file [None]')

pars.add_argument("-la", "--latitude", default=None,
                  type=float, help='latitude (deg) [None]')

pars.add_argument("-lon", "--longitude", default=None,
                  type=float, help='longitude (deg)[None]')

pars.add_argument('-d', "--date", default=None,
                  type=float, help='date (MJD)[None]')


pars.add_argument("-c", "--cutoff", default=10., type=float,
                  help='cut off rigidity (GeV) [10.]')

pars.add_argument("-g", "--geomlat", default=-0.1, type=float,
                  help='geomagnetic lat (rad) [-0.1]')

pars.add_argument("-mc", "--mcosima", default=None, type=int,
                  help='number of core to use for mcosima [None]')

pars.add_argument('-i', '--inclination', type=float, nargs='?',
                  default=0., help='Inclination of the orbit in degree [0.]')

pars.add_argument('-a', '--altitude', type=float, nargs='?',
                  default=550., help='Altitude of the orbit in km [550.]')

pars.add_argument('-s', '--solarmodulation', type=float, nargs='?',
                  default=650., help='solar modulation (550 min and 1100 max) [650.]')

pars.add_argument('-ir', '--irradiationtime', type=float, nargs='?',
                  default=None, help='Irradiation time for the step 2 of activations bck [time]')




# parse the arguments
args = pars.parse_args()


outputpath = args.Outputpath
geopath = args.geopath
latitude = args.latitude
longitude = args.longitude
date = args.date
orientationfile = args.orientationfile
time = args.time
inclination = args.inclination
altitude = args.altitude
solarmodulation = args.solarmodulation
IrradiationTime = args.irradiationtime
Rcutoff = args.cutoff
geomlat = args.geomlat
mc = args.mcosima
beamfilepath = args.beamfilepath
Lightcurvepath = args.Lightcurvepath
component = args.component



os.system("cp {1}/*.beam.dat {0}".format(outputpath, beamfilepath))


# split the time to each CPU
if mc is not None:

    time = time/mc
    print("**{0} CPU for mcosima so time per job will be {1}s**\n".format(mc, time))


# if irradiation time is not set take the same as exposure time
if IrradiationTime == None:
    print("**Irradiation time for step 2 of activation not specified : we take exposure time as reference **\n")
    IrradiationTime = args.time

###############################################################################################


if date and latitude and longitude and altitude is not None:

    print("**Compute the Rcutoff with the coordinate and date you entered**\n")

    # Convert MJD to datetime
    t_bin = Time(date, format='mjd').datetime
    # get the geomagnetic coordinates
    geocordinate = np.array(aacgmv2.get_aacgm_coord(
        latitude, longitude, altitude, t_bin))

    # compute the Rcutof
    """ Average Geomagnetic cutoff in GV
    for a dipole approximations
    Equation 4 Smart et al. 2005
    doi:10.1016/j.asr.2004.09.015
    """
    EarthRadius = R_earth.to('km').value
    R_E = R_earth.to('cm').value
    # g 01 term (in units of G) from IGRF-12 for 2020-25
    g10 = 29405 * 10**(-9) * 10**4  # G

    M = g10*R_E*300/10**9  # GV/cm2
    # formula took from : 10.1007/s10686-019-09624-0
    cutoff = (M/4*(1+altitude/EarthRadius)**(-2.0) *
              np.cos(np.deg2rad(geocordinate[0]))**4)

    print("#####################################################\n")
    print("Mean latitude in the orientation bin {0} deg".format(latitude))
    print("Mean longitude in the orientation bin {0} deg".format(longitude))
    print("Corresponding Geomagnetic latitude {0} deg".format(longitude))
    print("Rcutoff value is {0:.3f} GeV".format(cutoff))
    print("#####################################################\n")

else:
    cutoff = Rcutoff
    print("#####################################################\n")
    print("Rcutoff input value is {0:.3f} GV\n".format(cutoff))
    print("#####################################################\n")


# generate the .dat input spectrum file for MEGALib (by default altitude is 550km, solar mod 650MV , inclination 0 deg and cutoff 10~GV)


os.system("python3 CreateBackgroundSpectrumMEGAlib.py -c {0:.3f} -i {1} -a {2} -s {3} -o {4} -g {5} -f {6}".format(
    cutoff, inclination, altitude, solarmodulation, outputpath, geomlat, component))


os.chdir(outputpath)


###############################################################################################






# get the flux for the component in  #/cm^2/s

Infile = open(outputpath+"/"+component+"_Spec_{0}km_{1}deg_{2:.3f}cutoff_{3}solarmod.dat".format(
        altitude, inclination, cutoff, solarmodulation), "r")
lines = Infile.readlines()
flux = (float(lines[4].split()[3]))



print("#####################################################")

"""
#General source file
f = open("General.source","w")
f.write("Geometry {0} \n".format(geopath))
f.write("DetectorTimeConstant 0.000002 \n")
f.write("StoreSimulationInfo  init-only \n")
f.write("Include PrimaryElectron.source\n")
f.write("Include SecondaryElectron.source\n")
f.write("Include PrimaryPositron.source\n")
f.write("Include SecondaryPositron.source\n")
f.write("Include CosmicPhoton.source\n")
f.write("Include AlbedoPhoton.source\n")
f.write("Include NeutronAtmospheric_step1.source\n")
f.write("Include NeutronAtmospheric_step2.source\n")
f.write("Include NeutronAtmospheric_step3.source\n")
f.write("Include SecondaryProtons_step1.source\n")
f.write("Include SecondaryProtons_step2.source\n")
f.write("Include SecondaryProtons_step3.source\n")
f.write("Include PrimaryAlpha_step1.source\n")
f.write("Include PrimaryAlpha_step2.source\n")
f.write("Include PrimaryAlpha_step3.source\n")
f.write("Include PrimaryProtons_step1.source\n")
f.write("Include PrimaryProtons_step2.source\n")
f.write("Include PrimaryProtons_step3.source\n")
f.close()
"""

if component == "PrimaryElectrons" :

    # Primary electron
    f = open("PrimaryElectrons.source", "w")
    f.write("Geometry {0} \n".format(geopath))
    f.write("DetectorTimeConstant 0.000002 \n")
    f.write("StoreSimulationInfo  init-only \n")
    f.write("PhysicsListHD   qgsp-bic-hp\n")
    f.write("PhysicsListEM   LivermorePol\n\n")
    f.write("Run SpaceSim\n")
    f.write("SpaceSim.FileName  PrimaryElectron\n")
    f.write("SpaceSim.Time  {0}\n".format(time))
    if orientationfile is not None:
        f.write("SpaceSim.OrientationSky Galactic File NoLoop {0}\n\n".format(
            orientationfile))
    f.write("SpaceSim.Source PrimaryElectron\n")
    f.write("PrimaryElectron.ParticleType  3\n")
    f.write("PrimaryElectron.Beam   FarFieldFileZenithDependent CosmicElectronsMizuno.beam.dat\n")
    f.write("PrimaryElectron.Spectrum File PrimaryElectrons_Spec_{0}km_{1}deg_{2:.3f}cutoff_{3}solarmod.dat\n".format(
        altitude, inclination, cutoff, solarmodulation))
    f.write("PrimaryElectron.Flux		      {0}\n\n".format(flux))
    f.write("PrimaryElectron.LightCurve File false {0}\n".format(Lightcurvepath))
    f.close()


if component =="SecondaryElectrons" :
    # Secondary electron
    f = open("SecondaryElectrons.source", "w")
    f.write("Geometry {0} \n".format(geopath))
    f.write("DetectorTimeConstant 0.000002 \n")
    f.write("StoreSimulationInfo  init-only \n")
    f.write("PhysicsListHD   qgsp-bic-hp\n")
    f.write("PhysicsListEM   LivermorePol\n\n")
    f.write("Run SpaceSim\n")
    f.write("SpaceSim.FileName  SecondaryElectron\n")
    f.write("SpaceSim.Time  {0}\n".format(time))
    if orientationfile is not None:
        f.write("SpaceSim.OrientationSky Galactic File NoLoop {0}\n\n".format(
            orientationfile))
    f.write("SpaceSim.Source SecondaryElectrons\n")
    f.write("SecondaryElectrons.ParticleType  3\n")
    f.write("SecondaryElectrons.Beam   FarFieldFileZenithDependent AlbedoElectronsAlcarazMizuno.beam.dat\n")
    f.write("SecondaryElectrons.Spectrum File SecondaryElectrons_Spec_{0}km_{1}deg_{2:.3f}cutoff_{3}solarmod.dat\n".format(
        altitude, inclination, cutoff, solarmodulation))
    f.write("SecondaryElectrons.Flux		      {0}\n\n".format(flux))
    f.write("SecondaryElectron.LightCurve File false {0}\n".format(Lightcurvepath))
    f.close()


if component == "PrimaryPositrons" :

    # Primary Positron
    f = open("PrimaryPositrons.source", "w")
    f.write("Geometry {0} \n".format(geopath))
    f.write("DetectorTimeConstant 0.000002 \n")
    f.write("StoreSimulationInfo  init-only \n")
    f.write("PhysicsListHD   qgsp-bic-hp\n")
    f.write("PhysicsListEM   LivermorePol\n\n")
    f.write("Run SpaceSim\n")
    f.write("SpaceSim.FileName  PrimaryPositron\n")
    f.write("SpaceSim.Time  {0}\n".format(time))
    if orientationfile is not None:
        f.write("SpaceSim.OrientationSky Galactic File NoLoop {0}\n\n".format(
            orientationfile))
    f.write("SpaceSim.Source PrimaryPositrons\n")
    f.write("PrimaryPositrons.ParticleType  2\n")
    f.write("PrimaryPositrons.Beam   FarFieldFileZenithDependent CosmicPositronsMizuno.beam.dat\n")
    f.write("PrimaryPositrons.Spectrum File PrimaryPositrons_Spec_{0}km_{1}deg_{2:.3f}cutoff_{3}solarmod.dat\n".format(
        altitude, inclination, cutoff, solarmodulation))
    f.write("PrimaryPositrons.Flux		      {0}\n\n".format(flux))
    f.write("PrimaryPositron.LightCurve File false {0}\n".format(Lightcurvepath))
    f.close()

if component == "SecondaryPositrons" :
    
    # Secondary Positron
    f = open("SecondaryPositrons.source", "w")
    f.write("Geometry {0} \n".format(geopath))
    f.write("DetectorTimeConstant 0.000002 \n")
    f.write("StoreSimulationInfo  init-only \n")
    f.write("PhysicsListHD   qgsp-bic-hp\n")
    f.write("PhysicsListEM   LivermorePol\n\n")
    f.write("Run SpaceSim\n")
    f.write("SpaceSim.FileName  SecondaryPositron\n")
    f.write("SpaceSim.Time  {0}\n".format(time))
    if orientationfile is not None:
        f.write("SpaceSim.OrientationSky Galactic File NoLoop {0}\n\n".format(
            orientationfile))
    f.write("SpaceSim.Source SecondaryPositrons\n")
    f.write("SecondaryPositrons.ParticleType  2\n")
    f.write("SecondaryPositrons.Beam   FarFieldFileZenithDependent AlbedoPositronsAlcarazMizuno.beam.dat\n")
    f.write("SecondaryPositrons.Spectrum File SecondaryPositrons_Spec_{0}km_{1}deg_{2:.3f}cutoff_{3}solarmod.dat\n".format(
        altitude, inclination, cutoff, solarmodulation))
    f.write("SecondaryPositrons.Flux		      {0}\n\n".format(flux))
    f.write("SecondaryPositrons.LightCurve File false {0}\n".format(Lightcurvepath))
    f.close()

if component == "CosmicPhotons" :

    # Cosmic Photons
    f = open("CosmicPhotons.source", "w")
    f.write("Geometry {0} \n".format(geopath))
    f.write("DetectorTimeConstant 0.000002 \n")
    f.write("StoreSimulationInfo  init-only \n")
    f.write("PhysicsListEM LivermorePol\n\n")
    f.write("Run SpaceSim\n")
    f.write("SpaceSim.FileName  CosmicPhotons\n")
    f.write("SpaceSim.Time  {0}\n".format(time))
    if orientationfile is not None:
        f.write("SpaceSim.OrientationSky Galactic File NoLoop {0}\n\n".format(
            orientationfile))
    f.write("SpaceSim.Source CosmicPhotons\n")
    f.write("CosmicPhotons.ParticleType	    1\n")
    f.write(
        "CosmicPhotons.Beam	FarFieldFileZenithDependent CosmicPhotonsGruber.beam.dat\n")
    f.write("CosmicPhotons.Spectrum	File CosmicPhotons_Spec_{0}km_{1}deg_{2:.3f}cutoff_{3}solarmod.dat\n".format(
        altitude, inclination, cutoff, solarmodulation))
    f.write("CosmicPhotons.Flux		   {0}\n\n".format(flux))
    f.close()

if component == "AlbedoPhotons" :

    # Albedo Photons
    f = open("AlbedoPhotons.source", "w")
    f.write("Geometry {0} \n".format(geopath))
    f.write("DetectorTimeConstant 0.000002 \n")
    f.write("StoreSimulationInfo  init-only \n")
    f.write("PhysicsListEM LivermorePol\n")
    f.write("Run SpaceSim\n")
    f.write("SpaceSim.FileName  AlbedoPhotons\n")
    f.write("SpaceSim.Time  {0}\n".format(time))
    if orientationfile is not None:
        f.write("SpaceSim.OrientationSky Galactic File NoLoop {0}\n\n".format(
            orientationfile))
    f.write("SpaceSim.Source AlbedoPhotons\n")
    f.write("AlbedoPhotons.ParticleType	       1\n")
    f.write("AlbedoPhotons.Beam FarFieldFileZenithDependent AlbedoPhotonsTuerlerMizunoAbdo.beam.dat\n")
    f.write("AlbedoPhotons.Spectrum  File AlbedoPhotons_Spec_{0}km_{1}deg_{2:.3f}cutoff_{3}solarmod.dat\n".format(
        altitude, inclination, cutoff, solarmodulation))
    f.write("AlbedoPhotons.Flux	  {0}\n\n".format(flux))
    f.write("AlbedoPhotons.LightCurve File false {0}\n".format(Lightcurvepath))
    f.close()

if component == "AtmosphericNeutrons" :

    # Neutron Atm (step 1)
    f = open("AtmosphericNeutrons_step1.source", "w")
    f.write("Geometry {0} \n".format(geopath))
    f.write("DetectorTimeConstant 0.000002 \n")
    f.write("StoreSimulationInfo  init-only \n")
    f.write("PhysicsListHD   qgsp-bic-hp\n")
    f.write("PhysicsListEM   LivermorePol\n")
    f.write("DecayMode	 ActivationBuildup\n")
    f.write("Run SpaceSim\n")
    f.write("SpaceSim.FileName		    NeutronAtm\n\n")
    f.write("SpaceSim.Time  {0}\n".format(time))
    f.write("SpaceSim.IsotopeProductionFile     NeutronAtmIsotopes\n")
    if orientationfile is not None:
        f.write("SpaceSim.OrientationSky Galactic File NoLoop {0}\n\n".format(
            orientationfile))
    f.write("SpaceSim.Source NeutronAtm\n")
    f.write("NeutronAtm.ParticleType   6\n")
    f.write("NeutronAtm.Beam FarFieldFileZenithDependent AlbedoNeutronsKole.beam.dat\n")
    f.write("NeutronAtm.Spectrum File AtmosphericNeutrons_Spec_{0}km_{1}deg_{2:.3f}cutoff_{3}solarmod.dat\n".format(
        altitude, inclination, cutoff, solarmodulation))
    f.write("NeutronAtm.Flux   {0}\n\n".format(flux))
    f.write("NeutronAtm.LightCurve File false {0}\n".format(Lightcurvepath))
    f.close()

    # Neutron Atm (step 2)
    f = open("AtmosphericNeutrons_step2.source", "w")
    f.write("Geometry {0} \n".format(geopath))
    f.write("DetectorTimeConstant               0.000002\n")
    f.write("PhysicsListHD                      qgsp-bic-hp\n")
    f.write("PhysicsListEM                      LivermorePol\n")
    f.write("StoreSimulationInfo                all\n")
    f.write("Activator A\n")
    f.write("A.ActivationMode          ConstantIrradiation  {0}\n".format(
        IrradiationTime))
    f.write("A.ActivationFile          NeutronAtm_Activation.dat\n")

    if mc == None:
        f.write("A.IsotopeProductionFile   NeutronAtmIsotopes.dat\n")

    else:
        for i in range(1, mc+1):
            f.write(
                "A.IsotopeProductionFile   NeutronAtmIsotopes.p1.inc{0}.dat\n".format(i))

    f.close()

    # Neutron Atm (step 3)
    f = open("AtmosphericNeutrons_step3.source", "w")
    f.write("Geometry {0} \n".format(geopath))
    f.write("DetectorTimeConstant               0.000002\n")
    f.write("StoreSimulationInfo                all\n")
    f.write("PhysicsListHD                      qgsp-bic-hp\n")
    f.write("PhysicsListEM                      LivermorePol\n")
    f.write("DecayMode                          ActivationDelayedDecay\n")
    f.write("Run ActivationStep3\n")
    f.write("ActivationStep3.FileName                         NeutronAtm_Decay\n")
    f.write(
        "ActivationStep3.Time                             {0}\n".format(time))
    f.write(
        "ActivationStep3.ActivationSources                NeutronAtm_Activation.dat\n")
    f.close()

if component == "SecondaryProtonsUpward" or component == "SecondaryProtonsDownward"  :

    # Secondary Protons (step 1)
    f = open("SecondaryProtons_step1.source", "w")
    f.write("Geometry {0} \n".format(geopath))
    f.write("DetectorTimeConstant 0.000002 \n")
    f.write("StoreSimulationInfo  init-only \n")
    f.write("PhysicsListHD   qgsp-bic-hp\n")
    f.write("PhysicsListEM   LivermorePol\n")
    f.write("DecayMode	 ActivationBuildup\n")
    f.write("Run SpaceSim\n")
    f.write("SpaceSim.FileName		    SecondaryProtons\n\n")
    f.write("SpaceSim.Time  {0}\n".format(time))
    f.write("SpaceSim.IsotopeProductionFile     SecondaryProtonsIsotopes\n")
    if orientationfile is not None:
        f.write("SpaceSim.OrientationSky Galactic File NoLoop {0}\n\n".format(
            orientationfile))
    f.write("SpaceSim.Source SecondaryProtons\n")
    f.write("SecondaryProtons.ParticleType           4\n")
    f.write(
        "SecondaryProtons.Beam  FarFieldFileZenithDependent AlbedoProtonMizuno.beam.dat\n")
    f.write("SecondaryProtons.Spectrum File SecondaryProtonsUpward_Spec_{0}km_{1}deg_{2:.3f}cutoff_{3}solarmod.dat\n".format(
        altitude, inclination, cutoff, solarmodulation))
    f.write("SecondaryProtons.Flux                 {0}\n\n".format(flux))
    f.write("SecondaryProtons.LightCurve File false {0}\n".format(Lightcurvepath))
    f.close()

    # Secondary Protons (step 2)
    f = open("SecondaryProtons_step2.source", "w")
    f.write("Geometry {0} \n".format(geopath))
    f.write("DetectorTimeConstant               0.000002\n")
    f.write("PhysicsListHD                      qgsp-bic-hp\n")
    f.write("PhysicsListEM                      LivermorePol\n")
    f.write("StoreSimulationInfo                all\n")
    f.write("Activator A\n")
    f.write("A.ActivationMode          ConstantIrradiation  {0}\n".format(
        IrradiationTime))
    f.write("A.ActivationFile          SecondaryProtons_Activation.dat\n")

    if mc == None:
        f.write("A.IsotopeProductionFile   SecondaryProtonsIsotopes.dat\n")

    else:
        for i in range(1, mc+1):
            f.write(
                "A.IsotopeProductionFile   SecondaryProtonsIsotopes.p1.inc{0}.dat\n".format(i))

    f.close()

    # Secondary Protons (step 3)
    f = open("SecondaryProtons_step3.source", "w")
    f.write("Geometry {0} \n".format(geopath))
    f.write("DetectorTimeConstant               0.000002\n")
    f.write("StoreSimulationInfo                all\n")	
    f.write("PhysicsListHD                      qgsp-bic-hp\n")
    f.write("PhysicsListEM                      LivermorePol\n")
    f.write("DecayMode                          ActivationDelayedDecay\n")
    f.write("Run ActivationStep3\n")
    f.write("ActivationStep3.FileName                         SecondaryProtons_Decay\n")
    f.write(
        "ActivationStep3.Time                             {0}\n".format(time))
    f.write("ActivationStep3.ActivationSources                SecondaryProtons_Activation.dat\n")
    f.close()

if component == "PrimaryProtons" :

    # Primary Protons (step 1)
    f = open("PrimaryProtons_step1.source", "w")
    f.write("Geometry {0} \n".format(geopath))
    f.write("DetectorTimeConstant 0.000002 \n")
    f.write("StoreSimulationInfo  init-only \n")
    f.write("PhysicsListHD   qgsp-bic-hp\n")
    f.write("PhysicsListEM   LivermorePol\n")
    f.write("DecayMode	 ActivationBuildup\n")
    f.write("Run SpaceSim\n")
    f.write("SpaceSim.FileName		    PrimaryProtons\n\n")
    f.write("SpaceSim.Time  {0}\n".format(time))
    f.write("SpaceSim.IsotopeProductionFile     PrimaryProtonsIsotopes\n")
    if orientationfile is not None:
        f.write("SpaceSim.OrientationSky Galactic File NoLoop {0}\n\n".format(
            orientationfile))
    f.write("SpaceSim.Source PrimaryProtons\n")
    f.write("PrimaryProtons.ParticleType	     4\n")
    f.write(
        "PrimaryProtons.Beam FarFieldFileZenithDependent CosmicProtonsSpenvis.beam.dat\n")
    f.write("PrimaryProtons.Spectrum  	     File PrimaryProtons_Spec_{0}km_{1}deg_{2:.3f}cutoff_{3}solarmod.dat\n".format(
        altitude, inclination, cutoff, solarmodulation))
    f.write("PrimaryProtons.Flux		     {0}\n\n".format(flux))
    f.write("PrimaryProtons.LightCurve File false {0}\n".format(Lightcurvepath))
    f.close()

    # Primary Protons (step 2)
    f = open("PrimaryProtons_step2.source", "w")
    f.write("Geometry {0} \n".format(geopath))
    f.write("DetectorTimeConstant               0.000002\n")
    f.write("PhysicsListHD                      qgsp-bic-hp\n")
    f.write("PhysicsListEM                      LivermorePol\n")
    f.write("StoreSimulationInfo                all\n")
    f.write("Activator A\n")
    f.write("A.ActivationMode          ConstantIrradiation  {0}\n".format(
        IrradiationTime))
    f.write("A.ActivationFile          PrimaryProtons_Activation.dat\n")

    if mc == None:
        f.write("A.IsotopeProductionFile   PrimaryProtonsIsotopes.dat\n")

    else:
        for i in range(1, mc+1):
            f.write(
                "A.IsotopeProductionFile   PrimaryProtonsIsotopes.p1.inc{0}.dat\n".format(i))

    f.close()

    # Primary Protons (step 3)
    f = open("PrimaryProtons_step3.source", "w")
    f.write("Geometry {0} \n".format(geopath))
    f.write("DetectorTimeConstant               0.000002\n")
    f.write("StoreSimulationInfo                all\n")	
    f.write("PhysicsListHD                      qgsp-bic-hp\n")
    f.write("PhysicsListEM                      LivermorePol\n")
    f.write("DecayMode                          ActivationDelayedDecay\n")
    f.write("Run ActivationStep3\n")
    f.write("ActivationStep3.FileName                         PrimaryProtons_Decay\n")
    f.write(
        "ActivationStep3.Time                             {0}\n".format(time))
    f.write(
        "ActivationStep3.ActivationSources                PrimaryProtons_Activation.dat\n")
    f.close()


if component == "PrimaryAlphas" :

    # Primary Alpha (step 1)
    f = open("PrimaryAlphas_step1.source", "w")
    f.write("Geometry {0} \n".format(geopath))
    f.write("DetectorTimeConstant 0.000002 \n")
    f.write("StoreSimulationInfo  init-only \n")
    f.write("PhysicsListHD   qgsp-bic-hp\n")
    f.write("PhysicsListEM   LivermorePol\n")
    f.write("DecayMode	 ActivationBuildup\n")
    f.write("Run SpaceSim\n")
    f.write("SpaceSim.FileName		    PrimaryAlpha\n\n")
    f.write("SpaceSim.Time  {0}\n".format(time))
    f.write("SpaceSim.IsotopeProductionFile     PrimaryAlphaIsotopes\n")
    if orientationfile is not None:
        f.write("SpaceSim.OrientationSky Galactic File NoLoop {0}\n\n".format(
            orientationfile))
    f.write("SpaceSim.Source PrimaryAlpha\n")
    f.write("PrimaryAlpha.ParticleType	    21\n")
    f.write(
        "PrimaryAlpha.Beam FarFieldFileZenithDependent  CosmicAlphasSpenvis.beam.dat\n")
    f.write("PrimaryAlpha.Spectrum File PrimaryAlphas_Spec_{0}km_{1}deg_{2:.3f}cutoff_{3}solarmod.dat\n".format(
        altitude, inclination, cutoff, solarmodulation))
    f.write("PrimaryAlpha.Flux		    {0}\n\n".format(flux))
    f.write("PrimaryAlpha.LightCurve File false {0}\n".format(Lightcurvepath))
    f.close()

    # Primary Alpha (step 2)
    f = open("PrimaryAlphas_step2.source", "w")
    f.write("Geometry {0} \n".format(geopath))
    f.write("DetectorTimeConstant               0.000002\n")
    f.write("PhysicsListHD                      qgsp-bic-hp\n")
    f.write("PhysicsListEM                      LivermorePol\n")
    f.write("StoreSimulationInfo                all\n")
    f.write("Activator A\n")
    f.write("A.ActivationMode          ConstantIrradiation  {0}\n".format(
        IrradiationTime))
    f.write("A.ActivationFile          PrimaryAlpha_Activation.dat\n")

    if mc == None:
        f.write("A.IsotopeProductionFile   PrimaryAlphaIsotopes.dat\n")

    else:
        for i in range(1, mc+1):
            f.write(
                "A.IsotopeProductionFile   PrimaryAlphaIsotopes.p1.inc{0}.dat\n".format(i))

    f.close()

    # Primary Alpha (step 3)
    f = open("PrimaryAlphas_step3.source", "w")
    f.write("Geometry {0} \n".format(geopath))
    f.write("DetectorTimeConstant               0.000002\n")
    f.write("StoreSimulationInfo                all\n")	
    f.write("PhysicsListHD                      qgsp-bic-hp\n")
    f.write("PhysicsListEM                      LivermorePol\n")
    f.write("DecayMode                          ActivationDelayedDecay\n")
    f.write("Run ActivationStep3\n")
    f.write("ActivationStep3.FileName                         PrimaryAlpha_Decay\n")
    f.write(
        "ActivationStep3.Time                             {0}\n".format(time))
    f.write(
        "ActivationStep3.ActivationSources                PrimaryAlpha_Activation.dat\n")
    f.close()
