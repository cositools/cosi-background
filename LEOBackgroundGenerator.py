__author__ = 'Paolo Cumani, Jurgen Kiener, Vincent Tatischeff, Andreas Zoglauer | Updated by Savitri Gallego and Francesco Fenu'

import numpy as np
import pandas as pd
from astropy.constants import R_earth, m_p, m_n, c
from scipy.optimize import fsolve
from scipy.interpolate import interp1d
from scipy.integrate import quad
import aacgmv2

class LEOBackgroundGenerator:
    """
    Class to generate a background spectrum for a low Earth orbit (LEO).
    It uses equations/data from:
    - Albedo Neutrons: Kole et al. 2015
                         doi:10.1016/j.astropartphys.2014.10.002
                       Lingenfelter 1963
                         doi:10.1029/JZ068i020p05633
    - Cosmic Photons: Türler et al. 2010
                        doi:10.1051/0004-6361/200913072
                      Mizuno et al. 2004
                        http://stacks.iop.org/0004-637X/614/i=2/a=1113
                      Ackermann et al. 2015
                        doi:10.1088/0004-637X/799/1/86
    - Galactic Center/Disk: Fermi-LAT collaboration
                              https://fermi.gsfc.nasa.gov/ssc/data/analysis/software/aux/gll_iem_v06.fits
    - Primary Protons: Aguilar et al. 2015
                           doi:10.1103/PhysRevLett.114.171103
    - Secondary Protons: Mizuno et al. 2004
                           http://stacks.iop.org/0004-637X/614/i=2/a=1113
    - Primary Alphas:  Aguilar et al. 2015b
                         doi:10.1103/PhysRevLett.115.211101
    - Primary Electrons: Aguilar et al. 2014
                           doi:10.1103/PhysRevLett.113.121102
                         Mizuno et al. 2004
                           http://stacks.iop.org/0004-637X/614/i=2/a=1113
    - Primary Positrons: Aguilar et al. 2014
                           doi:10.1103/PhysRevLett.113.121102
                         Mizuno et al. 2004
                           http://stacks.iop.org/0004-637X/614/i=2/a=1113
    - Secondary Electrons: Mizuno et al. 2004
                           http://stacks.iop.org/0004-637X/614/i=2/a=1113
    - Secondary Positrons: Mizuno et al. 2004
                           http://stacks.iop.org/0004-637X/614/i=2/a=1113
    - Atmospheric Photons: Sazonov et al. 2007
                        doi:10.1111/j.1365-2966.2007.11746.x
                      Churazov et al. 2006
                        doi:10.1111/j.1365-2966.2008.12918.x
                      Türler et al. 2010
                        doi:10.1051/0004-6361/200913072
                      Mizuno et al. 2004
                        http://stacks.iop.org/0004-637X/614/i=2/a=1113
                      Abdo et al. 2009
                        doi:10.1103/PhysRevD.80.122004

    Parameters
    ----------
    altitude: float
        altitude of the orbit
    inclination: float
        inclination of the orbit
    solarmodulation: float
        solar modulation potential in MV at the time of the observation
    """

    def __init__(self, altitude, inclination, geomlat,GeoCutoff=None, solarmodulation=None ):
        self.Alt = altitude  # instrument altitude (km)
        self.magl = inclination  # orbit inclination (deg.)
        self.geomlat = geomlat  #geomagnetic latitude (rad)
        self.AvGeomagCutOff = GeoCutoff #cutoff rigidity (GV)
                
        

        """ solar modulation potential (MV): ~550 for solar minimum
                                            ~1100 for solar maximum
        """
        if solarmodulation is None:
            self.solmod = 500.
        else:
            self.solmod = solarmodulation

        #print("solar modulation = {0}".format(self.solmod))


        EarthRadius = R_earth.to('km').value

        AtmosphereHeight = 40  # km

        self.HorizonAngle = 90.0 + np.rad2deg(np.arccos(
                            (EarthRadius + AtmosphereHeight)
                            / (EarthRadius+self.Alt)))
        
        
        if GeoCutoff is None :
            self.AvGeomagCutOff = self.ComputeRcut(self.geomlat, self.Alt)
        

    def ComputeRcut(self,geomaglat, altitude):
        """ Average Geomagnetic cutoff in GV
            for a dipole approximations
            Equation 4 Smart et al. 2005
            doi:10.1016/j.asr.2004.09.015
            
         arguments : geomaglat [rad] 
         
         return : cutoff [GV]   
        """
        EarthRadius = R_earth.to('km').value
        R_E = R_earth.to('cm').value
        # g 01 term (in units of G) from IGRF-12 for 2020-25 : 29405
        # g 01 term extrapolated to 2027 from fitting the precedent value : 29147.79
        g10 = 29147.79 * 10**(-9) * 10**4  # G

        
        M = g10*R_E*300/10**9  # GV/cm2
        
        Rcut = (M/4*(1+altitude/EarthRadius)**(-2.0)
                             * np.cos(geomaglat)**4)
                           
        return Rcut                     



    def log_interp1d(self, xx, yy, fill='extrapolate', kind='linear'):
        """Functions for an interpolation in log-space
           https://stackoverflow.com/questions/29346292/
        """

        logx = np.log10(xx)
        logy = np.log10(yy)
        if fill != 'extrapolate':
            lin_interp = interp1d(logx, logy, kind=kind, fill_value=fill, bounds_error=False)
        else:
            lin_interp = interp1d(logx, logy, kind=kind, fill_value=fill)

        def log_interp(zz): return np.power(10.0, lin_interp(np.log10(zz)))
        return log_interp

    def LingenfelterNeutrons(self):
        """ Low energy neutrons spectrum at the top of the atmosphere
        as calculated in Lingenfelter 1963
        """
        filename = 'Data/Neutrons_Lingenfelter.dat'
        data = pd.read_table(filename, sep=',')

        data["Ener(MeV)"] = data["Ener(MeV)"]
        data["Flux(n/cm2MeVs)"] = data["Flux(n/cm2MeVs)"]

        self.LowENeutrons = data.copy()

    def AtmosphericNeutrons(self, E):
        """Atmospheric neutrons determinined after Kole et al. 2015
        Assumptions:
        * Angular distribution is flat out to the Earth-horizon
        * Downward component can be neglected
        """

        """ Solar activity calculated from the solar modulation
        as linear between minimum and maximum (page 10 Kole et al. 2015)
        """
        solac = (self.solmod - 250.0)/1109.0

        Pressure = 0.  # in hPa

        EnergyMeV = 0.001*np.copy(np.asarray(E, dtype=float))
        Flux = np.copy(np.asarray(E, dtype=float))

        a = 0.0003 + (7.0-5.0*solac)*0.001*(1-np.tanh(np.deg2rad(180-4.0*self.geomlat)))
        b = 0.0140 + (1.4-0.9*solac)*0.1*(1-np.tanh(np.deg2rad(180-3.5*self.geomlat)))
        c = 180 - 42*(1-np.tanh(np.deg2rad(180-5.5*self.geomlat)))
        d = -0.008 + (6.0-1.0*solac)*0.001*(1-np.tanh(np.deg2rad(180-4.4*self.geomlat)))

        Slope1 = -0.29 * np.exp(-Pressure/7.5) + 0.735
        Norm1 = (a*Pressure + b)*np.exp(-Pressure/c) + d
        Mask1 = EnergyMeV < 0.9

        Slope2 = -0.247 * np.exp(-Pressure/36.5) + 1.4
        Norm2 = Norm1*pow(0.9, -Slope1+Slope2)
        Mask2 = np.logical_and(EnergyMeV >= 0.9, EnergyMeV < 15)

        Slope3 = -0.40 * np.exp(-Pressure/40.0) + 0.9
        Norm3 = Norm2*pow(15, -Slope2+Slope3)
        Mask3 = np.logical_and(EnergyMeV >= 15, EnergyMeV < 70)

        Slope4 = -0.46 * np.exp(-Pressure/100.0) + 2.53
        Norm4 = Norm3*pow(70, -Slope3+Slope4)
        Mask4 = EnergyMeV >= 70

        Flux[Mask1] = Norm1 * pow(EnergyMeV[Mask1], -Slope1)
        Flux[Mask2] = Norm2 * pow(EnergyMeV[Mask2], -Slope2)
        Flux[Mask3] = Norm3 * pow(EnergyMeV[Mask3], -Slope3)
        Flux[Mask4] = Norm4 * pow(EnergyMeV[Mask4], -Slope4)

        try:
            self.LowENeutrons
        except AttributeError:
            self.LingenfelterNeutrons()

        data = self.LowENeutrons
        f = self.log_interp1d(data["Ener(MeV)"].loc[data['Flux(n/cm2MeVs)'] > 0.],
                              data["Flux(n/cm2MeVs)"].loc[data['Flux(n/cm2MeVs)'] > 0.])

        LowEnergyNeutron = self.LingenfelterNeutrons

        Scaler = (Norm1 * pow(0.008, -Slope1))/f(0.008)

        Flux[EnergyMeV < 0.008] = f(EnergyMeV[EnergyMeV < 0.008]) * Scaler

        # View angle of the atmosphere = 4 PI - 2 PI (1-cos(HorizonAngle))
        AngleFactor = 2*np.pi * (np.cos(np.deg2rad(self.HorizonAngle)) + 1)

        return Flux / (AngleFactor * 1000.0)  # Switch from n/MeV/cm2/s to n/keV/cm2/s/sr.

    def TuerlerCosmicPhotons(self, E):
        """Equation 5 from Tuerler et al. 2010
           Return a flux in ph /cm2 /s /keV /sr
        """
        Flux = np.copy(np.asarray(E, dtype=float))
        E = np.asarray(E, dtype=float)
        Flux[E > 1000] = 0.
        Flux[E <= 1000] = 0.109 / ((E[E <= 1000]/28)**1.4+(E[E <= 1000]/28)**2.88)

        return Flux
        
    def GruberCosmicPhotons(self, E):
        """Equation 1 from Gruber et al. 1999
           Return a flux in ph /cm2 /s /keV /sr
        """
        Flux = np.copy(np.asarray(E, dtype=float))
        E = np.asarray(E, dtype=float)
        Flux[E <= 60] = 7.877*(E[E <= 60]**(-1.29))*np.exp(-E[E <= 60]/41.13)
        Flux[E > 60] = 0.0259*(E[E > 60]**(-6.5))/(60.**(-5.5)) + 0.504*(E[E > 60]**(-2.58))/(60.**(-1.58)) + 0.0288*(E[E > 60]**(-2.05))/(60.**(-1.05))

        return Flux

    def MizunoCosmicPhotons(self, E):
        """Equation 18 from Mizuno et al. 2004
           Return a flux in ph /cm2 /s /keV /sr
        """

        Flux = np.copy(np.asarray(E, dtype=float))
        E = np.asarray(E, dtype=float)
        Flux[E < 800] = 0.
        Flux[E >= 800] = 40.*pow(E[E >= 800]/1000, -2.15)/(10**7)
        return Flux

    def AckermannCosmicPhotons(self, E):
        """Equation 1 from Ackermann et al. 2015
           Using the foreground model A
           Return a flux in ph /cm2 /s /keV /sr
        """
        I100 = 0.95*10**(-7)/1000
        gamma = 2.32
        Ecut = 279*10**6
        E = np.asarray(E, dtype=float)
        Flux = np.copy(np.asarray(E, dtype=float))

        Flux[E < 800] = 0.
        Flux[E >= 800] = I100 * (E[E >= 800]/(100*1000))**(-gamma)*np.exp(-E[E >= 800]/Ecut)

        return Flux

    def CosmicPhotons(self, E):

        E = np.asarray(E, dtype=float)
        Flux = np.copy(E)

        Eint = fsolve(lambda x: self.AckermannCosmicPhotons(x)
                      - self.MizunoCosmicPhotons(x), 1200)[0]

        mask = np.logical_and(E >= 890, E < Eint)

        Flux[E < 890] = self.TuerlerCosmicPhotons(E[E < 890])
        Flux[mask] = self.MizunoCosmicPhotons(E[mask])
        Flux[E >= Eint] = self.AckermannCosmicPhotons(E[E >= Eint])

        return Flux

    def GalacticCenter(self, E):
        """ Read Table created by LATBackground.py
        with data from the Fermi-LAT collaboration
        https://fermi.gsfc.nasa.gov/ssc/data/analysis/software/aux/gll_iem_v06.fits
        for the average Galactic center region (b+-1 deg, l+-2.5deg),
        Return a flux in ph /cm2 /s /keV /sr
        """
        filename = 'Data/LATBackground.dat'
        data = pd.read_table(filename, sep='\s+', header=0, comment='#')

        fGC = self.log_interp1d(data['Energy'], data['FluxGCAv'], fill="NaN")

        return fGC(E)

    def GalacticDisk(self, E):
        """ Read Table created by LATBackground.py
        with data from the Fermi-LAT collaboration
        https://fermi.gsfc.nasa.gov/ssc/data/analysis/software/aux/gll_iem_v06.fits
        for the average Galactic Disk region (b+-1 deg, l+-90 deg) excluding the
        Galactic center region (b+-1 deg, l+-2.5deg),
        Return a flux in ph /cm2 /s /keV /sr
        """
        filename = 'Data/LATBackground.dat'
        data = pd.read_table(filename, sep='\s+', header=0, comment='#')

        fDisk = self.log_interp1d(data['Energy'], data['FluxDiskAv'], fill="NaN")

        return fDisk(E)

    def ChurazovAlbedoPhotons(self, E):
        """ Equation 9 from Churazov et al. 2006, valid up to 1 MeV
            Compute the reflected cosmic X-ray background
            using Tuerler et al. 2010
           Return a flux in ph /cm2 /s /keV /sr
        """

        thetamax = 180 - self.HorizonAngle  # deg max polar angle wrt nadir
        omega = 2*np.pi*(1-np.cos(np.deg2rad(thetamax)))

        Flux = np.copy(np.asarray(E, dtype=float))
        E = np.asarray(E, dtype=float)

        first = 1.22/((E/28.5)**(-2.54)+(E/51.3)**1.57-0.37)
        second = (2.93+(E/3.08)**4)/(1+(E/3.08)**4)
        third = (0.123+(E/91.83)**3.44)/(1+(E/91.83)**3.44)

        Flux = self.GruberCosmicPhotons(E)*first*second*third

        return Flux

    def SazonovAlbedoPhotons(self, E, geomaglat, solarmod):
        """ Equation 7 and 1 from Sazonov et al. 2007,
            hard X-ray surface brightness of the Earth’s atmosphere
            the rigidity cut-off is computed at the Earth's surface
            the zenith angle is integrated to compute the angle-averaged surface brightness
           Return a flux in ph /cm2 /s /keV /sr
        """

        Flux = np.copy(np.asarray(E, dtype=float))
        E = np.asarray(E, dtype=float)

        thetamax = 180 - self.HorizonAngle  # deg max polar angle wrt nadir
        omega = 2*np.pi*(1-np.cos(np.deg2rad(thetamax)))
        cosomega = np.cos(np.deg2rad(thetamax))
        
        phi = solarmod / 1000  # GV
        
        # the rigidity cut-off used in the formula is computed at the geomagnetic latitude of the orbit position
        # the C intensity factor is ony integrated in zenith angle
        R_c = self.ComputeRcut(geomaglat, 40)
        part1 = (3./(5*np.pi))*1.47*0.0178*((((phi/2.8)**0.4) + ((phi/2.8)**1.5))**(-1))
        part2 = 1.3*(phi**0.25)*(1. + 2.5*(phi**0.4))
        f_E = 1./(((E/44.)**(-5)) + ((E/44.)**(1.4)))
        I_albedo = quad(lambda x: 2.*np.pi*np.cos(x)*(1. + np.cos(x))*np.sin(x), 0, np.deg2rad(thetamax))
        C_int = part1*I_albedo[0]*(1./np.sqrt(1. + ((R_c/part2)**2.)))
        Flux[E > 2000] = 0.
        Flux[E <= 2000] = C_int*f_E/omega

        return Flux

    def MizunoAlbedoPhotons(self, E):
        """ Equation 21 to 23 from Mizuno et al. 2004,
           Vertically upward
           Return a flux in ph /cm2 /s /keV /sr
        """

        Flux = np.copy(np.asarray(E, dtype=float))
        E = np.asarray(E, dtype=float)

        mask = np.logical_and(E >= 20000, E < 1000000.)
        maskLE = np.logical_and(E >= 1000, E < 20000.)

        Flux[E < 1000] = 0
        Flux[maskLE] = 1010.0*pow(E[maskLE]/1000, -1.34) / 10**7
        Flux[mask] = 7290.0*pow(E[mask]/1000, -2.0) / 10**7
        Flux[E >= 1000000.] = 29000*pow(E[E >= 1000000.]/1000, -2.2) / 10**7

        return Flux

    def AbdoAlbedoPhotons(self, E):
        """ From Abdo et al. 2010,
           Return a flux in ph /cm2 /s /keV /sr
        """

        Flux = np.copy(np.asarray(E, dtype=float))
        E = np.asarray(E, dtype=float)

        Flux[E < 100000] = 0
        Flux[E >= 100000] = 1.823e-8*pow(E[E >= 100000]/200000, -2.8)
        return Flux

    def AlbedoPhotons(self, E):
        """ Generate an albedo photon spectrum after
        Sazonov et al. 2007 & Churazov et al. 2006
        Gruber et al 1999
        Mizuno et al. 2004
        Abdo et al. 2010
        Sazonov is used to normalize Mizuno solar modulation dependence
        Abdo is normalized to Mizuno (not correct?)
        Return a flux in ph /cm2 /s /keV /sr
        """
        
        # Scaling from Mizuno et al. 2004
        Rcut_desired = self.ComputeRcut(self.geomlat, 40)
        Rcut_Mizuno = 4.5
        ScalerMizuno = pow(Rcut_desired/Rcut_Mizuno, -1.13)
        
        # Computing the scaling for Mizuno
        # we compute the Sazonov model at Phi = 500 MV and rescale the Mizuno flux for this value.
        # the Mizuno solar modulation is introduced as a scaling factor using the Sazonov dependence
        # the Abdo model is not scaled
        
        ChurazovSazonovValueAtPhi500MV_1850keV = self.ChurazovAlbedoPhotons(1850.) + self.SazonovAlbedoPhotons(1850., self.geomlat, 500)
        MizunoValue_1850keV = ScalerMizuno * self.MizunoAlbedoPhotons(1850.)
        ScalerMizuno_Phi500MV = ChurazovSazonovValueAtPhi500MV_1850keV/MizunoValue_1850keV
        ScalerPhi = self.SazonovAlbedoPhotons(1850., self.geomlat, self.solmod)/ChurazovSazonovValueAtPhi500MV_1850keV
        
        Flux = np.copy(np.asarray(E, dtype=float))
        E = np.asarray(E, dtype=float)
        
        Flux[E < 1850.] = self.ChurazovAlbedoPhotons(E[E < 1850.]) + self.SazonovAlbedoPhotons(E[E < 1850.], self.geomlat, self.solmod)
        
        mask = np.logical_and(E >= 1850, E < 200000.)
        maskabdo = E >= 200000.
        
        Flux[mask] = ScalerPhi * ScalerMizuno_Phi500MV * ScalerMizuno * self.MizunoAlbedoPhotons(E[mask])
        
        MizunoValue = ScalerPhi * ScalerMizuno_Phi500MV * ScalerMizuno * self.MizunoAlbedoPhotons(200000)
        AbdoValue = self.AbdoAlbedoPhotons(200000)
        ScalerAbdo = MizunoValue/AbdoValue
        Flux[maskabdo] = ScalerAbdo * self.AbdoAlbedoPhotons(E[maskabdo])
        
        return Flux

    def MizunoCutoffpl(self, f0, f1, a, ec, E):
        """Function describing a power-law with a cutoff
        """
        ec = ec * 1000
        Flux = np.copy(np.asarray(E, dtype=float))
        mask = np.logical_and(E >= 1, E < 100)
        maskHE = E >= 100

        Flux[E < 1] = 0.
        Flux[mask] = f0*pow(E[mask]/100, -1)
        Flux[maskHE] = f1*pow(E[maskHE]/1000, -a) * np.exp(
                                         -pow(E[maskHE]/ec, -a+1))

        return Flux

    def MizunoBrokenpl(self, f0, a, eb, b, E):
        """Function describing a power-law with a cutoff
        """
        eb = eb
        Flux = np.copy(np.asarray(E, dtype=float))
        mask = np.logical_and(E >= 1, E < 100)
        maskHE = np.logical_and(E >= 100, E < eb)

        Flux[E < 1] = 0.
        Flux[mask] = f0*pow(E[mask]/100, -1)
        Flux[maskHE] = f0*pow(E[maskHE]/100, -a)
        Flux[E >= eb] = f0*pow(eb/100, -a) * pow(E[E >= eb]/eb, -b)

        return Flux

    def SecondaryProtons(self, E):
        """ Equation 8 from Mizuno et al. 2004,
            the geomagnetic latitude intervals
            are translated in cutoff intervals
            considering the spacecraft orbit
           Return a flux in ph /cm2 /s /keV /sr
        """

        EnergyMeV = 0.001*np.copy(np.asarray(E, dtype=float))

        Rcut = self.AvGeomagCutOff
        
        if Rcut >= self.ComputeRcut(0.2, self.Alt) and Rcut <= self.ComputeRcut(0, self.Alt):
            FluxU = self.MizunoCutoffpl(0.136, 0.123, 0.155, 0.51, EnergyMeV)
            FluxD = self.MizunoCutoffpl(0.136, 0.123, 0.155, 0.51, EnergyMeV)
        elif Rcut >= self.ComputeRcut(0.3, self.Alt) and Rcut <= self.ComputeRcut(0.2, self.Alt):
            FluxU = self.MizunoBrokenpl(0.1, 0.87, 600, 2.53, EnergyMeV)
            FluxD = self.MizunoBrokenpl(0.1, 0.87, 600, 2.53, EnergyMeV)
        elif Rcut >= self.ComputeRcut(0.4, self.Alt) and Rcut <= self.ComputeRcut(0.3, self.Alt):
            FluxU = self.MizunoBrokenpl(0.1, 1.09, 600, 2.40, EnergyMeV)
            FluxD = self.MizunoBrokenpl(0.1, 1.09, 600, 2.40, EnergyMeV)
        elif Rcut >= self.ComputeRcut(0.5, self.Alt) and Rcut <= self.ComputeRcut(0.4, self.Alt):
            FluxU = self.MizunoBrokenpl(0.1, 1.19, 600, 2.54, EnergyMeV)
            FluxD = self.MizunoBrokenpl(0.1, 1.19, 600, 2.54, EnergyMeV)
        elif Rcut >= self.ComputeRcut(0.6, self.Alt) and Rcut <= self.ComputeRcut(0.5, self.Alt):
            FluxU = self.MizunoBrokenpl(0.1, 1.18, 400, 2.31, EnergyMeV)
            FluxD = self.MizunoBrokenpl(0.1, 1.18, 400, 2.31, EnergyMeV)
        elif Rcut >= self.ComputeRcut(0.7, self.Alt) and Rcut <= self.ComputeRcut(0.6, self.Alt):
            FluxD = self.MizunoBrokenpl(0.13, 1.1, 300, 2.25, EnergyMeV)
            FluxU = self.MizunoBrokenpl(0.13, 1.1, 300, 2.95, EnergyMeV)
        elif Rcut >= self.ComputeRcut(0.8, self.Alt) and Rcut <= self.ComputeRcut(0.7, self.Alt):
            FluxD = self.MizunoBrokenpl(0.2, 1.5, 400, 1.85, EnergyMeV)
            FluxU = self.MizunoBrokenpl(0.2, 1.5, 400, 4.16, EnergyMeV)
        elif Rcut >= self.ComputeRcut(0.9, self.Alt) and Rcut <= self.ComputeRcut(0.8, self.Alt):
            FluxD = self.MizunoCutoffpl(0.23, 0.017, 1.83, 0.177, EnergyMeV)
            FluxU = self.MizunoBrokenpl(0.23, 1.53, 400, 4.68, EnergyMeV)
        elif Rcut >= self.ComputeRcut(1.0, self.Alt) and Rcut <= self.ComputeRcut(0.9, self.Alt):
            FluxD = self.MizunoCutoffpl(0.44, 0.037, 1.98, 0.21, EnergyMeV)
            FluxU = self.MizunoBrokenpl(0.44, 2.25, 400, 3.09, EnergyMeV)

        return (FluxU+FluxD)/10**7, (FluxU)/10**7, (FluxD)/10**7

    def SecondaryProtonsUpward(self, E):
        return self.SecondaryProtons(E)[1]

    def SecondaryProtonsDownward(self, E):
        return self.SecondaryProtons(E)[2]

    def AguilarElectronPositron(self,date=2021):
        """ Read Table I from Aguilar et al. 2014 or Aguilar et al. 2018 or Aguilar et al. 2019 or Aguilar 2021 (rescaled tom 2027)
            Return a dataframe to be used by
            PrimaryElectrons and PrimaryPositrons
        """
        filename = 'Data/AguilarElectronPositron.dat'
        filename_2018 = 'Data/AguilarElectronPositron_2018.dat'
        filename_2019 = 'Data/AguilarPositron_2019.dat' #to do
        filename_2021 = 'Data/AMS21_Electrons_Positrons_rescaledTo27.dat'
        
        data_2014 = pd.read_table(filename, sep='\s+')
        data_2018 = pd.read_csv(filename_2018,sep=" ",skipinitialspace=True)
        data_2019 = pd.read_csv(filename_2019,sep=" ",skipinitialspace=True)
        data_2021 = pd.read_csv(filename_2021,sep="\s+",skipinitialspace=True)
        
        data = pd.DataFrame({})
        
        if date == 2014 :
            data["Fluxele"] = data_2014["Fluxele"]/10**10
            data['Fluxpos'] = data_2014['Fluxpos']/10**10
            data["EkeV"] = data_2014["EkeV"]
        
        elif date ==2018 :
            data["Fluxele"] = data_2018["Fluxele"]/10**10
            data['Fluxpos'] = data_2018['Fluxpos']/10**10
            data["EkeV"] = data_2018["EGeV"]*1e6 #GeV to keV
             
        elif date ==2021 :
            data["Fluxele"] = data_2021["Fluxele"]/10**10
            data['Fluxpos'] = data_2021['Fluxpos']/10**10
            data["EkeV"] = data_2021["EGeV"]*1e6 #GeV to keV

        self.PrimElecPosi = data.copy()

    def PrimaryElectrons(self, E):
        """ Table I from Aguilar et al. 2014,
            Reduction factor from equation 5 in Mizuno et al.2004
            Return a flux in ph /cm2 /s /keV /sr
        """

        try:
            self.PrimElecPosi
        except AttributeError:
            self.AguilarElectronPositron()

        data = self.PrimElecPosi

        EnergyGeV = 0.000001*np.asarray(E, dtype=float)
        E0 = 0.511/1000

        Rigidity = np.sqrt(EnergyGeV*EnergyGeV + 2*EnergyGeV*0.000511)

        f = self.log_interp1d(data["EkeV"].loc[data['Fluxele'] > 0.],
                              data["Fluxele"].loc[data['Fluxele'] > 0.])

        """ Solar modulation factor from Gleeson & Axford 1968"""
        solmodfac = ((EnergyGeV+E0)**2-E0**2)/(
                    (EnergyGeV+E0+self.solmod/1000)**2-E0**2)



        
        redfac = 1/(1+(Rigidity/self.AvGeomagCutOff)**-6.0)
            
        return f(E)*redfac   
        
        

    def PrimaryPositrons(self, E):
        """ Table I from Aguilar et al. 2014,
            Reduction factor from equation 5 in Mizuno et al.2004
            Return a flux in ph /cm2 /s /keV /sr
        """

        try:
            self.PrimElecPosi
        except AttributeError:
            self.AguilarElectronPositron()

        data = self.PrimElecPosi

        EnergyGeV = 0.000001*np.asarray(E, dtype=float)
        E0 = 0.511/1000

        Rigidity = np.sqrt(EnergyGeV*EnergyGeV + 2*EnergyGeV*0.000511)

        f = self.log_interp1d(data['EkeV'].loc[data['Fluxpos'] > 0.],
                              data['Fluxpos'].loc[data['Fluxpos'] > 0.])

        """ Solar modulation factor from Gleeson & Axford 1968"""
        solmodfac = ((EnergyGeV+E0)**2-E0**2)/(
                    (EnergyGeV+E0+self.solmod/1000)**2-E0**2)

        
        
        redfac = 1/(1+(Rigidity/self.AvGeomagCutOff)**-6.0)

        return f(E)*redfac
            
       
            

            
            
    def MizunoPl(self, f0, a, E):
        """Function describing a power-law
        """

        Flux = np.copy(np.asarray(E, dtype=float))
        mask = np.logical_and(E >= 1, E < 100)

        Flux[E < 1] = 0.
        Flux[mask] = f0*pow(E[mask]/100, -1)
        Flux[E >= 100] = f0 * pow(E[E >= 100]/100, -a)

        return Flux

    def MizunoPlhump(self, f0, a, f1, b, ec, E):
        """Function describing a power-law with hump
        """

        Flux = np.copy(np.asarray(E, dtype=float))
        mask = np.logical_and(E >= 1, E < 100)

        Flux[E < 1] = 0.
        Flux[mask] = f0*pow(E[mask]/100, -1)
        Flux[E >= 100] = (f0 * pow(E[E >= 100]/100, -a)
                          + f1 * pow(E[E >= 100]/1000, b)
                          * np.exp(-pow(E[E >= 100]/(ec*1000), b+1)))
        return Flux

    def SecondaryElectrons(self, E):
        """ Secondary electrons determinined after  section 3.4 of
            Mizuno et al. 2004
            Return a flux in ph /cm2 /s /keV /sr
        """
        EnergyMeV = 0.001*np.copy(np.asarray(E, dtype=float))

        Rcut = self.AvGeomagCutOff
       
        if Rcut >= self.ComputeRcut(0.3, self.Alt) and Rcut <= self.ComputeRcut(0.0, self.Alt):
            Flux = self.MizunoBrokenpl(0.3, 2.2, 3000, 4.0, EnergyMeV)
        elif Rcut >= self.ComputeRcut(0.6, self.Alt) and Rcut <= self.ComputeRcut(0.3, self.Alt):
            Flux = self.MizunoPl(0.3, 2.7, EnergyMeV)
        elif Rcut >= self.ComputeRcut(0.8, self.Alt) and Rcut <= self.ComputeRcut(0.6, self.Alt):
            Flux = self.MizunoPlhump(0.3, 3.3, 2/10000, 1.5, 2.3, EnergyMeV)
        elif Rcut >= self.ComputeRcut(0.9, self.Alt) and Rcut <= self.ComputeRcut(0.8, self.Alt):
            Flux = self.MizunoPlhump(0.3, 3.5, 1.6/1000, 2.0, 1.6, EnergyMeV)
        elif Rcut >= self.ComputeRcut(1.0, self.Alt) and Rcut <= self.ComputeRcut(0.9, self.Alt):
            Flux = self.MizunoPl(0.3, 2.5, EnergyMeV)

        return Flux/10**7

    def SecondaryPositrons(self, E):
        """ Secondary positrons determinined after section 3.4 of
            Mizuno et al. 2004
            Return a flux in ph /cm2 /s /keV /sr
        """
        EnergyMeV = 0.001*np.copy(np.asarray(E, dtype=float))

        Rcut = self.AvGeomagCutOff
        
        if Rcut >= self.ComputeRcut(0.3, self.Alt) and Rcut <= self.ComputeRcut(0.0, self.Alt):
            Flux = self.MizunoBrokenpl(0.3, 2.2, 3000, 4.0, EnergyMeV)
            ratio = 3.3
        elif Rcut >= self.ComputeRcut(0.6, self.Alt) and Rcut <= self.ComputeRcut(0.3, self.Alt):
            Flux = self.MizunoPl(0.3, 2.7, EnergyMeV)
            ratio = 1.66
        elif Rcut >= self.ComputeRcut(0.8, self.Alt) and Rcut <= self.ComputeRcut(0.6, self.Alt):
            Flux = self.MizunoPlhump(0.3, 3.3, 2/10000, 1.5, 2.3, EnergyMeV)
            ratio = 1.0
        elif Rcut >= self.ComputeRcut(0.9, self.Alt) and Rcut <= self.ComputeRcut(0.8, self.Alt):
            Flux = self.MizunoPlhump(0.3, 3.5, 1.6/1000, 2.0, 1.6, EnergyMeV)
            ratio = 1.0
        elif Rcut >= self.ComputeRcut(1.0, self.Alt) and Rcut <= self.ComputeRcut(0.9, self.Alt):
            Flux = self.MizunoPl(0.3, 2.5, EnergyMeV)
            ratio = 1.0

        return ratio*Flux/10**7

    def PrimaryProtons(self, E):
        """ Read Table from Aguilar et al. 2015,
            Rigidity in GV and Flux in /m2 /sr /s /GV
            Return a flux in ph /cm2 /s /keV /sr
        """
        filename = 'Data/AguilarProton.dat'
        data = pd.read_table(filename, sep='\s+')

        E0 = ((m_p * c**2).to('GeV')).value

        data["Flux"] = data["Flux"]*data['RigidityGV']
        data['RigidityGV'] = (np.sqrt(E0**2+data['RigidityGV']**2)-E0)*10**6
        data["Flux"] = data["Flux"]/(data['RigidityGV'])/10**4

        EnergyGeV = 0.000001*np.asarray(E, dtype=float)

        Rigidity = np.sqrt(EnergyGeV*EnergyGeV + 2*EnergyGeV*E0)

        f = self.log_interp1d(data['RigidityGV'].loc[data['Flux'] > 0.],
                              data['Flux'].loc[data['Flux'] > 0.])

        

        """ Solar modulation factor from Gleeson & Axford 1968"""
        solmodfac = ((EnergyGeV+E0)**2-E0**2)/(
                    (EnergyGeV+E0+self.solmod/1000)**2-E0**2)
                    
        
        """ Geomagnetic modulation factor from Mizuno et al. 2004"""
        redfac = 1/(1+(Rigidity/self.AvGeomagCutOff)**-12.0)            

        return f(E)*redfac
            
        
            
    def PrimaryProtons_HelMod(self, E):
        """ Read Table from HelMod public online tool,
            Rigidity in GV and Flux in /m2 /sr /s /GV
            Return a flux in ph /cm2 /s /keV /sr
        """
        filename = './Data/HelmodProton_Mar27_Jun27.dat'
        data = pd.read_table(filename, sep='\s+')
        
        E0 = ((m_p * c**2).to('GeV')).value

        data["Flux"] = data["Flux"]*data['RigidityGV']
        data['RigidityGV'] = (np.sqrt(E0**2+data['RigidityGV']**2)-E0)*10**6
        data["Flux"] = data["Flux"]/(data['RigidityGV'])/10**4

        EnergyGeV = 0.000001*np.asarray(E, dtype=float)

        Rigidity = np.sqrt(EnergyGeV*EnergyGeV + 2*EnergyGeV*E0)

        f = self.log_interp1d(data['RigidityGV'].loc[data['Flux'] > 0.],
                              data['Flux'].loc[data['Flux'] > 0.])

        """ Geomagnetic modulation factor from Mizuno et al. 2004"""
        redfac = 1/(1+(Rigidity/self.AvGeomagCutOff)**-12.0)

        return f(E)*redfac


    def PrimaryAlphas(self, E):
        """ Read Table from Aguilar et al. 2015b,
            Rigidity in GV and Flux in /m2 /sr /s /GV
            Return a flux in ph /cm2 /s /keV /sr
        """
        filename = 'Data/AguilarAlphas.dat'
        data = pd.read_table(filename, sep='\s+')

        E0 = 2*((m_p * c**2 + m_n * c**2).to('GeV')).value

        data["Flux"] = data["Flux"]*data['RigidityGV']
        data['RigidityGV'] = 4*(np.sqrt(E0**2+(data['RigidityGV']/2)**2)-E0)*10**6

        data["Flux"] = data["Flux"]/(data['RigidityGV'])/10**4

        EnergyGeV = 0.000001*np.asarray(E, dtype=float)

        Rigidity = np.sqrt(EnergyGeV*EnergyGeV + 2*EnergyGeV*E0)/2.

        f = self.log_interp1d(data['RigidityGV'].loc[data['Flux'] > 0.],
                              data['Flux'].loc[data['Flux'] > 0.])

        

        """ Solar modulation factor from Gleeson & Axford 1968"""
        solmodfac = ((EnergyGeV+E0)**2-E0**2)/(
                    (EnergyGeV+E0+2*self.solmod/1000)**2-E0**2)
                
        """ Geomagnetic modulation factor from Mizuno et al. 2004"""
        redfac = 1/(1+(Rigidity/self.AvGeomagCutOff)**-12.0)
                
        return f(E)*redfac

        
    def PrimaryAlphas_HelMod(self, E):
        """ Read Table from HelMod public online tool
            Rigidity in GV and Flux in /m2 /sr /s /GV
            Return a flux in ph /cm2 /s /keV /sr
        """
        filename = './Data/HelmodHelium_Mar27_Jun27.dat'
        data = pd.read_table(filename, sep='\s+')

        E0 = 2*((m_p * c**2 + m_n * c**2).to('GeV')).value

        data["Flux"] = data["Flux"]*data['RigidityGV']
        data['RigidityGV'] = 4*(np.sqrt(E0**2+(data['RigidityGV']/2)**2)-E0)*10**6

        data["Flux"] = data["Flux"]/(data['RigidityGV'])/10**4

        EnergyGeV = 0.000001*np.asarray(E, dtype=float)

        Rigidity = np.sqrt(EnergyGeV*EnergyGeV + 2*EnergyGeV*E0)/2.

        f = self.log_interp1d(data['RigidityGV'].loc[data['Flux'] > 0.],
                              data['Flux'].loc[data['Flux'] > 0.])

        """ Geomagnetic modulation factor from Mizuno et al. 2004"""
        redfac = 1/(1+(Rigidity/self.AvGeomagCutOff)**-12.0)

        return f(E)*redfac
