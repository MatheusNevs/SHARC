# -*- coding: utf-8 -*-
"""
Auxiliary functions for ITU-R P.2001-4 Implementation
Adapted from the reference Python implementation
"""

import numpy as np
from typing import Tuple
import os
import warnings


class P2001Aux:
    """
    Auxiliary class containing all helper functions for P.2001-4 calculations
    """
    
    def __init__(self, digital_maps_path: str = None):
        """
        Initialize with constants and digital maps
        
        Parameters
        ----------
        digital_maps_path : str, optional
            Path to P2001.npz file. If None, looks in Dataset folder.
            If not found, uses typical values as fallback.
        """
        self.c0 = 2.998e8  # Speed of light (m/s)
        self.Re = 6371.0   # Average Earth radius (km)
        
        # Try to load digital maps
        self.DigitalMaps = {}
        self.maps_loaded = False
        
        if digital_maps_path is None:
            # Default path: sharc/propagation/Dataset/P2001.npz
            base_dir = os.path.dirname(__file__)
            digital_maps_path = os.path.join(base_dir, 'Dataset', 'P2001.npz')
        
        if os.path.exists(digital_maps_path):
            try:
                with np.load(digital_maps_path) as data:
                    for key in data.files:
                        self.DigitalMaps[key] = data[key]
                self.maps_loaded = True
                print(f"✅ Loaded P.2001 digital maps from: {digital_maps_path}")
            except Exception as e:
                warnings.warn(f"Failed to load digital maps: {e}. Using typical values.")
        else:
            warnings.warn(f"Digital maps not found at: {digital_maps_path}. Using typical values.")
    
    def interp2(self, map_name: str, lon: float, lat: float, long_spacing: float, lat_spacing: float, default_value: float = None) -> float:
        """
        Interpolate value from digital map at given coordinates
        
        Parameters
        ----------
        map_name : str
            Name of map in DigitalMaps (e.g., 'DN_Median', 'surfwv_50_fixed')
        lon : float
            Longitude in degrees (-180 to +180)
        lat : float
            Latitude in degrees (-90 to +90)
        long_spacing: float
            Longitude spacing in degrees
        lat_spacing: float
            Latitude spacing in degrees
        default_value : float, optional
            Value to return if maps not loaded or interpolation fails
            
        Returns
        -------
        float
            Interpolated value from map, or default_value if not available
        """
        if not self.maps_loaded or map_name not in self.DigitalMaps:
            if default_value is None:
                warnings.warn(f"Map '{map_name}' not available. Using fallback value.")
                return 0.0
            return default_value
        
        try:
            # Digital maps may use 1.5° or 1.125° spacing
            # Latitude: 90 to -90 (index 0 to 120)
            # Longitude: 0 to 360 (index 0 to 239), but input is -180 to +180
            
            # Convert longitude from -180..+180 to 0..360
            lon_360 = lon if lon >= 0 else lon + 360
            
            # Map coordinates to indices
            lat_idx = (90 - lat) / lat_spacing
            lon_idx = lon_360 / long_spacing
            
            # Get map data
            map_data = self.DigitalMaps[map_name]
            
            # Bilinear interpolation
            i0 = int(np.floor(lat_idx))
            i1 = min(i0 + 1, map_data.shape[0] - 1)
            j0 = int(np.floor(lon_idx)) % map_data.shape[1]
            j1 = (j0 + 1) % map_data.shape[1]
            
            # Fractional parts
            di = lat_idx - i0
            dj = lon_idx - int(np.floor(lon_idx))
            
            # Bilinear interpolation
            val = ((1 - di) * (1 - dj) * map_data[i0, j0] +
                   (1 - di) * dj * map_data[i0, j1] +
                   di * (1 - dj) * map_data[i1, j0] +
                   di * dj * map_data[i1, j1])
            
            return float(val)
            
        except Exception as e:
            warnings.warn(f"Interpolation failed for {map_name}: {e}")
            return default_value if default_value is not None else 0.0
        
    def bt_loss(self, d: np.ndarray, h: np.ndarray, z: np.ndarray, GHz: float, 
                Tpc: float, Phire: float, Phirn: float, Phite: float, Phitn: float,
                Hrg: float, Htg: float, Grx: float, Gtx: float, FlagVP: int) -> float:
        """
        Calculate basic transmission loss according to ITU-R P.2001-4
        
        Parameters
        ----------
        d : np.ndarray
            Distance from transmitter of i-th profile point (km)
        h : np.ndarray
            Height of i-th profile point (m amsl)
        z : np.ndarray
            Zone code (1=Sea, 3=Coastal Land, 4=Inland)
        GHz : float
            Frequency (GHz)
        Tpc : float
            Percentage of average year for which predicted loss is not exceeded
        Phire, Phirn : float
            Receiver longitude (deg, +E), latitude (deg, +N)
        Phite, Phitn : float
            Transmitter longitude (deg, +E), latitude (deg, +N)
        Hrg, Htg : float
            Receiving/Transmitting antenna height above ground (m)
        Grx, Gtx : float
            Receiving/Transmitting antenna gain (dBi)
        FlagVP : int
            Polarization: 1=vertical, 0=horizontal
            
        Returns
        -------
        float
            Basic transmission loss (dB)
        """
        
        # 3.1 Limited percentage time
        Tpcp = Tpc + 0.00001 * (50 - Tpc) / 50  # Eq (3.1.1)
        Tpcq = 100 - Tpcp  # Eq (3.1.2)
        
        # Validate inputs
        if not np.all(np.diff(d) >= 0):
            raise ValueError("Distance array d must be in ascending order")
        
        if d[0] > 0.0:
            raise ValueError(f"First distance d[0]={d[0]} must be zero")
            
        if len(d) <= 10:
            raise ValueError("Path profile must have more than 10 points")
            
        if not np.all((z == 1) | (z == 3) | (z == 4)):
            raise ValueError("Zone codes must be 1, 3, or 4")
            
        if not (0 < Tpc < 100):
            raise ValueError("Time percentage must be in range (0, 100)")
            
        if Htg <= 0 or Hrg <= 0:
            raise ValueError("Antenna heights must be positive")
            
        if FlagVP not in (0, 1):
            raise ValueError("Polarization must be 0 (horizontal) or 1 (vertical)")
        
        # 3.2 Path length, intermediate points, and fraction over sea
        dt = d[-1]  # Eq (3.2.1)
        
        # Calculate mid-point
        dpnt = 0.5 * dt
        Phime, Phimn, _, _ = self.great_circle_path(
            Phire, Phite, Phirn, Phitn, self.Re, dpnt
        )
        Phime1 = Phime if Phime >= 0 else Phime1 = Phime + 360
        
        # Ground height at mid-point (3.2.2)
        n = len(d)
        if n % 2 == 1:  # odd
            mp = int(0.5 * (n + 1) - 1)
            Hmid = h[mp]
        else:  # even
            mp = int(0.5 * n - 1)
            Hmid = 0.5 * (h[mp] + h[mp + 1])
        
        # Fraction of path over sea
        omega = self.path_fraction(d, z, 1)
        
        # 3.3 Antenna altitudes and path inclinations
        Hts = Htg + h[0]  # Tx height amsl
        Hrs = Hrg + h[-1]  # Rx height amsl
        
        Hhi = max(Hts, Hrs)  # Higher antenna
        Hlo = min(Hts, Hrs)  # Lower antenna
        
        Sp = (Hhi - Hlo) / dt  # Path inclination (3.3.3)
        
        # 3.6 Wavelength
        Wave = 1e-9 * self.c0 / GHz  # Eq (3.6.1)
        
        # 3.5 Effective Earth-radius geometry
        # Get refractivity gradient from digital maps or use typical value
        Nd1km50 = self.interp2('DN_Median', Phime1, Phimn, 1.5, 1.5, default_value=-40.0)
        Reff50 = 157.0 * self.Re / (157.0 + Nd1km50)  # Eq (3.5.1)
        
        # For p% time
        Nd1kmp = Nd1km50 + self.interp2('DN_SupSlope', Phime1, Phimn, 1.5, 1.5) * np.log10(0.02 * Tpcp) if Tpcp < 50 \
            else Nd1km50 - self.interp2('DN_SubSlope', Phime1, Phimn, 1.5, 1.5) * np.log10(0.02 * Tpcp)
        Cp = (157.0 + Nd1kmp) / (157.0 * self.Re)  # Eq (3.5.2)
        Reffp = 1.0 / Cp if Cp > 1e-6 else 1e6  # Eq (3.5.3)
        
        Thetae = dt / Reff50  # Eq (3.5.4)
        
        # 3.7-3.8 Path classification and effective heights
        (theta_t, theta_r, theta_tpos, theta_rpos, dlt, dlr, lt, lr,
         hstip, hsrip, hstipa, hsripa, htea, hrea, mses, hm, 
         hst, hsr, htep, hrep, FlagLos50) = self.smooth_earth_heights(
            d, h, Hts, Hrs, Reff50, Wave
        )
        
        # ====================================================================
        # SUB-MODEL 1: Line-of-sight and diffraction
        # ====================================================================
        
        # Diffraction loss
        Ld_pol, _, _, _, _, _, _, _ = self.dl_p(
            d, h, Hts, Hrs, htep, hrep, GHz, omega, Reffp, Cp
        )
        Ld = Ld_pol[FlagVP]
        
        # Free-space loss
        dfs = np.sqrt(dt**2 + ((Hts - Hrs) / 1000.0)**2)
        Lbfs = self.tl_free_space(GHz, dfs)  # Eq (3.11.1)
        
        # Gaseous absorption on surface paths (3.10)
        Aosur, Awsur, Awrsur, _, _, _, _ = self.gaseous_abs_surface(
            Phime, Phimn, Hmid, Hts, Hrs, dt, GHz
        )
        Agsur = Aosur + Awsur  # Eq (3.10.1)
        
        # Multipath activity (Attachment B.2)
        Nd65m1 = self.interp2('dndz_01', Phime1, Phimn, 1.5, 1.5)
        Q0ca = self.multi_path_activity(
            GHz, dt, Hts, Hrs, dlt, dlr, h[lt], h[lr], Hlo,
            theta_t, theta_r, Sp, Nd65m1, Phimn, FlagLos50
        )
        
        # Precipitation fade preliminary calculations (Attachment C.2) - Simplified
        Q0ra = 2.0  # Typical percentage for rain (simplified)
        Fwvr = 0.5  # Simplified water vapor factor
        A1 = 0.0    # Simplified: no additional attenuation for now

        # a, b, c, dr, Q0ra, Fwvr, kmod, alpha_mod, Gm, Pm, flagrain = self.precipitation_fade_initial(GHz, Tpcq, Phimn, Phime, Hlo, Hhi, dt, FlagVP)
        
        # Sub-model 1 basic transmission loss (4.1.4)
        Lbm1 = Lbfs + Ld + A1 + Fwvr * (Awrsur - Awsur) + Agsur
        
        # ====================================================================
        # SUB-MODEL 2: Anomalous propagation (ducting/layer-reflection)
        # ====================================================================
        
        Lba, _, _, _, _, _, _, _ = self.tl_anomalous_reflection(
            GHz, d, z, Hts, Hrs, htea, hrea, hm, theta_t, theta_r,
            dlt, dlr, Phimn, omega, Reff50, Tpcp, Tpcq
        )
        Lbm2 = Lba + Agsur  # (4.2.1)
        
        # ====================================================================
        # SUB-MODEL 3: Troposcatter propagation
        # ====================================================================
        
        # Tropospheric path segments (3.9)
        Dtcv, Drcv, Phicve, Phicvn, Hcv, Phitcve, Phitcvn, Phircve, Phircvn = \
            self.tropospheric_path(
                dt, Hts, Hrs, Thetae, theta_tpos, theta_rpos, Reff50,
                Phire, Phite, Phirn, Phitn, self.Re
            )
        
        # Troposcatter basic transmission loss (Attachment E)
        Lbs, _, _ = self.tl_troposcatter(
            GHz, dt, theta_t, theta_r, Thetae, Phicvn, Phicve,
            Phitn, Phite, Phirn, Phire, Gtx, Grx, Reff50, Tpcp
        )
        
        # Limit troposcatter loss
        Lbs = max(Lbs, Lbfs)

        # a, b, c, dr, Q0ra, Fwvrtx, kmod, alpha_mod, Gm, Pm, flagrain = precipitation_fade_initial(GHz, Tpcq, Phitcvn, Phitcve, Hts, Hcv, Dtcv, FlagVP)

        
        # Gaseous absorption for troposcatter (Attachment F.3)
        Aos, Aws, Awrs, _, _, _, _, _, _, _, _ = self.gaseous_abs_tropo(
            Phite, Phitn, Phire, Phirn, h[0], h[-1],
            theta_tpos, theta_rpos, Dtcv, Drcv, GHz
        )
        Ags = Aos + Aws
        
        # Simplified attenuation for troposcatter paths
        A2 = 0.0  # Simplified
        Fwvrtx = 0.5  # Simplified
        Fwvrrx = 0.5  # Simplified
        
        # Sub-model 3 basic transmission loss (4.3.8)
        Lbm3 = Lbs + A2 + 0.5*(Fwvrtx + Fwvrrx)*(Awrs - Aws) + Ags
        
        # ====================================================================
        # SUB-MODEL 4: Sporadic-E propagation
        # ====================================================================
        
        Lbm4, _, _, _, _, _, _, _, _, _, _, _, _, _, _ = self.tl_sporadic_e(
            GHz, dt, theta_t, theta_r, Phimn, Phime, Phitn, Phite,
            Phirn, Phire, dlt, dlr, Reff50, self.Re, Tpcp
        )
        
        # ====================================================================
        # COMBINING SUB-MODEL RESULTS (Section 5)
        # ====================================================================
        
        # 5.1 Combining sub-models 1 and 2
        Lm = min(Lbm1, Lbm2)
        Lbm12 = Lm - 10.0*np.log10(10**(-0.1*(Lbm1 - Lm)) + 10**(-0.1*(Lbm2 - Lm)))
        
        # 5.2 Combining sub-models 1+2, 3, and 4
        Lm = min([Lbm12, Lbm3, Lbm4])
        Lb = Lm - 5.0*np.log10(10**(-0.2*(Lbm12 - Lm)) + 10**(-0.2*(Lbm3 - Lm)) + 10**(-0.2*(Lbm4 - Lm)))
        
        return Lb
    
    def great_circle_path(self, Phire: float, Phite: float, Phirn: float, 
                          Phitn: float, Re: float, dpnt: float) -> Tuple:
        """
        Calculate great-circle path intermediate point
        
        Returns
        -------
        Phipnte : float
            Longitude of intermediate point (deg)
        Phipntn : float
            Latitude of intermediate point (deg)
        Bt2r : float
            Bearing from Tx to Rx (deg)
        dgc : float
            Great-circle path length (km)
        """
        Dlon = Phire - Phite
        
        r = self.sind(Phitn) * self.sind(Phirn) + \
            self.cosd(Phitn) * self.cosd(Phirn) * self.cosd(Dlon)
        
        Phid = np.arccos(np.clip(r, -1, 1))
        dgc = Phid * Re
        
        x1 = self.sind(Phirn) - r * self.sind(Phitn)
        y1 = self.cosd(Phitn) * self.cosd(Phirn) * self.sind(Dlon)
        
        if abs(x1) < 1e-9 and abs(y1) < 1e-9:
            Bt2r = Phire
        else:
            Bt2r = self.atan2d(y1, x1)
        
        Phipnt = dpnt / Re
        
        s = self.sind(Phitn) * np.cos(Phipnt) + \
            self.cosd(Phitn) * np.sin(Phipnt) * self.cosd(Bt2r)
        
        Phipntn = np.arcsin(np.clip(s, -1, 1)) * 180.0 / np.pi
        
        x2 = np.cos(Phipnt) - s * self.sind(Phitn)
        y2 = self.cosd(Phitn) * np.sin(Phipnt) * self.sind(Bt2r)
        
        if x2 < 1e-9 and y2 < 1e-9:
            Phipnte = Bt2r
        else:
            Phipnte = Phite + self.atan2d(y2, x2)
        
        return Phipnte, Phipntn, Bt2r, dgc
    
    def path_fraction(self, d: np.ndarray, zone: np.ndarray, zone_r: int) -> float:
        """
        Calculate fraction of path in a given zone
        
        Parameters
        ----------
        d : np.ndarray
            Distances (km)
        zone : np.ndarray
            Zone codes
        zone_r : int
            Reference zone
            
        Returns
        -------
        float
            Fraction of path in zone_r
        """
        dm = 0.0
        start, stop = self.find_intervals(zone == zone_r)
        
        n = len(start)
        nmax = len(d)
        
        for i in range(n):
            delta = 0.0
            if stop[i] < nmax - 1:
                delta += (d[stop[i] + 1] - d[stop[i]]) / 2.0
            if start[i] > 0:
                delta += (d[start[i]] - d[start[i] - 1]) / 2.0
            
            dm = max(d[stop[i]] - d[start[i]] + delta, dm)
        
        omega = dm / (d[-1] - d[0]) if d[-1] > d[0] else 0.0
        return omega
    
    def find_intervals(self, series: np.ndarray) -> Tuple:
        """Find intervals of consecutive True values"""
        series_int = series.astype(int)
        
        if np.max(series_int) == 1:
            k1 = np.where(np.diff(np.append(0, series_int)) == 1)[0]
            k2 = np.where(np.diff(np.append(series_int, 0)) == -1)[0]
        else:
            k1 = np.array([])
            k2 = np.array([])
        
        return k1, k2
    
    def smooth_earth_heights(self, d: np.ndarray, h: np.ndarray, hts: float,
                            hrs: float, ae: float, lam: float) -> Tuple:
        """
        Calculate smooth-Earth effective antenna heights
        
        Returns a tuple with many terrain parameters needed for diffraction calculations
        """
        n = len(d)
        dtot = d[-1]
        
        # Find horizon angles (simplified)
        ii = np.arange(1, n - 1)
        
        theta = (h[ii] - hts) / d[ii] - 500 * d[ii] / ae
        theta_tim = np.max(theta) if len(theta) > 0 else 0
        
        theta_tr = (hrs - hts) / dtot - 500 * dtot / ae
        
        if theta_tim < theta_tr:
            FlagLos50 = 1
            nu = (h[ii] + 500 * d[ii] * (dtot - d[ii]) / ae - 
                  (hts * (dtot - d[ii]) + hrs * d[ii]) / dtot) * \
                 np.sqrt(0.002 * dtot / (lam * d[ii] * (dtot - d[ii])))
            
            numax = np.max(nu) if len(nu) > 0 else 0
            kindex = np.where(nu == numax)[0]
            lt = kindex[-1] + 1 if len(kindex) > 0 else 1
            dlt = d[lt]
            dlr = dtot - dlt
            lr = lt
            
            theta_t = theta_tr
            theta_r = -theta_tr - 1000 * dtot / ae
        else:
            FlagLos50 = 0
            theta_ti = (h[ii] - hts) / d[ii] - 500 * d[ii] / ae
            kindex = np.where(theta_ti == theta_tim)[0]
            lt = kindex[-1] + 1 if len(kindex) > 0 else 1
            dlt = d[lt]
            theta_t = theta_tim
            
            theta_ri = (h[ii] - hrs) / (dtot - d[ii]) - 500 * (dtot - d[ii]) / ae
            theta_rim = np.max(theta_ri) if len(theta_ri) > 0 else 0
            kindex = np.where(theta_ri == theta_rim)[0]
            lr = kindex[-1] + 1 if len(kindex) > 0 else n - 1
            dlr = dtot - d[lr]
            theta_r = theta_rim
        
        theta_tpos = max(theta_t, 0)
        theta_rpos = max(theta_r, 0)
        
        # Calculate smooth surface heights (simplified)
        hstip = h[0]
        hsrip = h[-1]
        hstipa = min(hstip, h[0])
        hsripa = min(hsrip, h[-1])
        
        mses = (hsripa - hstipa) / dtot
        
        htea = hts - hstipa
        hrea = hrs - hsripa
        
        # Path roughness (simplified)
        if lt < lr and lr < n:
            hm = np.max(h[lt:lr+1] - (hstipa + mses * d[lt:lr+1]))
        else:
            hm = 0
        
        hst = h[0]
        hsr = h[-1]
        htep = hts - hst
        hrep = hrs - hsr
        
        return (theta_t, theta_r, theta_tpos, theta_rpos, dlt, dlr, lt, lr,
                hstip, hsrip, hstipa, hsripa, htea, hrea, mses, hm,
                hst, hsr, htep, hrep, FlagLos50)
    
    def dl_p(self, d: np.ndarray, h: np.ndarray, hts: float, hrs: float,
             hte: float, hre: float, f: float, omega: float, 
             ap: float, Cp: float) -> Tuple:
        """
        Calculate diffraction loss for p% time
        
        Returns tuple with diffraction losses for both polarizations
        """
        dtot = d[-1]
        
        # Spherical-Earth diffraction
        Ldsph = self.dl_se(dtot, hte, hre, ap, f, omega)
        
        # Bullington diffraction (actual path)
        Ldba, Ldbka, FlagLospa = self.dl_bull_actual(d, h, hts, hrs, Cp, f)
        
        # Bullington diffraction (smooth path)
        Ldbs, Ldbks, FlagLosps = self.dl_bull_smooth(d, h, hte, hre, ap, f)
        
        # Total diffraction loss
        Ld = np.zeros(2)
        Ld[0] = Ldba + max(Ldsph[0] - Ldbs, 0)
        Ld[1] = Ldba + max(Ldsph[1] - Ldbs, 0)
        
        return Ld, Ldsph, Ldba, Ldbs, Ldbka, Ldbks, FlagLospa, FlagLosps
    
    def dl_se(self, d: float, hte: float, hre: float, ap: float, 
              f: float, omega: float) -> np.ndarray:
        """Spherical-Earth diffraction loss"""
        lam = 1e-9 * self.c0 / f
        
        dlos = np.sqrt(2 * ap) * (np.sqrt(0.001 * hte) + np.sqrt(0.001 * hre))
        
        if d >= dlos:
            return self.dl_se_ft(d, hte, hre, ap, f, omega)
        
        # Below LoS distance
        c = (hte - hre) / (hte + hre)
        m = 250 * d * d / (ap * (hte + hre))
        
        b = 2 * np.sqrt((m + 1) / (3 * m)) * \
            np.cos(np.pi / 3.0 + (1.0 / 3.0) * np.arccos(
                3.0 * c / 2.0 * np.sqrt(3.0 * m / ((m + 1.0)**3))))
        
        dse1 = d / 2.0 * (1.0 + b)
        dse2 = d - dse1
        
        hse = ((hte - 500 * dse1 * dse1 / ap) * dse2 + 
               (hre - 500 * dse2 * dse2 / ap) * dse1) / d
        
        hreq = 17.456 * np.sqrt(dse1 * dse2 * lam / d)
        
        if hse > hreq:
            return np.zeros(2)
        
        aem = 500 * (d / (np.sqrt(hte) + np.sqrt(hre)))**2
        Ldft = self.dl_se_ft(d, hte, hre, aem, f, omega)
        
        Ldsph = np.zeros(2)
        Ldsph[0] = 0 if Ldft[0] < 0 else (1 - hse / hreq) * Ldft[0]
        Ldsph[1] = 0 if Ldft[1] < 0 else (1 - hse / hreq) * Ldft[1]
        
        return Ldsph
    
    def dl_se_ft(self, d: float, hte: float, hre: float, adft: float,
                 f: float, omega: float) -> np.ndarray:
        """First-term spherical-Earth diffraction"""
        # Over land
        Ldft_land = self.dl_se_ft_inner(22, 0.003, d, hte, hre, adft, f)
        # Over sea
        Ldft_sea = self.dl_se_ft_inner(80, 5, d, hte, hre, adft, f)
        
        return omega * Ldft_sea + (1 - omega) * Ldft_land
    
    def dl_se_ft_inner(self, epsr: float, sigma: float, d: float,
                       hte: float, hre: float, adft: float, f: float) -> np.ndarray:
        """Inner calculation for first-term diffraction"""
        K = np.zeros(2)
        K[0] = 0.036 * (adft * f)**(-1.0/3.0) * \
               ((epsr - 1)**2 + (18 * sigma / f)**2)**(-1.0/4.0)
        K[1] = K[0] * (epsr**2 + (18.0 * sigma / f)**2)**(1.0/2.0)
        
        beta_dft = (1 + 1.6 * K**2 + 0.67 * K**4) / \
                   (1 + 4.5 * K**2 + 1.53 * K**4)
        
        X = 21.88 * beta_dft * (f / adft**2)**(1.0/3.0) * d
        Yt = 0.9575 * beta_dft * (f**2.0 / adft)**(1.0/3.0) * hte
        Yr = 0.9575 * beta_dft * (f**2.0 / adft)**(1.0/3.0) * hre
        
        Fx = np.zeros(2)
        GYt = np.zeros(2)
        GYr = np.zeros(2)
        
        for ii in range(2):
            if X[ii] >= 1.6:
                Fx[ii] = 11 + 10 * np.log10(X[ii]) - 17.6 * X[ii]
            else:
                Fx[ii] = -20 * np.log10(X[ii]) - 5.6488 * X[ii]**1.425
        
        Bt = beta_dft * Yt
        Br = beta_dft * Yr
        
        for ii in range(2):
            if Bt[ii] > 2:
                GYt[ii] = 17.6 * (Bt[ii] - 1.1)**0.5 - \
                         5 * np.log10(Bt[ii] - 1.1) - 8
            else:
                GYt[ii] = 20 * np.log10(Bt[ii] + 0.1 * Bt[ii]**3)
            
            if Br[ii] > 2:
                GYr[ii] = 17.6 * (Br[ii] - 1.1)**0.5 - \
                         5 * np.log10(Br[ii] - 1.1) - 8
            else:
                GYr[ii] = 20 * np.log10(Br[ii] + 0.1 * Br[ii]**3)
            
            GYr[ii] = max(GYr[ii], 2 + 20 * np.log10(K[ii]))
            GYt[ii] = max(GYt[ii], 2 + 20 * np.log10(K[ii]))
        
        return -Fx - GYt - GYr
    
    def dl_bull_actual(self, d: np.ndarray, h: np.ndarray, hts: float,
                       hrs: float, Cp: float, f: float) -> Tuple:
        """Bullington diffraction for actual path"""
        lam = 1e-9 * self.c0 / f
        dtot = d[-1]
        
        di = d[1:-1]
        hi = h[1:-1]
        
        if len(di) == 0:
            return 0, 0, 1
        
        Stim = np.max((hi + 500 * Cp * di * (dtot - di) - hts) / di)
        Str = (hrs - hts) / dtot
        
        if Stim < Str:
            FlagLospa = 1
            nu = (hi + 500 * Cp * di * (dtot - di) - 
                  (hts * (dtot - di) + hrs * di) / dtot) * \
                 np.sqrt(0.002 * dtot / (lam * di * (dtot - di)))
            numax = np.max(nu)
            Ldbka = self.dl_knife_edge(numax)
        else:
            FlagLospa = 0
            Srim = np.max((hi + 500 * Cp * di * (dtot - di) - hrs) / (dtot - di))
            dbp = (hrs - hts + Srim * dtot) / (Stim + Srim)
            nub = (hts + Stim * dbp - (hts * (dtot - dbp) + hrs * dbp) / dtot) * \
                  np.sqrt(0.002 * dtot / (lam * dbp * (dtot - dbp)))
            Ldbka = self.dl_knife_edge(nub)
        
        Ldba = Ldbka + (1 - np.exp(-Ldbka / 6.0)) * (10 + 0.02 * dtot)
        return Ldba, Ldbka, FlagLospa
    
    def dl_bull_smooth(self, d: np.ndarray, h: np.ndarray, htep: float,
                       hrep: float, ap: float, f: float) -> Tuple:
        """Bullington diffraction for smooth path"""
        lam = 1e-9 * self.c0 / f
        dtot = d[-1]
        
        di = d[1:-1]
        
        if len(di) == 0:
            return 0, 0, 1
        
        Stim = np.max(500 * (dtot - di) / ap - htep / di)
        Str = (hrep - htep) / dtot
        
        if Stim < Str:
            FlagLosps = 1
            nu = (500 * di * (dtot - di) / ap - 
                  (htep * (dtot - di) + hrep * di) / dtot) * \
                 np.sqrt(0.002 * dtot / (lam * di * (dtot - di)))
            numax = np.max(nu)
            Ldbks = self.dl_knife_edge(numax)
        else:
            FlagLosps = 0
            Srim = np.max(500 * di / ap - hrep / (dtot - di))
            dbp = (hrep - htep + Srim * dtot) / (Stim + Srim)
            nub = (htep + Stim * dbp - 
                   (htep * (dtot - dbp) + hrep * dbp) / dtot) * \
                  np.sqrt(0.002 * dtot / (lam * dbp * (dtot - dbp)))
            Ldbks = self.dl_knife_edge(nub)
        
        Ldbs = Ldbks + (1 - np.exp(-Ldbks / 6.0)) * (10 + 0.02 * dtot)
        return Ldbs, Ldbks, FlagLosps
    
    def dl_knife_edge(self, nu: float) -> float:
        """Knife-edge diffraction loss"""
        if nu > -0.78:
            return 6.9 + 20 * np.log10(np.sqrt((nu - 0.1)**2 + 1) + nu - 0.1)
        return 0
    
    def tl_free_space(self, f: float, d: float) -> float:
        """Free-space basic transmission loss"""
        return 92.4 + 20 * np.log10(f) + 20 * np.log10(d)
    
    # Trigonometric helpers
    def sind(self, x: float) -> float:
        return np.sin(np.deg2rad(x))
    
    def cosd(self, x: float) -> float:
        return np.cos(np.deg2rad(x))
    
    def atan2d(self, y: float, x: float) -> float:
        return np.rad2deg(np.arctan2(y, x))
    
    # ============================================================================
    # SUB-MODEL 2: ANOMALOUS PROPAGATION (DUCTING/LAYER-REFLECTION)
    # ============================================================================
    
    def tl_anomalous_reflection(self, f: float, d: np.ndarray, z: np.ndarray, 
                                hts: float, hrs: float, htea: float, hrea: float,
                                hm: float, thetat: float, thetar: float, 
                                dlt: float, dlr: float, phimn: float, omega: float,
                                ae: float, p: float, q: float) -> Tuple:
        """
        Basic transmission loss associated with anomalous propagation (Attachment D)
        
        Parameters
        ----------
        f : float
            Frequency (GHz)
        d : np.ndarray
            Vector of distances (km)
        z : np.ndarray
            Zone codes
        hts, hrs : float
            Tx/Rx antenna heights (m amsl)
        htea, hrea : float
            Effective antenna heights (m amsl)
        hm : float
            Path roughness parameter (m)
        thetat, thetar : float
            Horizon elevation angles (mrad)
        dlt, dlr : float
            Horizon distances (km)
        phimn : float
            Mid-point latitude (deg)
        omega : float
            Fraction of path over sea
        ae : float
            Effective Earth radius (km)
        p, q : float
            Percentage of time (p: not exceeded, q: exceeded)
            
        Returns
        -------
        Tuple containing Lba, Aat, Aad, Aac, dct, dcr, dtm, dlm
        """
        dt = d[-1] - d[0]
        
        # D.1 Characterize radio-climatic zones
        zoner = 34  # Coastal + inland
        dtm = self.longest_cont_dist(d, z, zoner)
        
        zoner = 4  # Inland only
        dlm = self.longest_cont_dist(d, z, zoner)
        
        # D.2 Point incidence of ducting
        tau = 1 - np.exp(-0.000412 * dlm**2.41)  # Eq (D.2.1)
        
        mu1 = min((10**(-dtm / (16 - 6.6*tau)) + 10**(-(2.48 + 1.77*tau)))**0.2, 1)  # Eq (D.2.2)
        
        if abs(phimn) <= 70:
            mu4 = 10**((-0.935 + 0.0176*abs(phimn)) * np.log10(mu1))  # Eq (D.2.3)
        else:
            mu4 = 10**(0.3 * np.log10(mu1))
        
        if abs(phimn) <= 70:
            b0 = mu1 * mu4 * 10**(-0.015*abs(phimn) + 1.67)  # Eq (D.2.4)
        else:
            b0 = 4.17 * mu1 * mu4
        
        # D.3 Site-shielding losses
        gtr = 0.1 * dlt  # Eq (D.3.1)
        grr = 0.1 * dlr
        
        thetast = thetat - gtr  # Eq (D.3.2)
        thetasr = thetar - grr
        
        # Eq (D.3.3), (D.3.4)
        if thetast > 0:
            Ast = 20*np.log10(1 + 0.361*thetast*np.sqrt(f*dlt)) + 0.264*thetast*f**(1/3)
        else:
            Ast = 0
            
        if thetasr > 0:
            Asr = 20*np.log10(1 + 0.361*thetasr*np.sqrt(f*dlr)) + 0.264*thetasr*f**(1/3)
        else:
            Asr = 0
        
        # D.4 Over-sea surface duct coupling corrections
        dct, dcr = self.distance_to_sea(d, z)
        
        # Eq (D.4.2), (D.4.3)
        if omega >= 0.75 and dct <= dlt and dct <= 5:
            Act = -3*np.exp(-0.25*dct**2) * (1 + np.tanh(0.07*(50 - hts)))
        else:
            Act = 0
            
        if omega >= 0.75 and dcr <= dlr and dcr <= 5:
            Acr = -3*np.exp(-0.25*dcr**2) * (1 + np.tanh(0.07*(50 - hrs)))
        else:
            Acr = 0
        
        # D.5 Total coupling loss
        if f < 0.5:
            Alf = (45.375 - 137.0*f + 92.5*f**2) * omega  # Eq (D.5.2)
        else:
            Alf = 0
        
        Aac = 102.45 + 20*np.log10(f*(dlt + dlr)) + Alf + Ast + Asr + Act + Acr  # Eq (D.5.1)
        
        # D.6 Angular-distance dependent loss
        gammad = 5e-5 * ae * f**(1/3)  # Eq (D.6.1)
        
        theta_at = min(thetat, gtr)  # Eq (D.6.2)
        theta_ar = min(thetar, grr)
        
        theta_a = 1000*dt/ae + theta_at + theta_ar  # Eq (D.6.3)
        
        Aad = gammad * theta_a if theta_a > 0 else 0  # Eq (D.6.4)
        
        # D.7 Distance and time-dependent loss
        dar = min(dt - dlt - dlr, 40)  # Eq (D.7.1)
        
        mu3 = np.exp(-4.6e-5*(hm - 10)*(43 + 6*dar)) if hm > 10 else 1  # Eq (D.7.2)
        
        alpha = -0.6 - 3.5e-9*dt**3.1 * tau  # Eq (D.7.3)
        alpha = max(alpha, -3.4)
        
        mu2 = min((500*dt**2 / (ae*(np.sqrt(htea) + np.sqrt(hrea))**2))**alpha, 1)  # Eq (D.7.4)
        
        bduct = b0 * mu2 * mu3  # Eq (D.7.5)
        
        Gamma = 1.076*np.exp(-1e-6*dt**1.13 * (9.51 - 4.8*np.log10(bduct) + 
                0.198*(np.log10(bduct))**2)) / ((2.0058 - np.log10(bduct))**1.012)  # Eq (D.7.6)
        
        Aat = -12 + (1.2 + 0.0037*dt)*np.log10(p/bduct) + 12*(p/bduct)**Gamma + 50.0/q  # Eq (D.7.7)
        
        # D.8 Basic transmission loss
        Lba = Aac + Aad + Aat  # Eq (D.8.1)
        
        return Lba, Aat, Aad, Aac, dct, dcr, dtm, dlm
    
    def distance_to_sea(self, d: np.ndarray, zone: np.ndarray) -> Tuple[float, float]:
        """
        Distance from each terminal to the sea (Attachment D)
        
        Returns
        -------
        dct : float
            Coast distance from transmitter (km)
        dcr : float
            Coast distance from receiver (km)
        """
        kk = np.where(zone == 1)[0]
        
        if len(kk) == 0:  # Complete path is inland
            dct = d[-1]
            dcr = d[-1]
        else:
            nt = kk[0]
            nr = kk[-1]
            
            if nt == 0:  # Tx is over sea
                dct = 0
            else:
                dct = (d[nt] + d[nt-1]) / 2.0 - d[0]
            
            if nr == len(zone) - 1:  # Rx is over sea
                dcr = 0
            else:
                dcr = d[-1] - (d[nr] + d[nr+1]) / 2.0
        
        return dct, dcr
    
    def longest_cont_dist(self, d: np.ndarray, zone: np.ndarray, zone_r: int) -> float:
        """
        Longest continuous path belonging to zone_r
        
        Parameters
        ----------
        d : np.ndarray
            Vector of distances (km)
        zone : np.ndarray
            Vector of zones
        zone_r : int
            Reference zone (1=Sea, 3=Coastal, 4=Inland, 34=Coastal+Inland)
            
        Returns
        -------
        float
            Longest continuous section (km)
        """
        dm = 0.0
        
        if zone_r == 34:
            start, stop = self.find_intervals((zone == 3) | (zone == 4))
        else:
            start, stop = self.find_intervals(zone == zone_r)
        
        n = len(start)
        nmax = len(d)
        
        for i in range(n):
            delta = 0.0
            if stop[i] < nmax - 1:
                delta += (d[stop[i]+1] - d[stop[i]]) / 2.0
            if start[i] > 0:
                delta += (d[start[i]] - d[start[i]-1]) / 2.0
            
            dm = max(d[stop[i]] - d[start[i]] + delta, dm)
        
        return dm
    
    # ============================================================================
    # SUB-MODEL 3: TROPOSCATTER PROPAGATION
    # ============================================================================
    
    def tropospheric_path(self, dt: float, hts: float, hrs: float, theta_e: float,
                         theta_tpos: float, theta_rpos: float, ae: float,
                         phi_re: float, phi_te: float, phi_rn: float, 
                         phi_tn: float, Re: float) -> Tuple:
        """
        Tropospheric path segments (Section 3.9)
        
        Returns
        -------
        Tuple containing d_tcv, d_rcv, phi_cve, phi_cvn, h_cv, 
                        phi_tcve, phi_tcvn, phi_rcve, phi_rcvn
        """
        # Eq (3.9.1a)
        d_tcv = (dt*np.tan(0.001*theta_rpos + 0.5*theta_e) - 0.001*(hts - hrs)) / \
                (np.tan(0.001*theta_tpos + 0.5*theta_e) + np.tan(0.001*theta_rpos + 0.5*theta_e))
        
        d_tcv = np.clip(d_tcv, 0, dt)
        d_rcv = dt - d_tcv  # Eq (3.9.1b)
        
        # Common volume position
        phi_cve, phi_cvn, _, _ = self.great_circle_path(
            phi_re, phi_te, phi_rn, phi_tn, Re, d_tcv
        )
        
        # Height of common volume (Eq 3.9.2)
        h_cv = hts + 1000*d_tcv*np.tan(0.001*theta_tpos) + 1000*d_tcv**2/(2*ae)
        
        # Midpoint Tx to CV
        d_pnt = 0.5 * d_tcv
        phi_tcve, phi_tcvn, _, _ = self.great_circle_path(
            phi_re, phi_te, phi_rn, phi_tn, Re, d_pnt
        )
        
        # Midpoint CV to Rx
        d_pnt = dt - 0.5*d_rcv
        phi_rcve, phi_rcvn, _, _ = self.great_circle_path(
            phi_re, phi_te, phi_rn, phi_tn, Re, d_pnt
        )
        
        return d_tcv, d_rcv, phi_cve, phi_cvn, h_cv, phi_tcve, phi_tcvn, phi_rcve, phi_rcvn
    
    def tl_troposcatter(self, f: float, dt: float, thetat: float, thetar: float,
                       thetae: float, phicvn: float, phicve: float,
                       phitn: float, phite: float, phirn: float, phire: float,
                       Gt: float, Gr: float, ae: float, p: float) -> Tuple:
        """
        Troposcatter basic transmission loss (Attachment E)
        
        Returns
        -------
        Tuple containing Lbs, theta, climzone
        """
        # E.2 Climatic classification from digital maps or use typical zone 4
        # Use mid-point coordinates
        phi_mn = 0.5 * (phitn + phirn)
        phi_me = 0.5 * (phite + phire)
        climzone = int(self.interp2('TropoClim', phi_me, phi_mn, 1.5, 1.5, default_value=4))
        
        # Table E.1 parameters
        M_arr = [116, 129.6, 119.73, 109.3, 128.5, 119.73, 123.2]
        gamma_arr = [0.27, 0.33, 0.27, 0.32, 0.27, 0.27, 0.27]
        eq_arr = [7, 8, 6, 9, 10, 6, 6]
        
        M = M_arr[climzone]
        gamma = gamma_arr[climzone]
        eq = eq_arr[climzone]
        
        # E.3 Troposcatter basic transmission loss
        theta = 1000*thetae + thetat + thetar  # mrad, Eq (E.1)
        
        H = 0.25e-3 * theta * dt  # Eq (E.3)
        htrop = 0.125e-6 * theta**2 * ae  # Eq (E.4)
        
        LN = 20*np.log10(5 + gamma*H) + 4.34*gamma*htrop  # Eq (E.2)
        
        ds = 0.001 * theta * ae  # Eq (E.5)
        
        # Calculate Y90 based on equation number
        if eq == 6:
            Y90 = -2.2 - (8.1 - 0.23*min(f, 4))*np.exp(-0.137*htrop)
        elif eq == 7:
            Y90 = -9.5 - 3*np.exp(-0.137*htrop)
        elif eq == 8:
            if ds < 100:
                Y90 = -8.2
            elif ds >= 1000:
                Y90 = -3.4
            else:
                Y90 = 1.006e-8*ds**3 - 2.569e-5*ds**2 + 0.02242*ds - 10.2
        elif eq == 9:
            if ds < 100:
                Y90 = -10.845
            elif ds >= 465:
                Y90 = -8.4
            else:
                Y90 = -4.5e-7*ds**3 + 4.45e-4*ds**2 - 0.122*ds - 2.645
        elif eq == 10:
            if ds < 100:
                Y90 = -11.5
            elif ds >= 550:
                Y90 = -4
            else:
                Y90 = -8.519e-8*ds**3 + 7.444e-5*ds**2 - 4.18e-4*ds - 12.1
        else:
            Y90 = 0
        
        # Conversion factor (E.11)
        if p >= 50:
            C = 1.26 * (-np.log10((100.0 - p)/50.0))**0.63
        else:
            C = -1.26 * (-np.log10(p/50.0))**0.63
        
        Yp = C * Y90  # Eq (E.12)
        
        theta = max(theta, 1e-6)
        
        Ldist = max(10*np.log10(dt) + 30*np.log10(theta) + LN,
                   20*np.log10(dt) + 0.573*theta + 20)  # Eq (E.13), (E.14)
        
        Lfreq = 25*np.log10(f) - 2.5*(np.log10(0.5*f))**2  # Eq (E.14)
        
        Lcoup = 0.07*np.exp(0.055*(Gt + Gr))  # Eq (E.15)
        
        Lbs = M + Lfreq + Ldist + Lcoup - Yp  # Eq (E.16)
        
        return Lbs, theta, climzone
    
    def gaseous_abs_tropo(self, phi_te: float, phi_tn: float, phi_re: float,
                         phi_rn: float, h1: float, hn: float,
                         thetatpos: float, thetarpos: float,
                         dtcv: float, drcv: float, f: float) -> Tuple:
        """
        Gaseous absorption for troposcatter path (Attachment F.3)
        
        Returns
        -------
        Tuple containing Aos, Aws, Awrs, Aotcv, Awtcv, Awrtcv, 
                        Aorcv, Awrcv, Awrrcv, Wvsurtx, Wvsurrx
        """
        # Simplified - use typical water vapor density
        rho_sur_tx = 7.5  # g/m^3 (typical)
        rho_sur_rx = 7.5
        
        # Tx to CV path
        Aotcv, Awtcv, Awrtcv = self.gaseous_abs_tropo_t2cv(
            rho_sur_tx, h1, thetatpos, dtcv, f
        )
        
        # Rx to CV path
        Aorcv, Awrcv, Awrrcv = self.gaseous_abs_tropo_t2cv(
            rho_sur_rx, hn, thetarpos, drcv, f
        )
        
        # Total for complete path
        Aos = Aotcv + Aorcv
        Aws = Awtcv + Awrcv
        Awrs = Awrtcv + Awrrcv
        
        return Aos, Aws, Awrs, Aotcv, Awtcv, Awrtcv, Aorcv, Awrcv, Awrrcv, rho_sur_tx, rho_sur_rx
    
    def gaseous_abs_tropo_t2cv(self, rho_sur: float, h_sur: float,
                               theta_el: float, dcv: float, f: float) -> Tuple:
        """
        Gaseous absorption for terminal to common volume (Attachment F.4)
        
        Returns
        -------
        Tuple containing Ao, Aw, Awr
        """
        gamma_o, gamma_w = self.specific_sea_level_attenuation(f, rho_sur, h_sur)
        
        rho_surr = self.water_vapour_density_rain(rho_sur, h_sur)
        _, gamma_wr = self.specific_sea_level_attenuation(f, rho_surr, h_sur)
        
        d0 = 5.0 / (0.65*np.sin(0.001*theta_el) + 0.35*np.sqrt((np.sin(0.001*theta_el))**2 + 0.00304))
        dw = 2.0 / (0.65*np.sin(0.001*theta_el) + 0.35*np.sqrt((np.sin(0.001*theta_el))**2 + 0.00122))
        
        deo = d0 * (1 - np.exp(-dcv/d0)) * np.exp(-h_sur/5000.0)
        dew = dw * (1 - np.exp(-dcv/dw)) * np.exp(-h_sur/2000.0)
        
        Ao = gamma_o * deo
        Aw = gamma_w * dew
        Awr = gamma_wr * dew
        
        return Ao, Aw, Awr
    
    def gaseous_abs_surface(self, phi_me: float, phi_mn: float, h_mid: float,
                           hts: float, hrs: float, dt: float, f: float) -> Tuple:
        """
        Gaseous absorption on surface paths (Attachment F.2)
        
        Returns
        -------
        Tuple containing Aosur, Awsur, Awrsur, gamma_o, gamma_w, gamma_wr, rho_sur
        """
        # Get water vapor density from digital maps or use typical value
        rho_sur = self.interp2('surfwv_50_fixed', phi_me, phi_mn, 1.5, 1.5, default_value=7.5)
        h_sur = h_mid
        
        gamma_o, gamma_w = self.specific_sea_level_attenuation(f, rho_sur, h_sur)
        
        rho_surr = self.water_vapour_density_rain(rho_sur, h_sur)
        _, gamma_wr = self.specific_sea_level_attenuation(f, rho_surr, h_sur)
        
        h_rho = 0.5 * (hts + hrs)
        
        Aosur = gamma_o * dt * np.exp(-h_rho/5000.0)
        Awsur = gamma_w * dt * np.exp(-h_rho/2000.0)
        Awrsur = gamma_wr * dt * np.exp(-h_rho/2000.0)
        
        return Aosur, Awsur, Awrsur, gamma_o, gamma_w, gamma_wr, rho_sur
    
    def specific_sea_level_attenuation(self, f: float, rho_sur: float, 
                                      h_sur: float) -> Tuple[float, float]:
        """
        Specific sea-level attenuations (Attachment F.6)
        
        Returns
        -------
        gamma_o : float
            Oxygen attenuation (dB/km)
        gamma_w : float
            Water vapor attenuation (dB/km)
        """
        rho_sea = rho_sur * np.exp(h_sur/2000)
        eta = 0.955 + 0.006*rho_sea
        
        gamma_o = (7.2/(f**2 + 0.34) + 0.62/((54 - f)**1.16 + 0.83)) * f**2 * 1e-3
        
        gamma_w = (0.046 + 0.0019*rho_sea + 
                  3.98*eta/((f - 22.235)**2 + 9.42*eta**2) * 
                  (1 + ((f - 22.0)/(f + 22.0))**2)) * f**2 * rho_sea * 1e-4
        
        return gamma_o, gamma_w
    
    def water_vapour_density_rain(self, rho_sur: float, h_sur: float) -> float:
        """
        Atmospheric water-vapour density in rain (Attachment F.5)
        """
        if h_sur <= 2600:
            return rho_sur + 0.4 + 0.0003*h_sur
        else:
            return rho_sur + 5*np.exp(-h_sur/1800.0)
    
    # ============================================================================
    # SUB-MODEL 4: SPORADIC-E PROPAGATION
    # ============================================================================
    
    def tl_sporadic_e(self, f: float, dt: float, thetat: float, thetar: float,
                     phimn: float, phime: float, phitn: float, phite: float,
                     phirn: float, phire: float, dlt: float, dlr: float,
                     ae: float, Re: float, p: float) -> Tuple:
        """
        Sporadic-E transmission loss (Attachment G)
        
        Returns
        -------
        Tuple containing Lbe, Lbes1, Lbes2, Lp1t, Lp2t, Lp1r, Lp2r,
                        Gamma1, Gamma2, FoEs1hop, FoEs2hop, 
                        Phi1qe, Phi1qn, Phi3qe, Phi3qn
        """
        # G.2 Derivation of FoEs from digital maps or use typical values
        # Interpolate FoEs based on time percentage
        if p <= 1:
            FoEs1hop = self.interp2('FoEs01', phime, phimn, 1.5, 1.5, default_value=6.5)
            FoEs2hop = FoEs1hop * 0.9
        elif p <= 10:
            FoEs01 = self.interp2('FoEs01', phime, phimn, 1.5, 1.5, default_value=6.5)
            FoEs10 = self.interp2('FoEs10', phime, phimn, 1.5, 1.5, default_value=5.0)
            # Linear interpolation between 1% and 10%
            w = (p - 1) / 9.0
            FoEs1hop = FoEs01 * (1 - w) + FoEs10 * w
            FoEs2hop = FoEs1hop * 0.9
        elif p <= 50:
            FoEs10 = self.interp2('FoEs10', phime, phimn, 1.5, 1.5, default_value=5.0)
            FoEs50 = self.interp2('FoEs50', phime, phimn, 1.5, 1.5, default_value=3.5)
            # Linear interpolation between 10% and 50%
            w = (p - 10) / 40.0
            FoEs1hop = FoEs10 * (1 - w) + FoEs50 * w
            FoEs2hop = FoEs1hop * 0.9
        else:
            FoEs1hop = self.interp2('FoEs50', phime, phimn, 1.5, 1.5, default_value=3.5)
            FoEs2hop = FoEs1hop * 0.9
        
        # G.2 1-hop propagation
        Gamma1 = (40.0/(1 + dt/130.0 + (dt/250.0)**2) + 0.2*(dt/2600.0)**2) * \
                 (1000.0*f/FoEs1hop)**2 + np.exp((dt - 1660)/280.0)
        
        hes = 120.0  # km
        l1 = 2 * (ae**2 + (ae + hes)**2 - 2*ae*(ae + hes)*np.cos(dt/(2*ae)))**0.5
        
        Lbfs1 = self.tl_free_space(f, l1)
        
        alpha1 = dt / (2*ae)
        epsr1 = 0.5*np.pi - np.arctan(ae*np.sin(alpha1)/(hes + ae*(1 - np.cos(alpha1)))) - alpha1
        
        delta1t = 0.001*thetat - epsr1
        delta1r = 0.001*thetar - epsr1
        
        nu1t = (1 if delta1t >= 0 else -1) * 3.651 * \
               np.sqrt(1000*f*dlt*(1 - np.cos(delta1t))/np.cos(0.001*thetat))
        nu1r = (1 if delta1r >= 0 else -1) * 3.651 * \
               np.sqrt(1000*f*dlr*(1 - np.cos(delta1r))/np.cos(0.001*thetar))
        
        Lp1t = self.dl_knife_edge(nu1t)
        Lp1r = self.dl_knife_edge(nu1r)
        
        Lbes1 = Lbfs1 + Gamma1 + Lp1t + Lp1r
        
        # G.3 2-hop propagation
        Phi1qe, Phi1qn, _, _ = self.great_circle_path(phire, phite, phirn, phitn, Re, 0.25*dt)
        Phi3qe, Phi3qn, _, _ = self.great_circle_path(phire, phite, phirn, phitn, Re, 0.75*dt)
        
        Gamma2 = (40.0/(1 + (dt/260.0) + (dt/500.0)**2) + 0.2*(dt/5200.0)**2) * \
                 (1000.0*f/FoEs2hop)**2 + np.exp((dt - 3220.0)/560.0)
        
        l2 = 4 * (ae**2 + (ae + hes)**2 - 2*ae*(ae + hes)*np.cos(dt/(4.0*ae)))**0.5
        
        Lbfs2 = self.tl_free_space(f, l2)
        
        alpha2 = dt / (4.0*ae)
        epsr2 = 0.5*np.pi - np.arctan(ae*np.sin(alpha2)/(hes + ae*(1 - np.cos(alpha2)))) - alpha2
        
        delta2t = 0.001*thetat - epsr2
        delta2r = 0.001*thetar - epsr2
        
        nu2t = (1 if delta2t >= 0 else -1) * 3.651 * \
               np.sqrt(1000*f*dlt*(1 - np.cos(delta2t))/np.cos(0.001*thetat))
        nu2r = (1 if delta2r >= 0 else -1) * 3.651 * \
               np.sqrt(1000*f*dlr*(1 - np.cos(delta2r))/np.cos(0.001*thetar))
        
        Lp2t = self.dl_knife_edge(nu2t)
        Lp2r = self.dl_knife_edge(nu2r)
        
        Lbes2 = Lbfs2 + Gamma2 + Lp2t + Lp2r
        
        # G.4 Basic transmission loss
        if Lbes1 < Lbes2 - 20:
            Lbe = Lbes1
        elif Lbes2 < Lbes1 - 20:
            Lbe = Lbes2
        else:
            Lbe = -10*np.log10(10**(-0.1*Lbes1) + 10**(-0.1*Lbes2))
        
        return Lbe, Lbes1, Lbes2, Lp1t, Lp2t, Lp1r, Lp2r, Gamma1, Gamma2, \
               FoEs1hop, FoEs2hop, Phi1qe, Phi1qn, Phi3qe, Phi3qn
    
    # ============================================================================
    # MULTIPATH AND PRECIPITATION FUNCTIONS
    # ============================================================================
    
    def multi_path_activity(self, f: float, dt: float, hts: float, hrs: float,
                           dlt: float, dlr: float, hlt: float, hlr: float,
                           hlo: float, thetat: float, thetar: float, epsp: float,
                           Nd65m1: float, phimn: float, FlagLos50: int) -> float:
        """
        Multipath fading calculation (Attachment B.2)
        
        Returns
        -------
        Q0ca : float
            Notional zero-fade annual percentage time
        """
        K = 10**(-(4.6 + 0.0027*Nd65m1))
        
        if FlagLos50 == 1:
            Q0ca = self.zero_fade_annual_time(dt, epsp, hlo, f, K, phimn)
        else:
            Q0cat = self.zero_fade_annual_time(dlt, abs(thetat), min(hts, hlt), f, K, phimn)
            Q0car = self.zero_fade_annual_time(dlr, abs(thetar), min(hrs, hlr), f, K, phimn)
            Q0ca = max(Q0cat, Q0car)
        
        return Q0ca
    
    def zero_fade_annual_time(self, dca: float, epsca: float, hca: float,
                             f: float, K: float, phimn: float) -> float:
        """
        Notional zero-fade annual percentage time (Attachment B.3)
        """
        qw = K * dca**3.1 * (1 + epsca)**(-1.29) * f**0.8 * 10**(-0.00089*hca)
        
        if abs(phimn) <= 45:
            Cg = min(10.5 - 5.6*np.log10(1.1 + (abs(self.cosd(2*phimn)))**0.7) - 
                    2.7*np.log10(dca) + 1.7*np.log10(1 + epsca), 10.8)
        else:
            Cg = min(10.5 - 5.6*np.log10(1.1 - (abs(self.cosd(2*phimn)))**0.7) - 
                    2.7*np.log10(dca) + 1.7*np.log10(1 + epsca), 10.8)
        
        Q0ca = qw * 10**(-0.1*Cg)
        return Q0ca
    
