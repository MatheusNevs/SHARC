# -*- coding: utf-8 -*-
"""
ITU-R P.2001-4 Propagation Model Implementation for SHARC
Based on the reference implementation from ITU-R P.2001-4
"""

import numpy as np
from sharc.propagation.propagation import Propagation
from sharc.propagation.p2001_aux import P2001Aux
from sharc.parameters.parameters import Parameters
from sharc.station_manager import StationManager


class PropagationP2001(Propagation):
    """
    Implementation of ITU-R P.2001-4 propagation model.
    
    This model covers frequencies from 30 MHz to 50 GHz and is most accurate 
    for distances from 3 km to at least 1000 km. It includes multiple propagation
    mechanisms:
    - Sub-model 1: Line-of-sight and diffraction
    - Sub-model 2: Anomalous propagation (ducting/layer-reflection)
    - Sub-model 3: Troposcatter
    - Sub-model 4: Sporadic-E
    """

    def __init__(self, random_number_gen: np.random.RandomState, param_p2001):
        super().__init__(random_number_gen)
        self.param_p2001 = param_p2001
        self.aux = P2001Aux()
        # P.2001 can handle earth-space paths depending on configuration
        self.is_earth_space_model = False

    def get_loss(
        self,
        params: Parameters,
        frequency: float,
        station_a: StationManager,
        station_b: StationManager,
        station_a_gains=None,
        station_b_gains=None,
    ) -> np.array:
        """
        Calculate path loss using ITU-R P.2001-4 model.

        Parameters
        ----------
        params : Parameters
            Simulation parameters
        frequency : float
            Frequency in MHz
        station_a : StationManager
            Transmitter station
        station_b : StationManager
            Receiver station
        station_a_gains : np.array, optional
            Transmitter antenna gains
        station_b_gains : np.array, optional
            Receiver antenna gains

        Returns
        -------
        np.array
            Path loss matrix (dB) with shape (num_stations_a, num_stations_b)
        """
        # Convert frequency from MHz to GHz
        freq_ghz = frequency / 1000.0

        num_stations_a = station_a.num_stations
        num_stations_b = station_b.num_stations

        # Initialize path loss matrix
        path_loss = np.zeros((num_stations_a, num_stations_b))

        # Get percentage time from parameters (default to 50% if not specified)
        Tpc = getattr(self.param_p2001, 'time_percentage', 50.0)

        # Get polarization (0 = horizontal, 1 = vertical)
        FlagVP = getattr(self.param_p2001, 'polarization', 0)

        # Loop through all station pairs
        for idx_a in range(num_stations_a):
            for idx_b in range(num_stations_b):
                # Get station positions
                lat_a = station_a.get_station(idx_a).x
                lon_a = station_a.get_station(idx_a).y
                h_a = station_a.get_station(idx_a).height

                lat_b = station_b.get_station(idx_b).x
                lon_b = station_b.get_station(idx_b).y
                h_b = station_b.get_station(idx_b).height

                # Get antenna heights above ground
                htg = getattr(station_a.get_station(idx_a), 'height_above_ground', h_a)
                hrg = getattr(station_b.get_station(idx_b), 'height_above_ground', h_b)

                # Get antenna gains if provided
                Gtx = station_a_gains[idx_a, idx_b] if station_a_gains is not None else 0.0
                Grx = station_b_gains[idx_a, idx_b] if station_b_gains is not None else 0.0

                # Get or generate path profile
                d, h, z = self._get_path_profile(
                    lon_a, lat_a, h_a,
                    lon_b, lat_b, h_b,
                    params
                )

                # Calculate basic transmission loss using P.2001-4
                Lb = self.aux.bt_loss(
                    d=d,
                    h=h,
                    z=z,
                    GHz=freq_ghz,
                    Tpc=Tpc,
                    Phire=lon_b,
                    Phirn=lat_b,
                    Phite=lon_a,
                    Phitn=lat_a,
                    Hrg=hrg,
                    Htg=htg,
                    Grx=Grx,
                    Gtx=Gtx,
                    FlagVP=FlagVP
                )

                path_loss[idx_a, idx_b] = Lb

        return path_loss

    def _get_path_profile(self, lon_a, lat_a, h_a, lon_b, lat_b, h_b, params):
        """
        Get or generate path profile between two stations.

        Parameters
        ----------
        lon_a, lat_a, h_a : float
            Transmitter longitude, latitude, and height (m amsl)
        lon_b, lat_b, h_b : float
            Receiver longitude, latitude, and height (m amsl)
        params : Parameters
            Simulation parameters

        Returns
        -------
        d : np.array
            Distances from transmitter (km)
        h : np.array
            Heights above mean sea level (m)
        z : np.array
            Zone codes (1=Sea, 3=Coastal, 4=Inland)
        """
        # Calculate great-circle distance
        Re = 6371.0  # Earth radius in km
        
        # Simple great-circle distance calculation
        lat_a_rad = np.deg2rad(lat_a)
        lat_b_rad = np.deg2rad(lat_b)
        dlon_rad = np.deg2rad(lon_b - lon_a)
        
        a = np.sin((lat_b_rad - lat_a_rad) / 2) ** 2 + \
            np.cos(lat_a_rad) * np.cos(lat_b_rad) * np.sin(dlon_rad / 2) ** 2
        c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
        distance_km = Re * c

        # Generate simplified path profile
        # In a full implementation, this should use actual terrain data
        num_points = max(int(distance_km) + 1, 11)
        
        d = np.linspace(0, distance_km, num_points)
        
        # Simple linear interpolation for heights
        h = np.linspace(h_a, h_b, num_points)
        
        # Default zone to inland (4)
        # In a full implementation, this should use actual land/sea data
        z = np.ones(num_points, dtype=int) * 4

        return d, h, z
