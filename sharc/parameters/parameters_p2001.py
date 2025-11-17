# -*- coding: utf-8 -*-
"""Parameters definitions for ITU-R P.2001 propagation model

The ITU-R P.2001 model is a general purpose terrestrial propagation prediction method
for point-to-point services in the frequency range 30 MHz to 50 GHz.
"""
from dataclasses import dataclass
from sharc.parameters.parameters_base import ParametersBase


@dataclass
class ParametersP2001(ParametersBase):
    """Dataclass containing the P.2001 propagation model parameters
    
    Attributes
    ----------
    frequency : float
        Frequency in MHz (30 MHz to 50 GHz)
    atmospheric_pressure : float
        Total air pressure in hPa (default: 1013.25)
    air_temperature : float
        Temperature in Kelvin (default: 288.15)
    time_percentage : float
        Time percentage for which the predicted basic transmission loss 
        is not exceeded (0.01% to 50%, default: 50%)
    transmitter_height : float
        Height of the transmitter antenna center above ground level in meters
    receiver_height : float
        Height of the receiver antenna center above ground level in meters
    polarization : str
        Antenna polarization: "horizontal" or "vertical" (default: "horizontal")
    tx_lat : float
        Latitude of transmitter in degrees
    tx_lon : float
        Longitude of transmitter in degrees
    rx_lat : float
        Latitude of receiver in degrees
    rx_lon : float
        Longitude of receiver in degrees
    zone : str
        Radio-climatic zone ("COASTAL", "INLAND", "SEA" - default: "INLAND")
    delta_N : float
        Average radio-refractivity lapse-rate through the lowest 1 km of the atmosphere
        (N-units/km). If not provided, it will be calculated from maps.
    N0 : float
        Sea-level surface refractivity (N-units). If not provided, 
        it will be calculated from maps.
    path_center_lat : float
        Latitude of path center in degrees (optional)
    path_center_lon : float
        Longitude of path center in degrees (optional)
    percentage_sea : float
        Percentage of the path over sea (0 to 100%, default: 0%)
    clutter_loss : bool
        Determine whether clutter loss is added (default: True)
    """
    
    # Frequency parameters
    frequency: float = 2000.0  # MHz
    
    # Atmospheric parameters
    atmospheric_pressure: float = 1013.25  # hPa
    air_temperature: float = 288.15  # Kelvin (15°C)
    
    # Time percentage
    time_percentage: float = 50.0  # percent
    
    # Antenna heights
    transmitter_height: float = 10.0  # meters
    receiver_height: float = 10.0  # meters
    
    # Polarization
    polarization: str = "horizontal"  # or "vertical"
    
    # Geographic coordinates
    tx_lat: float = -23.55028  # degrees
    tx_lon: float = -46.63396  # degrees
    rx_lat: float = -23.17889  # degrees
    rx_lon: float = -46.87611  # degrees
    
    # Radio-climatic zone
    zone: str = "INLAND"  # COASTAL, INLAND, SEA
    
    # Refractivity parameters (optional - can be calculated from ITU maps)
    delta_N: float = 40.0  # N-units/km
    N0: float = 325.0  # N-units
    
    # Path characteristics
    path_center_lat: float = None  # degrees (optional)
    path_center_lon: float = None  # degrees (optional)
    percentage_sea: float = 0.0  # percent
    
    # Additional options
    clutter_loss: bool = True
    
    def __post_init__(self):
        """Validate parameters after initialization"""
        self._validate_parameters()
    
    def _validate_parameters(self):
        """Validate P.2001 parameters according to ITU-R P.2001 recommendation"""
        
        # Validate frequency range
        if not (30 <= self.frequency <= 50000):
            raise ValueError(
                f"ParametersP2001: Frequency must be between 30 MHz and 50 GHz. "
                f"Got {self.frequency} MHz"
            )
        
        # Validate time percentage
        if not (0.01 <= self.time_percentage <= 50):
            raise ValueError(
                f"ParametersP2001: Time percentage must be between 0.01% and 50%. "
                f"Got {self.time_percentage}%"
            )
        
        # Validate polarization
        if self.polarization.lower() not in ["horizontal", "vertical"]:
            raise ValueError(
                f"ParametersP2001: Polarization must be 'horizontal' or 'vertical'. "
                f"Got '{self.polarization}'"
            )
        
        # Validate zone
        valid_zones = ["COASTAL", "INLAND", "SEA"]
        if self.zone.upper() not in valid_zones:
            raise ValueError(
                f"ParametersP2001: Zone must be one of {valid_zones}. "
                f"Got '{self.zone}'"
            )
        
        # Validate percentage sea
        if not (0 <= self.percentage_sea <= 100):
            raise ValueError(
                f"ParametersP2001: Percentage of path over sea must be between 0% and 100%. "
                f"Got {self.percentage_sea}%"
            )
        
        # Validate heights (must be positive)
        if self.transmitter_height < 0 or self.receiver_height < 0:
            raise ValueError(
                "ParametersP2001: Antenna heights must be positive values"
            )
        
        # Normalize string parameters
        self.polarization = self.polarization.lower()
        self.zone = self.zone.upper()
    
    def load_from_parameters(self, param: ParametersBase):
        """Used to load parameters of P.2001 from IMT or system parameters

        Parameters
        ----------
        param : ParametersBase
            IMT or system parameters
        """
        # Load common parameters if they exist
        if hasattr(param, 'tx_lat'):
            self.tx_lat = param.tx_lat
        if hasattr(param, 'tx_lon'):
            self.tx_lon = param.tx_lon
        if hasattr(param, 'rx_lat'):
            self.rx_lat = param.rx_lat
        if hasattr(param, 'rx_lon'):
            self.rx_lon = param.rx_lon
        if hasattr(param, 'atmospheric_pressure'):
            self.atmospheric_pressure = param.atmospheric_pressure
        if hasattr(param, 'air_temperature'):
            self.air_temperature = param.air_temperature
        
        # Validate after loading
        self._validate_parameters()
