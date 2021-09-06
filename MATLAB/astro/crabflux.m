function F = crabflux(E)
%CRABFLUX Flux of Crab nebula at given energy (>=20keV).
%
%Description
%E is energy in keV.
%F is photon flux in photons per second per cm^2 per keV.
%Reference:
%1. Jourdain, E. and Roques, J. P.,  The high-energy emission of the Crab
%nebula from 20keV to 6 MeV with INTEGRAL/SPI. ApJ, 2009.
%2. Massaro, E., et al., Fine phase resolved spectroscopy of the X-ray
%emission of the Crab pulsar (PSR B0531+21) observed with BeppoSAX.

a = 1.79;
b = 0.134;
A = 3.87; % pts / cm^2 / s / keV
E0 = 20;
F = A * (E).^(-(a + b*log10(E/E0)));
return