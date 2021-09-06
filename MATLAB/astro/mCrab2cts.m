function cts = mCrab2cts(mCrab,AreaEff,keVEff)
%MCRAB2CTS Convert mCrab to counts/s according to the detector parameters
%and the energy spectrum.
%
%INPUTS
% mCrab   is flux of the source in milli-Crab.
% AreaEff is effective area of the detector in cm^2.
% keVEff  is effective photon energy in keV.
cts = 1e-3*mCrab*AreaEff*15/keVEff;
return