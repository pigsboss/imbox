function mCrab = cts2mCrab(cts,AreaEff,keVEff)
% mCrab = cts/AreaEff/keV2erg(keVEff);
mCrab = cts*1000/AreaEff/15*keVEff;
return