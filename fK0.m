
function mK0 = fK0(dmu, dsig2, dphi, dk0)

mK0 = abs(dmu) + dk0*sqrt(dsig2/(1-dphi^2));