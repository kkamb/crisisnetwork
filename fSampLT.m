%%  fSampLT samples latent threshold conditional on time-varying parameter
%%  based on Ox Code by Nakajima
%%    y_t = X_t*b_t + e_t,    e_t ~ N(0, G2)

function [vdo, ca] = fSampLT(my, amX, mp, amG2inv, vd, vdc)
%vd - vda (na*1 alpha threshold) or vdb (nb*1 beta threshold)
%vdc - vdac alpha prior for U or vdbc beta prior for U
%amG2inv - nk*nk*(ns - # of lags) or am0minv 

ns = size(mp,1);
nk = size(mp,2);
vdo = vd;
ca = 0;

for i=1:nk
    vdn = vdo;
    vdn(i) = rand() * vdc(i);
    mdn = mp.*bsxfun(@gt,abs(mp),vdn');
    mdo= mp.*bsxfun(@gt,abs(mp),vdo');
    dln = 0;
    dlo=0;
    for t=1:ns
    	ve = my(t,:) - mdn(t,:) * amX(:,:,t)';
        dln = dln+ve * amG2inv(:,:,t) * ve';
    	ve = my(t,:) - mdo(t,:) * amX(:,:,t)';
        dlo = dlo + ve * amG2inv(:,:,t) * ve';
    end
    if rand(1, 1) < exp(-0.5 * (dln - dlo))
    	vdo = vdn;
        ca=ca+1;
    end
end