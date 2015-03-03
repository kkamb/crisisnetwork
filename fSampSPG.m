%%  fSampSPG samples (Sig, Phi, gamma)
%%  conditional on time-varying parameter
%%
%%  [NOTE]
%%    y_t = sqrt(gam)*exp(h_t/2)eps_t,  eps_t ~ N(0,1)
%%    h_{t+1} = Phi*h_t + eta_t,
%%    eta_t ~ N(0, Sig),  eta_0 ~ N(0, Sig/(I-Phi^2))
%%
function [mSign, mPhin, vgamn] = fSampSPG(my, mh, mPhi, vgam, dnu, dV0, da0, db0, dg0, dG0)
    ns = size(my,1);
    nk = size(my,2);
    mSign = eye(nk); mPhin = mPhi; vgamn  = vgam;
    
    mh0=mh(1:ns-1,:); mh1=mh(2:end,:);
    vG = dG0+sum(my.^2./exp(mh));
    vsum=sum(mh(2:ns-1,:).^2);
    vphii=sum(mh1.*mh0)./vsum;
    
    for i=1:nk
        
        %%--- Sampling Sig ---*/
        phi=mPhi(i,i);
        e=mh1(:,i)- phi*mh0(:,i);
        q=dV0+e'*e+(1-phi^2)*mh(1,i)^2;
        mSign(i,i) = 1 / gamrnd(dnu/2, 2/q);
    
        %%---- sampling gamma------%%/
	    vgamn(i) = 1 / gamrnd((dg0+ns)/2, 2./vG(i));
    
        %%--- sampling Phi from truncated normal posterior---*/
        vsigi=mSign(i,i)/vsum(i);
        A=normcdf([0 1]',vphii(i),sqrt(vsigi));
        dphin=norminv(A(1)+rand(1)*(A(2)-A(1)),vphii(i),sqrt(vsigi));
        dfrac = betapdf((dphin+1)/2, da0, db0)/ betapdf((phi+1)/2, da0, db0) ...
				* sqrt(1 - dphin^2) / sqrt(1 - phi^2);
        if (rand(1, 1) < dfrac)
            mPhin(i,i) = dphin;
        end
        
    end