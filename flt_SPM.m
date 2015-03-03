%%  fSampSPM samples (Sig, Phi, mu)
%%  conditional on time-varying parameter and truncated variables vd and dk0
%%
%%   B_{t+1} = mu + Phi*(B_t - mu) + eta_t,
%%   eta_t ~ N(0, Sig),  eta_0 ~ N(0, Sig/(I-Phi^2))
%%

function [mSign, mPhin, vmun] = flt_SPM(mp, mSig, mPhi, vmu, dnu, dV0, da0, db0, dm0, ds0, vd, dk0)
    ns = size(mp,1);
    np = size(mp,2);
    mSign = diag(mSig);
    mPhin = diag(mPhi);
    vmun=vmu;
    a=dnu-ns-1;
    G=ds0; g=dm0;
    for i=1:np
        z=mp(:,i);
        phi=mPhin(i); v=mSign(i); mu=vmun(i);
        
        %sample phi
        fldo=0; fl=0;
        xlag=z(1:(ns-1))-mu;
        x=z(2:end)-mu;
        B=xlag'*xlag; 
        Gphi=v/B; gphi=xlag'*x/B;
        Gphi=sqrt(Gphi); A=normcdf([0 1]',gphi,Gphi); 
        while (fl==0 && fldo<100)
            phistar=norminv(A(1)+rand(1)*(A(2)-A(1)),gphi,Gphi);
            dup=fK0(mu,v,phistar,dk0);
            fldo=fldo+1;
            if vd(i)<dup
                fl=1;
            end
            aphi=betapdf((phistar+1)/2,da0,db0)/betapdf((phi+1)/2,da0,db0)*sqrt(1-phistar^2)/sqrt(1-phi^2);
            dfrac=aphi*fK0(mu,v,phi,dk0)/dup;
            if rand(1)<dfrac
                    phi=phistar;
            end
        end
        %sample mu
        z0=z(1);zlag=z(1:(ns-1));z=z(2:end);
        fldo=0; fl=0;
        s=1/(1/G+(1-phi^2)/v+ns*(1-phi)^2/v);
        mustar=s*(g/G+(1-phi^2)*z0/v+sum(z-phi*zlag)*(1-phi)/v);
        while fldo<100 && fl==0						
            mustar = mustar+sqrt(s)*randn(1,1);
            dup = fK0(mustar,v,phi,dk0);
            if vd(i)<dup
                fl=1;
                dfrac = fK0(mu,v,phi,dk0) / dup;
                if rand(1) < dfrac
                    mu=mustar;
                end
            end
            fldo=fldo+1;
        end
        
        %sample sigma
        fldo=0; fl=0;
        e=(z-mu)-phi*(zlag-mu);
        q=(1-phi^2)*(z0-mu)^2+e'*e; 
        while fldo<100 && fl==0		
            vstar=1/gamrnd((a+ns+1)/2,(a*dV0+q)/2);
            dup = fK0(mu,vstar,phi,dk0);
            if vd(i)<dup
                fl=1;
                dfrac = fK0(mu,v,phi,dk0)/ dup;
                if rand(1) < dfrac
                    v=vstar;
                end
            end
            fldo=fldo+1;
        end
        mPhin(i)=phi;
        mSign(i)=v;
        vmun(i)=mu;
    end
    mSign=diag(mSign);
    mPhin=diag(mPhin);