%%  fSampLTBv samples beta with latent threshold in multivariate model
%%  by singe-move sampler
%%  based on Nakajima's Ox Code
%%
%%  [model]
%%      y_t = X_t*b_t + e_t       e_t ~ N(0, G2)
%%      b_t = B_t*s_t             s_t ~ I(B_t>=d)
%%      B_{t+1} = mu + Phi*(B_t - mu) + eta_t,
%%      eta_t ~ N(0, Sig), eta_0 ~ N(0, Sig/(I-Phi^2))
%%
%%      y_t:     nk*1 vector
%%      X_t:     nk*np matrix
%%      b_t:     np*1 vector
%%
%%  [input]
%%      my:     response (ns*nk vector)
%%      amX:    independent variable (nk*np*ns array)
%%  [output]
%%      mba, ca

function [mba, ca] = fSampLTBv(my, amX, mPhi, vmu, mSig, amG2inv, vd, mbo)

%%--- set variables ---%%

ns = size(my, 1);    % # of time periods
nk = size(my, 2);    % # of series
np = size(amX(:,:,1), 2);    % # of state
ca=0;
mba=mbo;
mSigi=inv(mSig);
mSigp=mSigi*(eye(np)+mPhi.^2);
vpm1=(eye(np)-mPhi)*vmu;
vpm2=(eye(np)-2*mPhi+mPhi.^2)*vmu;
mSig0i=(eye(np)-mPhi.^2)*mSigi;

for i = 1 : ns
    if i==1
        mSigh= finvm(amX(:,:,i)'*amG2inv(:,:,i)*amX(:,:,i)+mSig0i+mSigi+mPhi.^2);
        vbh=mSigh*(amX(:,:,i)'*amG2inv(:,:,i)*my(i,:)'+mSig0i*vmu+mSigi*mPhi*(mba(i+1,:)'-vpm1));
    elseif i<ns
        mSigh = finvm(amX(:,:,i)'*amG2inv(:,:,i)*amX(:,:,i)+mSigp);
        vbh=mSigh*(amX(:,:,i)'*amG2inv(:,:,i)*my(i,:)'+mSigi*(mPhi*(vba_1+mba(i+1,:)')+vpm2));
    else
        mSigh = finvm(amX(:,:,i)'*amG2inv(:,:,i)*amX(:,:,i)+mSigi);
        vbh=mSigh*(amX(:,:,i)'*amG2inv(:,:,i)*my(i,:)'+mSigi*(mPhi*vba_1+vpm1));
    end
    vbn = vbh + chol(mSigh) * randn(np, 1);  %candidate
    vbo = mba(i,:)';
    mSighi = finvm(mSigh);
    dhn = -0.5 * log(det(mSigh));
    dho = -0.5 * log(det(mSigh));
    dhn = dhn -0.5 * (vbn - vbh)'* mSighi * (vbn - vbh);
    dho = dho -0.5 * (vbo - vbh)'* mSighi * (vbo - vbh);
    if sum(abs(vbn)<vd)~=0
        mxh=bsxfun(@times,amX(:,:,i),(abs(vbn)>=vd)');
        if i==1
            mSigh = finvm( mxh'*amG2inv(:,:,i)*mxh + mSig0i + mSigi * mPhi.^2 );
            vbh = mSigh*( mxh'*amG2inv(:,:,i)*my(i,:)' + mSig0i*vmu + mSigi*mPhi*(mba(i+1,:)'-vpm1));
        elseif i<ns
            mSigh = finvm(mxh'*amG2inv(:,:,i)*mxh+mSigp);
            vbh = mSigh*( mxh'*amG2inv(:,:,i)*my(i,:)'+ mSigi*(mPhi*(vba_1+mba(i+1,:)')+vpm2));
        else
            mSigh = finvm(mxh'*amG2inv(:,:,i)*mxh+mSigi);
            vbh = mSigh*(mxh'*amG2inv(:,:,i)*my(i,:)'+ mSigi*( mPhi*vba_1 + vpm1 ) );
        end
        dln =-0.5 * log(det(mSigh))-0.5 * (vbn-vbh)'*inv(mSigh) * (vbn-vbh);
    else
        dln = dhn;
    end
    if sum(abs(vbo)<vd)~=0	
        mxh=bsxfun(@times,amX(:,:,i),(abs(vbo)>=vd)');
        if i==1
            mSigh = finvm( mxh'*amG2inv(:,:,i)*mxh + mSig0i + mSigi * mPhi.^2 );
            vbh = mSigh*( mxh'*amG2inv(:,:,i)*my(i,:)' + mSig0i*vmu + mSigi*mPhi*(mba(i+1,:)'-vpm1));
        elseif i<ns
            mSigh = finvm(mxh'*amG2inv(:,:,i)*mxh+mSigp);
            vbh = mSigh*( mxh'*amG2inv(:,:,i)*my(i,:)'+ mSigi*(mPhi*(vba_1+mba(i+1,:)')+vpm2));
        else
            mSigh = finvm(mxh'*amG2inv(:,:,i)*mxh+mSigi);
            vbh = mSigh*(mxh'*amG2inv(:,:,i)*my(i,:)'+ mSigi*( mPhi*vba_1 + vpm1 ) );
        end
	dlo =-0.5 * log(det(mSigh))-0.5 * (vbo-vbh)'*inv(mSigh) * (vbo-vbh);
    else
        dlo = dho;
    end
    dfrac = exp(dln - dhn - dlo + dho);
    if rand(1, 1) < dfrac
        mba(i,:) = vbn';
        vba_1 = vbn;
        ca=ca+1;
    else
        vba_1 = mba(i,:)';
    end
end