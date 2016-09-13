%code largely replicates Nakajima's Ox-code for latent-threshold multivariate time-varying AR model

fcrisis;
nsim=20000;
nburn=nsim*.1;
nstore=1000;

%function lt_mcmc2(nsim)
global m_my  m_nl m_ns m_nk m_fli m_iseed m_flLTa m_flLTb ...
    dvb0 dVb0 dva0 dVa0 dvh0 dVh0 da0 db0 dm0 ds0 m_k;
tic;

%%--- set default options ---%%
if isempty(m_fli) == 1
    m_fli = 0;
end
if isempty(m_iseed) == 1
    m_iseed = 1;
end
if isempty(m_flLTb) == 1
    m_flLTb = 1;
end
if isempty(m_flLTa) == 1
    m_flLTa = 1;
end
rng('default'); rng(m_iseed);

%%--- set variables ---%%
ns = m_ns;  % # of time periods
nk = m_nk;  % # of series
nl = m_nl;  % # of lags
nb = nk * (nk*nl + m_fli);  % # of coefficients in beta
na = nk * (nk-1) / 2;       % # of parameters in a

if m_fli == 1
    vym = zeros(1, nk);
else
    vym = mean(m_my);
end
m_my = m_my - ones(ns, 1) * vym;  %demeaning the series if no intercept

myh = zeros(ns, nk); 
mya = zeros(ns, nk); 
amX = zeros(nk, nb, ns); %the nk*nb design matrix for the betas over time (ns)
amXh = zeros(nk, na, ns); %the nk*na design matrix for alphas over time
amG2 = zeros(nk, nk, ns); amG2inv = zeros(nk,nk,ns); 
amOm=zeros(nk,nk,ns); amOminv=zeros(nk,nk,ns); 
for i = nl+1 : ns %creates the nk*nk design matrix for betas for each time period
    amX(:, :, i) = fXt(m_my(i-nl:i-1, :), m_fli);
    if m_flLTb==1
        amOminv(:,:,i)=eye(nk);
    else
        amOm(:,:,i)=eye(nk);
    end
end


mb = zeros(ns, nb); %matrix for holding betas
ma = zeros(ns, na); %matrix or holding alphas
mh = zeros(ns, nk); %matrix for holding latent volatility parameter (h). since exchangeable time series, one per series
mbd=mb; mad=ma;

mSigb = eye(nb) * 0.01; %AR(1) process variance for beta
mSiga = eye(na) * 0.01; %AR(1) process variance for alphas
mSigh = eye(nk) * 0.01;

%for AR case
vmub=zeros(nb,1); %level mean for betas
vmua=zeros(na,1); %level mean for alphas
mPhib=eye(nb)*.95; %underlying AR(1) process diagonal beta matrix
mPhia=eye(na)*.95; %underlying AR(1) process diagonal alpha matrix
mPhih=eye(nk)*.95;
vgam=ones(1,nk);

vdb=zeros(nb,1); %threshold for b
vda=zeros(na,1); %threshold for a

%%--- prior ---%%
dvb0 = 20;          % sigb_i^2 ~ IG(vb0/2, Vb0/2) 
dVb0 = .02;
dva0=4;             %siga_i^2 ~ IG(va0/2, Va0/2)
dVa0=.02;
dvh0=4;
dVh0=.02;

dm0=0;             %mu ~ N(m0, s0^2)
ds0=1;             
da0=20;             %(phi+1)/2 ~ Beta(a0,b0)
db0=1.5;
dg0=6;
dG0=.06;

dnub = dvb0 + ns - nl;
dnua = dva0 + ns - nl;
dnuh= dvh0 + ns - nl;

dk0=3;                    % latent threshold prior;
    
%%--- set sampling option ---%%
nburn = 0.1 * nsim;         % burn-in period
ARpmt=3*(nb+na);        % # of parameter to store
vARsamp=zeros(ARpmt,1);
ARsamp=zeros(nsim, ARpmt);  % sample box for AR parameters
volpmt=3*nk;
vvolsamp=zeros(volpmt,1);
volsamp=zeros(nsim,volpmt);
Tpmt=nb+na;
vTsamp=zeros(Tpmt,1);
Tsamp=zeros(nsim,Tpmt); % sample box for threshold parameters
betasamp=zeros(nstore,ns,nb);
alphasamp=zeros(nstore,ns,na); 
hsamp=zeros(nstore,ns,nk);

vac=zeros(4,1);
nK=floor(m_ns/30);

%%--- MCMC sampling ---%%

fprintf('\nIteration:\n');

%%------------- S A M P L I N G   S T A R T --------------%%

for m_k = -nburn : nsim
fprintf('\nIteration: %i', m_k);

  %%--- sampling beta ---%%

    if m_flLTb==1   %threshold case
        [mb(nl+1:end, :), ca] = fSampLTBv(m_my(nl+1:end,:), amX(:,:,nl+1:end), ...
                                      mPhib, vmub, mSigb, amOminv(:,:,nl+1:end), vdb, mb(nl+1:end,:));
        if m_k>0
            vac(1) = vac(1) + ca;
        end
    else  %no-threshold case
        vW=(eye(nb)-mPhib)*vmub;
        mSb0=mSigb/(eye(nb)-mPhib.^2); %b/A = b*inv(A)
        mb(nl+1:end, :) = lt_ssmooth(m_my(nl+1:end,:), amX(:,:,nl+1:end), ...
                              amOm(:,:,nl+1:end), vW, mPhib, mSigb, vmub, mSb0)';
    end
           
  %%--- sampling (Sigb, phib, mub) ---%%
  [mSigb, mPhib, vmub] = flt_SPM(mb(nl+1:end,:),mSigb, mPhib, vmub, dnub, dVb0, da0, db0, dm0, ds0, vdb, dk0);
          
  %%--- sampling threshold of beta -----%%
    if m_flLTb==1   %threshold case
        vdbc=abs(vmub)' + dk0 *sqrt(diag(mSigb/(eye(nb)-mPhib.^2))');
        [vdb, ca] = fSampLT(m_my(nl+1:end,:), amX(:,:,nl+1:end), mb(nl+1:end,:), amOminv(:,:,nl+1:end), vdb, vdbc);
        if m_k>0
            vac(2) = vac(2)+ca;
        end
    end
    mbd=abs(mb).*bsxfun(@gt,abs(mb),vdb'); %comparing values to threshold values, setting to zero if less than
    
  %%--- sampling a ---%%
    
    for i = nl+1 : ns
       myh(i, :) = m_my(i, :) - mbd(i, :) * amX(:, :, i)';
       amXh(:, :, i) = fXh(myh(i, :), nk, na);
       if m_flLTa==1
           amG2inv(:,:,i)=diag(exp(-mh(i,:))./vgam);
       else
           amG2(:, :, i) = diag(exp(mh(i, :)).*vgam);
       end
    end
    
    if m_flLTa==1   %threshold case
        [ma(nl+1:end, :), ca] = fSampLTBv(myh(nl+1:end,:), amXh(:,:,nl+1:end), ...
                                      mPhia, vmua, mSiga, amG2inv(:,:,nl+1:end), vda, ma(nl+1:end,:));
        if m_k>0
            vac(3) = vac(3) + ca;
        end
    else  %no-threshold case
        vW=(eye(na)-mPhia)*vmua;
        mSa0=mSiga/(eye(na)-mPhia.^2);
        ma(nl+1:end, :) = lt_ssmooth(myh(nl+1:end,:), amXh(:,:,nl+1:end), ...
                              amG2(:,:,nl+1:end), vW, mPhia, mSiga, vmua, mSa0)';
    end
  
   %%--- sampling (Siga, phia, mua) ---%%
   [mSiga, mPhia, vmua] = flt_SPM(ma(nl+1:end,:),mSiga, mPhia, vmua, dnua, dVa0, da0, db0, dm0, ds0, vda, dk0);
  
   %%--- sampling threshold of alpha -----%%
    if m_flLTa==1   %threshold case
        vdac=abs(vmua)' + dk0 *sqrt(diag(mSiga/(eye(na)-mPhia.^2))');
        [vda, ca] = fSampLT(myh(nl+1:end,:), amXh(:,:,nl+1:end), ...
                            ma(nl+1:end,:), amG2inv(:,:,nl+1:end), vda, vdac);
        if m_k>0
            vac(4) = vac(4)+ca;
        end
    end
    mad=abs(ma).*bsxfun(@ge,abs(ma),vda');  %comparing values to threshold values, setting to zero if less than
    
   %%--- sampling (h, Sigh, phih, muh ---%%
    for i = nl+1 : ns
        mya(i, :) = myh(i, :) * fAt(mad(i, :), nk)';
    end
    mSh0=mSigh*inv(eye(nk)-mPhih.^2);
    for i = 1 : nk
        mh(nl+1:end, i) = lt_svsamp(mya(nl+1:end,i)/sqrt(vgam(i)), mh(nl+1:end,i), ...
                                 mPhih(i,i), mSigh(i,i), 0, mSh0(i,i), nK);
    end

    
    %%--- sampling (Sigh, phih, muh) ---%%
  [mSigh, mPhih, vgam] = fSampSPG(mya(nl+1:end,:),mh(nl+1:end,:), mPhih, ...
                                  vgam, dnuh, dVh0, da0, db0, dg0, dG0);

    for i = nl+1 : ns
        mA=fAt(mad(i,:),nk);
        mAinv=inv(mA);
        if m_flLTb==1
            amOminv(:,:,i)=mA'*diag(exp(-mh(i,:))./vgam)*mA;
        else
            amOm(:,:,i)=mAinv*diag(exp(mh(i,:)).*vgam)*mAinv';
        end
    end


    if m_k > 0
        vARsamp=[vmub; diag(mPhib); diag(sqrt(mSigb)); vmua; diag(mPhia); diag(sqrt(mSiga))];
        ARsamp(m_k, :) = vARsamp';
        vvolsamp=[vgam'; diag(mPhih); diag(mSigh)];
        volsamp(m_k,:)=vvolsamp';
        vTsamp=[vdb;vda];
        Tsamp(m_k,:)=vTsamp;
        if (nsim-m_k)<(nstore+1)
            betasamp(nsim-m_k,:,:)=mbd;
            alphasamp(nsim-m_k,:,:)=mad;
            hsamp(nsim-m_k,:,:)=mh;
        end
    end
    if mod(m_k,1000) == 0       % print counter
        fprintf('%i \n', m_k);
    end

end
