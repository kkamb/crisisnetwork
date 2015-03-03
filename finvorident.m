function omatrix = finvorident(imatrix,isize)
    drD=rcond(imatrix);
    if isnan(drD) || (drD<eps*10^2)
        omatrix=eye(isize)*10;
    else
        omatrix=inv(imatrix);
    end