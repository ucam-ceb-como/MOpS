function [D,DgofD] = kernelest_w(m0,dn,sigma,npts,D,w)
    
    N     = length(dn);
    if isempty(D)
        D = (linspace(min(dn)*.25,max(dn)*2,npts))';
    else
        npts = length(D);
    end
    c1    = 1/(sigma*sqrt(2*pi));
    c2    = -1/(2*sigma^2);
    mu    = log(dn);
    lnD   = log(D);
    fofD  = zeros(npts,1);
    
    for i = 1:N
        fofD = fofD+(1./D).*exp(c2*(lnD-mu(i)).^2)*w(i);
    end
    
    fofD = fofD/sum(w);
    
    fofD  = fofD*c1;
    gofD  = m0*fofD;
    DgofD = gofD.*D;
    
    