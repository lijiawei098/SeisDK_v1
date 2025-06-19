function [res] = Cal_fmd(mag,mbin,mc)
    if isempty(mc)
        mincr = [min(round(mag/mbin)*mbin):mbin:max(round(mag/mbin)*mbin)];
    else
        mincr = [mc:mbin:max(round(mag/mbin)*mbin)];
    end
    nbm = length(mincr);
    cumnbmag = zeros(nbm,1);
    nbmag = zeros(nbm,1);
    
    for n = 1:nbm
        cumnbmag(n,1) = length(find(mag > mincr(n)-mbin/2));
    end
    
    cumnbmagtmp = [cumnbmag;0];
    nbmag = abs(diff(cumnbmagtmp));
	    
    res.cumFMD  = cumnbmag';
    res.ncumFMD = nbmag';
    res.mi = mincr;
end