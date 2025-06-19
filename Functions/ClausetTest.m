function  [stat,h,dobs] = ClausetTest(orgmag,a,b,Mc,dm,nsim,pval)
    %{
    % INPUTS
    orgmag = binned magnitude
    a: avalue (can be left empty)
    b: b value (can be left empty)
    Mc: magnitude of completeness
    dm: magnitude bin size
    Mmax: maximum magnitude
    nsim: number of simulations (reasonable value : >=1000)
    pval: minimum fraction of synthetic cases with ks distance larger in the
    synthetic catalogs than in the real catalog
    
    OUTPUTS:
    stat: fraction of synthetic cases with ks distance larger in the
    synthetic catalogs than in the real catalog
    h :  1 or 0 depending on whether GR law can be rejected or accepted
    dobs: ks distance between the emprical and best fit GR law
    %}
    
    %%
    Mmax = max(orgmag);
    orgmag(orgmag<Mc) = [];
    [C,count] = Count_Unique(orgmag);
    orgcdf = cumsum(count)/sum(count);
    
    % estimate of a and b
    if isempty(a) || isempty(b)
        [beta,a] = Tinti_aval(orgmag,Mc,dm,200);
        b = beta/log(10);
    end
    
    NC = BinMags([Mc:dm:Mmax+dm]',0,dm);
    [~,~,binid] = histcounts(NC,[C;Mmax]);
    
    indo = find(binid~=0);
    Norgcdf = orgcdf(binid(indo));
    
    theocount = 10.^(a-b*(NC(indo)-dm/2)) - 10.^(a-b*(NC(indo)+dm/2));
    theocdf = cumsum(theocount)/sum(theocount);
    
    dobs = max(abs(Norgcdf-theocdf));
    
    Nevt = length(orgmag);
    for i = 1:nsim
        
        [mag] = GR_truncated_simulator_ver2(Nevt,Mmax,Mc-dm/2,b);
        mag = BinMags(mag',0,dm);
        
        [C,count] = Count_Unique(mag);
        
        simcdf = cumsum(count)/sum(count);
        [~,~,binid] = histcounts(NC,[C;Mmax+dm]);
        
        inds = find(binid~=0);
        Nsimcdf = simcdf(binid(inds));
        
        [beta,aval,~] = Tinti_aval(mag,Mc,dm,0);
        
        simb = beta/log(10);
        
        theocount = 10.^(aval-simb*(NC(inds)-dm/2)) - 10.^(aval-simb*(NC(inds)+dm/2));
        theocdf = cumsum(theocount)/sum(theocount);
        
        dsim(i,1) =  max(abs(Nsimcdf-theocdf));
        %[beta1(i),aval1(i),LL_sim(i)] = Tinti_aval(mag,Mc,dm,0);
        check = 1;
    end
    
    % dsim(dsim==0)=[];
    stat = sum(dsim>dobs)/length(dsim);
    h = stat<pval;
    
    check = 1;
end