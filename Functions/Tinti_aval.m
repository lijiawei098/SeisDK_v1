function [b,aval,LL] = Tinti_aval(mag,M0,dm,bs)
    % INPUTS:
    % mag: data
    % M0: cut-off magnitude
    % dm: magnitude bin
    % bs: flag for bootstrap
    
    
    mag(mag<M0) = [];
    
%     mmag = sum(mag-M0)/length(mag); % average magnitude
%     beta = 1/dm*log(1+dm/mmag);
%     b = beta/log(10);
    b = log10( exp(1) ) / ( mean( mag( mag > M0 - dm/2 ) ) - (M0-dm/2) );
    beta = b*log(10);
    
    
    [umag,iaf]=unique(sort(mag));
    [umag,ial]=unique(sort(mag),'last');
    count=ial-iaf+1;
    
    aval=log10( sum(count) / sum(10.^(-b*(umag-dm/2))) / (1-10^(-b*dm)) );
    LL = sum((log(1 - exp(-beta *dm)) - beta * (mag - M0)));
    
    check = 1;
    
    %% bootstrap for error
    % if bs==1
    %     Nbs=1000;
    %     betabs=zeros(Nbs,1);
    %     avalabs=zeros(Nbs,1);
    %     for i=1:Nbs
    %         
    %         rndind = randi(length(mag),length(mag),1);
    %         nmag=mag(rndind);
    %         mmag=sum(nmag-M0)/length(mag);
    %         betabs(i)=1/dm*log(1+dm/mmag);
    %         b=betabs(i)/log(10);
    %         
    %         [umag,iaf]=unique(sort(nmag));
    %         [umag,ial]=unique(sort(nmag),'last');
    %         count=ial-iaf+1;
    %         
    %         avalabs(i)=log10(sum(count)/sum(10.^(-b*(umag-dm/2)))...
    %             /(1-10^(-b*dm)));
    %     end
    %     [betaquants]=quantile(betabs,[0.025,0.975]);
    %     [avalquants]=quantile(avalabs,[0.025,0.975]);
    %     
    %     beta=[beta,betaquants];
    %     aval=[aval,avalquants];
    %     check=1;
    % end
    
    check=1;
    
    %[beta,aval]=Tinti(mag,M0,dm)
end