function [udat,count,iaf,ial,dat] = Count_Unique(dat)
    % INPUTS:
    % dat: data point
    %
    % OUTPUTS:
    % udat: data sorted in ascending order after eliminating duplicates
    % count: count corresponding to udat
    % iaf: index vector
    % ial: index vector
    
    dat = sort(dat) ; % in ascending order
    [udat, iaf] = unique(dat) ;
    [~, ial] = unique(dat,'last') ;
    count = ial-iaf+1;

end