function result = TestStatistic(Mag,mag,loop,Method)
    % Mag:  Mc to the end
    % mag:  Mc to Mt
    % loop: Sequence number of event being testing
    % Mag = Mag'; mag = mag';
    switch Method
        case 'SS'  % sum-sum test statistic
            result = sum(Mag(loop:end)) / sum(Mag);
        case 'SRS' % sum-robust-sum test statistic
            result = sum(Mag(loop:end)) / sum(mag);
        case 'MS'  % max-sum test statistic
            result = Mag(loop) / sum(Mag);
        case 'MRS' % max-robust-sum test statistic
            result = Mag(loop) / sum(mag);
        case 'D'   % Dixon test statistic
            result = max(Mag) / max(mag);
        case 'DK'  % dragon-king test statistic
            disp('waiting to be develop') %%%%%%
    end
end