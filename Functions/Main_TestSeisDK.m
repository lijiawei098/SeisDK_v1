function [ Mc, mt, DataTable ] = Main_TestSeisDK( CAT, Method, Statistics, alpha )
% Input:
% MAG:        magnitude
% Method:     Inward, Outward
% Statistics: SS(block), SRS(block), MS(Inward/Outward), MRS(Inward/Outward), D, DK
% alpha:      Significance level
%
% Output:

    MAG = CAT(:,9);
    
    % Initial parameters
    DK = nan; mag = nan; 
    Mt = nan; mt  = [];
    Mc = [];  mc = [];
    B = [];
    p = nan;

    % Parameters
    setting.mag_lim = [floor(min(MAG)*2)/2, ceil(max(MAG)*2)/2];           % Low and up limit of magnitude
    % Rank magnitude in ascending order
    MAG = sort(MAG);
    CAT = sortrows(CAT, 9);
    

    % -------------- Framework I. Estimate mc, mt, β, and r ------------- %
    disp(['% =========================== Framework I. Estimate mc, mt, β,' ...
        ' and r ============================= %'])
    % Estimate mc, mt, beta
    [ mc, mt, b ] = Cal_mc_mt_New( MAG );   
    if mc < min(MAG)
        mc = min(MAG);
    end
    if ~any(MAG==mc)
        mc_temp = MAG(MAG-mc > 0);
        mc = mc_temp(1);
    end
    if ~any(MAG==mt)
        mt_temp = MAG(MAG-mt > 0);
        mt = mt_temp(1);
    end
    Mc = mc; 
    setting.mc = mc;  setting.mt = mt;  setting.b = b;
            

    % ------------ Framework II. Dragon-king earthquake test ------------ %
    disp(['% =========================== Framework II. Dragon-king' ...
        ' earthquake test =========================== %'])
    mag = MAG(MAG>=mc & MAG < mt);
    Mag = MAG(MAG>=mc);
    Data = Mag;                                                            % Data: [Mag, Statistic, S_critical, p, DKflag, b-value]
    Data(:,9:14) = CAT(MAG>=mc,1:6);
    OC_ind = find(Mag>=mt);                                                % Index of outlier candidates
    
    switch Method
        case 'Inward' % Inward test
            for loop = size(Mag,1):-1:1
                loop;
                Mag_temp = Mag;
                % beta value for excluding outlier tested
                % [setting.b,~] = Tinti_aval(Mag_temp(1:loop-1),mc,0.1,0); 
                [setting.b,~] = Tinti_aval(mag,mc,0.1,0);
                % calculate test statistic; Methods: SS, SRS, MS, MRS, D, DK
                statistic = TestStatistic(Mag_temp,mag,loop,Statistics);
                setting.Num = length(mag);
                Data(loop,2) = statistic;
                % calculate critical value
                [Sc,testStatistics] = Cal_CriticalValue(loop, Mag_temp, alpha, setting, Statistics); % alpha = 0.01
                Data(loop,3) = Sc;
                p = length(testStatistics(testStatistics>statistic))/length(testStatistics);
                Data(loop,4) = p;
                if statistic > Sc
                    statistic > Sc
                    Data(loop,5) = 1;
                    Mt = Mag(loop);
                    Data(loop,6) = Tinti_aval(Mag_temp(1:loop-1),mc,0.1,0);
                    DK = Mag(loop:end);
                else
                    Data(loop,6) = Tinti_aval(Mag_temp(1:loop),mc,0.1,0);
                    break                    
                end                
            end
        case 'Outward' % Outward test
            for loop = OC_ind(1):OC_ind(end)
                % beta value for excluding outlier tested
                [beta,~] = Tinti_aval(Mag_temp(1:loop-1),mc,0.1,0); 
                % calculate test statistic; Methods: SS, SRS, MS, MRS, D, DK
                statistic = TestStatistic(Mag_temp,mag,loop,Statistics);
                % calculate critical value
                [Sc,BK_test] = Cal_CriticalValue(loop, Mag_temp, alpha, path, setting, Statistics); % alpha = 0.01
                if statistic <= Sc && loop ~= OC_ind(end)
                    continue
                elseif statistic <= Sc && loop == OC_ind(end)
                    disp('There is no DK!');
                    return
                else
                    setting.Mt = Mag(loop);
                    setting.beta = beta;
                    DK = Mag(loop:end);
                    BK = BK_test;
                    ST = statistic;
                    disp(['Testing DK: NumDK = ',num2str(length(DK)),', k = ', num2str(loop+1),', Mt = ', num2str(setting.Mt),' ...']);
                    break
                end
            end
    end
       

    % -------- Framework III. Dragon-king earthquake validation --------- %
    disp(['% =========================== Framework III. Dragon-king' ...
        ' earthquake validation ==================== %'])
    % Numerical number of DK in 
    ind_NumDK = find(Mag >= Mt);
    if isempty(ind_NumDK)
        disp('There is no DK!');
    elseif length(ind_NumDK) > 1
        statistic = TestStatistic(Mag,mag,ind_NumDK(1),'SRS');
        [setting.b,~] = Tinti_aval(Mag(1:ind_NumDK(1)-1),mc,0.1,0);
        [~,BK] = Cal_CriticalValue(ind_NumDK(1), Mag, alpha, setting, 'SRS');
        Data(ind_NumDK(1),7) = length(BK(BK>statistic))/length(BK);
        Data(ind_NumDK(1),8) = alpha;
    elseif length(ind_NumDK) == 1
        statistic = TestStatistic(Mag,mag,ind_NumDK,Statistics);
        [setting.b,~] = Tinti_aval(Mag(1:ind_NumDK-1),mc,0.1,0);
        [~,BK] = Cal_CriticalValue(ind_NumDK, Mag, alpha, setting, Statistics);
        Data(ind_NumDK,7) = length(BK(BK>statistic))/length(BK);
        Data(ind_NumDK,8) = alpha;
    end 
    DataNames = {'Mag','Statistic','S_critical','p-value','DKflag',...
        'b_value','p-value validation','alpha',...
        'year','month','day','hour','minute','second'};
    DataTable = array2table(Data, 'VariableNames', DataNames);
end