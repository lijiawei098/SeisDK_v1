clc; clear; close;
% ================================================= %
% Parameters
TAG     = 'Haicheng';      % Haicheng; Tangshan; Bohai
PART    = 'Foreshocks';    % Foreshocks; Aftershocks; All
alpha   = 0.05;            % Significance level 0.05
Method  = 'Inward';
fthbin  = 0.001;           % Density bin
fth_end = 1;               % End num of fth
% ================================================= %

% Path and functions
PATH = '.\Seis_DK\';
addpath(fullfile(PATH,'Functions/'));

% Adjust parameter
if strcmp(TAG, 'Haicheng')
    Target_event = [ 1975,02,04,00,00,00,7.4 ];
elseif strcmp(TAG, 'Tangshan')
    Target_event = [ 1976,07,28,00,00,00,7.9 ];
elseif strcmp(TAG, 'Bohai')
    Target_event = [ 1969,07,18,00,00,00,7.4 ];
end

% Density threshold boundary
fth = [fthbin:fthbin:fth_end]';

for LOOP = 1:size(fth,1)
    % ------------------------ Prepare catalog -------------------------- %
    % Load data
    Data_PathName = fullfile(PATH,'Outputs','Data',TAG,num2str(LOOP));
    load(Data_PathName);
    % cat(cat(:,1)<1900,:) = [];
    % Prepare catalog
    if strcmp(PART, 'Foreshocks')
        CAT = cat(datenum(cat(:,1:6)) <= datenum(Target_event(1:6)), :);
    elseif strcmp(PART, 'Aftershocks')
        CAT = cat(datenum(cat(:,1:6)) >= datenum(Target_event(1:6)), :);
    elseif strcmp(PART, 'All')
        CAT = cat;
    end
    
    % ------------------------------ DK TESTs --------------------------- %
    % 对震级值进行排序并获取排名
    [MAG_sorted, rank_order] = sort(CAT(:,9), 'descend'); % 排序震级值
    rank_values = 1:length(MAG_sorted); % 生成排名值
    % DK tests
    [ Mc, mt, Data ] = Main_TestSeisDK( CAT, Method,'MRS', alpha );
    Data = flipud(Data);
    
    % ------------------------------- Plot ------------------------------ %
    figure; hold on;
    % Plot incompleteness part
    ind_Mc = find(MAG_sorted == Mc);
    scatter(MAG_sorted(ind_Mc(end)+1:end), rank_values(ind_Mc(end)+1:end), 100, 'filled', ...
        'MarkerFaceColor', [230, 230, 230] ./ 255, ...
        'MarkerEdgeColor', [200, 200, 200] ./ 255, 'LineWidth',0.3);
    plot(MAG_sorted(ind_Mc(end):end), rank_values(ind_Mc(end):end), ...
        '-', 'LineWidth', 1.5,'Color',[200, 200, 200] ./ 255);
    xline(MAG_sorted(ind_Mc(end))-0.1/2, '-', 'LineWidth', 1.5, 'Color', [100, 100, 100] ./ 255);
    % Plot power-law part
    ind_mt = find(MAG_sorted == mt);
    scatter(MAG_sorted(ind_mt(end)+1:ind_Mc(end)), rank_values(ind_mt(end)+1:ind_Mc(end)), 100, 'filled', ...
        'MarkerFaceColor', [180, 180, 180] ./ 255, ...
        'MarkerEdgeColor', 'k', 'LineWidth',0.3);
    plot(MAG_sorted(ind_mt(end):ind_Mc(end)), rank_values(ind_mt(end):ind_Mc(end)), ...
        '-', 'LineWidth', 1.5,'Color',[180, 180, 180] ./ 255);
    xline(MAG_sorted(ind_mt(end))-0.1/2, '--', 'LineWidth', 1.5, 'Color', [0, 0, 0] ./ 255);
    % Plot outlier candidate part
    scatter(MAG_sorted(1:ind_mt(end)), rank_values(1:ind_mt(end)), 100, 'filled', ...
        'MarkerFaceColor', 'k', ...
        'MarkerEdgeColor', 'k', 'LineWidth',0.3);
    plot(MAG_sorted(1:ind_mt(end)), rank_values(1:ind_mt(end)), ...
        '-', 'LineWidth', 1.5,'Color','k');
    % Plot DK part
    b = Data{1,6};
    if any(Data{:,5} == 1)                                                 % DK event
        ind_DK = find(Data{:,5}==1);
        b = Data{ind_DK(end)+1,6};
        Mt = MAG_sorted(ind_DK(end));
        xline(Mt-0.1/2, '-', 'LineWidth', 1.5, 'Color', [0, 0, 0] ./ 255);
        scatter(MAG_sorted(1:ind_DK(end)), rank_values(1:ind_DK(end)), 100, 'filled', ...
            'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'LineWidth',0.3);
    end
    
    % Plot GR
    % 根据GR关系计算对应的震级直线
    rank_mc = rank_values(MAG_sorted==Mc);
    a = b*Mc+log10(rank_mc(end));
    predicted_M = (a - log10(rank_values)) / b;  % 根据GR关系计算震级
    % 绘制基于GR关系计算出的震级直线
    plot(MAG_sorted(1:ind_Mc(end)), 10.^(a-MAG_sorted(1:ind_Mc(end))*b), '--', 'LineWidth', 1.5, 'Color', [255, 69, 0] / 255);  % 使用虚线表示预测直线
    
    % Plot Target event
    IND = find( datenum(Data{:,9:14}) == datenum(Target_event(1:6)));
    scatter(MAG_sorted(IND), rank_values(IND), 400, 'filled', ...
            'h', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k', 'LineWidth',0.3);

    % Legend
    if ~any(Data{:,5} == 1)                                                      % No DK event
        legend('Incomplete', '', ['$\mathit{m}_\mathrm{c}=' sprintf('%.1f', Mc) '$'], ...
           'Gutenberg-Richter', '', ['$\mathit{m}_\mathrm{t}=' sprintf('%.1f', mt) '$'], ...
           'DK candidate', '', ['$\mathrm{rank}(\mathit{m}) = 10^{' sprintf('%.2f', a) ' - ' sprintf('%.2f', b) '\mathit{m}}$'], ...
           [TAG ' $\mathit{m}_\mathrm{L} = ' num2str(MAG_sorted(IND)) ' \ \mathrm{earthquake}$' ' (' sprintf('%.4f', Data{1,4}),')'], ...
           'Interpreter', 'latex','Location','SouthOutside','NumColumns',2);
    elseif any(Data{:,5} == 1)                                                   % DK event
        legend('Incomplete', '', ['$\mathit{m}_\mathrm{c}=' sprintf('%.1f', Mc) '$'], ...
           'Gutenberg-Richter', '', ['$\mathit{m}_\mathrm{t}=' sprintf('%.1f', mt) '$'], ...
           'DK candidate', '', ['$\mathit{m}_\mathrm{DK}=' sprintf('%.1f', Mt) '$'], 'DK', ...
           ['$\mathrm{rank}(\mathit{m}) = 10^{' sprintf('%.2f', a) ' - ' sprintf('%.2f', b) '\mathit{m}}$'], ...
           [TAG ' $\mathit{m}_\mathrm{L} = ' num2str(MAG_sorted(IND)) ' \ \mathrm{earthquake}$' ' (' sprintf('%.4f', Data{1,4}),')'], ...
           'Interpreter', 'latex','Location','SouthOutside','NumColumns',2);
    end
    
    
    hold off; grid on; box on;
    set(gcf, 'Position', [100,100,800,800]);
    set(gca,'YScale','log')
    current_ylim = ylim;
    ylim([1, current_ylim(2)])
    % Set the filtered even ticks and labels
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);
    % Optional: label the axes
    xlabel('Magnitude', 'FontName', 'Times New Roman', 'FontSize', 25);
    ylabel('Rank', 'FontName', 'Times New Roman', 'FontSize', 25);
    title([TAG,' earthquake sequence: ',PART,'; ','\it f_{\rmth} = \rm',num2str(fth(LOOP))], ...
        'FontName', 'Times New Roman', 'FontSize', 20)
        
    % Save figure
    Fig_PathName = fullfile(PATH,'Outputs','Figures',TAG,PART,num2str(LOOP));
    print('-djpeg', '-painters', '-r300', Fig_PathName);
    close
    % Saving data
    Data_PathName = fullfile(PATH,'Outputs','Data',[TAG,'_DK'],PART,num2str(LOOP));
    save(Data_PathName,"Data","Mc");
end
rmpath(fullfile(PATH,'Functions/'));