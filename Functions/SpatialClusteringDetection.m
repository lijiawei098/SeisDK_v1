function [data] = SpatialClusteringDetection(CAT, SpaceLim, PATH)
% Input parameters:
% CAT - Catalog with columns: yr, mo, da, hr, mu, se, lat, lon, magn
% SpaceLim - Spatial bounds vector: [lon1, lon2, lat1, lat2]
% [~] = SpatialClusteringDetection_New( CAT, [116, 124, 38, 42] )

    % CN border
    fid = fopen(fullfile(PATH,'Data','CN-border.dat'));
    CN_temp = textscan(fid, '%s %s', 'collectoutput',1);
    fclose(fid);
    china = str2double(CN_temp{1,1}); clear fid CN_temp 

    % ---------------------------- Parameters --------------------------- %
    % Target earthquakes
    HC = [ 1975,02,04,00,00,00,7.4 ];
    TS = [ 1976,07,28,00,00,00,7.9 ];
    BH = [ 1969,07,18,00,00,00,7.4];
    % Density bin
    fthbin = 0.001; % 0.001
    % Define grid resolution
    Sbin = 0.01; % 0.01
    % End num of fth
    fth_end = 1;
    
    % Extract latitude, longitude and magnitude data
    CAT(CAT(:, 8) < SpaceLim(1) | CAT(:, 8) > SpaceLim(2) |  ...
        CAT(:, 7) < SpaceLim(3) | CAT(:, 7) > SpaceLim(4),:) = [];
    lats = CAT(:, 7);     lons = CAT(:, 8);
    mags = CAT(:, 9);
    
    % Broader spatial boundary constraints
    lon1 = SpaceLim(1)-1;    lon2 = SpaceLim(2)+1;
    lat1 = SpaceLim(3)-1;    lat2 = SpaceLim(4)+1;
    
    % ------------------------- Calculate density ----------------------- %
    % 使用Sbin来形成X和Y网格
    [X, Y] = meshgrid(linspace(lon1, lon2, (lon2-lon1)/Sbin), ...
                      linspace(lat1, lat2, (lat2-lat1)/Sbin));
    % 使用核密度估计（KDE）来估计空间密度
    density = kde2(lons, lats, mags, X, Y);
    % contour(X, Y, density, 2000)
    % colorbar

    % -------- Calculate relative change of EQ number in cluster -------- %
    % Maximum number of the cluster
    N_max = size(X,1)*size(X,2);
    % Density threshold boundary
    fth = [fthbin:fthbin:fth_end]';
    % Initialize data
    data = zeros(size(fth,1),N_max);
    % Deal with each Density threshold boundary
    for N = 250%1:size(fth,1)
        disp(['Current N:',num2str(N),'/',num2str(size(fth,1))])
        % Initial
        C = []; indx = []; region = []; NUM_EQ = [];
        % 找到fth等值线的边界
        hold on; [C, ~] = contour(X, Y, density, [fth(N) fth(N)]); hold off
        if isempty(C)
            continue
        end
        C = C';       [indx,~] = find(C(:,1) == fth(N));
        % Deal with each cluster
        for n = 1:size(indx,1)-1
            cat = [];
            region = [C(indx(n)+1:indx(n+1)-1,1),C(indx(n)+1:indx(n+1)-1,2)];
            region(end+1,:) = region(1,:);
            REGION.(sprintf('n%d', n)) = region;
            cat = CAT(logical(inPoly([CAT(:,8),CAT(:,7)],[region(:,1),region(:,2)])),:);
            NUM_EQ(n) = size(cat,1);
            % Saving Haicheng cat
            if any(ismember([cat(:, 1:6),cat(:, 9)], HC, 'rows'))
                Data_PathName = fullfile(PATH,'Outputs','Data','Haicheng',[num2str(N),'.mat']);
                save(Data_PathName,"cat","region");
            end
            % Saving Tangshan cat
            if any(ismember([cat(:, 1:6),cat(:, 9)], TS, 'rows'))
                Data_PathName = fullfile(PATH,'Outputs','Data','Tangshan',[num2str(N),'.mat']);
                save(Data_PathName,"cat","region");
            end
            % Saving Bohai cat
            if any(ismember([cat(:, 1:6),cat(:, 9)], BH, 'rows'))
                Data_PathName = fullfile(PATH,'Outputs','Data','Bohai',[num2str(N),'.mat']);
                save(Data_PathName,"cat","region");
            end
        end
        NUM_EQ = sort(NUM_EQ,'descend');
        data(N,1:size(NUM_EQ,2)) = NUM_EQ;
    end
    
    Data_PathName = fullfile(PATH,'Outputs','Data','Sigma.mat');
    save(Data_PathName,"data");
    % load(Data_PathName)

    % Pruning data   
    data(:,sum(data, 1)==0) = [];
    fth(sum(data, 2)==0,:) = [];
    data(sum(data, 2)==0,:) = [];
    
    % data 是 M 行 N 列的矩阵
    [M, N] = size(data); % 获取矩阵的大小
    % 初始化一个新的矩阵，用于存储差分结果
    diff_data = zeros(M-1, N); % 差分矩阵大小为 (M-1) x N
    % 计算每一行与上一行的差分
    data(data == 0) = 1e-10;
    for m = 2:M-1
        diff_data(m, :) = (abs( data(m, :) - data(m+1, :) ) + abs( data(m, :) - data(m-1, :) ))./data(m, :);
    end

    % ---------- Plot my relative change of EQ number in cluster -------- %
%     for m = 2:M-1
%         delta_MaxNum(m, 1) = (abs(data(m-1,1)-data(m,1)) + abs(data(m+1,1)-data(m,1)))/data(m,1);
%     end
    for NUM = 1:15
        for m = 1:M-1
            delta_MaxNum(m, 1) = sum(diff_data(m,1:NUM));
        end
        % Plot
        figure; plot(fth(2:M-1,1),delta_MaxNum(2:M-1, 1),'wo','MarkerFaceColor','r')
        hold on; plot(fth(2:M-1,1),delta_MaxNum(2:M-1, 1),'k-'); hold off
        set(gca,'YScale','log'); grid on; box on;
        % xticks([0:fthbin*100:0.14]); xlim([0,0.1])
        set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);
        xlabel('Density threshold \it f_{\rmth}','FontSize',30, 'FontName','Times New Roman')
        ylabel(['\Phi_{\itN_{\rm\Phi} \rm = ',num2str(NUM),'}'],'FontSize',30, 'FontName','Times New Roman')
        set(gcf, 'Position', [100,100,1500,500]);
        Fig_PathName = fullfile(PATH,'Outputs','Figures','cluster',['Delta_fth',num2str(NUM)]);
        print('-djpeg', '-painters', '-r300', Fig_PathName);
        close
    end

    % ------- Plot CDF of relative change of EQ number in cluster ------- %
    for M = 1:size(fth,1)-1

        if data(M,1) == 0
            continue
        end
        figure;
        
        % Plot density
        subplot(2,2,[1,3])
        contourf(X, Y, density, [0,fth(M)]);
        hold on; plot(lons,lats,'r.'); 
        contour(X, Y, density, [fth(M) fth(M)], 'LineWidth', 2, 'LineColor', 'k');
        % Plot CN border
        plot(china(:,1),china(:,2),'k','LineWidth',1)
        hold off
        xlim([SpaceLim(1),SpaceLim(2)]); ylim([SpaceLim(3),SpaceLim(4)]);
        title(['\it f_{\rmth} = \rm',num2str(fth(M))],'FontSize',15,'FontName', 'Times New Roman')
        % Generate tick labels for even numbers only
        xticks = get(gca, 'XTick');     
        yticks = get(gca, 'YTick');
        % Filter out only even ticks
        even_xticks = xticks(mod(xticks, 1) == 0);     
        even_yticks = yticks(mod(yticks, 1) == 0);
        % Format the tick labels to include the degree symbol
        xticklabels = arrayfun(@(x) sprintf('%d^{\\circ}', x), even_xticks, 'UniformOutput', false);
        yticklabels = arrayfun(@(y) sprintf('%d^{\\circ}', y), even_yticks, 'UniformOutput', false);
        % Set the filtered even ticks and labels
        set(gca, 'XTick', even_xticks, 'XTickLabel', xticklabels, 'FontName', 'Times New Roman', 'FontSize', 18);
        set(gca, 'YTick', even_yticks, 'YTickLabel', yticklabels, 'FontName', 'Times New Roman', 'FontSize', 18);
        % Optional: label the axes
        xlabel('Longitude', 'FontName', 'Times New Roman', 'FontSize', 25);
        ylabel('Latitude', 'FontName', 'Times New Roman', 'FontSize', 25);
        
        % Plot \sigma
        % 初始化矩阵，用于存储每列的累计分布函数 (CDF)
        Diff_Data_temp = [];
        Diff_Data_temp = diff_data(M,data(M,:)>10^-5);
        N = length(Diff_Data_temp);
        Diff_Data_temp(Diff_Data_temp == 0) = 10^-5;
        subplot(2,2,2)
        plot(Diff_Data_temp,1:size(Diff_Data_temp,2),'ko','MarkerFaceColor','k')
        hold on; % plot(Diff_Data_temp,1:size(Diff_Data_temp,2),'k-'); hold off; 
        grid on; box on;
        set(gca,'YScale','log');set(gca,'XScale','log');
        set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);
        xlabel('\sigma_{(\itc\rm)}','FontSize',25, 'FontName', 'Times New Roman');
        ylabel('Rank','FontSize',25, 'FontName', 'Times New Roman');
        
        % Plot number of event
        subplot(2,2,4)
        plot(data(M,1:N),1:size(data(M,1:N),2),'k^','MarkerFaceColor','k')
        hold on; plot(data(M,1:N),1:size(data(M,1:N),2),'k-'); hold off; 
        grid on; box on;
        set(gca,'YScale','log');set(gca,'XScale','log')
        set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);
        xlabel('Number of events in cluster','FontSize',25, 'FontName', 'Times New Roman');
        ylabel('Rank','FontSize',25, 'FontName', 'Times New Roman');
                
        set(gcf, 'Position', [100,100,1500,800]);
        Fig_PathName = fullfile(PATH,'Outputs','Figures','cluster',['fth_',num2str(M)]);
        print('-djpeg', '-painters', '-r300', Fig_PathName);
        close
    end
    
end


