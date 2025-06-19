function [f] = kde2(x, y, m, xi, yi)
% 自适应带宽的二维核密度估计 (带震级权重, 地理单位 km)
% x, y - 地震经纬度 (列向量, 单位: 度)
% m    - 地震震级 (列向量)
% xi, yi - 网格点经纬度（矩阵）
% f - 每个网格点的密度估计值（单位: 1/km^2）

    alpha = 1;               % productivity exponent
    h0 = 0.01;                 % 全局固定带宽，单位：km
    epsilon = 1e-10;

    % Step 1: 计算局部密度（固定带宽）并确定自适应带宽
    local_density = kde2_fixed(x, y, m, x, y, h0);
    adaptive_bandwidth = h0 ./ sqrt(local_density + epsilon);  % Abramson 带宽

    % Step 2: 核密度估计（带震级 + 自适应带宽）
    [nrow, ncol] = size(xi);
    f = zeros(nrow, ncol);
    Snorm = 0;
    N = length(x);

    for i = 1:N
        % 计算样本点 i 到所有网格点的地理距离（单位：km）
        dist = sqrt((xi - x(i)).^2 + (yi - y(i)).^2);
        % dist_km = haversine(x(i), y(i), xi, yi);
        h_i = adaptive_bandwidth(i);  % 单位：km
        weight = exp(alpha * m(i));
        kernel = weight .* exp(-dist.^2 / (2 * h_i^2)) ./ (2 * pi * h_i^2);
        f = f + kernel;
        Snorm = Snorm + weight;
        % disp([num2str(i), '/',num2str(N)])
    end

    f = f / (Snorm);  % 归一化，使密度积分为 1
end
