function [f] = kde2_fixed(x, y, m, xi, yi, bandwidth)
% 使用固定带宽的核密度估计（地理距离 + 震级权重）
% x, y - 样本点经纬度
% m    - 样本点震级
% xi, yi - 网格经纬度（矩阵）
% bandwidth - 固定带宽（单位：km）

    [nrow, ncol] = size(xi);
    f = zeros(nrow, ncol);
    alpha = 1;
    Snorm = 0;
    N = length(x);

    for i = 1:N
        dist = sqrt((xi - x(i)).^2 + (yi - y(i)).^2);% haversine(x(i), y(i), xi, yi);
        weight = exp(alpha * m(i));
        kernel = weight * exp(-dist.^2 / (2 * bandwidth^2)) / (2 * pi * bandwidth^2);
        f = f + kernel;
        Snorm = Snorm + weight;
        % disp([num2str(i), '/',num2str(N)])
    end

    f = f / (Snorm);
end
