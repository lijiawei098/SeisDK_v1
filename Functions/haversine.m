function d = haversine(lon1, lat1, lon2, lat2)
% Haversine 公式计算两点之间地球大圆距离 (单位: km)
% 输入单位: 度
    R = 6371;  % 地球半径 (km)
    lon1 = deg2rad(lon1);
    lat1 = deg2rad(lat1);
    lon2 = deg2rad(lon2);
    lat2 = deg2rad(lat2);

    dlon = lon2 - lon1;
    dlat = lat2 - lat1;

    a = sin(dlat/2).^2 + cos(lat1).*cos(lat2).*sin(dlon/2).^2;
    c = 2 * atan2(sqrt(a), sqrt(1 - a));
    d = R * c;
end
