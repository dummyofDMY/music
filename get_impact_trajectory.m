function [a0, a1, a2, t] = get_impact_trajectory(x0, xt, v0)
    % GET_IMPACT_TRAJECTORY 计算敲击轨迹参数，采用二次轨迹，N为轨迹维数
    % INPUT:
    % x0       起始位置(1xN)
    % xt       顶点位置(1xN)
    % v0       起始速度(1xN)
    % OUTPUT:
    % [a0, a1, a2, t]    二次多项式的0、1、2次项系数，每个都是1xN，t为运动时间
    
    a0 = x0;
    a1 = v0;
    t = 4 * norm(xt - x0) ./ norm(v0);
    a2 = -v0 / t;  % 这个地方要记得，二次项系数是半个加速度
end

