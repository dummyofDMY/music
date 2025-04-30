function [a0, a1, a2, a3] = get_transfer_trajectory(x0, xf, v0, vf, t)
    % GET_TRANSFER_TRAJECTORY 计算两个撞击轨迹间的转移轨迹参数，采用三次规划，N为轨迹维数
    % INPUT:
    % x0       起始位置(1xN)
    % xf       终止位置(1xN)
    % v0       起始速度(1xN)
    % vf       终止速度(1xN)
    % t        走完轨迹所用时间
    % OUTPUT:
    % [a0, a1, a2, a3]    三次多项式的0、1、2、3次项系数，每个都是1xN

    a0 = x0;
    a1 = v0;
    h = xf - x0;
    a2 = (3 * h - (2 * v0 + vf) * t) / t^2;
    a3 = (-2 * h + (v0 + vf) * t) / t^3;
end
