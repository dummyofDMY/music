function theta = subproblem1(Xi, r, p, q)
    % SUBPROBLEM1 求解子问题一
    % INPUT:
    % Xi       转轴旋量
    % r        转轴上一点
    % p        起始点
    % q        终止点
    % OUTPUT:
    % res      求得的结果，已经规范化了

    omega = Xi(4:6)';  % 转轴为列向量
    omega = omega / norm(omega);  % 单位化转轴

    % 定义输入都是行向量
    u = (p - r)';
    v = (q - r)';
    du = omega * omega' * u;
    dv = omega * omega' * v;
    u_ = u - du;
    v_ = v - dv;  % 计算投影
    
    tem1 = norm(u) - norm(v);
    tem2 = omega' * u - omega' * v;
    if norm(u) - norm(v) > 1e-4 || omega' * u - omega' * v > 1e-4
%         disp('no solution for subproblem1')
        theta = NaN;
    else
        if norm(u_) < 1e-4 || norm(v_) < 1e-4
            disp('WARNING: subproblem1 singular!!!')
            theta = 0;
        else
            theta = atan2((omega' * cross(u_, v_)), (u_' * v_));  % 计算 theta
            theta = leagalize_theta(theta);  % 归一化 theta
        end
    end
end