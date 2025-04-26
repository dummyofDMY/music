function theta = subproblem3(Xi, r, p, q, sigma)
    % 定义输入都是行向量
    u = (p - r)';
    v = (q - r)';
    omega = Xi(4:6)';  % 转轴为列向量
    omega = omega / norm(omega);  % 单位化转轴
    u_ = u - omega * omega' * u;
    v_ = v - omega * omega' * v;  % 计算投影

    c = (u_'*u_ + v_'*v_ - sigma^2) / (2 * norm(u_) * norm(v_));
    if norm(u_) < 1e-6 || norm(v_) < 1e-6
        disp('WARNING: subproblem3 singular!!!')
        if c < -1 || c > 1
            theta = zeros(2, 1);
            theta(:, 1) = NaN;
        else
            theta = zeros(2, 1);
        end
    else
        theta0 = atan2((omega' * cross(u_, v_)), (u_'*v_));
        if -1 - c > 1e-6 || c - 1 > 1e-6
    %         disp('no solution for subproblem3');
            theta = zeros(2, 1);
            theta(:, 1) = NaN;
        else
            c = max(min(c, 1), -1);
            dtheta = acos(c);
            theta = zeros(2, 1);  % 初始化 theta 向量
            theta(1, 1) = leagalize_theta(theta0 + dtheta);  % 第一个解
            theta(2, 1) = leagalize_theta(theta0 - dtheta);  % 第二个解
        end
    end
end