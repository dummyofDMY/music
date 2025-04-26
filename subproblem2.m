function theta = subproblem2(Xi1, Xi2, r, p, q)
    omega1 = Xi1(4:6)';  % 转轴为列向量
    omega2 = Xi2(4:6)';  % 转轴为列向量
    omega1 = omega1 / norm(omega1);  % 单位化转轴
    omega2 = omega2 / norm(omega2);  % 单位化转轴
    c = cross(omega1, omega2);

    % 定义输入都是行向量
    u = (p - r)';
    v = (q - r)';
    
    alpha = ((omega1' * omega2) * omega2' * u - omega1' * v) / ((omega1' * omega2)^2 - 1);
    beta = ((omega1' * omega2) * omega1' * v - omega2' * u) / ((omega1' * omega2)^2 - 1);
    gamma2 = (u'*u - alpha^2 - beta^2 - 2*alpha*beta*(omega1'*omega2)) / norm(cross(omega1, omega2))^2;
    theta = zeros(2, 2);
    if gamma2 < 0
%         disp('no solution for subproblem2');
        theta(:, :) = NaN;
    else
        if gamma2 < 1e-6 && (alpha < 1e-6 || beta < 1e-6)
            disp('WARNING: subproblem2 singular!!!')
        end
        gamma = gamma2^0.5;
        z1 = alpha * omega1 + beta * omega2 + gamma * c;  % 计算 z 向量
        z2 = alpha * omega1 + beta * omega2 - gamma * c;  % 计算 z 向量
        theta(1, 1) = subproblem1(Xi1, r, z1' + r, q);  % 计算 theta1
        theta(1, 2) = subproblem1(Xi2, r, p, z1' + r);  % 计算 theta2
        theta(2, 1) = subproblem1(Xi1, r, z2' + r, q);  % 计算 theta3
        theta(2, 2) = subproblem1(Xi2, r, p, z2' + r);  % 计算 theta4
    end
end