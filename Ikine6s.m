function theta = Ikine6s(T)
    % IKINE6S 这是六轴机器人逆运动学函数，适用于六轴机器人
    % INPUT:
    % T                机器人末端姿态矩阵，4x4矩阵    
    % OUTPUT:
    % theta           关节转角矩阵，8x6，单位rad，已经规范化了

    L1 = 491; L2 = 450; L3 = 450; L4 = 84;
    q1 = [0, 0, 0]; q2 = [0, 0, L1];
    q3 = [0, 0, L1+L2]; q4 = [0, 0, L1+L2];
    q5 = [0, 0, L1+L2+L3]; q6 = [0, 0, L1+L2+L3];
    w1 = [0, 0, 1]; w2 = [0, 1, 0];
    w3 = [0, 1, 0]; w4 = [0, 0, 1];
    w5 = [0, 1, 0]; w6 = [0, 0, 1];
    
    g0 = [-1 0 0 0; 
        0 -1 0 0; 
        0 0 1 L1+L2+L3+L4; 
        0 0 0 1];
    
    Xi1 = zeros(1, 6); Xi1(1:3) = cross(q1, w1); Xi1(4:6) = w1;
    Xi2 = zeros(1, 6); Xi2(1:3) = cross(q2, w2); Xi2(4:6) = w2;
    Xi3 = zeros(1, 6); Xi3(1:3) = cross(q3, w3); Xi3(4:6) = w3;
    Xi4 = zeros(1, 6); Xi4(1:3) = cross(q4, w4); Xi4(4:6) = w4;
    Xi5 = zeros(1, 6); Xi5(1:3) = cross(q5, w5); Xi5(4:6) = w5;
    Xi6 = zeros(1, 6); Xi6(1:3) = cross(q6, w6); Xi6(4:6) = w6;

    q1_ = [q1, 1]'; q2_ = [q2, 1]'; q3_ = [q3, 1]';
    q4_ = [q4, 1]'; q5_ = [q5, 1]'; q6_ = [q6, 1]';

    % 初始化解
    theta_solve = zeros(8, 6);
    theta_solve(:, :) = NaN;
    g = squeeze(T) / g0;
    % 求解theta3
    gq5_q2 = g * q5_ - q2_;
    gq5_q2 = gq5_q2(1:3, 1)';
    sigma1 = norm(gq5_q2);
    solve1 = subproblem3(Xi3, q3, q5, q2, sigma1);
    theta_solve(1:4, 3) = solve1(1,1);
    theta_solve(5:8, 3) = solve1(2,1);

    for n = 1:2
        if isnan(solve1(n, 1))
            continue
        end
        % 求解theta1、theta2
        e3 = screwToTransformationMatrix(Xi3, solve1(n, 1));
        e3_q5 = e3 * q5_;
        e3_q5 = e3_q5(1:3, 1)';
        gq5 = g * q5_;
        gq5 = gq5(1:3, 1)';
        solve2 = subproblem2(Xi1, Xi2, q2, e3_q5, gq5);
        theta_solve((n-1)*4+1:(n-1)*4+2, 1) = solve2(1, 1);
        theta_solve((n-1)*4+3:(n-1)*4+4, 1) = solve2(2, 1);
        theta_solve((n-1)*4+1:(n-1)*4+2, 2) = solve2(1, 2);
        theta_solve((n-1)*4+3:(n-1)*4+4, 2) = solve2(2, 2);
        for k = 1:2
            if any(isnan(solve2(k, :)))
                continue
            end
            % 求解theta4、theta5
            e123 = screwToTransformationMatrix(Xi3, solve1(n, 1));
            e123 = screwToTransformationMatrix(Xi2, solve2(k, 2)) * e123;
            e123 = screwToTransformationMatrix(Xi1, solve2(k, 1)) * e123;
            en1n2n3_gq1 = e123 \ g * q1_;
            en1n2n3_gq1 = en1n2n3_gq1(1:3, 1)';
            solve3 = subproblem2(Xi4, Xi5, q5, q1, en1n2n3_gq1);
            theta_solve((n-1)*4+(k-1)*2+1, 4) = solve3(1, 1);
            theta_solve((n-1)*4+(k-1)*2+2, 4) = solve3(2, 1);
            theta_solve((n-1)*4+(k-1)*2+1, 5) = solve3(1, 2);
            theta_solve((n-1)*4+(k-1)*2+2, 5) = solve3(2, 2);
            for m = 1:2
                if any(isnan(solve3(m, :)))
                    continue
                end
                % 求解theta6
                tem = [1, 0, 0];
                tem_ = [tem, 1]';
                e12345 = screwToTransformationMatrix(Xi5, solve3(m, 2));
                e12345 = screwToTransformationMatrix(Xi4, solve3(m, 1)) * e12345;
                e12345 = screwToTransformationMatrix(Xi3, solve1(n, 1)) * e12345;
                e12345 = screwToTransformationMatrix(Xi2, solve2(k, 2)) * e12345;
                e12345 = screwToTransformationMatrix(Xi1, solve2(k, 1)) * e12345;
                en1n2n3n4n5_gtem = e12345 \ g * tem_;
                en1n2n3n4n5_gtem = en1n2n3n4n5_gtem(1:3)';
                theta_solve((n-1)*4+(k-1)*2+m, 6) = subproblem1(Xi6, q6, tem, en1n2n3n4n5_gtem);
            end
        end
    end
    theta = theta_solve;
end