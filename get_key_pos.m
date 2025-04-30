function [gsts] = get_key_pos(l)
    % GET_KEY_POS 优化算法得出最优关键点的位姿
    % INPUT:
    % l             木琴两个音条中轴线的距离(标量)
    % OUTPUT:
    % gsts          关键点位姿(Nx4x4)

    % 初始猜测值 (需为单位四元数)
    typical_q = [0.030982429914554,-0.100066112070343,0.994442974947455,-0.010489605777688];
    typical_p = [5.122186101582196e+02,3.133945035755696,4.610435346606103e+02];
    x0 = typical_p; % 示例: 表示绕y轴旋转90度
    for i = 1:7
        x0 = [x0, typical_q];
    end
    
    % 设置非线性约束
    nonlcon = @constraint;
    
    % 创建优化选项
    options = optimoptions('fmincon', 'Algorithm', 'sqp', ...
        'Display', 'iter', 'MaxIterations', 400, 'MaxFunctionEvaluations', 1500);
    
    % 调用fmincon求解
    [x, fval] = fmincon(@obj_fun, x0, [], [], [], [], [], [], nonlcon, options);
    
    % 输出结果
    disp(['最优解(单位四元数): ', num2str(x)]);
    disp(['最小化的目标函数值: ', num2str(fval)]);

    q = zeros(N, 4);
    for i = 1:N
        q(i, :) = x(1, (4 * i):(4 * i + 3));
    end
    gsts = zeros(N, 4, 4);
    for i = 1:N
        R = quaternion_to_rotation_matrix(q(i, :));
        now_p = p + [0, (i - (N + 1) / 2) * l, 0];
        gsts(i, :, :) = [R, now_p';
                         [0, 0, 0], 1];
    end

    save('best_pos', "gsts");
end

function f = obj_fun(x)
    % OBJ_FUN 评价函数，衡量一组最优关键点位姿的优秀程度
    % INPUT:
    % x         待优化参数(1x(3+4N))，内容为[p, q1, q2, ..., qN]
    %           其中p(1x3)为关键点笛卡尔坐标均值，
    %           q为各关键点姿态的四元数(1x4)
    % OUTPUT:
    % f         输出评价指标（越小越好）

    l = 20;  % 木琴两个音条中轴线的距离(标量)
    p = x(1, 1:3);
    N = (length(x) - 3) / 4;
    q = zeros(N, 4);
    for i = 1:N
        q(i, :) = x(1, (4 * i):(4 * i + 3));
    end
    gsts = zeros(N, 4, 4);
    for i = 1:N
        R = quaternion_to_rotation_matrix(q(i, :));
        now_p = p + [0, (i - (N + 1) / 2) * l, 0];
        gsts(i, :, :) = [R, now_p';
                         [0, 0, 0], 1];
    end
    tem0 = squeeze(gsts(1, :, :));
    thetas = get_nearest_theta(gsts);
    tem = rad2deg(thetas);
    if any(isnan(gsts))
        f = inf;
        return;
    end
    % 计算几个特殊位姿在关节坐标系下的速度
    % 注意这个速度是竖直向下的
    key_dtheta = zeros(N, 6);
    v = 1000;
    for i = 1:N
        dg = zeros(4, 4);
        dg(3, 4) = -v;
        V_hat = dg / squeeze(gsts(i, :, :));
        V = anti_se3_hat(V_hat);
        J = Jacobian(thetas(i, :));
        key_dtheta(i, :) = J \ V';
    end
    f = 0;
    td = 0.5;
    for i = 1:N
        for j = (i + 1):N
            th1 = thetas(i, :); th2 = thetas(j, :);
            w1 = -key_dtheta(i, :); w2 = key_dtheta(j, :);
            s = td^2 * (5 * (w1 + w2).^2 - 2 * w1 .* w2) + ...
            18 * td * (w1 + w2) .* (th1 - th2) + ...
            18 * (th1 - th2).^2;
            f = f + max(s);
        end
    end
    % 把约束放这里作为罚函数
    angle_limit = [-170, -120, -170, -170, -120, -360;
                   170, 120, 170, 170, 120, 360];
    angle_limit = deg2rad(angle_limit);
    e = 0;
    for i = 1:N
        u = thetas(i, :) - angle_limit(2, :);
        u = max(u, 0);
        l = angle_limit(1, :) - thetas(i, :);
        l = max(l, 0);
        e = e + sum([u, l]);
    end
    f = f + exp(e) - 1;
end

function [c, ceq] = constraint(x)
    % CONSTRAINT 对待优化变量的约束
    % INPUT:
    % x         待优化参数(1x(3+4N))，内容为[p, q1, q2, ..., qN]
    %           其中p(1x3)为关键点笛卡尔坐标均值，
    %           q为各关键点姿态的四元数(1x4)
    % OUTPUT:
    % c         不等式约束(c <= 0)
    % ceq       等式约束(ceq = 0)

    p = x(1, 1:3);
    N = (length(x) - 3) / 4;
    q = zeros(N, 4);
    for i = 1:N
        q(i, :) = x(1, (4 * i):(4 * i + 3));
    end
    c = [];
    ceq = zeros(1, N);
    for i = 1:N
        ceq(1, i) = norm(q(i, :));
    end
end

function best_thetas = get_nearest_theta(gsts)
    % GET_NEAREST_THETA 在每个关键位姿中选取一组使得位姿两两间关节空间距离之和最小
    % INPUT:
    % gsts          各关键点的位姿(Nx4x4)
    % OUTPUT:
    % best_thetas   计算出的两两间距离最小一组theta(Nx6)
    
    % 先算出关键位姿的8个逆解
    N = size(gsts);
    N = N(1);
    key_theta = zeros(N, 8, 6);
    for i = 1:N
        key_theta(i, :, :) = Ikine6s(gsts(i, :, :));
    end

    % % 关节转角限制
    % angle_limit = [-170, -120, -170, -170, -120, -360;
    %                170, 120, 170, 170, 120, 360];
    % angle_limit = deg2rad(angle_limit);
    
    matches = NaN(8, N);
    matches(:, 1) = 1:8;
    score_mat = inf(8, N);
    score_mat(:, 1) = 0;
    for i = 2:N
        for j = 1:8
            if isnan(matches(j, i - 1))
                continue
            end
            last_th = squeeze(key_theta(i - 1, matches(j, i - 1), :));
            best_id = 7;
            best_score = inf;
            for k = 1:8
                now_theta = squeeze(key_theta(i, k, :));
                % 这里不作为强制指标检查了，否则梯度丢失
                % if any(now_theta' > angle_limit(2, :)) || ...
                %    any(now_theta' < angle_limit(1, :))
                %     continue;
                % end
                now_score = sum(abs(now_theta - last_th));
                if now_score < best_score
                    best_id = k;
                    best_score = now_score;
                end
            end
            matches(j, i) = best_id;
            score_mat(j, i) = best_score;
        end
    end
    best_line = 0;
    line_score = inf;
    for i = 1:8
        s = sum(score_mat(i, :));
        if s < line_score
            best_line = i;
            line_score = s;
        end
    end
    if best_line == 0
        % disp('WTF???')
        best_line = 7;
    end
    % if best_line > 8 || best_line < 1
    %     disp('WTF???')
    % end
    best_thetas = zeros(N, 6);
    for i = 1:N
        best_thetas(i, :) = key_theta(i, matches(best_line, i), :);
    end
    
end
