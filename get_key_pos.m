function [res] = get_key_pos(l)
    % GET_KEY_POS 优化算法得出最优关键点的位姿
    % INPUT:
    % l             木琴两个音条中轴线的距离(标量)
    % OUTPUT:
    % gsts          关键点位姿(Nx4x4)

    % 初始猜测值 (需为单位四元数)
    typical_q = [0.382683432365090	0	0.923879532511287	0];
    typical_p = [5.122186101582196e+02,3.133945035755696,4.610435346606103e+02];
    x0 = typical_p; % 示例: 表示绕y轴旋转90度
    gsts = zeros(7, 4, 4);
    R = quaternion_to_rotation_matrix(typical_q);
    for i = 1:7
        gsts(i, 1:3, 1:3) = R;
        gsts(i, 1:3, 4) = typical_p + [0, (i - 4) * l, 0];
        gsts(i, 4, 4) = 1;
    end
    
    
    % 设置非线性约束
    nonlcon = @constraint_xyz;

    % 创建优化选项
    options = optimoptions('fmincon', 'Display', 'iter', ...
         'MaxIterations', 400, 'MaxFunctionEvaluations', 1500);

    fun = @(x) obj_fun_xyz(x, l);
    % 调用fmincon求解
    [x, fval] = fmincon(fun, x0, [], [], [], [], [], [], nonlcon, options);
    res = zeros(7, 4, 4);
    for i = 1:7
        res(i, 1:3, 1:3) = R;
        res(i, 1:3, 4) = x + [0, (i - 4) * l, 0];
        res(i, 4, 4) = 1;
    end
    % 角度优化
    options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'interior-point', ...
          'AlwaysHonorConstraints', 'bounds');
    init_thetas = get_nearest_theta(res);
    x0 = typical_q;
    nonlcon = @constraint;
    for i = 1:7
        thetas = [init_thetas(1:(i - 1), :); init_thetas((i + 1):7, :)];
        xyz = squeeze(res(i, 1:3, 4));
        fun = @(x) obj_fun(x, thetas, xyz);
        
        % 调用fmincon求解
        [x, fval] = fmincon(fun, x0, [], [], [], [], [], [], nonlcon, options);

        % % 输出结果
        % disp(['最优解(单位四元数): ', num2str(x)]);
        % disp(['最小化的目标函数值: ', num2str(fval)]);

        R = quaternion_to_rotation_matrix(x);
        gst = [R, xyz';
               zeros(1, 3), 1];
        my_theta = Ikine6s(gst);
        scores = zeros(1, 6);

        angle_limit = [-170, -120, -170, -170, -120, -360;
                       170, 120, 170, 170, 120, 360];
        angle_limit = deg2rad(angle_limit);
        for j = 1:6
            nqq = 5 * my_theta(j, :) * my_theta(j, :)';
            db_q0qj = 2 * my_theta(j, :) * sum(thetas, 1)';
            scores(1, j) = nqq - db_q0qj;
            % % 罚函数，防止超出角度限制范围
            % punishment1 = squeeze(my_theta(j, :)) - squeeze(angle_limit(2, :));
            % punishment1(punishment1 < 0) = 0;
            % punishment1 = exp(punishment1) - 1;
            % punishment2 = squeeze(angle_limit(2, :)) - squeeze(my_theta(j, :));
            % punishment2(punishment2 < 0) = 0;
            % punishment2 = exp(punishment2) - 1;
            % scores(1, j) = scores(1, j) + sum(punishment1 + punishment2);
            % % 防止接近奇异
            % J = Jacobian(squeeze(my_theta(j, :)));
            % punishment3 = exp(det(J));
            % scores(1, j) = scores(1, j) + punishment3;
        end
        [~, idx] = min(scores);
        init_thetas(i, :) = my_theta(i, :);
        res(i, :, :) = gst;
    end

    save('best_pos', "res");
end

function f = obj_fun_xyz(x, l)
    % OBJ_FUN_XYZ 评价函数，衡量一组最优关键点位姿的优秀程度
    % INPUT:
    % x         待优化参数(1x(3+4N))，内容为[w, x, y, z]，是待优化旋转的四元数
    % thetas    其余的不优化的坐标的关节空间坐标((N-1)x6)
    % xyz       该关键点的xyz坐标(1x3)
    % OUTPUT:
    % f         输出评价指标（越小越好）
    typical_q = [0.382683432365090	0	0.923879532511287	0];
    R = quaternion_to_rotation_matrix(typical_q);
    gsts = zeros(7, 4, 4);
    for i = 1:7
        gsts(i, 1:3, 1:3) = R;
        gsts(i, 1:3, 4) = x + [0, (i - 4) * l, 0];
        gsts(i, 4, 4) = 1;
    end
    key_theta = get_nearest_theta(gsts);
    scores = 0;

    angle_limit = [-170, -120, -170, -170, -120, -360;
                   170, 120, 170, 170, 120, 360];
    angle_limit = deg2rad(angle_limit);
    for i = 1:7
        for j = (i+1):7
            scores = scores + norm(key_theta(i, :) - key_theta(j, :));
        end
        % % 罚函数，防止超出角度限制范围
        % punishment1 = squeeze(key_theta(i, :)) - squeeze(angle_limit(2, :));
        % punishment1(punishment1 < 0) = 0;
        % punishment1 = exp(punishment1) - 1;
        % punishment2 = squeeze(angle_limit(2, :)) - squeeze(key_theta(i, :));
        % punishment2(punishment2 < 0) = 0;
        % punishment2 = exp(punishment2);
        % scores = scores + sum(punishment1 + punishment2);
        % % 防止接近奇异
        % J = Jacobian(squeeze(key_theta(i, :)));
        % sigma = svd(J);
        % punishment3 = exp(abs(sigma(end)));
        % scores = scores + punishment3;
    end
    f = scores;
    is_valid = isfinite(f) & isreal(f);
    if ~is_valid
        disp('not valid!');
        disp(x);
        disp(f);
    end
end

function f = obj_fun(x, thetas, xyz)
    % OBJ_FUN 评价函数，衡量一组最优关键点位姿的优秀程度
    % INPUT:
    % x         待优化参数(1x(3+4N))，内容为[w, x, y, z]，是待优化旋转的四元数
    % thetas    其余的不优化的坐标的关节空间坐标((N-1)x6)
    % xyz       该关键点的xyz坐标(1x3)
    % OUTPUT:
    % f         输出评价指标（越小越好）
    R = quaternion_to_rotation_matrix(x);
    gst = [R, xyz';
           zeros(1, 3), 1];
    my_theta = Ikine6s(squeeze(gst));
    if any(isnan(my_theta))
        disp('nan!!!');
        disp(gst);
        disp(my_theta);
    end
    scores = zeros(1, 6);

    angle_limit = [-170, -120, -170, -170, -120, -360;
                   170, 120, 170, 170, 120, 360];
    angle_limit = deg2rad(angle_limit);
    for i = 1:6
        nqq = 5 * my_theta(i, :) * my_theta(i, :)';
        db_q0qi = 2 * my_theta(i, :) * sum(thetas, 1)';
        scores(1, i) = nqq - db_q0qi;
        % % 罚函数，防止超出角度限制范围
        % punishment1 = squeeze(my_theta(i, :)) - squeeze(angle_limit(2, :));
        % punishment1(punishment1 < 0) = 0;
        % punishment1 = exp(punishment1);
        % punishment2 = squeeze(angle_limit(2, :)) - squeeze(my_theta(i, :));
        % punishment2(punishment2 < 0) = 0;
        % punishment2 = exp(punishment2);
        % scores(1, i) = scores(1, i) + sum(punishment1 + punishment2);
        % % 防止接近奇异
        % J = Jacobian(squeeze(my_theta(i, :)));
        % sigma = svd(J);
        % punishment3 = exp(abs(sigma(end)));
        % scores(1, i) = scores(1, i) + punishment3;
    end
    f = min(scores);
    is_valid = isfinite(f) & isreal(f);
    if ~is_valid
        disp('not valid!');
        disp(x);
        disp(thetas);
        disp(xyz);
        disp(f);
    end
end

function [c, ceq] = constraint(x)
    % CONSTRAINT 对待优化变量的约束
    % INPUT:
    % x         待优化参数(1x(3+4N))，内容为[w, x, y, z]
    %           其中p(1x3)为关键点笛卡尔坐标均值，
    %           q为各关键点姿态的四元数(1x4)
    % OUTPUT:
    % c         不等式约束(c <= 0)
    % ceq       等式约束(ceq = 0)

    c = 1 - 2*x(2)^2 - 2*x(3)^2;
    ceq = norm(x) - 1;
end

function [c, ceq] = constraint_xyz(x)
    % CONSTRAINT 对待优化变量的约束
    % INPUT:
    % x         待优化参数(1x(3+4N))，内容为[w, x, y, z]
    %           其中p(1x3)为关键点笛卡尔坐标均值，
    %           q为各关键点姿态的四元数(1x4)
    % OUTPUT:
    % c         不等式约束(c <= 0)
    % ceq       等式约束(ceq = 0)

    % c = 1 - 2*x(2)^2 - 2*x(3)^2;
    % ceq = norm(x) - 1;
    r = norm(x);

    c = [500 - r, r - 800, -200 - x(1, 3)];
    ceq = [];
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
        key_theta(i, :, :) = Ikine6s(squeeze(gsts(i, :, :)));
    end

    % 关节转角限制
    angle_limit = [-170, -120, -170, -170, -120, -360;
                   170, 120, 170, 170, 120, 360];
    angle_limit = deg2rad(angle_limit);
    
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
            best_id = nan;
            best_score = inf;
            for k = 1:8
                now_theta = squeeze(key_theta(i, k, :));
                if any(now_theta' > angle_limit(2, :)) || ...
                   any(now_theta' < angle_limit(1, :))
                    continue;
                end
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
    disp('theta match result:')
    disp(matches);
    best_line = 0;
    line_score = inf;
    for i = 1:8
        s = sum(score_mat(i, :));
        if s < line_score
            best_line = i;
            line_score = s;
        end
    end
    best_thetas = nan(N, 6);
    if best_line == 0
        return
    end
    disp('best line:')
    disp(best_line);
    for i = 1:N
        best_thetas(i, :) = key_theta(i, matches(best_line, i), :);
    end
    
end
