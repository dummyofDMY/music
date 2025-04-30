function [t, xyz, theta] = get_trajectory(music, gst, v, h, saft_gst)
    % GET_TRAJIECTORY 根据乐谱计算整体轨迹
    % INPUT:
    % music     乐谱(2xM，M为关键点数)，记录哪个时间点应该敲击哪个位置，第一行是时间，第二行是敲击点下标
    %           要求空出第1s（此时认为在等待位置）
    % gst       进入敲击轨迹的位姿(Nx4x4，N为敲击点数)
    % v        进入敲击轨迹时的速度(标量)
    % h         敲击轨迹长度(标量)
    % OUTPUT:
    % t         每个坐标对应的时间点(1xlen，len为采样点数)
    % xyz       最终轨迹在笛卡尔空间的坐标(3xlen)
    % theta     最终轨迹在关节空间的坐标(6xlen)
    key_t = music(1, :);
    key_p_id = music(2, :);
    dt = 0.001;
    len = ceil((key_t(end) + 1) / dt);  % 包括了演奏前的等待（1s）和演奏完回归安全位置（1s）
    t = 0:(len);  % 在时间划分上统一认为分界点属于较晚的一段轨迹
    t = t * dt;
    xyz = zeros(3, len);
    theta = zeros(6, len);

    M = size(music);
    M = M(2);
    N = size(gst);
    N = N(1);
    
    % 计算8个位姿中应该选哪一个
    key_thetas = get_nearest_theta(gst);  % (Nx6)
    
    % 先处理第一段等待演奏开始的部分
    saft_xyz = saft_gst(1:3, 4)';  % 行向量
    saft_thetas = get_saft_theta(saft_gst, key_thetas);  % 行向量
    wait_len = ceil(1 / dt);
    xyz(:, 1:wait_len) = repmat(saft_xyz', 1, wait_len);
    theta(:, 1:wait_len) = repmat(saft_thetas', 1, wait_len);

    % 计算在敲击轨迹上要用多少时间
    t_im = 4 * h / v;

    % 计算几个特殊位姿在关节坐标系下的速度
    % 注意这个速度是竖直向下的
    key_dtheta = zeros(N, 6);
    for i = 1:N
        dg = zeros(4, 4);
        dg(3, 4) = -v;
        V_hat = dg / squeeze(gst(i, :, :));
        V = anti_se3_hat(V_hat);
        J = Jacobian(key_thetas(i, :));
        key_dtheta(i, :) = J \ V';
    end

    % 关节转角限制
    angle_limit = [-170, -120, -170, -170, -120, -360;
                   170, 120, 170, 170, 120, 360];
    angle_limit = deg2rad(angle_limit);

    % 逐段计算轨迹
    for i = 1:M
        if i == 1
            start_t = 1;
            x0 = saft_thetas;
            v0 = zeros(1, 6);
        else
            start_t = key_t(i - 1);
            x0 = key_thetas(key_p_id(i - 1), :);
            v0 = -key_dtheta(key_p_id(i - 1), :);
        end
        end_t = key_t(i);
        mid_t = end_t - t_im;
        start_id = ceil(start_t / dt) + 1;
        end_id = ceil(end_t / dt) + 1;
        mid_id = ceil(mid_t / dt) + 1;
        % 计算转移轨迹
        xf = key_thetas(key_p_id(i), :);
        vf = key_dtheta(key_p_id(i), :);
        [a0, a1, a2, a3] = get_transfer_trajectory(x0, xf, v0, vf, mid_t - start_t);
        
        A = [a0', a1', a2', a3'];
        t_line = start_id:mid_id - 1;
        t_line = (t_line - 1) * dt;
        t_line = t_line - start_t;
        T = [ones(1, mid_id - start_id); t_line; t_line.^2; t_line.^3];
        tem = size(A);
        tem = tem(1);
        if tem ~= 6
            disp('WTF?')
        end
        theta(:, start_id:(mid_id - 1)) = A * T;
        for j = start_id:(mid_id - 1)
            now_gst = Fkine(theta(:, j)');
            xyz(:, j) = now_gst(1:3, 4);
        end
        % visualize(xyz(1, start_id:(mid_id - 1)), xyz(2, start_id:(mid_id - 1)), xyz(3, start_id:(mid_id - 1)))
        
        % 计算敲击轨迹
        x0 = squeeze(gst(key_p_id(i), 1:3, 4));
        v0 = [0, 0, -1] * v;
        xt = x0;
        xt(3) = xt(3) - h;
        [a0, a1, a2, ~] = get_impact_trajectory(x0, xt, v0);
        A = [a0', a1', a2'];
        t_line = mid_id:(end_id - 1);
        t_line = (t_line - 1) * dt;
        t_line = t_line - mid_t;
        T = [ones(1, end_id - mid_id); t_line; t_line.^2];
        xyz(:, mid_id:(end_id - 1)) = A * T;

        % 观察撞击曲线算错没
        z_list = squeeze(xyz(3, mid_id:(end_id - 1)));
        z_list_len = length(z_list);
        plot(1:z_list_len, z_list);
        grid on;

        now_gst = squeeze(gst(key_p_id(i), :, :));
        for j = mid_id:(end_id - 1)
            now_gst(1:3, 4) = xyz(:, j);
            solve = Ikine6s(now_gst);
            best_solve_id = 0;
            best_solve_score = inf;
            for k = 1:8
                ths = solve(k, :);
                if any(ths > angle_limit(2, :)) || ...
                   any(ths < angle_limit(1, :))
                    continue
                end
                last_th = theta(:, j - 1)';
                d = norm(ths - last_th);
                if d < best_solve_score
                    best_solve_id = k;
                    best_solve_score = d;
                end
            end
            if best_solve_score > 0.1
                disp('WARNING: thetas leap in trajectory solving! Gap is:');
                disp(best_solve_score);
            end
            theta(:, j) = solve(best_solve_id, :);
        end
    end
    
    % 生成最后一段回归原点的轨迹
    start_t = key_t(end);
    end_t = start_t + 1;
    start_id = ceil(start_t / dt) + 1;
    end_id = ceil(end_t / dt) + 1;
    x0 = key_thetas(key_p_id(end), :);
    v0 = -key_dtheta(key_p_id(end), :);
    xf = saft_thetas;
    vf = zeros(1, 6);
    [a0, a1, a2, a3] = get_transfer_trajectory(x0, xf, v0, vf, end_t - start_t);
    A = [a0', a1', a2', a3'];
    t_line = start_id:end_id - 1;
    t_line = (t_line - 1) * dt;
    t_line = t_line - start_t;
    T = [ones(1, end_id - start_id); t_line; t_line.^2; t_line.^3];
    theta(:, start_id:end_id - 1) = A * T;
    for j = start_id:(end_id - 1)
        now_gst = Fkine(theta(:, j)');
        xyz(:, j) = now_gst(1:3, 4);
    end
end

function visualize(x, y, z)
    % VISUALIZE 可视化三维轨迹
    % INPUT:
    % x         x坐标
    % y         y坐标
    % z         z坐标

    figure();
    plot3(x, y, z);
    grid on;
    hold off;
end

function best_thetas = get_saft_theta(gsts, key_thetas)
    % GET_SAFT_THETA 在安全点处选一个关节坐标使得与关键点间关节空间距离之和最小
    % INPUT:
    % gsts          安全点的位姿(4x4)
    % key_thetas    关键点关节坐标(Nx6)
    % OUTPUT:
    % best_thetas   计算出的两两间距离最小一组theta(Nx6)

    % 关节转角限制
    angle_limit = [-170, -120, -170, -170, -120, -360;
                   170, 120, 170, 170, 120, 360];
    angle_limit = deg2rad(angle_limit);

    saft_thetas = Ikine6s(gsts);  % (8x6)
    best_id = 0;
    best_score = inf;
    for i = 1:8
        now_score = 0;
        now_th = saft_thetas(i, :);
        if any(now_th > angle_limit(2, :)) || ...
           any(now_th < angle_limit(1, :))
            continue;
        end
        for j = 1:7
            d = norm(now_th - key_thetas(j, :));
            now_score = now_score + d;
        end
        if now_score < best_score
            best_id = i;
            best_score = now_score;
        end
    end
    best_thetas = saft_thetas(best_id, :);
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
    best_thetas = zeros(N, 6);
    disp('best line:')
    disp(best_line);
    for i = 1:N
        best_thetas(i, :) = key_theta(i, matches(best_line, i), :);
    end
    
end



% function best_thetas = search_nearest_theta(gsts)
%     % SEARCH_NEAREST_THETA 在每个关键位姿中选取一组使得位姿两两间关节空间距离之和最小
%     % INPUT:
%     % gsts          各关键点的位姿(Nx4x4)
%     % OUTPUT:
%     % best_thetas   计算出的两两间距离最小一组theta(Nx6)
% 
%     % 这个问题可以归类为一个带约束的线性整数规划问题
%     % 先算出关键位姿的8个逆解
%     N = size(gsts);
%     N = N(1);
%     key_theta = zeros(N, 8, 6);
%     angle_limit = [-170, -120, -170, -170, -120, -360;
%                    170, 120, 170, 170, 120, 360];
%     angle_limit = deg2rad(angle_limit);
%     for i = 1:N
%         key_theta(i, :, :) = Ikine6s(gsts(i, :, :));
%         for j = 1:8
%             % 排除无法转到的角度
%             if any(squeeze(key_theta(i, j, :))' > angle_limit(2, :)) || ...
%                any(squeeze(key_theta(i, j, :))' < angle_limit(1, :))
%                 key_theta(i, j, :) = nan;
%             end
%         end
%     end
%     best_score = inf;
%     best_id_vec = zeros(1, N);
%     now_id = 0;
%     id_vec = ones(1, N);  % 存储各个位姿当前迭代取的是哪个theta
%     id_vec(8) = 0;
%     % 初始化进度条
%     while any(id_vec ~= 8)
%         % 更新表示循环进度的id_vec
%         id_vec(8) = id_vec(8) + 1;
%         now_id = now_id + 1;
%         for i = 8:-1:1
%             if id_vec(i) >= 8
%                 id_vec(i) = 1;
%                 id_vec(i - 1) = id_vec(i - 1) + 1;
%             end
%         end
%         now_score = 0;
%         for i = 1:N
%             for j = (i + 1):N
%                 th1 = squeeze(key_theta(i, id_vec(i), :));
%                 th2 = squeeze(key_theta(j, id_vec(j), :));
%                 if any(isnan([th1; th2]))
%                     now_score = inf;
%                 else
%                     now_score = now_score + norm(th1 - th2);
%                 end
%             end
%         end
%         if now_score < best_score
%             best_score = now_score;
%             best_id_vec = id_vec;
%         end
%     end
%     best_thetas = zeros(N, 6);
%     for i = 1:N
%         best_thetas(i, :) = key_theta(i, best_id_vec(i), :);
%     end
%     disp('best key thetas selection done');
%     disp(best_thetas);
% end
% 
% function best_thetas = get_nearest_theta(gsts)
%     % GET_NEAREST_THETA 在每个关键位姿中选取一组使得位姿两两间关节空间距离之和最小
%     % INPUT:
%     % gsts          各关键点的位姿(Nx4x4)
%     % OUTPUT:
%     % best_thetas   计算出的两两间距离最小一组theta(Nx6)
% 
%     % 这个问题可以归类为一个带约束的线性整数规划问题
%     % 先算出关键位姿的8个逆解
%     key_theta = zeros(N, 8, 6);
%     for i = 1:N
%         key_theta(i, :, :) = Ikine6s(gsts(i, :, :));
%     end
% 
%     % 绷不住了，写完这个函数才发现由于约束太强，这个整数规划类似于暴力搜索
%     % 懒得重写了，暴力就暴力，毁灭吧
%     % 仔细想想如果用分支定界法的话还不如暴力搜索（
%     % 之前以为只有N*8个变量，实际有N^8个变量……
% 
%     % 决策变量定义(以只有三个关键位姿为例),一共会有N^8个决策变量，它们是：
%     % x111, x112, x113, ..., x118, x121, x122, ..., x888
%     % 其中xijk表示第一、二、三个位姿分别选第i、j、k个解
%     % 定义目标函数向量f
%     f = zeros(N * 8, 1);
%     id_vec = ones(1, N);  % 存储各个位姿当前迭代取的是哪个theta
%     id_vec(8) = 0;
%     f_id = 0;  % 表示当前迭代在算f的哪一个元素
%     while any(id_vec ~= 8)
%         % 更新表示循环进度的id_vec
%         id_vec(8) = id_vec(8) + 1;
%         for i = 8:-1:1
%             if id_vec(i) >= 8
%                 id_vec(i) = 1;
%                 id_vec(i - 1) = id_vec(i - 1) + 1;
%             end
%         end
%         f_id = f_id + 1;
%         % for i = 1:N
%         %     f_id = f_id + (id_vec(i) - 1) * 8^(8 - i) + 1;
%         % end
%         % 计算各关键点两两间在关节空间中的距离
%         for i = 1:N
%             for j = (i + 1):N
%                 f(f_id) = f(f_id) + norm(key_theta(i, id_vec(i), :) - key_theta(j, id_vec(j), :));
%             end
%         end
%     end
% 
%     % 规定哪些是整数变量
%     intcon = 1:(N * 8);  % 全都是
%     % 规定不等式约束矩阵和向量
%     A = [];
%     b = [];
%     % 规定等式约束
%     Aeq = ones(1, N * 8);
%     beq = 1;
%     % 规定上下界
%     lb = zeros(N * 8, 1);
%     ub = ones(N * 8, 1);
% 
%     % 求解
%     disp('start sloving best theta')
%     [x, ~] = intlinprog(f, intcon, A, b, Aeq, beq, lb, ub);
%     disp('solving done');
%     best_thetas = zeros(N, 6);
%     for i = 1:N
%         best_thetas(i, :) = key_theta(i, x(i), :);
%     end
% end