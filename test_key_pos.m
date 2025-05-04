a = 1/2^0.5;
R = [-a, 0, a;
     0, 1, 0;
     -a, 0, -a];
q = rotation_matrix_to_quaternion(R);

gsts = get_key_pos(20);
% 计算8个位姿中应该选哪一个
key_thetas = get_nearest_theta(gsts);  % (Nx6)
disp(squeeze(gsts(4, :, :)));
disp(gsts)

gsts = permute(gsts, [2 3 1]);  % 变成 4×4×7

figure();
hold on;
colors = lines(7);
scale = 20;

for i = 1:7
    T = gsts(:, :, i);  % 每个 4x4 齐次矩阵
    pt = T(1:3, 4);     % 取出平移部分
    R = T(1:3, 1:3);    % 取出旋转部分
    
    scatter3(pt(1), pt(2), pt(3), 50, 'filled', 'MarkerFaceColor', colors(i, :));
    
    % 绘制Z轴方向
    vec = R * [0; 0; 1];
    quiver3(pt(1), pt(2), pt(3), vec(1)*scale, vec(2)*scale, vec(3)*scale, ...
        'Color', colors(i, :), 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
end

xlabel('X'); ylabel('Y'); zlabel('Z');
title('End-effector Frames');
grid on;
% axis equal;
view(3);
xlim padded; ylim padded; zlim padded;
hold off;


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
