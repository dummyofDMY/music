% music = [1.1, 1.2;
%          1, 7];

% % Do Re Mi检验
% t = [1.500, 2.000, 2.500, 3.000, 3.500, 4.000, 4.500, 5.000, 5.500, 6.000, 6.500, 7.000, 7.500, 8.000, 8.500
% ];
% scale = 1:15;
% % scale = [scale, 1];
% music = [t; scale];

% 小星星简谱
t = [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5,...
    9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.5, 14.0, 14.5, 15.0, 15.5,...
    16.0, 16.5, 17.5, 18.0, 18.5, 19.0, 19.5, 20.0, 20.5, 21.5, 22.0, 22.5,...
    23.0, 23.5, 24.0, 24.5];
scale = [4, 4, 8, 8, 9, 9, 8, 7, 7, 6, 6, 5, 5, 4, 8, 8, 7, 7, 6, 6, 5, 8,...
    8, 7, 7, 6, 6, 5, 4, 4, 8, 8, 9, 9, 8, 7, 7, 6, 6, 5, 5, 4];
music = [t; scale];

% gst = zeros(15, 4, 4);
% th0 = [-16.551, 30.381, 107.605, 0.379, 38.273, 151.667];
% th0 = deg2rad(th0);
% the = [17.783, 29.758, 108.697, -3.087, 38.362, 188.659];
% the = deg2rad(the);
% gst(1, :, :) = Fkine(th0);
% gst(15, :, :) = Fkine(the);
% xyz0 = gst(1, 1:3, 4);
% xyze = gst(15, 1:3, 4);
% R = squeeze(gst(1, 1:3, 1:3));
% for i = 1:13
%     xyz = (xyz0 * (14 - i) + xyze * i) / 14;
%     gst((i + 1), :, :) = [R, xyz';
%                           zeros(1, 3), 1];
% end
% thetas = get_nearest_theta(gst);
% writematrix(thetas, 'key_theta_test.txt', 'Delimiter', 'space');

v0 = 500;
h = 10;
% key_pt_Js = load("key_theta_test.txt");
key_pt_Js = load('key_thetas.mat');
key_pt_Js = key_pt_Js.key_thetas;
key_pt_Js = key_pt_Js(end:-1:1, :);
[pose_num, ~] = size(key_pt_Js);
gst = zeros(pose_num, 4, 4);
for i = 1:pose_num
    gst(i, :, :) = Fkine(squeeze(key_pt_Js(i, :)));
end
saft_gst = gst(round((1 + pose_num) / 2), :, :);
saft_gst = squeeze(saft_gst);

[t, xyz, theta] = get_trajectory(music, gst, v0, h, saft_gst);
plot3(xyz(1, :), xyz(2, :), xyz(3, :))
grid on;

theta = rad2deg(theta);
theta = round(theta, 4);
writematrix(theta', 'pt_list.txt', 'Delimiter', 'space');

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