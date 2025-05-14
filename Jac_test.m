% 检验贾克比矩阵计算的正确性
clear;
M = 10000;
thetas = zeros(6, M);
v = [0.1; 0.1; 0; 0; 0; 0];
dt = 0.01;

for i = 2:M
    thetas(:, i) = thetas(:, i - 1) + v * dt;
end

xyz = zeros(3, M);
for i = 1:M
    gst = Fkine(thetas(:, i)');
    xyz(:, i) = gst(1:3, 4);
end

velo_gt = zeros(3, M - 1);
for i = 2:M
    velo_gt(:, i - 1) = (xyz(:, i) - xyz(:, i - 1)) / dt;
end

velo_est = zeros(3, M - 1);
for i = 2:M
    V = Jacobian(thetas(:, i)') * v;
    V_hat = se3_hat(V');
    p_ = xyz(:, i);
    V_vec = V_hat * [p_; 1];
    velo_est(:, i - 1) = V_vec(1:3);
end

t = 0:(M - 1);
t = t * dt;

figure();
for i = 1:6
    plot(t, rad2deg(thetas(i, :)));
    hold on;
end
legend('1', '2', '3', '4', '5', '6');
title('theta')
grid on;
hold off

figure();
subplot(3, 2, 2);
plot(t(2:end), velo_est(1, :))
hold on;
plot(t(2:end), velo_gt(1, :))
hold off;
legend('est', 'gt');
title('v_x');

subplot(3, 2, 4);
plot(t(2:end), velo_est(2, :))
hold on;
plot(t(2:end), velo_gt(2, :))
hold off;
legend('est', 'gt');
title('v_y');

subplot(3, 2, 6);
plot(t(2:end), velo_est(3, :))
hold on;
plot(t(2:end), velo_gt(3, :))
hold off;
legend('est', 'gt');
title('v_z');

subplot(3, 2, 1);
plot(t, xyz(1, :))
title('x');

subplot(3, 2, 3);
plot(t, xyz(2, :))
title('y');

subplot(3, 2, 5);
plot(t, xyz(3, :))
title('z');
