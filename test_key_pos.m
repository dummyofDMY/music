gsts = get_key_pos(20);
% 计算8个位姿中应该选哪一个
key_thetas = get_nearest_theta(gsts);  % (Nx6)
disp(key_thetas);
figure();
for i = 1:N
    pt = [gsts(i, 1, 4), gsts(i, 2, 4), gsts(i, 3, 4)];
    scatter3(pt(1, 1), pt(1, 2), pt(1, 3));
    hold on;
    vec = squeeze(gsts(i, 1:3, 1:3)) * [0; 0; 1];
    quiver3(pt(1, 1), pt(1, 2), pt(1, 3), vec(1, 1), vec(2, 1), vec(3, 1));
    hold on;
end
hold off;
