function screw_vec = anti_se3_hat(screw_mat)
    % SE2_HAT 将旋量从矩阵转化为向量
    % INPUT:
    % screw    旋量矩阵(4x4)
    % OUTPUT:
    % res      求得的结果(1x6)
    v = screw_mat(1:3, 4)';
    w_hat = screw_mat(1:3, 1:3);
    w = [w_hat(3, 2), w_hat(1, 3), w_hat(2, 1)];
    screw_vec = [v, w];
end

