function Ad_mat = Ad(matrix)
    % AD 计算矩阵的伴随矩阵
    % INPUT:
    % matrix   待求伴随的变换矩阵
    % OUTPUT:
    % Ad_mat   求得的结果

    R = matrix(1:3, 1:3);
    p = matrix(1:3, 4);
    p_hat = [0, -p(3), p(2);
             p(3), 0, -p(1);
             -p(2), p(1), 0];
    Ad_mat = [R, p_hat * R;
              zeros(3), R];
end