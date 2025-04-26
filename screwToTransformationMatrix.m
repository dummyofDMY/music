function T = screwToTransformationMatrix(S, theta)
    % SCREWTOTRANSFORMATIONMATRIX 将旋量转换为齐次变换矩阵
    % Input:
    %   S: 6x1 旋量 [v; omega] (v: 线速度, omega: 旋转轴)
    %   theta: 旋转角度（弧度）
    % Output:
    %   T: 4x4 齐次变换矩阵
    
    v = S(1:3);       % 提取线速度分量
    omega = S(4:6);   % 提取旋转轴
    
    % 计算旋转矩阵 R（使用 Rodrigues 公式）
    omega_skew = [0, -omega(3), omega(2);
                  omega(3), 0, -omega(1);
                  -omega(2), omega(1), 0];  % 叉积矩阵 [omega]_×
    
    R = eye(3) + sin(theta) * omega_skew + (1 - cos(theta)) * (omega_skew^2);
    
    % 计算平移向量 p
    I = eye(3);
    omega_cross_v = cross(omega', v');
    p = (I - R) * omega_cross_v + omega' * omega * v' * theta;
    
    % 构建齐次变换矩阵
    T = [R, p;
         0, 0, 0, 1];
end