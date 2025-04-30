function res = Jacobian(thetas)
    % JACOBIAN 计算当前姿态的雅可比矩阵
    % INPUT:
    % thetas   各关节当前转角(1xN)
    % OUTPUT:
    % res      求得的结果(6xN)

    % 初始化旋量、初始位姿
    L1 = 491; L2 = 450; L3 = 450; L4 = 84;
    q1 = [0, 0, 0]; q2 = [0, 0, L1];
    q3 = [0, 0, L1+L2]; q4 = [0, 0, L1+L2];
    q5 = [0, 0, L1+L2+L3]; q6 = [0, 0, L1+L2+L3];
    w1 = [0, 0, 1]; w2 = [0, 1, 0];
    w3 = [0, 1, 0]; w4 = [0, 0, 1];
    w5 = [0, 1, 0]; w6 = [0, 0, 1];
    
    g0 = [-1 0 0 0; 
        0 -1 0 0; 
        0 0 1 L1+L2+L3+L4; 
        0 0 0 1];
    
    Xi1 = zeros(1, 6); Xi1(1:3) = cross(q1, w1); Xi1(4:6) = w1;
    Xi2 = zeros(1, 6); Xi2(1:3) = cross(q2, w2); Xi2(4:6) = w2;
    Xi3 = zeros(1, 6); Xi3(1:3) = cross(q3, w3); Xi3(4:6) = w3;
    Xi4 = zeros(1, 6); Xi4(1:3) = cross(q4, w4); Xi4(4:6) = w4;
    Xi5 = zeros(1, 6); Xi5(1:3) = cross(q5, w5); Xi5(4:6) = w5;
    Xi6 = zeros(1, 6); Xi6(1:3) = cross(q6, w6); Xi6(4:6) = w6;
    Xi = [Xi1; Xi2; Xi3; Xi4; Xi5; Xi6];

    % 计算雅可比
    N = length(thetas);
    res = zeros(6, N);
    for i = 1:N
        e_i_1 = eye(4);
        for j = 1:(i-1)
            e_i_1  = e_i_1 * screwToTransformationMatrix(Xi(j, :), thetas(j));
        end
        Ad_mat = Ad(e_i_1);
        res(:, i) = Ad_mat * Xi(i, :)';
        
    end
end

