% AUTHOR :dummy
%
% ABSTRACT： 这是计算机器人正解函数，通用函数，需调用Transformation函数
% 
% INPUT： theta   机器人关节位移，1xN向量，单位m和rad
%
% OUTPUT: g_st    机器人末端位姿位姿， 4X4矩阵
% 
function g_st = Fkine(theta)
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
    g_st = g0;
    for i = flip(1:size(Xi,1))
        g_st = screwToTransformationMatrix(Xi(i,:),theta(i)) * g_st; % 齐次变换矩阵
    end
end


    
        
    
    
            
            
    





    