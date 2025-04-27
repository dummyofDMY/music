% AUTHOR :dummy
%
% ABSTRACT： 将旋量从向量转化为矩阵
% 
% INPUT： screw    旋量向量(1x6)
%
% OUTPUT: res      求得的结果
% 

function res = se3_hat(screw)
    w = screw(4:6);
    R = [0, -w(3), w(2);
         w(3), 0, -w(1);
         -w(2), w(1), 0];
    res = [R, screw(1:3)';
           zeros(1, 3), 0];
end

