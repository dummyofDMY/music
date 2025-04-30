function R = quaternion_to_rotation_matrix(q)
    % QUATERNION_TO_ROTATION_MATRIX 将四元数转换为旋转矩阵
    % INPUT:
    % q             四元数[w, x, y, z](1xN)
    % OUTPUT:
    % R             算得的旋转矩阵(3x3)
    w = q(1); x = q(2); y = q(3); z = q(4);
    R = [1 - 2*y^2 - 2*z^2, 2*x*y - 2*w*z, 2*x*z + 2*w*y;
         2*x*y + 2*w*z, 1 - 2*x^2 - 2*z^2, 2*y*z - 2*w*x;
         2*x*z - 2*w*y, 2*y*z + 2*w*x, 1 - 2*x^2 - 2*y^2];
end

