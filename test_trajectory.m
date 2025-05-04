music = [1.1, 1.2;
         1, 7];

gst = zeros(7, 4, 4);
th0 = [-16.551, 30.381, 107.605, 0.379, 38.273, 151.667];
th0 = deg2rad(th0);
the = [17.783, 29.758, 108.697, -3.087, 38.362, 188.659];
the = deg2rad(the);
gst(1, :, :) = Fkine(th0);
gst(7, :, :) = Fkine(the);
xyz0 = gst(1, 1:3, 4);
xyze = gst(7, 1:3, 4);
R = squeeze(gst(1, 1:3, 1:3));
for i = 1:5
    xyz = (xyz0 * (6 - i) + xyze * i) / 6;
    gst((i + 1), :, :) = [R, xyz';
                          zeros(1, 3), 1];
end
v0 = 1000;
h = 10;
R = squeeze(gst(4, 1:3, 1:3));
q = rotation_matrix_to_quaternion(R);
disp(q);
saft_gst = squeeze(gst(4, :, :));
saft_gst(2, 4) = saft_gst(2, 4) - 50;
saft_gst(3, 4) = saft_gst(3, 4) + 50;
[t, xyz, theta] = get_trajectory(music, gst, v0, h, saft_gst);
plot3(xyz(1, :), xyz(2, :), xyz(3, :))
grid on;

theta = rad2deg(theta);
writematrix(theta', 'pt_list.txt', 'Delimiter', 'space');