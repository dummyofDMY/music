function theta = leagalize_theta(theta)
    % LEAGALIZE_THETA 将角度限制在[-pi, pi]之间
    % INPUT:
    % theta     弧度制角度
    % OUTPUT:
    % theta     规范化后的弧度制角度

    if any(isnan(theta))
        return
    else
        % 将 theta 归一化到 [-pi, pi] 范围内
        for i = 1:length(theta)
            while theta(i) > pi
                theta(i) = theta(i) - 2 * pi;
            end
            while theta(i) < -pi
                theta(i) = theta(i) + 2 * pi;
            end
        end
    end
end