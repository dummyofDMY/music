function theta = leagalize_theta(theta)
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