function x = logQ(scalar, q)

x = (scalar == 0) .* 2 .* atan2(vecnorm(q(:,2:4), 2, 2), q(:,1)) .* normRow(q(:,2:4)) ...
    + (scalar == 4) .* 2 .* atan2(vecnorm(q(:,1:3), 2, 2), q(:,4)) .* normRow(q(:,1:3));  

end

