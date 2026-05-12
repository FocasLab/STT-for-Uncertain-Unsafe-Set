function u = control(tube,rad,robot)
     %% control input as in Theorem 4.1 or Algorithm 1
    e = norm(robot - tube)/rad;
    k = 10;
    u = -k * (robot - tube) * log((1+e)/(1-e));
    u = real(u);
end

