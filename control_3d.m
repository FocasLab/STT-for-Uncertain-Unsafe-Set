function u = control_3d(tube, rad, robot,velo,rho_L,rho_U,rho_d,rho_d_dot)
    %% control input as in Theorem 4.1 or Algorithm 1
    k_2=15;
    e = norm(robot - tube)/rad;
    k_1 =4;
    v_d = real(-k_1 * (robot - tube) * log((1+e)/(1-e)));
    for i=1:length(robot)   
    x=velo(i)-v_d(i);
    e1=(x-0.5*(rho_U(i)+rho_L(i)))/(0.5*(rho_d(i)));
    eps=log((1+e1)/(1-e1));
    zeta=4/(rho_d(i)*(1-e1^2));
    u(i)=-(k_2*zeta*eps-0.5*rho_d_dot(i)*e1);
    end
    u = real(u);
end

