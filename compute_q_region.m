function [ Qmasked, X1, X2] = compute_q_region(mu, s, Rs, eps)
    % this file helps to get contour for sub-level set of \hatq
    % the finer the region more accurate is the contour but the higehr the time, if the uncertainity level is not change with time we can
    % precompute this
    % x1 = linspace(-50, 340, 500);
    % x2 = linspace(-70, 430, 500);
    x1 = linspace(-15, 18, 1300);
    x2 = linspace(-15, 18, 1300);
    [X1, X2] = meshgrid(x1, x2);
    Q = zeros(size(X1));
    % Compute q(x) at each grid point
    for i = 1:length(x1)
        for j = 1:length(x2)
            tube = [X1(j,i); X2(j,i)];
            lambda = (norm(mu - tube) / s)^2;
            p_i = ncx2cdf((Rs/s)^2, 2, lambda);
            q = 1 - p_i;
            Q(j,i) = q;
        end
    end
 if Rs<0.34
    % Mask region
    Qmasked = Q;
    Qmasked(Q >min(1.4*(1-eps),1)) = NaN;
 else
     Qmasked = Q;
    Qmasked(Q >min(1.1*(1-eps),1)) = NaN;
 end
     % Qmasked(Q <0.4) = NaN;
end