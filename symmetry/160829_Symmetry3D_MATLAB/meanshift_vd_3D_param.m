function [K,N] = meanshift_vd_3D_param(R)
    V = R(:,1:3);
    D = R(:,4);

    s123 = max(std(V));
    s4 = std(D);
    K = @(x,x0) exp(-( s123*sum((x(1:3)-x0(1:3)).^2) + s4*(x(4)-x0(4))^2)); % kernel
    N = @(x,x0) ( norm(x(1:3)-x0(1:3)) < s123 && abs(x(4)-x0(4)) < s4 ); % neighborhood
end