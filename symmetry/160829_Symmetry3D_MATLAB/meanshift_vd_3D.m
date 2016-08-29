function m = meanshift_vd_3D(R,m,K,N)

    den = 0;
    for i = 1:size(R,1)
        den = den+N(R(i,:),m)*K(R(i,:),m);
    end
    num = zeros(1,4);
    for i = 1:size(R,1)
        num = num+N(R(i,:),m)*K(R(i,:),m)*R(i,:);
    end
    m = num/den;
    m(1:3) = m(1:3)/norm(m(1:3)); % project to sphere

end