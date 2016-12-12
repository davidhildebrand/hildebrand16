function [x,y,z] = symplanecoord(P,sp,V)
    mp = mean(sp);
    MP = repmat(mp,[size(P,1) 1]);

    vy = sp(2,:)-sp(1,:);
    vy = vy/norm(vy);
    vx = V;
    vz = cross(vx,vy);

    VV = repmat(vx,[size(P,1) 1]);
    x = sum((P-MP).*VV,2);
    VV = repmat(vy,[size(P,1) 1]);
    y = sum((P-MP).*VV,2);
    VV = repmat(vz,[size(P,1) 1]);
    z = sum((P-MP).*VV,2);
end