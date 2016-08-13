function [pperp, displ, avgprjdst,sp] = symmetrytest_pairs(ms,B)
    dotB = @(ms2,ms3) B(1)+B(2)*ms2+B(3)*ms3;
    [sp,y] = symmetryplanepoints(ms,dotB);
    
    % vector perpendicular to plane
    p12 = sp(2,:)-sp(1,:); p12 = p12/norm(p12);
    p14 = sp(4,:)-sp(1,:); p14 = p14/norm(p14);
    pperp = cross(p12,p14);
    
    pointonplane = sp(1,:);
    displ = dot(pointonplane,pperp); % distance(plane,origin)

    avgprjdst = mean(abs(ms(:,1)-y));
end