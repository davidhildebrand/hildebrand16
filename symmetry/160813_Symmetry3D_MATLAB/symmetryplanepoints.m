function [sp,y] = symmetryplanepoints(ms,dotB)

    y = dotB(ms(:,2),ms(:,3)); % x component of projection of symmetry points into plane

    mn2 = min(ms(:,2)); mx2 = max(ms(:,2));
    mn3 = min(ms(:,3)); mx3 = max(ms(:,3));
    sp1 = [dotB(mn2,mn3) mn2 mn3];
    sp2 = [dotB(mx2,mn3) mx2 mn3];
    sp3 = [dotB(mx2,mx3) mx2 mx3];
    sp4 = [dotB(mn2,mx3) mn2 mx3];
    sp = [sp1; sp2; sp3; sp4]; % boundaries points to draw plane

end