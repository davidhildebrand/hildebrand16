function [rmin, rmax] = getrange(D,iskels)
    rmins = zeros(length(iskels),3);
    rmaxs = rmins;
    parfor i = 1:length(iskels)
        rows = find(D(:,1) == iskels(i));

        DS = D(rows,:); % skel id, node id, parent id, x, y, z

        rmin = [Inf Inf Inf];
        rmax = -rmin;
        for j = 1:length(rows) % loop through nodes
            rowp = find(DS(:,3) == DS(j,2)); % row of parent
            if size(rowp,1) == 1
                p = DS(rowp,4:6); % (x,y,z)
                rmin = min(rmin,p);
                rmax = max(rmax,p);
            end
        end
        rmins(i,:) = rmin;
        rmaxs(i,:) = rmax;
    end
    rmin = min(rmins);
    rmax = max(rmaxs);
end