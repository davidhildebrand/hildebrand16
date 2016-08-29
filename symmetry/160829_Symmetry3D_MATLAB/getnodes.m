function [P,TP] = getnodes(D,iskelsindex)
    rows = find(D(:,1) == iskelsindex);
    DS = D(rows,:); % skel id, node id, parent id, x, y, z
    P = []; % points
    TP = []; % tangents
    for j = 1:length(rows) % loop through nodes
        rowp = find(DS(:,3) == DS(j,2)); % row of parent
        if size(rowp,1) == 1 % only one parent
            p = DS(j,4:6); % (x,y,z)
            pp = DS(rowp,4:6);
            P = [P; pp]; % parent
            TP = [TP; p-pp]; % direction from parent to child
        end
    end
end