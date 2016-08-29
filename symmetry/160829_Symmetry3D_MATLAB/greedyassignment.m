function [assignment,cost] = greedyassignment(C)
    assignment = zeros(1,size(C,1));
    cost = 0;
    for i = 1:size(C,1)-1
        c = Inf;
        ic = 0;
        for j = i+1:size(C,2)
            if C(i,j) < c
                c = C(i,j);
                ic = j;
            end
        end
        if c < Inf
            assignment(i) = ic;
            assignment(ic) = i;
            C(ic,:) = Inf;
            cost = cost+2*c; % (i,ic) and (ic, i)
        end
%         disp(C)
%         pause
    end
end

