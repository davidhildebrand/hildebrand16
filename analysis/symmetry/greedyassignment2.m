function [assignment,cost,costs] = greedyassignment2(C)
    % greedy on the matching cost

    for i = 1:size(C,1)
        for j = 1:i
            C(i,j) = Inf;
        end
    end
    
    cost = 0;
    costs = inf(1,ceil(size(C,1)/2));
    assignment = zeros(1,size(C,1));
    count = 0;
    while min(min(C)) < Inf
        [m,rowOfMinAtEachCol] = min(C);
        [~,colOfMin] = min(m);
        rowOfMin = rowOfMinAtEachCol(colOfMin);
        
        cost = cost+2*C(rowOfMin,colOfMin);
        
        count = count+1;
        costs(count) = C(rowOfMin,colOfMin);
        
        assignment(rowOfMin) = colOfMin;
        assignment(colOfMin) = rowOfMin;
        
        C(rowOfMin,:) = Inf;
        C(colOfMin,:) = Inf;
        C(:,rowOfMin) = Inf;
        C(:,colOfMin) = Inf;
    end
end

