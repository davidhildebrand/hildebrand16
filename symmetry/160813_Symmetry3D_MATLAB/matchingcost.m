function c = matchingcost(P,Q,sp,V)
    SP = reflectedpoints(P,sp,V);
    s = SP(1:50:size(SP,1),:);
    t = Q(1:50:size(Q,1),:);
    [c, is, it] = dtw(s,t);
    if isempty(is)
        c = Inf;
    end
end