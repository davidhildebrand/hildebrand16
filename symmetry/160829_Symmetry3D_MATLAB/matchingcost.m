function c = matchingcost(P,Q,sp,V,sampfreq)
    SP = reflectedpoints(P,sp,V);
    s = SP(1:sampfreq:size(SP,1),:);
    t = Q(1:sampfreq:size(Q,1),:);
    [c, is, it] = dtw(s,t);
    if isempty(is)
        c = Inf;
    end
end