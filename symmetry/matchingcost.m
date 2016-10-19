function c = matchingcost(P,Q,sp,V,sampfreq)
    if nargin < 5 || ~isempty(sampfreq)
        % if no sampling frequency specified, default to 50
        sampfreq = 50;
    end
    
    SP = reflectedpoints(P,sp,V);
    s = SP(1:sampfreq:size(SP,1),:);
    t = Q(1:sampfreq:size(Q,1),:);
    [c, is, ~] = dtw(s,t);
    if isempty(is)
        c = Inf;
    end
end