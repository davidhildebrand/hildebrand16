function SP = reflectedpoints(P,sp,V)
    mp = mean(sp);
    MP = repmat(mp,[size(P,1) 1]);
    VV = repmat(V,[size(P,1) 1]);
    DT = 2*sum((P-MP).*VV,2);
    SP = P-repmat(DT,[1 3]).*VV;
end