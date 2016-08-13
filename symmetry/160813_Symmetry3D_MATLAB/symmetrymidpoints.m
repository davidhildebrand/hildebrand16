function [ms,Sr,Lr] = symmetrymidpoints(P,TP,Q,TQ,samples,symthr,rg)
    if size(P,1) < size(Q,1)
        Sr = P; Lr = Q; % shorter / longer
        TSr = TP; TLr = TQ;
    else
        Sr = Q; Lr = P;
        TSr = TQ; TLr = TP;
    end

    n = samples;% 50; % n points on first skeleton
    ip = round(linspace(1,size(Sr,1),n));
    ms = [];
    dst = 0.01*rg;
    for j = 1:n
        p = Sr(ip(j),:);
        tp = TSr(ip(j),:)/norm(TSr(ip(j),:));
        indices = find(abs(Lr(:,2)-p(2)) < dst & abs(Lr(:,3)-p(3)) < dst);
        mmaxsym = [];
        maxsym = -Inf;
        for l = 1:length(indices)
            q = Lr(indices(l),:);
            tq = TLr(indices(l),:)/norm(TLr(indices(l),:));
            [m,s] = pairwisesymmetry(p,tp,q,tq);
            if s > symthr && s > maxsym
                maxsym = s;
                mmaxsym = m;
            end
        end
        if ~isempty(mmaxsym)
            ms = [ms; mmaxsym];
        end
    end
end