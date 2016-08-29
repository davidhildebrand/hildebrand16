function [ms,Sr,Lr] = symmetrymidpoints(P,TP,Q,TQ,sampfreq,symthr,rg)
    if size(P,1) < size(Q,1)
        Sr = P; Lr = Q; % shorter / longer
        TSr = TP; TLr = TQ;
    else
        Sr = Q; Lr = P;
        TSr = TQ; TLr = TP;
    end

    ip = 1:sampfreq:size(Sr,1);
    n = size(ip,2);
    %samples = 50;
    %ip = round(linspace(1,size(Sr,1),samples));
    %n = samples;
    
    fprintf('DEBUG: size Sr    1 = %d\n',size(Sr,1));
    fprintf('DEBUG: size ip    2 = %d\n',size(ip,2));
    %fprintf('DEBUG: size ipold 2 = %d\n',size(ipold,2));
    fprintf('DEBUG: n = %d\n',n);
    
    %n = samples;% 50; % n points on first skeleton
    %ip = round(linspace(1,size(Sr,1),n));
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
    size(ms)
%     fprintf('DEBUG: size ms = %d\n',size(ms));
end