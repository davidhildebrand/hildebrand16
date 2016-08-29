function gm = go_evmf(V,k,s)
    % Performs Mean Shift Global Optimization on the Hypercylinder
    % Reference 1:
    % T. Kobayashi and N. Otsu, "Von Mises-Fisher Mean Shift for Clustering on a Hypersphere,"
    % Pattern Recognition (ICPR), 2010 20th International Conference on, Istanbul, 2010, pp. 2130-2133.
    % Reference 2:
    % https://en.wikipedia.org/wiki/Mean_shift

    [n,d] = size(V);
    KV = @(x,x0) exp(k*dot(x,x0)); % kernel for 'direction'
    KD = @(x,x0) exp(-((x-x0)/s)^2); % kernel for 'displacement'
    ND = @(x,x0) (norm(x-x0) < 2*s); % neighborhood for 'displacement'

    nglobmaxcand = 6; % number of global maxima candidates
    % if using 'parfor', could be a multiple of the number of cores
    
    C = zeros(nglobmaxcand,d); % convergence points
    
    randpermidx = randperm(n);
    idx = randpermidx(1:nglobmaxcand);
    
    parfor i = 1:nglobmaxcand
        m = V(idx(i),:); % starting point

        previousm = [-m(1:d-1) 0];
        count = 0;
        while count < 100 && dot(previousm(1:d-1),m(1:d-1)) < 0.999999 || abs(previousm(d)-m(d)) > 0.0001
            count = count+1;
            previousm = m;

            M = zeros(1,d);
            nD = 0; % normalization of 'displacement' dimension
            for j = 1:n
                M(1:d-1) = M(1:d-1)+KV(V(j,1:d-1),m(1:d-1))*V(j,1:d-1);
                kD = ND(V(j,d),m(d))*KD(V(j,d),m(d));
                M(d) = M(d)+kD*V(j,d);
                nD = nD+kD;
            end
            m = [M(1:d-1)/norm(M(1:d-1)) M(d)/nD];
        end

        C(i,:) = m;
    end

    c = C(1,:); % convergence points (one per cluster)
    l = zeros(nglobmaxcand,1); % cluster labels (one per selected sample)
    for i = 1:nglobmaxcand
        didbreak = 0;
        for j = 1:size(c,1)
            if dot(C(i,1:d-1),c(j,1:d-1)) > 0.999 && abs(C(i,d)-c(j,d)) < 0.01
                l(i) = j;
                didbreak = 1;
                break;
            end
        end
        if ~didbreak
            c = [c; C(i,:)];
        end
    end
    
    % pick cluster with the most votes
    nc = size(c,1);
    if nc == 1
        gm = c;
    else
        count = zeros(nc,1);
        for i = 1:nc
            count(i) = sum(l == i);
        end
        [~,im] = max(count);
        gm = c(im,:);
    end
end