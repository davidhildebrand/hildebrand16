function symmetryshow_pair(P,Q,sp,V,rmin,rmax,i,iname,k,kname,sampfreq)

    if nargin < 11 || isempty(sampfreq)
        % if no sampling frequency specified, default to 50
        sampfreq = 50;
    end

    plot3(P(:,1),P(:,2),P(:,3),'.r'), hold on
    plot3(Q(:,1),Q(:,2),Q(:,3),'.g'),

    fill3(sp(:,1),sp(:,2),sp(:,3),'k')
    alpha(0.1)

    SP = reflectedpoints(P,sp,V);
    plot3(SP(:,1),SP(:,2),SP(:,3),'.r','MarkerSize',0.5)
    
    s = SP(1:sampfreq:size(SP,1),:);
    t = Q(1:sampfreq:size(Q,1),:);
    [aligncost, is, it] = dtw(s,t,'PenalizeUnmatched');
    if isempty(is)
        aligncost = Inf;
    end
    if ~isempty(is)
        for pair = 1:1:length(is)
            plot3([s(is(pair),1) t(it(pair),1)], [s(is(pair),2) t(it(pair),2)], [s(is(pair),3) t(it(pair),3)],'-k')
        end
    end
    title(sprintf('%s (%d) vs %s (%d) has cost %f', ...
        iname, i, kname, k, aligncost))
    hold off
    grid on, axis equal
%     d = 0.1*(rmax-rmin);
%     axis([rmin(1)-d(1) rmax(1)+d(1) rmin(2)-d(2) rmax(2)+d(2) rmin(3)-d(3) rmax(3)+d(3)])
%     view(0,90)
%     pause
end