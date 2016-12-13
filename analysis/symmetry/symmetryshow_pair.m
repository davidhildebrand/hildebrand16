function symmetryshow_pair(P,Q,sp,V,rmin,rmax,i,iname,k,kname,sampfreq)
    if nargin < 11 || isempty(sampfreq)
        % if no sampling frequency specified, default to 50
        sampfreq = 50;
    end

    plot3d = 1;
    plot2d = 0;
    
    SP = reflectedpoints(P,sp,V);
    s = SP(1:sampfreq:size(SP,1),:);
    t = Q(1:sampfreq:size(Q,1),:);
    n = size(t,1);
    [aligncost, is, it] = dtw(s,t,'PenalizeUnmatched');
    %[aligncost, is, it] = dtw(s,t);
    if isempty(is)
        aligncost = Inf;
    end
    % change ordering of one skeleton for plotting if necessary
    % (note this is done in dtw function if needed)
    if norm(s(1,:)-t(n,:)) < norm(s(1,:)-t(1,:))
        tp = flipud(t);
    else
        tp = t;
    end
    title(sprintf('%s (%d) vs %s (%d) has cost %f', ...
        iname, i, kname, k, aligncost))
    
    Pcol = [228/255 102/255 38/255];
    Qcol = [8/255 132/255 255/255];
    Lcol = [220/255 220/255 220/255];
    marksize = 3;
    
    if plot3d
        hold on
        %plot3(P(:,1),P(:,2),P(:,3),'.','Color',Pcol,'MarkerSize',marksize)
        plot3(SP(:,1),SP(:,2),SP(:,3),'.','Color',Pcol,'MarkerSize',marksize)
        plot3(Q(:,1),Q(:,2),Q(:,3),'.','Color',Qcol,'MarkerSize',marksize)
        fill3(sp(:,1),sp(:,2),sp(:,3),'k'), alpha(0.1) % plane
        % connect matches with a line
        if ~isempty(is)
            for pair = 1:1:length(is)
                plot3([s(is(pair),1) tp(it(pair),1)],...
                    [s(is(pair),2) tp(it(pair),2)],...
                    [s(is(pair),3) tp(it(pair),3)],...
                    '-','Color',Lcol)
            end
        end
        xlabel('x'), ylabel('z')
        hold off
        axis equal
        axis tight
        grid on
        axis on
        ax = gca;
        %ax.Projection = 'perspective';
        view(180,0)
        
        %dt = [s; t0; t];
        %mdt = min(dt);
        %Mdt = max(dt);
        % axis([mdt(1) Mdt(1) mdt(2) Mdt(2) mdt(3) Mdt(3)])

        % xlabel('x'), ylabel('y'), zlabel('z')
        % title(sprintf('dtw cost: %.02f',aligncost))
        % view(-200,60)
    end
    
    [SPX,SPY,SPZ] = symplanecoord(SP,sp,V);
    [sX,sY,sZ] = symplanecoord(s,sp,V);
    [QX,QY,QZ] = symplanecoord(Q,sp,V);
    [tpX,tpY,tpZ] = symplanecoord(tp,sp,V);
    
    if plot2d
        % subplot(1,4,3); grid off, axis equal, hold on
        subplot(1,2,1); grid off, axis equal, hold on
        xlim([min(vertcat(SPX(:),QX(:))) max(vertcat(SPX(:),QX(:)))])
        ylim([min(vertcat(SPZ(:),QZ(:))) max(vertcat(SPZ(:),QZ(:)))])
        plot(SPX,SPZ,'.','Color',Pcol,'MarkerSize',marksize), title('top')
        plot(QX,QZ,'.','Color',Qcol,'MarkerSize',marksize), title('top')
        if ~isempty(is)
            for pair = 1:1:length(is)
                plot([sX(is(pair)) tpX(it(pair))],...
                    [sZ(is(pair)) tpZ(it(pair))],...
                    '-','Color',Lcol)
            end
        end
        hsp3 = get(gca,'Position');
        xlabel('x'), ylabel('z')
        hold off
        
        %subplot(1,4,4); grid off, axis equal, hold on
        subplot(1,2,2); grid off, axis equal, hold on
        xlim([min(vertcat(SPY(:),QY(:))) max(vertcat(SPY(:),QY(:)))])
        ylim([min(vertcat(SPZ(:),QZ(:))) max(vertcat(SPZ(:),QZ(:)))])
        plot(SPY,SPZ,'.','Color',Pcol,'MarkerSize',marksize), title('side')
        plot(QY,QZ,'.','Color',Qcol,'MarkerSize',marksize), title('side')
        if ~isempty(is)
            for pair = 1:1:length(is)
                plot([sY(is(pair)) tpY(it(pair))],...
                    [sZ(is(pair)) tpZ(it(pair))],...
                    '-','Color',Lcol)
            end
        end
        hsp4 = get(gca,'Position');
        set(gca,'Position',[hsp4(1:2) hsp3(3:4)])
        xlabel('y'), ylabel('z')
    
        % d = 0.1*(rmax-rmin);
        % axis([rmin(1)-d(1) rmax(1)+d(1) rmin(2)-d(2) rmax(2)+d(2) rmin(3)-d(3) rmax(3)+d(3)])
        % view(0,90)
    end
end