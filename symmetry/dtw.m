function [aligncost, is, it, flipt] = dtw(s,t,varargin)

m = size(s,1);
n = size(t,1);

% sort skeletons so that both directions match
if norm(s(1,:)-t(n,:)) < norm(s(1,:)-t(1,:))
    t = flipud(t);
    flipt = 1;
else
    flipt = 0;
end
% can sort each by z, but this will not respect treeline
% [~,idx] = sort(s(:,3)); s = s(idx,:);
% [~,idx] = sort(t(:,3)); t = t(idx,:);

% plot3(s(:,1),s(:,2),s(:,3)), hold on
% plot3(s(1,1),s(1,2),s(1,3),'o')
% plot3(t(1,1),t(1,2),t(1,3),'o')
% plot3(t(:,1),t(:,2),t(:,3)), hold off
% pause

% cgap = Inf;
% for i = 1:m
%     for j = 1:n
%         d = norm(s(i,:)-t(j,:));
%         if d < cgap
%             cgap = d;
%         end
%     end
% end
cgap = 0;

C = zeros(m+1,n+1);
PR = C; % row predecessor
PC = C; % col predecessor

C(1,1) = 0;
for i = 2:m+1
    C(i,1) = C(i-1,1)+cgap;
    PR(i,1) = i-1;
    PC(i,1) = 1;
end
for j = 2:n+1
    C(1,j) = C(1,j-1)+cgap;
    PR(1,j) = 1;
    PC(1,j) = j-1;
end

for i = 2:m+1
%     if mod(i,floor(m/10)) == 1
%         fprintf('.')
%     end
    for j = 2:n+1
        cij = norm(s(i-1,:)-t(j-1,:));
        [mn,imn] = min([C(i-1,j)+cgap C(i-1,j-1) C(i,j-1)+cgap]);
        C(i,j) = cij+mn;
        if imn == 1
            PR(i,j) = i-1;
            PC(i,j) = j;
        elseif imn == 2
            PR(i,j) = i-1;
            PC(i,j) = j-1;
        else
            PR(i,j) = i;
            PC(i,j) = j-1;
        end
    end
end

% imshow(C/max(max(C)))
% pause

% col = n+1;
% row = m+1;
% if C(m+1,n) < C(m,n+1) % sequence s 'ends earlier'
%     while col > 1 && C(m+1,col) > C(m+1,col-1)
%         col = col-1;
%     end
% else % sequence t 'ends earlier'
%     while row > 1 && C(row,n+1) > C(row-1,n+1)
%         row = row-1;
%     end
% end

col = n+1;
while col > 1 && C(m+1,col) > C(m+1,col-1)
    col = col-1;
end
r1 = m+1; c1 = col;
row = m+1;
while row > 1 && C(row,n+1) > C(row-1,n+1)
    row = row-1;
end
r2 = row; c2 = n+1;
if C(r1,c1) < C(r2,c2)
    if C(r1,c1) > 0
        row = r1; col = c1;
    else
        row = r2; col = c2;
    end
else
    if C(r2,c2) > 0
        row = r2; col = c2;
    else
        row = r1; col = c1;
    end
end

is = [];
it = [];
P = repmat(C/max(max(C)),[1 1 3]);
P(row,col,[1 3]) = 0;
P(row,col,2) = 1;
aligncost = C(row,col);

if row > 1 && col > 1
    is = [is row-1];
    it = [it col-1];
end

if ~isempty(is)
    while 1
        r = PR(row,col);
        c = PC(row,col);
        row = r;
        col = c;
        P(row,col,[1 3]) = 0;
        P(row,col,2) = 1;
        if row == 1 || col == 1
            break
        else
            is = [is row-1];
            it = [it col-1]; 
        end
    end
    aligncost = (aligncost-C(row,col))/length(is);
else
    aligncost = Inf;
end

is = fliplr(is);
it = fliplr(it);

% penalize unmatched portions
if ~isempty(varargin) && strcmp(varargin{1},'PenalizeUnmatched') && aligncost < Inf
    lengthMatchedS = 0;
    for i = is(1)+1:is(end)
        lengthMatchedS = lengthMatchedS+norm(s(i,:)-s(i-1,:));
    end
    lengthMatchedT = 0;
    for i = it(1)+1:it(end)
        lengthMatchedT = lengthMatchedT+norm(t(i,:)-t(i-1,:));
    end
    lengthS = 0;
    for i = 2:length(s)
        lengthS = lengthS+norm(s(i,:)-s(i-1,:));
    end
    lengthT = 0;
    for i = 2:length(t)
        lengthT = lengthT+norm(t(i,:)-t(i-1,:));
    end
    % propMatchedS = lengthMatchedS/lengthS;
    % propMatchedT = lengthMatchedT/lengthT;

    penaltyUnmatchedS = lengthS/lengthMatchedS;
    penaltyUnmatchedT = lengthT/lengthMatchedT;

    aligncost = aligncost*penaltyUnmatchedS*penaltyUnmatchedT;
end

% if ~isempty(is)
%     figure
%     imshow(P)
%     figure
%     plot3(s(:,1),s(:,2),s(:,3),'r.'), hold on
%     plot3(t(:,1),t(:,2),t(:,3),'g.')
%     for i = 1:1:length(is)
%         plot3([s(is(i),1) t(it(i),1)], [s(is(i),2) t(it(i),2)], [s(is(i),3) t(it(i),3)], '-k')
%     end
%     hold off
%     legend('s','t')
%     title(sprintf('%.02f', aligncost))
%     xlabel('x'), ylabel('y'), zlabel('z')
%     axis equal
%     pause
% end

end