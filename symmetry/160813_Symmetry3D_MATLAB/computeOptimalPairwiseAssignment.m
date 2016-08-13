% --------------------------------------------------
fprintf('loading data\n');
% --------------------------------------------------

% parameters
load('/home/mc457/files/CellBiology/IDAC/Marcelo/Hildebrand/sp.mat'); % symmetry plane points
load('/home/mc457/files/CellBiology/IDAC/Marcelo/Hildebrand/V.mat'); % symmetry plane perpendicular vector

% skeletons
D = importdata('/home/mc457/files/CellBiology/IDAC/Marcelo/Hildebrand/160803t1831_130201zf142_160515SWiFT_ANNOTsymmetry_20000lengththresh_PHYScoord.txt');
S = importdata('/home/mc457/files/CellBiology/IDAC/Marcelo/Hildebrand/160803t1831_130201zf142_160515SWiFT_ANNOTsymmetry_20umLengthThresh_UniqSkels.txt');
iskels = zeros(1,length(S));
if iscell(S)
    for i = 1:length(S)
        ss = strsplit(S{i},' ');
        iskels(i) = str2double(ss{1});
    end
elseif isnumeric(S)
    for i = 1:length(S)
        iskels(i) = S(i);
    end
end

% --------------------------------------------------
fprintf('getting range\n');
% --------------------------------------------------

if ~exist('rg','var')
    [rmin, rmax] = getrange(D,iskels);
    rg = max(rmax-rmin);
end

% --------------------------------------------------
fprintf('compute pairwise matching cost\n');
% --------------------------------------------------

% setup pair indexing
nskels = length(iskels);
npairs = nskels*(nskels-1)/2;
ik = zeros(npairs,2);
count = 0;
for i = 1:nskels-1
    for k = i+1:nskels
        count = count+1;
        ik(count,:) = [i k];
    end
end

% compute pairwise matching costs
% all (distinct) skeleton pairs, not only selected
c = zeros(npairs,1);
tic
parfor pairindex = 1:npairs
    i = ik(pairindex,1);
    k = ik(pairindex,2);
    [P,~] = getnodes(D,iskels(i));
    [Q,~] = getnodes(D,iskels(k));
    c(pairindex) = matchingcost(P,Q,sp,V);
end
toc

% convert to matrix for pairwise assignment
C = zeros(nskels,nskels);
for pairindex = 1:npairs
    i = ik(pairindex,1);
    k = ik(pairindex,2);
    C(i,k) = c(pairindex);
    C(k,i) = c(pairindex);
end
for i = 1:size(C,1)
    C(i,i) = Inf;
end

% --------------------------------------------------
fprintf('pairwise assignment\n');
% --------------------------------------------------
%%
[assignment,cost] = munkres(C);
% [assignment,cost] = greedyassignment(C);

disp('assignment')
disp([1:length(assignment); assignment]')
%%
% show final pairs

scsz = get(0,'ScreenSize'); % scsz = [left botton width height]
figure('Position',[scsz(3)/4 scsz(4)/4 scsz(3)/2 scsz(4)/2])
for i = 1:length(assignment)
    [P,~] = getnodes(D,iskels(i));
    if assignment(i) > 0
        [Q,~] = getnodes(D,iskels(assignment(i)));
        
%         SP = reflectedpoints(P,sp,pperp);
%         plot3(P(:,1),P(:,2),P(:,3),'r.'), hold on
%         plot3(Q(:,1),Q(:,2),Q(:,3),'g.')
%         fill3(sp(:,1),sp(:,2),sp(:,3),'k'), alpha(0.1), hold off
%         grid on, axis equal
%         xlabel('x'), ylabel('y'), zlabel('z')
%         pause

        symmetrytest_all(P,Q,sp,V,rmin,rmax,i,assignment(i));
    end
end
close all