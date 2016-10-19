clear
clc

%% settings
DataPath = 'D:\Dropbox (Personal)\MATLAB\Data\';
DataFile = '161017t1551_130201zf142_160515SWiFT_ProjOrLngstLtL_ANNOTsymmetry_IGNblacklistsymblack_1umLenThresh_PHYScoord_rootToNaN.txt';
SubsetFile = '161018t1840_130201zf142_160515SWiFT_SUBSETspinalbackfillsIDENT.txt';

DateString = datestr(now,30);
DateString = strrep(DateString(3:length(DateString)-2),'T','t');
Prefix = strcat(DateString,'_');

%% --------------------------------------------------
fprintf('loading data\n');
% --------------------------------------------------

% plane parameters
% symmetry plane points
load(strcat(DataPath,'161018t2106_sp_from161017t1551exp161018t1840subsetICP.mat')); 
% symmetry plane perpendicular vector
load(strcat(DataPath,'161018t2106_V_from161017t1551exp161018t1840subsetICP.mat')); 

% skeletons
D = importdata(strcat(DataPath,filesep,DataFile));
% subset
S = importdata(strcat(DataPath,filesep,SubsetFile));

iskels = zeros(1,length(S));
iskelnames = cell(1,length(S));
if iscell(S)
    for i = 1:length(S)
        ss = strsplit(S{i},' ');
        iskels(i) = str2double(ss{1});
        iskelnames{i} = ss{3};
    end
elseif isnumeric(S)
    for i = 1:length(S)
        iskels(i) = S(i);
        % set name to skelID if no name in file
        iskelnames{i} = S(i);
    end
end

%% --------------------------------------------------
fprintf('getting range\n');
% --------------------------------------------------

if ~exist('rg','var')
    [rmin, rmax] = getrange(D,iskels);
    rg = max(rmax-rmin);
end

%% --------------------------------------------------
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

% find which side root falls in
nroots = sum(sum(isnan(D(:,3))));
roots = zeros(nroots,3);
rootside = zeros(nroots,2); % skeleton id, side (-1 or 1)
count = 0;
for i = 1:size(D,1)
    if isnan(D(i,3))
        count = count+1;
        roots(count,:) = D(i,4:6);
        rootside(count) = D(i,1);
    end
end
mp = mean(sp);
MP = repmat(mp,[size(roots,1) 1]);
VV = repmat(V,[size(roots,1) 1]);
DP = sum((roots-MP).*VV,2);
rootside(:,2) = sign(DP);
figure
plot3(roots(DP < 0,1),roots(DP < 0,2),roots(DP < 0,3),'.r'), hold on
plot3(roots(DP > 0,1),roots(DP > 0,2),roots(DP > 0,3),'.g'), hold on
fill3(sp(:,1),sp(:,2),sp(:,3),'k'), alpha(0.1), hold off
legend('left','right'), title('roots')

% compute pairwise matching costs
% set sampling frequency (one sample per sampfreq nodes)
sampfreq = round(2*1000/60);
% all (distinct) skeleton pairs, not only selected
c = inf(npairs,1);
tic
parfor pairindex = 1:npairs
% for pairindex = 1:npairs
%     fprintf('%d / %d\n', pairindex, npairs);
    i = ik(pairindex,1);
    k = ik(pairindex,2);
    
    % select by opposite sides
    sideI = rootside(rootside(:,1) == iskels(i),2);
    sideK = rootside(rootside(:,1) == iskels(k),2);
    oppositesides = 0;
    if ~isempty(sideI) && ~isempty(sideK)
        oppositesides = sideI(1)*sideK(1) < 0;
    end
    
    if oppositesides
        [P,~] = getnodes(D,iskels(i));
        [Q,~] = getnodes(D,iskels(k));
        %c(pairindex) = matchingcost(P,Q,sp,V);
        %c(pairindex) = matchingcost(P,Q,sp,V,sampfreq);
        %c(pairindex) = matchingcost_overhangpenalty(P,Q,sp,V);
        c(pairindex) = matchingcost_overhangpenalty(P,Q,sp,V,sampfreq);
    else
        c(pairindex) = Inf;
    end
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

%% --------------------------------------------------
fprintf('pairwise assignment\n');
% --------------------------------------------------

[assignment_mr,cost_mr] = munkres(C);
[assignment_gd,cost_gd,costs_gd] = greedyassignment2(C);

assignment = assignment_gd;

% disp('assignment')
disp([1:length(assignment); assignment]')
for i = 1:length(assignment)
    if assignment(i) == 0
        fprintf('%s (%d) matched NONE\n',iskelnames{i},iskels(i))
        continue
    end
    fprintf('%s (%d) matched to %s (%d)\n',...
        iskelnames{i},iskels(i),...
        iskelnames{assignment(i)},iskels(assignment(i)));
end

% check matching by names
matches = NaN(size(assignment));
for i = 1:length(assignment)
    if assignment(i) == 0
        continue
    end
    if ~isempty(strfind(iskelnames{i},'_R'))
        targ = iskelnames{i};
        test = strrep(iskelnames{assignment(i)},'_L','_R');
    end
    if ~isempty(strfind(iskelnames{i},'_L'))
        targ = iskelnames{i};
        test = strrep(iskelnames{assignment(i)},'_R','_L');
    end
    if strcmp(targ,test)
        matches(i) = 1;
        fprintf('%s (%d) and %s (%d) MATCHED as expected\n',...
            iskelnames{i},iskels(i),...
            iskelnames{assignment(i)},iskels(assignment(i)));
    else
        matches(i) = 0;
        fprintf('%s (%d) and %s (%d) did NOT match expected\n',...
            iskelnames{i},iskels(i),...
            iskelnames{assignment(i)},iskels(assignment(i)));
    end
end
clear targ test
unexpected = sort(assignment(matches==0)); %iskelnames(assignment(matches==0))

%% --------------------------------------------------
fprintf('show cost matrix\n');
% --------------------------------------------------
HeatMap(C,'RowLabels',iskelnames(:),'ColumnLabels',iskelnames(:),'DisplayRange',min(C(C<Inf))*5,'Symmetric',false)
CnoInf=C;
CnoInf(C==Inf)=9999999999999;
clustergram(CnoInf,'RowLabels',iskelnames(:),'ColumnLabels',iskelnames(:),'DisplayRange',min(C(C<Inf))*5,'Symmetric',false)

%% save
% save(strcat(DataPath,filesep,Prefix,'C_160908T1518planeSUBSETsbackfillsICPadj_dtwSampFreq10.mat'),'C','-v7.3');
save(strcat(DataPath,filesep,Prefix,'compOptPairAssgn_161018t1840subset_161018t2106plane_OHpenalty_dtwFreq',num2str(sampfreq),'.mat'),'-v7.3');
save(strcat(DataPath,filesep,Prefix,'assignment_161018t1840subset_161018t2106plane_OHpenalty_dtwFreq',num2str(sampfreq),'.mat'),...
    'assignment_mr','assignment_gd','iskelnames','iskels','cost_mr',...
    'cost_gd','costs_gd','C','-v7.3')

%% show unexpected matches
scsz = get(0,'ScreenSize');
% scsz = [left botton width height]
figure('Position',[scsz(3)/4 scsz(4)/4 scsz(3)/2 scsz(4)/2])
for i = 1:length(unexpected)
    [P,~] = getnodes(D,iskels(unexpected(i)));
    if unexpected(i) > 0
        [Q,~] = getnodes(D,iskels(assignment(unexpected(i))));
        symmetryshow_pair(P,Q,sp,V,rmin,rmax,...
            iskels(unexpected(i)),iskelnames{unexpected(i)},...
            iskels(assignment(unexpected(i))),iskelnames{assignment(unexpected(i))},...
            sampfreq);
        pause
%         saveas(gcf,sprintf('%spair_%d_%d.png',DataPath,iskels(i),iskels(assignment(i))));
    end
end

%% show all final pairs
% scsz = get(0,'ScreenSize'); % scsz = [left botton width height]
% figure('Position',[scsz(3)/4 scsz(4)/4 scsz(3)/2 scsz(4)/2])
% for i = 1:length(assignment)
%     [P,~] = getnodes(D,iskels(i));
%     if assignment(i) > 0
%         [Q,~] = getnodes(D,iskels(assignment(i)));
%         
% %         SP = reflectedpoints(P,sp,pperp);
% %         plot3(P(:,1),P(:,2),P(:,3),'r.'), hold on
% %         plot3(Q(:,1),Q(:,2),Q(:,3),'g.')
% %         fill3(sp(:,1),sp(:,2),sp(:,3),'k'), alpha(0.1), hold off
% %         grid on, axis equal
% %         xlabel('x'), ylabel('y'), zlabel('z')
% %         pause
% 
%         symmetrytest_all(P,Q,sp,V,rmin,rmax,i,assignment(i));
%         pause
% %         saveas(gcf,sprintf('Pair_%d_%d.png',iskels(i),iskels(assignment(i))));
%     end
% end
% close all


