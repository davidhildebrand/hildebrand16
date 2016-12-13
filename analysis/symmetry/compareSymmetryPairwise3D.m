clear
clc

%% settings
DataPath = 'D:\Dropbox (Personal)\MATLAB\Data\';
DataFile = '161107t1546_130201zf142_160515SWiFT_ProjOrLngstLtL_ANNOTsymmetry_IGNblacklistsymblack_1umLenThresh_PHYScoord_rootToNaN.txt';
%SubsetFile = '161107t1203_130201zf142_160515SWiFT_SUBSETspinalbackfillsIDENTnucMLFandMauthner.txt';
SubsetFile = '161107t1202_130201zf142_160515SWiFT_SUBSETspinalbackfillsIDENT.txt';

DateString = datestr(now,30);
DateString = strrep(DateString(3:length(DateString)-2),'T','t');
Prefix = strcat(DateString,'_');

%% --------------------------------------------------
fprintf('loading data\n');
% --------------------------------------------------

% load plane parameters (perpendicular vector and points)
load(strcat(DataPath,filesep,'161107t1647_plane_161107t1546data_161107t1201subsetMLF_ICPnosubsamp.mat'));

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

% --------------------------------------------------
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
figure(20)
plot3(roots(DP < 0,1),roots(DP < 0,2),roots(DP < 0,3),'.r'), hold on
plot3(roots(DP > 0,1),roots(DP > 0,2),roots(DP > 0,3),'.g'), hold on
fill3(sp(:,1),sp(:,2),sp(:,3),'k'), alpha(0.1), hold off
legend('left','right'), title('roots')

% compute pairwise matching costs
% set sampling frequency (one sample per sampfreq nodes)
sampfreq = round(1000/60);
%sampfreq = 50;
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
asgnm_mr = [1:length(assignment_mr); assignment_mr]';
asgnm_gd = [1:length(assignment_gd); assignment_gd]';

assignment = assignment_mr;
asgnm = asgnm_mr;

disp(asgnm)
fprintf('\n');

for i = 1:length(assignment)
    if assignment(i) == 0
        fprintf('%s (%d) matched NONE\n',iskelnames{i},iskels(i))
        continue
    end
    fprintf('%s (%d) matched to %s (%d)\n',...
        iskelnames{i},iskels(i),...
        iskelnames{assignment(i)},iskels(assignment(i)));
end

fprintf('\n');

% check matching by names
matches = NaN(size(assignment));
for i = 1:length(assignment)
    if assignment(i) == 0
        continue
    end
    if ~isempty(strfind(iskelnames{i},'_R'))
        targ = regexprep(iskelnames{i},'_\d','');
        test = strrep(regexprep(iskelnames{assignment(i)},'_\d',''),'_L','_R');
    end
    if ~isempty(strfind(iskelnames{i},'_L'))
        targ = regexprep(iskelnames{i},'_\d','');
        test = strrep(regexprep(iskelnames{assignment(i)},'_\d',''),'_R','_L');
    end
    if strcmp(targ,test)
        matches(i) = 1;
        %fprintf('%s (%d) and %s (%d) MATCHED as expected\n',...
        %    iskelnames{i},iskels(i),...
        %    iskelnames{assignment(i)},iskels(assignment(i)));
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

hmcolmap = flipud(double(parula));

% construct cost matrix with right vs left 
Clen = length(C);
Cord = C(round(Clen/2)+1:end,1:round(Clen/2));
NordL = iskelnames(1:round(Clen/2));
NordR = iskelnames(round(Clen/2)+1:end);

[assgnmCord_mr,costCord_mr] = munkres(Cord);
[assgnmCord_gd,costCord_gd,costsCord_gd] = greedyassignment2(Cord);
asgnmCord_mr = [1:length(assgnmCord_mr); assgnmCord_mr]';
asgnmCord_gd = [1:length(assgnmCord_gd); assgnmCord_gd]';
asgnmCord = asgnmCord_mr;

% figure(30); clf;
% heatmapcust(Cord,NordL(:),NordR(:),[],'ColorBar',1,'GridLines','-',...
%     'TickAngle',270,'ShowAllTicks',1,'UseLogColormap',false,...
%     'Colormap',hmcolmap,'MaxColorValue',7000,'MinColorValue',0);
% axis square
% ax = gca;
% ax.TickLength = [0 0];
% xlim = get(ax,'XLim');
% ylim = get(ax,'YLim');
% % add yellow outlines  where assignment was expected match
% % add magenta outlines where assignment was unexpected
% asgrect = asgnm(length(asgnm)/2+1:end,:);
% for a=1:length(asgrect)
%     if asgrect(a,2) == asgrect(a,1)-length(asgrect)
%         rectangle('Position',horzcat([(asgrect(a,2)-1) (asgrect(a,1)-length(asgrect)-1)]+0.5,[1 1]),...
%             'EdgeColor','y','LineWidth',1)
%         fprintf('expected\n')
%     else
%         rectangle('Position',horzcat([(asgrect(a,2)-1) (asgrect(a,1)-length(asgrect)-1)]+0.5,[1 1]),...
%             'EdgeColor','m','LineWidth',1)
%         fprintf('unexpected\n')
%     end
% end
% ix = find(Cord(:)==max(Cord(:)));
% [rx,cx] = ind2sub(size(Cord),ix);
% Cord(rx,cx)

% reorganize cost matrix to have left on vertical axis 
Cordreorg = flipud(rot90(Cord,1));
[assgnmCordreorg_mr,costCordreorg_mr] = munkres(Cordreorg);
[assgnmCordreorg_gd,costCordreorg_gd,costsCordreorg_gd] = greedyassignment2(Cordreorg);
asgnmCordreorg_mr = [1:length(assgnmCordreorg_mr); assgnmCordreorg_mr]';
asgnmCordreorg_gd = [1:length(assgnmCordreorg_gd); assgnmCordreorg_gd]';
asgnmCordreorg = asgnmCordreorg_mr;
figure(31); clf;
heatmapcust(Cordreorg,NordR(:),NordL(:),[],'ColorBar',1,'GridLines','-',...
    'TickAngle',270,'ShowAllTicks',1,'UseLogColormap',false,...
    'Colormap',hmcolmap,'MaxColorValue',6000,'MinColorValue',2000);
axis square
ax = gca;
ax.TickLength = [0 0];
% add yellow outlines  where assignment was expected match
% add magenta outlines where assignment was unexpected
asgrect = asgnmCordreorg; %asgnm(round(Clen/2)+1:end,:);
for c=1:length(asgrect)
    rectpos = horzcat([(asgrect(c,1)-1) (asgrect(c,2)-1)]+0.5,[1 1]);
    textpos = [(asgrect(c,1)) (asgrect(c,2)+0.25)];
    if asgrect(c,2) == asgrect(c,1)
        rectangle('Position',rectpos,'EdgeColor','k','LineWidth',1)
        fprintf('expected\n')
        text('Position',textpos,'String','*','Color','k',...
            'FontSize',20,'FontUnits','normalized',...
            'HorizontalAlignment','center','VerticalAlignment','middle');
    else
        rectangle('Position',rectpos,'EdgeColor','k','LineWidth',1)
        fprintf('unexpected\n')
        text('Position',textpos,'String','*','Color','r',...
            'FontSize',20,'FontUnits','normalized',...
            'HorizontalAlignment','center','VerticalAlignment','middle');
    end
end
% saveas(gcf,strcat(DataPath,filesep,Prefix,'DTWcostMatrix_AnteriorPosterior'),'epsc');
% saveas(gcf,strcat(DataPath,filesep,Prefix,'DTWcostMatrix_AnteriorPosterior'),'svg');
figure(32); clf;
heatmapcust(Cordreorg,NordR(:),NordL(:),[],'ColorBar',1,'GridLines','-',...
    'TickAngle',270,'ShowAllTicks',1,'UseLogColormap',false,...
    'Colormap',hmcolmap,'MaxColorValue',6000,'MinColorValue',2000);
axis square
ax = gca;
ax.TickLength = [0 0];
% ix = find(Cordreorg(:)==max(Cordreorg(:)));
% [rx,cx] = ind2sub(size(Cordreorg),ix);
% Cordreorg(rx,cx)
% NordL(rx)
% NordR(cx)

% reorganize cost matrix: sorted by diagonal values, left vertical axis 
Cordrsortd = flipud(rot90(Cord,1));
[~,srtidx] = sort(diag(Cordrsortd));
Cordrsortd = Cordrsortd(srtidx,:);
Cordrsortd = Cordrsortd(:,srtidx);
[assgnmCordrsortd_mr,costCordrsortd_mr] = munkres(Cordrsortd);
[assgnmCordrsortd_gd,costCordrsortd_gd,costsCordrsortd_gd] = greedyassignment2(Cordrsortd);
asgnmCordrsortd_mr = [1:length(assgnmCordrsortd_mr); assgnmCordrsortd_mr]';
asgnmCordrsortd_gd = [1:length(assgnmCordrsortd_gd); assgnmCordrsortd_gd]';
asgnmCordrsortd = asgnmCordrsortd_mr;
figure(33); clf;
heatmapcust(Cordrsortd,NordR(srtidx),NordL(srtidx),[],'ColorBar',1,'GridLines','-',...
    'TickAngle',270,'ShowAllTicks',1,'UseLogColormap',false,...
    'Colormap',hmcolmap,'MaxColorValue',6000,'MinColorValue',2000);
axis square
ax = gca;
ax.TickLength = [0 0];
% add yellow outlines  where assignment was expected match
% add magenta outlines where assignment was unexpected
asgrect = asgnmCordrsortd; %asgnm(round(Clen/2)+1:end,:);
for c=1:length(asgrect)
    rectpos = horzcat([(asgrect(c,1)-1) (asgrect(c,2)-1)]+0.5,[1 1]);
    textpos = [(asgrect(c,1)) (asgrect(c,2)+0.25)];
    if asgrect(c,2) == asgrect(c,1)
        rectangle('Position',rectpos,'EdgeColor','k','LineWidth',1)
        fprintf('expected\n')
        text('Position',textpos,'String','*','Color','k',...
            'FontSize',20,'FontUnits','normalized',...
            'HorizontalAlignment','center','VerticalAlignment','middle');
    else
        rectangle('Position',rectpos,'EdgeColor','k','LineWidth',1)
        fprintf('unexpected\n')
        text('Position',textpos,'String','*','Color','r',...
            'FontSize',20,'FontUnits','normalized',...
            'HorizontalAlignment','center','VerticalAlignment','middle');
    end
end
% saveas(gcf,strcat(DataPath,filesep,Prefix,'DTWcostMatrix_DiagonalSorted'),'epsc');
% saveas(gcf,strcat(DataPath,filesep,Prefix,'DTWcostMatrix_DiagonalSorted'),'svg');
figure(34); clf;
heatmapcust(Cordrsortd,NordR(srtidx),NordL(srtidx),[],'ColorBar',1,'GridLines','-',...
    'TickAngle',270,'ShowAllTicks',1,'UseLogColormap',false,...
    'Colormap',hmcolmap,'MaxColorValue',6000,'MinColorValue',2000);
axis square
ax = gca;
ax.TickLength = [0 0];
ix = find(Cordrsortd(:)==max(Cordrsortd(:)));
[rx,cx] = ind2sub(size(Cordrsortd),ix);
Cordrsortd(rx,cx)
ix = find(Cordrsortd(:)==min(Cordrsortd(:)));
[rx,cx] = ind2sub(size(Cordrsortd),ix);
Cordrsortd(rx,cx)

% generate cost matrix organized by hierarchical clustering
% T = clusterdata(Cordreorg,length(Cordreorg)/2);
% CordreorgnoInf=Cordreorg;
% CordreorgnoInf(Cordreorg==Inf)=9999999999999;
% cg = clustergram(CordreorgnoInf,'RowLabels',NordL(:),'ColumnLabels',NordR(:)','DisplayRange',min(C(:))*3,'Symmetric',false)

% HeatMap(Cord,'RowLabels',NordL(:),'ColumnLabels',NordR(:),'DisplayRange',min(Cord(:))*5,'Symmetric',false);

% HeatMap(C,'RowLabels',iskelnames(:),'ColumnLabels',iskelnames(:),'DisplayRange',min(C(:))*3,'Symmetric',false)
% CnoInf=C;
% CnoInf(C==Inf)=9999999999999;
% figure(31)
% clustergram(CnoInf,'RowLabels',iskelnames(:),'ColumnLabels',iskelnames(:),'DisplayRange',min(C(:))*3,'Symmetric',false)

% logC = log(C);
% logCnoInf=log(CnoInf);
% HeatMap(logC,'RowLabels',iskelnames(:),'ColumnLabels',iskelnames(:),'Symmetric',false)
% clustergram(logCnoInf,'RowLabels',iskelnames(:),'ColumnLabels',iskelnames(:),'DisplayRange',min(logCnoInf(:))*1.2,'Symmetric',false)

%% save
% save(strcat(DataPath,filesep,Prefix,'C_160908T1518planeSUBSETsbackfillsICPadj_dtwSampFreq10.mat'),'C','-v7.3');
save(strcat(DataPath,filesep,Prefix,'compOptPairAssgn_161107t1546data_161107t1203subsetNucMLFMauth_161107t1647planeMLF_OHpenalty_dtwFreq',num2str(sampfreq),'.mat'),'-v7.3');
save(strcat(DataPath,filesep,Prefix,'assignment_161107t1546data_161107t1203subsetNucMLFMauth_161107t1647planeMLF_OHpenalty_dtwFreq',num2str(sampfreq),'.mat'),...
    'assignment_mr','assignment_gd','iskelnames','iskels','cost_mr',...
    'cost_gd','costs_gd','C','asgnm','asgnm_mr','asgnm_gd','-v7.3')

%% show unexpected matches
% fig40 = figure(40);
% scsz = get(0,'ScreenSize'); % scsz = [left botton width height]
% fig40.Position = [0 0 scsz(3) scsz(4)];
% for i = 1:length(unexpected)
%     clf
%     [P,~] = getnodes(D,iskels(unexpected(i)));
%     if unexpected(i) > 0
%         [Q,~] = getnodes(D,iskels(assignment(unexpected(i))));
%         symmetryshow_pair(P,Q,sp,V,rmin,rmax,...
%             iskels(unexpected(i)),iskelnames{unexpected(i)},...
%             iskels(assignment(unexpected(i))),iskelnames{assignment(unexpected(i))},...
%             sampfreq);
%         saveas(gcf,sprintf('%s%s161107t1656_pairs%spair__orng_%s__blue_%s.png',DataPath,...
%             filesep,filesep,iskelnames{unexpected(i)},iskelnames{assignment(unexpected(i))}));
%         saveas(gcf,sprintf('%s%s161107t1656_pairs%spair__orng_%s__blue_%s',DataPath,...
%             filesep,filesep,iskelnames{unexpected(i)},iskelnames{assignment(unexpected(i))}),'epsc');
%         pause(0.1)
%     end
% end
 
%% show all assignment pairs
% fig41 = figure(41);
% scsz = get(0,'ScreenSize'); % scsz = [left botton width height]
% fig41.Position = [0 0 scsz(3) scsz(4)];
% for i = 1:length(assignment)
%     [P,~] = getnodes(D,iskels(i));
%     if assignment(i) > 0
%         disp(i)
%         if i ~= 25
%             continue
%         end
%         clf
%         [Q,~] = getnodes(D,iskels(assignment(i)));
%         symmetryshow_pair(P,Q,sp,V,rmin,rmax,...
%             iskels(i),iskelnames{i},...
%             iskels(assignment(i)),iskelnames{assignment(i)},...
%             sampfreq);
%         %saveas(gcf,sprintf('%s%s161107t1656_pairs%spair__orng_%s__blue_%s.png',DataPath,...
%         %    filesep,filesep,iskelnames{i},iskelnames{assignment(i)}));
%         %saveas(gcf,sprintf('%s%s161107t1656_pairs%spair__orng_%s__blue_%s',DataPath,...
%         %    filesep,filesep,iskelnames{unexpected(i)},iskelnames{assignment(unexpected(i))}),'epsc');
%         pause(0.1)
%     end
% end

%% show all pairs
% fig42 = figure(42);
% scsz = get(0,'ScreenSize'); % scsz = [left botton width height]
% fig42.Position = [0 0 scsz(3) scsz(4)];
% for i = 1:size(Cord,1)
%     for j = 1:size(Cord,2)
%         clf
%         k = j+size(Cord,1);
%         [P,~] = getnodes(D,iskels(i));
%         [Q,~] = getnodes(D,iskels(k));
%         symmetryshow_pair(P,Q,sp,V,rmin,rmax,...
%             iskels(i),iskelnames{i},...
%             iskels(k),iskelnames{k},...
%             sampfreq);
%         saveas(gcf,sprintf('%s%s161109t1700_pairs%spair__orng_%s__blue_%s.png',DataPath,...
%             filesep,filesep,iskelnames{i},iskelnames{k}));
%         saveas(gcf,sprintf('%s%s161109t1700_pairs%spair__orng_%s__blue_%s',DataPath,...
%             filesep,filesep,iskelnames{i},iskelnames{k}),'epsc');
%         pause(0.1)
%     end
% end