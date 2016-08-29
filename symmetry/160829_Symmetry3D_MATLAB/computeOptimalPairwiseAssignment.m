%% --------------------------------------------------
clear
clc

%% --------------------------------------------------
fprintf('loading data\n');
% --------------------------------------------------

DataPath = 'D:\Dropbox (Personal)\MATLAB\Data\';
% set a date prefix for file saving
DateString = datestr(now,30);
DateString = DateString(3:length(DateString)-2);
Prefix = strcat(DateString,'_');

% load plane parameters
load(strcat(DataPath,filesep,'160812T1317_sp_300samples.mat')); % symmetry plane points
load(strcat(DataPath,filesep,'160812T1317_V_300samples.mat')); % symmetry plane perpendicular vector

% parse skeletons
D = importdata(strcat(DataPath,filesep,'160816t1040_130201zf142_160515SWiFT_ProjOrLngstLtL_ANNOTsymmetry_IGNblacklistsymblack_50nmLenThresh_PHYScoord_rootToNaN.txt'));
S = importdata(strcat(DataPath,filesep,'160803t1831_130201zf142_160515SWiFT_SUBSETspinalbackfills.txt'));
% S = importdata(strcat(DataPath,filesep,'160803t1831_130201zf142_160515SWiFT_SUBSETspinalbackfills.txt'));

iskels = zeros(1,length(S));
skelnames = cell(1,length(S));
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

%% compute pairwise matching costs
% all (distinct) skeleton pairs, not only selected
c = inf(npairs,1);
% set sampling frequency (one sample per sampfreq nodes)
sampfreq = 10;
tic
parfor pairindex = 1:npairs
    i = ik(pairindex,1);
    k = ik(pairindex,2);
    [P,~] = getnodes(D,iskels(i));
    [Q,~] = getnodes(D,iskels(k));
    %fprintf('comparing %d and %d\n',iskels(i),iskels(k));
%     if max(P(:,3))-min(P(:,3)) > 5*10^5 && max(Q(:,3))-min(Q(:,3)) > 5*10^5
        c(pairindex) = matchingcost(P,Q,sp,V,sampfreq);
%     end
end
toc

%% convert to matrix for pairwise assignment
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

save(strcat(DataPath,filesep,Prefix,'C.mat'),'C','-v7.3');

% --------------------------------------------------
fprintf('pairwise assignment\n');
% --------------------------------------------------

[assignment,cost] = munkres(C);
% [assignment,cost] = greedyassignment(C);

% for i = 1:size(C,1)
%     C(i,i) = 0;
% end
disp('assignment')
disp([1:length(assignment); assignment]')
for i = 1:length(assignment)
    fprintf('%s matched to %s\n',iskelnames{i},iskelnames{assignment(i)})
end

HeatMap(C,'RowLabels',iskelnames(:),'ColumnLabels',iskelnames(:),'DisplayRange',max(C(C<Inf))*0.5,'Symmetric',false)
CnoInf=C;
CnoInf(C==Inf)=9999999999999;
clustergram(CnoInf,'RowLabels',iskelnames(:),'ColumnLabels',iskelnames(:),'DisplayRange',max(C(C<Inf))*0.5,'Symmetric',false)

%%
% show final pairs
%%
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
        saveas(gcf,strcat(DataPath,filesep,Prefix,'pair_',...
            strrep(iskelnames{i},'_',''),'_',...
            strrep(iskelnames{assignment(i)},'_',''),'.png'));
    end
end
close all

% i = 667;
% j = 1056;
% [P,~] = getnodes(D,iskels(i));
% [Q,~] = getnodes(D,iskels(j));
% 
% SP = reflectedpoints(P,sp,V);
% plot3(P(:,1),P(:,2),P(:,3),'r.','MarkerSize',1), hold on
% plot3(P(1,1),P(1,2),P(1,3),'ro')
% plot3(Q(1,1),Q(1,2),Q(1,3),'go')
% plot3(Q(:,1),Q(:,2),Q(:,3),'g.','MarkerSize',1)
% fill3(sp(:,1),sp(:,2),sp(:,3),'k'), alpha(0.1), hold off
% grid on, axis equal
% xlabel('x'), ylabel('y'), zlabel('z')
