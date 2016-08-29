%%
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

D = importdata(strcat(DataPath,filesep,'160825t1936_130201zf142_160515SWiFT_ProjOrLngstLtL_ANNOTsymmetry_IGNblacklistsymblack_20umLenThresh_PHYScoord_rootToNaN.txt'));
S = importdata(strcat(DataPath,filesep,'160825t1936_130201zf142_160515SWiFT_ProjOrLngstLtL_ANNOTsymmetry_IGNblacklistsymblack_20umLenThresh_UniqSkels.txt'));
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

% normalize coordinate space
mn = mean(D(:,4:6));
sd = std(D(:,4:6));
D(:,4:6) = D(:,4:6)-repmat(mn,[size(D,1) 1]);
D(:,4:6) = D(:,4:6)./repmat(sd,[size(D,1) 1]);

% plot3(D(:,4),D(:,5),D(:,6),'.')
% grid on, axis equal

% --------------------------------------------------
fprintf('getting range\n');
% --------------------------------------------------

[rmin, rmax] = getrange(D,iskels);
rg = max(rmax-rmin);

% --------------------------------------------------
fprintf('computing global plane of symmetry\n');
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
load(strcat(DataPath,filesep,'160825T2144_sp_160825t1936physCoord_sampfreq20.mat'),'sp'); % symmetry plane points
load(strcat(DataPath,filesep,'160825T2144_V_160825t1936physCoord_sampfreq20.mat'),'V'); % symmetry plane perpendicular vector

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
roots = roots.*repmat(sd,[size(roots,1) 1]);
roots = roots+repmat(mn,[size(roots,1) 1]);
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
sp0 = sp; sp0 = sp0-repmat(mn,[4 1]); sp0 = sp0./repmat(sd,[4 1]);

% get pairwise plane parameters
planeparam = nan(npairs,4); % for every 'reasonable' pair i,k: pperp (3d), displ (1d)
tic
parfor pairindex = 1:npairs
    % for pairindex = 1:npairs
    
    sampfreq = 20;
    i = ik(pairindex,1);
    k = ik(pairindex,2);
    
    % select by opposite sides
    sideI = rootside(rootside(:,1) == iskels(i),2);
    sideK = rootside(rootside(:,1) == iskels(k),2);
    oppositesides = 0;
    if ~isempty(sideI) && ~isempty(sideK)
        oppositesides = sideI(1)*sideK(1) < 0;
    end
    
    [P,TP] = getnodes(D,iskels(i));
    [Q,TQ] = getnodes(D,iskels(k));
    if size(P,1) > 0 && size(Q,1) > 0 && oppositesides
        [ms,~,~] = symmetrymidpoints(P,TP,Q,TQ,sampfreq,0.99,rg);
        if size(ms,1) > 5
            B = regress(ms(:,1),[ones(size(ms,1),1) ms(:,2:3)]);
            %B = robustfit(ms(:,2:3),ms(:,1));
            [pperp, displ, avgprjdst, sp] = symmetrytest_pairs(ms,B);
            
            c = matchingcost(P,Q,sp,pperp,sampfreq);
            %c = matchingcost(P,Q,sp,pperp);
            
            %plot3(P(:,1),P(:,2),P(:,3),'r.'), hold on
            %if oppositesides
            %    plot3(Q(:,1),Q(:,2),Q(:,3),'g.')
            %else
            %    plot3(Q(:,1),Q(:,2),Q(:,3),'r.')
            %end
            %fill3(sp0(:,1),sp0(:,2),sp0(:,3),'k'), alpha(0.1), hold off
            %grid on, axis equal
            %xlabel('x'), ylabel('y'), zlabel('z')
            %pause
            
            if dot(pperp,[1 0 0]) > 0.75 && c < 0.2
                planeparam(pairindex,:) = [pperp displ];
            end
        end
    end
end
toc

% load skeletons for visualization purposes
% and setup array for meanshift optimization
PS = [];
pperps = [];
displs = [];
loadedskel = zeros(1,nskels);
for pairindex = 1:npairs
    if ~isnan(planeparam(pairindex,1));
        i = ik(pairindex,1);
        k = ik(pairindex,2);
        pp = planeparam(pairindex,:);
        pperps = [pperps; pp(1:3)];
        displs = [displs; pp(4)];
        if ~loadedskel(i)
            [P,~] = getnodes(D,iskels(i));
            PS = [PS; P(1:50:size(P,1),:)];
            loadedskel(i) = 1;
        end
        if ~loadedskel(k)
            [P,~] = getnodes(D,iskels(k));
            PS = [PS; P(1:50:size(P,1),:)];
            loadedskel(k) = 1;
        end
    end
end

% plot perp. vector and displacement
[x,y,z] = sphere;
figure, subplot(1,2,1)
surf(x,y,z), view(76,12), hold on
plot3(pperps(:,1),pperps(:,2),pperps(:,3),'.k'), hold off
axis equal, grid on, xlabel('x'), ylabel('y'), zlabel('z'), title('perp vectors')
subplot(1,2,2)
hist(displs), title('displacements')

[p1,p2,p3,p4,V] = meanshiftplane(pperps,displs,PS);

% plot mean shift results, with selected skeletons
figure, subplot(1,1,1)
sp = [p1; p2; p3; p4];
plot3(PS(:,1),PS(:,2),PS(:,3),'.b'), hold on
v1 = p1; v2 = p1+norm(p2-p1)/2*V;
plot3([v1(1) v2(1)],[v1(2) v2(2)],[v1(3) v2(3)],'-')
fill3(sp(:,1),sp(:,2),sp(:,3),'k'), hold off
grid on, axis equal
% d = 0.1*rg;
% axis([rmin(1)-d rmax(1)+d rmin(2)-d rmax(2)+d rmin(3)-d rmax(3)+d])
alpha(0.1)
xlabel('x'), ylabel('y'), zlabel('z'), title('optimal plane')

% original scale
figure
PS2 = PS.*repmat(sd,[size(PS,1) 1]);
PS2 = PS2+repmat(mn,[size(PS2,1) 1]);
plot3(PS2(:,1),PS2(:,2),PS2(:,3),'.b'), hold on

sp2 = sp.*repmat(sd,[size(sp,1) 1]);
sp2 = sp2+repmat(mn,[size(sp2,1) 1]);
fill3(sp2(:,1),sp2(:,2),sp2(:,3),'k'), alpha(0.1), hold off
grid on, axis equal
xlabel('x'), ylabel('y'), zlabel('z'), title('optimal plane, original scale')

sp = sp2;

save(strcat(DataPath,filesep,Prefix,'sp_LRrevised.mat'),'sp','sp2','-v7.3'); % symmetry plane points
save(strcat(DataPath,filesep,Prefix,'V_LRrevised.mat'),'V','-v7.3'); % symmetry plane perpendicular vector
save(strcat(DataPath,filesep,Prefix,'compGlobSymPln_LRrevised.mat'),'-v7.3')
