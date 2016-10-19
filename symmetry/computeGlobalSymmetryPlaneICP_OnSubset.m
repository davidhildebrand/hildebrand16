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
fprintf('initialize\n');
% --------------------------------------------------

% skeletons
D = importdata(strcat(DataPath,filesep,DataFile));
% subset
S = importdata(strcat(DataPath,filesep,SubsetFile));

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

P = [];
for i = 1:length(iskels)
    [PI,~] = getnodes(D,iskels(i));
    P = [P; PI];
end
D = P;

%% subsample
P = D;%(1:10:end,:);

%% --------------------------------------------------
fprintf('reflect\n');
% --------------------------------------------------

% z = P(:,3);
% P = P(z > 120000 & z < 500000,:);

m = mean(P);

mn = min(P);
mx = max(P);
rg = (mx-mn)/2;

sp = [m(1) m(2)-rg(2) m(3)-rg(3);...
      m(1) m(2)+rg(2) m(3)-rg(3);...
      m(1) m(2)+rg(2) m(3)+rg(3);...
      m(1) m(2)-rg(2) m(3)+rg(3)];
V = [1 0 0];

Q = reflectedpoints(P,sp,V);

%% --------------------------------------------------
fprintf('register\n');
% --------------------------------------------------

fixed = pointCloud(P);
moving = pointCloud(Q);

%[tform,ptCloudAligned,rmse] = pcregrigid(moving, fixed, 'MaxIterations', 100);
[tform,ptCloudAligned,rmse] = pcregrigid(moving, fixed);

R = ptCloudAligned.Location;

scsz = get(0,'ScreenSize'); % scsz = [left botton width height]
figure('Position',[scsz(3)/4 scsz(4)/4 scsz(3)/2 scsz(4)/2])

ax = zeros(1,3);
ax(1) = subplot(1,3,1);
plot3(P(:,1),P(:,2),P(:,3),'r.','MarkerSize',1),hold on
plot3(Q(:,1),Q(:,2),Q(:,3),'g.','MarkerSize',1)
fill3(sp(:,1),sp(:,2),sp(:,3),'k'), alpha(0.1), hold off
grid on, axis equal, axis tight, axis vis3d
xlabel('x'), ylabel('y'), zlabel('z')
view(0,0)
a = gca;
a.Projection = 'perspective';

ax(2) = subplot(1,3,2);
plot3(P(1:10:end,1),P(1:10:end,2),P(1:10:end,3),'r.','MarkerSize',1),hold on
plot3(R(1:10:end,1),R(1:10:end,2),R(1:10:end,3),'g.','MarkerSize',1)
MS = 0.5*(P+R);
plot3(MS(1:10:end,1),MS(1:10:end,2),MS(1:10:end,3),'.k','MarkerSize',1)
hold off
grid on, axis equal, axis tight, axis vis3d
xlabel('x'), ylabel('y'), zlabel('z')
view(0,0)
a = gca;
a.Projection = 'perspective';
title(sprintf('rmse = %f',rmse))

%% --------------------------------------------------
fprintf('fit plane\n');
% --------------------------------------------------

B = regress(MS(:,1),[ones(size(MS,1),1) MS(:,2:3)]);
% B = robustfit(MS(:,2:3),MS(:,1));
[pperp, displ, avgprjdst, sp] = symmetrytest_pairs(MS,B);
V = pperp;

ax(3) = subplot(1,3,3);
plot3(P(:,1),P(:,2),P(:,3),'r.','MarkerSize',1),hold on
fill3(sp(:,1),sp(:,2),sp(:,3),'k'), alpha(0.1), hold off
hold off
grid on, axis equal, axis tight, axis vis3d
xlabel('x'), ylabel('y'), zlabel('z')
view(0,0)
a = gca;
a.Projection = 'perspective';

linkprop(ax,{'CameraPosition','CameraViewAngle'});

%% save
% symmetry plane perpendicular vector
save(strcat(DataPath,filesep,Prefix,'V_from161017t1551exp161018t1840subsetICP.mat'),'V');
% symmetry plane points
save(strcat(DataPath,filesep,Prefix,'sp_from161017t1551exp161018t1840subsetICP.mat'),'sp');