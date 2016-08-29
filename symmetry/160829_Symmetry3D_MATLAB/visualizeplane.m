% clear
% clc

% --------------------------------------------------
fprintf('loading data\n');
% --------------------------------------------------

% skeletons
D = importdata('/home/mc457/files/CellBiology/IDAC/Marcelo/Hildebrand/160808t1608_130201zf142_160515SWiFT_ANNOTsymmetry_IGNblacklistsymblack_LongestLeafToLeaf_20umLenThresh_PHYScoord_rootToNaN.txt');
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
fprintf('plot data and planes\n');
% --------------------------------------------------

scsz = get(0,'ScreenSize'); % scsz = [left botton width height]
figure('Position',[scsz(3)/4 scsz(4)/4 scsz(3)/2 scsz(4)/2])

P = D(1:100:end,4:6);

plot3(P(:,1),P(:,2),P(:,3),'.'), hold on

load('/home/mc457/files/CellBiology/IDAC/Marcelo/Hildebrand/sp2.mat'); % symmetry plane points
load('/home/mc457/files/CellBiology/IDAC/Marcelo/Hildebrand/V2.mat'); % symmetry plane perpendicular vector
sp2 = sp; V2 = V;

fill3(sp2(:,1),sp2(:,2),sp2(:,3),'r'), alpha(0.1)

load('/home/mc457/files/CellBiology/IDAC/Marcelo/Hildebrand/sp.mat'); % symmetry plane points
load('/home/mc457/files/CellBiology/IDAC/Marcelo/Hildebrand/V.mat'); % symmetry plane perpendicular vector

fill3(sp(:,1),sp(:,2),sp(:,3),'g'), alpha(0.1), hold off

grid on, axis equal
xlabel('x'), ylabel('y'), zlabel('z'), title('optimal plane')
legend('data','plane before','plane after')
saveas(gcf,sprintf('~/Desktop/Planes.png'));

% --------------------------------------------------
fprintf('plot projections\n');
% --------------------------------------------------

mp = mean(sp);
MP = repmat(mp,[size(P,1) 1]);

vy = sp(2,:)-sp(1,:);
vy = vy/norm(vy);
vx = V2;
vz = cross(vx,vy);

VV = repmat(vx,[size(P,1) 1]);
PX = sum((P-MP).*VV,2);
VV = repmat(vy,[size(P,1) 1]);
PY = sum((P-MP).*VV,2);
VV = repmat(vz,[size(P,1) 1]);
PZ = sum((P-MP).*VV,2);

figure('Position',[scsz(3)/4 scsz(4)/4 scsz(3)/2 scsz(4)/2])
subplot(2,3,1), plot(PX,PY,'.','MarkerSize',1), title('top (before)')
subplot(2,3,2), plot(PX,PZ,'.','MarkerSize',1), title('front')
subplot(2,3,3), plot(PY,PZ,'.','MarkerSize',1), title('side')

vx = V;
vz = cross(vx,vy);

VV = repmat(vx,[size(P,1) 1]);
PX = sum((P-MP).*VV,2);
VV = repmat(vy,[size(P,1) 1]);
PY = sum((P-MP).*VV,2);
VV = repmat(vz,[size(P,1) 1]);
PZ = sum((P-MP).*VV,2);

subplot(2,3,4), plot(PX,PY,'.','MarkerSize',1), title('top (after)')
subplot(2,3,5), plot(PX,PZ,'.','MarkerSize',1), title('front')
subplot(2,3,6), plot(PY,PZ,'.','MarkerSize',1), title('side')
saveas(gcf,sprintf('~/Desktop/Projections.png'));