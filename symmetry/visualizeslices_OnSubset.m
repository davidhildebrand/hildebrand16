clear
clc

%% settings
DataPath = 'D:\Dropbox (Personal)\MATLAB\Data\';
DataFile = '161017t1551_130201zf142_160515SWiFT_ProjOrLngstLtL_ANNOTsymmetry_IGNblacklistsymblack_1umLenThresh_PHYScoord_rootToNaN.txt';
SubsetFile = '161020t1703_130201zf142_160515SWiFT_SUBSETspinalbackfillsIDENTnoRoM1R.txt';

DateString = datestr(now,30);
DateString = strrep(DateString(3:length(DateString)-2),'T','t');
Prefix = strcat(DateString,'_');

% z step size
%step = 100;
step = round(2*5000/60); % visualize every 5um

%% --------------------------------------------------
fprintf('loading data\n');
% --------------------------------------------------

% load plane parameters (perpendicular vector and points)
load(strcat(DataPath,filesep,'161020t1707_plane_from161017t1551exp161020t1700subsetICPnosubsamp.mat')); 

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
fprintf('plot data and planes\n');
% --------------------------------------------------

scsz = get(0,'ScreenSize'); % scsz = [left botton width height]
figure('Position',[scsz(3)/4 scsz(4)/4 scsz(3)/2 scsz(4)/2])

P = D(1:step:end,4:6);

plot3(P(:,1),P(:,2),P(:,3),'.'), hold on

fill3(sp(:,1),sp(:,2),sp(:,3),'g'), alpha(0.1), hold off

grid on, axis equal
xlabel('x'), ylabel('y'), zlabel('z'), title('optimal plane')
% saveas(gcf,sprintf('~/Desktop/Planes.png'));

%% --------------------------------------------------
fprintf('plot projections\n');
% --------------------------------------------------

[PX,PY,PZ] = symplanecoord(P,sp,V);

subplot(1,3,1), plot(PX,PY,'.','MarkerSize',1), title('top'), axis equal
subplot(1,3,2), plot(PX,PZ,'.','MarkerSize',1), title('front'), axis equal
subplot(1,3,3), plot(PY,PZ,'.','MarkerSize',1), title('side'), axis equal
% saveas(gcf,sprintf('~/Desktop/Projections.png'));

zMin = min(PZ);
zMax = max(PZ);
xMin = min(PX);
xMax = max(PX);
yMin = min(PY);
yMax = max(PY);

%% unique assignment pairs, and color scheme

load(strcat(DataPath,filesep,'161024t1132_assignment_161020t1701subset_161020t1707plane_OHpenalty_dtwFreq17.mat'),'asgnm_gd');
asgnm = asgnm_gd;
for i = size(asgnm,1):-1:1
    if asgnm(i,2) == 0
        asgnm(i,:) = [];
    end
end

asgnm2 = asgnm(1,:);
for i = 2:size(asgnm,1)
    t = 0;
    for j = 1:i-1
        if asgnm(i,1) == asgnm(j,2)
            t = 1;
            break;
        end
    end
    if t == 0
        asgnm2 = [asgnm2; asgnm(i,:)];
    end
end
asgnm = asgnm2;

npairs = size(asgnm,1);
hues = (randperm(npairs)-1)/npairs*2/3;
hsvs = [hues' ones(npairs,1) ones(npairs,1)];
rgbs = zeros(size(hsvs));

%% record skeleton pairs (for computational speed)
sqpairs = cell(npairs,2);
for i = 1:npairs
    fprintf('%f\n',i/npairs);
    rgbs(i,:) = hsv2rgb(hsvs(i,:));
    [sqA,~] = getnodes(D,iskels(asgnm(i,1)));
    [sqB,~] = getnodes(D,iskels(asgnm(i,2)));
    sqpairs{i,1} = sqA;%sqA(1:step:end,:);
    sqpairs{i,2} = sqB;%sqB(1:step:end,:);
end
%save('/home/mc457/Desktop/sqpairs.mat','sqpairs');
% load('/home/mc457/files/CellBiology/IDAC/Marcelo/Hildebrand/sqpairs.mat');

%% display slices in figure

fig40 = figure(40);
scsz = get(0,'ScreenSize'); % scsz = [left botton width height]
fig40.Position = [scsz(3)/4 scsz(4)/4 scsz(3)/2 scsz(4)/2];
allX = []; allY = []; allZ = [];
for pair = 1:npairs
    xyz = sqpairs{pair,1};
    xyz = [xyz; sqpairs{pair,2}];
    [x,y,z] = symplanecoord(xyz,sp,V);
    allX = [allX; x];
    allY = [allY; y];
    allZ = [allZ; z];
end
xmin = min(allX); xmax = max(allX);
ymin = min(allY); ymax = max(allY);
zmin = min(allZ); zmax = max(allZ);

nSlices = floor((zmax-zmin)/60); %500;

for i = 1:nSlices
    % entire set, 3D
    z0 = zmin+(i-1)/nSlices*(zmax-zmin);
    z1 = zmin+i/nSlices*(zmax-zmin);
    idx = find(PZ <= z0 | PZ >= z1);
    subplot(1,3,1)
    plot3(PX(idx),PY(idx),PZ(idx),'.r','MarkerSize',1), hold on
    idx = find(PZ > z0 & PZ < z1);
    plot3(PX(idx),PY(idx),PZ(idx),'.g','MarkerSize',1), hold off
    axis([xMin xMax yMin yMax zMin zMax])
    view(0,0)
    
    % subset, 3D
    subplot(1,3,2)
    
    xyz = sqpairs{1,1};
    xyz = [xyz; sqpairs{1,2}];
    [x,y,z] = symplanecoord(xyz,sp,V);
    idx = find(z <= z0 | z >= z1);
    plot3(x(idx),y(idx),z(idx),'.','MarkerSize',1,'color',rgbs(1,:)), hold on
    idx = find(z > z0 & z < z1);
    plot3(x(idx),y(idx),z(idx),'.k','MarkerSize',1)
    for pair = 2:npairs
        xyz = sqpairs{pair,1};
        xyz = [xyz; sqpairs{pair,2}];
        [x,y,z] = symplanecoord(xyz,sp,V);
        idx = find(z <= z0 | z >= z1);
        plot3(x(idx),y(idx),z(idx),'.','MarkerSize',1,'color',rgbs(pair,:)), hold on
        idx = find(z > z0 & z < z1);
        plot3(x(idx),y(idx),z(idx),'.k','MarkerSize',1)
    end
    hold off
    axis([xMin xMax yMin yMax zMin zMax])
    view(0,0)

    s3 = subplot(1,3,3);
    cla(s3)
    hold on
    for pair = 1:npairs
        xyz = sqpairs{pair,1};
        xyz = [xyz; sqpairs{pair,2}];
        [x,y,z] = symplanecoord(xyz,sp,V);
        plot(x,y,'.','MarkerSize',1,'color',hsv2rgb([1 0.25 1].*hsvs(pair,:)))
    end
    for pair = 1:npairs
        sqA = sqpairs{pair,1};
        sqB = sqpairs{pair,2};
        if ~isempty(sqA)
            [x,y,z] = symplanecoord(sqA,sp,V);
            idx = find(z > z0 & z < z1);
            if ~isempty(idx)
                plot(mean(x(idx)),mean(y(idx)),'.','MarkerSize',10,'color',rgbs(pair,:))
            end
        end
        if ~isempty(sqB)
            [x,y,z] = symplanecoord(sqB,sp,V);
            idx = find(z > z0 & z < z1);
            if ~isempty(idx)
                plot(mean(x(idx)),mean(y(idx)),'.','MarkerSize',10,'color',rgbs(pair,:))
            end
        end
    end
    axis([xmin xmax ymin ymax])
    hold off
%     pause
%     saveas(gcf,sprintf('~/Desktop/Slices/Slice%05d.png',i));
    pause(0.01)
end
% close all


%% display slices on image

I = zeros(1000,1000,3);
figure(50)
points = nan(nSlices,npairs,4);
for i = 1:nSlices
    disp(i/nSlices)
    I = zeros(1000,1000,3);
    z0 = zmin+(i-1)/nSlices*(zmax-zmin);
    z1 = zmin+i/nSlices*(zmax-zmin);
    for pair = 1:npairs
        sqA = sqpairs{pair,1};
        sqB = sqpairs{pair,2};
        rr = 10;
        if ~isempty(sqA)
            [x,y,z] = symplanecoord(sqA,sp,V);
            idx = find(z > z0 & z < z1);
            if ~isempty(idx)
                points(i,pair,1) = mean(x(idx)); points(i,pair,2) = mean(y(idx));
                x = round((mean(x(idx))-xmin)/(xmax-xmin)*999+1);
                y = round((mean(y(idx))-ymin)/(ymax-ymin)*999+1);
                for r = rr-1:rr+1
                    for a = linspace(0,2*pi,3*pi*r)
                        if x > r+1 && x < 1000-r && y > r+1 && y < 1000-r
                            I(round(x+r*cos(a)),round(y+r*sin(a)),:) = reshape(rgbs(pair,:),[1 1 3]);
                        end
                    end
                end
            end
        end
        if ~isempty(sqB)
            [x,y,z] = symplanecoord(sqB,sp,V);
            idx = find(z > z0 & z < z1);
            if ~isempty(idx)
%                 xmem = x; ymem = y;
                points(i,pair,3) = mean(x(idx)); points(i,pair,4) = mean(y(idx));
                x = round((mean(x(idx))-xmin)/(xmax-xmin)*999+1);
                y = round((mean(y(idx))-ymin)/(ymax-ymin)*999+1);
                for r = rr-1:rr+1
                    for a = linspace(0,2*pi,3*pi*r)
                        if x > r+1 && x < 1000-r && y > r+1 && y < 1000-r
                            I(round(x+r*cos(a)),round(y+r*sin(a)),:) = reshape(rgbs(pair,:),[1 1 3]);
                        end
                    end
                end
%                 x = -xmem; y = ymem;
%                 x = round((mean(x(idx))-xmin)/(xmax-xmin)*999+1);
%                 y = round((mean(y(idx))-ymin)/(ymax-ymin)*999+1);
%                 for r = rr-1:rr+1
%                     for a = linspace(0,2*pi,pi*r)
%                         if x > r+1 && x < 1000-r && y > r+1 && y < 1000-r
%                             I(round(x+r*cos(a)),round(y+r*sin(a)),:) = reshape(rgbs(pair,:),[1 1 3]);
%                         end
%                     end
%                 end
            end
        end
    end
    I(500:501,:,:) = 0.5;
    imshow(imrotate(I,90))
    pause(0.01)
%     imwrite(imrotate(I,90),sprintf('~/Desktop/Slices/Slice%05d.png',i));
end
%save('~/Desktop/points.mat','points')

%% video from frames

%  writerObj = VideoWriter('~/Desktop/Slices.avi');
%  writerObj.FrameRate = 12;
%  open(writerObj);
% for i = 2:215%nSlices
%     disp(i/nSlices)
%     I = imread(sprintf('~/Desktop/Slices/Slice%05d.png',i));
%     writeVideo(writerObj, im2frame(I));
% end
% close(writerObj);

%% analyze slices
%load('~/Desktop/points.mat')
doplot = 1;
if doplot
    fig60 = figure(60);
    scsz = get(0,'ScreenSize'); % scsz = [left botton width height]
    fig60.Position = [scsz(3)/4 scsz(4)/4 scsz(3)/2 scsz(4)/2];
end
frameindex = 0;
stack = [];
stackang = [];
stackdst = [];

% stack = NaN(nSlices,npairs*npairs);
for i = 1:nSlices
    fprintf('%d\n',i);
    t = squeeze(points(i,:,:));
    if sum(isnan(t)) == 0
        clf
        I = zeros(400,400,3);
        
        if doplot
            subplot(1,4,1)
            hold on
        end
        orgpoints = zeros(npairs,2);
        refpoints = zeros(npairs,2);
        for pair = 1:npairs
            x1 = points(i,pair,1);
            y1 = points(i,pair,2);
            x2 = points(i,pair,3);
            y2 = points(i,pair,4);
            if doplot
                plot([x1 -x2],[y1 y2],'o')
            end
            orgpoints(pair,:) = [x1 y1];
            refpoints(pair,:) = [-x2 y2];
            
            x1 = (x1-xmin)/(xmax-xmin)*400;
            y1 = 400-(y1-ymin)/(ymax-ymin)*400;
            I = insertShape(I,'circle',[x1 y1 5],'LineWidth',1,'Color','red');
            x2 = (-x2-xmin)/(xmax-xmin)*400;
            y2 = 400-(y2-ymin)/(ymax-ymin)*400;
            I = insertShape(I,'circle',[x2 y2 5],'LineWidth',1,'Color','green');
        end
        if doplot
            hold off
            axis([xmin xmax ymin ymax])
            title('pairs')
        end
        
        npoints = npairs;
        RD = zeros(npoints,npoints); % relative displacement
        mx = 0;
        for ii = 1:npoints
            for jj = ii:npoints
                p = orgpoints(ii,:);
                q = orgpoints(jj,:);
                r1 = q-p;
                p = refpoints(ii,:);
                q = refpoints(jj,:);
                r2 = q-p;
                n1 = norm(r1);
                n2 = norm(r2);
                nn = abs(n2-n1);
                if n1 > 0 && n2 > 0
                    RD(ii,jj) = (1-dot(r1/n1,r2/n2))/2;
                else
                    RD(ii,jj) = 0;
                end
                if jj > ii
                    RD(jj,ii) = nn;
                    if nn > mx
                        mx = nn;
                    end
                end
            end
        end
        for ii = 1:npoints-1
            for jj = ii+1:npoints
                RD(jj,ii) = RD(jj,ii)/(xmax-xmin);
            end
        end
        
        if doplot
            subplot(1,4,2)
            plot(orgpoints(:,1),orgpoints(:,2),'o')
            axis([xmin xmax ymin ymax])
            title('one side, original')

            subplot(1,4,3)
            plot(refpoints(:,1),refpoints(:,2),'o')
            axis([xmin xmax ymin ymax])
            title('other side, reflected')

            subplot(1,4,4)
            imshow(RD)
            title('correlation')
        end
        
        frameindex = frameindex+1;
%         saveas(gcf,sprintf('~/Desktop/Slices/Slice%05d.png',frameindex));
        I(:,200:201,:) = 0.5;
        J = imresize(RD,[400 400],'nearest');
        J = repmat(J,[1 1 3]);
        W = 0.25*ones(size(I,1),5,3);
        RGB = cat(2,cat(2,I,W),J);
%         imshow(RGB)
        pause(0.01)
%         imwrite(RGB,sprintf('~/Desktop/Slices/Slice%05d.png',frameindex))
%         return
        stack = [stack; reshape(RD,1,[])];
        stackang = [stackang; nonzeros(triu(RD)')'];
        stackdst = [stackdst; nonzeros(tril(RD)')'];
    else
        fprintf('t had NaN value\n');
    end
end
% if doplot
%     close all
% end
fig70 = figure(70);
imshow(flipud(imresize(stack,10,'nearest')))
fig80 = figure(80);
subplot(2,1,1), imshow(flipud(imresize(stackang,10,'nearest')))
subplot(2,1,2), imshow(flipud(imresize(stackdst,10,'nearest')))
