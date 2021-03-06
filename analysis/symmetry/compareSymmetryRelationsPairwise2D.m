clear
clc

%% settings
DataPath = 'D:\Dropbox (Personal)\MATLAB\Data\';
DataFile = '161107t1546_130201zf142_160515SWiFT_ProjOrLngstLtL_ANNOTsymmetry_IGNblacklistsymblack_1umLenThresh_PHYScoord.txt';
% subset of 6 pairs
SubsetFile = '161107t1203_130201zf142_160515SWiFT_SUBSETspinalbackfillsIDENTnucMLFandMauthner.txt';
% subset of 22 pairs
%SubsetFile = '161107t1202_130201zf142_160515SWiFT_SUBSETspinalbackfillsIDENT.txt';

DateString = datestr(now,30);
DateString = strrep(DateString(3:length(DateString)-2),'T','t');
Prefix = strcat(DateString,'_');

% z step size
step = round(1*1000/60); % visualize every 1um

%% --------------------------------------------------
fprintf('loading data\n');
% --------------------------------------------------

% load plane parameters (perpendicular vector and points)
load(strcat(DataPath,filesep,'161107t1647_plane_161107t1546data_161107t1201subsetMLF_ICPnosubsamp.mat'));
load(strcat(DataPath,filesep,'161121t1410_assignment_161107t1546data_161107t1203subsetNucMLFMauth_161107t1647planeMLF_OHpenalty_dtwFreq17.mat'),'asgnm_mr','C');
%load(strcat(DataPath,filesep,'161112t0241_assignment_161107t1546data_161107t1202subsetIDENT_161107t1647planeMLF_OHpenalty_dtwFreq17.mat'),'asgnm_mr','C');

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

Colors = importdata(strcat(DataPath,filesep,SubsetFile));
IskelsColors = [];
for i = 1:length(Colors)
    ss = strsplit(Colors{i},' ');
    IskelsColors = [IskelsColors; [str2double(ss{1}) hex2rgb(ss{2}(2:end))]];
end
Colors = zeros(length(iskels),3);
for i = 1:length(iskels)
    j = find(IskelsColors(:,1) == iskels(i));
    Colors(i,:) = IskelsColors(j,2:4);
end

P = D(1:step:end,4:6);

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

doplot = 0;
if doplot
    fig15 = figure(15);
    scsz = get(0,'ScreenSize'); % scsz = [left botton width height]
    fig15.Position = [scsz(3)/4 scsz(4)/4 scsz(3)/2 scsz(4)/2];
    plot3(P(:,1),P(:,2),P(:,3),'.'), hold on
    fill3(sp(:,1),sp(:,2),sp(:,3),'g'), alpha(0.1), hold off
    grid on, axis equal
    xlabel('x'), ylabel('y'), zlabel('z'), title('plane')
    % saveas(gcf,sprintf('~/Desktop/Planes.png'));
end

%% --------------------------------------------------
fprintf('plot projections\n');
% --------------------------------------------------

[PX,PY,PZ] = symplanecoord(P,sp,V);

doplot = 0;
if doplot
    fig16 = figure(16);
    scsz = get(0,'ScreenSize'); % scsz = [left botton width height]
    fig16.Position = [scsz(3)/4 scsz(4)/4 scsz(3)/2 scsz(4)/2];
    subplot(1,3,1), plot(PX,PY,'.','MarkerSize',1), title('top'), axis equal
    subplot(1,3,2), plot(PX,PZ,'.','MarkerSize',1), title('front'), axis equal
    subplot(1,3,3), plot(PY,PZ,'.','MarkerSize',1), title('side'), axis equal
    % saveas(gcf,sprintf('~/Desktop/Projections.png'));
end

zMin = min(PZ);
zMax = max(PZ);
xMin = min(PX);
xMax = max(PX);
yMin = min(PY);
yMax = max(PY);

%% unique assignment pairs

asgnm = asgnm_mr;
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

% shuffle assignment in two closest dtw pairs
paircosts = zeros(size(asgnm,1),1);
for i = 1:size(asgnm,1)
    paircosts(i) = C(asgnm(i,1),asgnm(i,2));
end
[~,im1] = min(paircosts);
paircosts(im1) = Inf;
[~,im2] = min(paircosts);

disp(asgnm)
assignment = asgnm(:,2);
for i = 1:length(assignment)
    if assignment(i) == 0
        fprintf('%s (%d) matched NONE\n',iskelnames{i},iskels(i))
        continue
    end
    fprintf('%s (%d) matched to %s (%d)\n',...
        iskelnames{i},iskels(i),...
        iskelnames{assignment(i)},iskels(assignment(i)));
end
fprintf('\n')
asgnm_shuf = asgnm;
asgnm_shuf([im1 im2],2) = flipud(asgnm([im1 im2],2));
npairs_shuf = size(asgnm_shuf,1);
disp(asgnm_shuf)
assignment_shuf = asgnm_shuf(:,2);
for i = 1:length(assignment_shuf)
    if assignment_shuf(i) == 0
        fprintf('%s (%d) matched NONE\n',iskelnames{i},iskels(i))
        continue
    end
    fprintf('%s (%d) matched to %s (%d)\n',...
        iskelnames{i},iskels(i),...
        iskelnames{assignment_shuf(i)},iskels(assignment_shuf(i)));
end

% generate colors in case not pulling them directly from subset file
hues = (randperm(npairs)-1)/npairs;%*2/3;
hsvs = [hues' ones(npairs,1) ones(npairs,1)];
rgbs = zeros(size(hsvs));
rgbs_shuf = zeros(size(hsvs));

%% record skeleton pairs (for computational speed)
sqpairs = cell(npairs,2);
for i = 1:npairs
    fprintf('%f\n',i/npairs);
    % rgbs(i,:) = hsv2rgb(hsvs(i,:));
    rgbs(i,:) = Colors(asgnm(i,1),:);
    [sqA,~] = getnodes(D,iskels(asgnm(i,1)));
    [sqB,~] = getnodes(D,iskels(asgnm(i,2)));
    sqpairs{i,1} = sqA;
    sqpairs{i,2} = sqB;
end

sqpairs_shuf = cell(npairs_shuf,2);
for i_shuf = 1:npairs_shuf
    fprintf('%f\n',i_shuf/npairs_shuf);
    rgbs_shuf(i_shuf,:) = Colors(asgnm_shuf(i_shuf,1),:);
    [sqA_shuf,~] = getnodes(D,iskels(asgnm_shuf(i_shuf,1)));
    [sqB_shuf,~] = getnodes(D,iskels(asgnm_shuf(i_shuf,2)));
    sqpairs_shuf{i_shuf,1} = sqA_shuf;
    sqpairs_shuf{i_shuf,2} = sqB_shuf;
end

%% compute/set visualization parameters

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
axmax = max(-xmin,xmax);
aymax = max(-ymin,ymax);
azmax = max(-zmin,zmax);

allX_shuf = []; allY_shuf = []; allZ_shuf = [];
for pair_shuf = 1:npairs_shuf
    xyz_shuf = sqpairs_shuf{pair_shuf,1};
    xyz_shuf = [xyz_shuf; sqpairs_shuf{pair_shuf,2}];
    [x_shuf,y_shuf,z_shuf] = symplanecoord(xyz_shuf,sp,V);
    allX_shuf = [allX_shuf; x_shuf];
    allY_shuf = [allY_shuf; y_shuf];
    allZ_shuf = [allZ_shuf; z_shuf];
end
xmin_shuf = min(allX_shuf); xmax_shuf = max(allX_shuf);
ymin_shuf = min(allY_shuf); ymax_shuf = max(allY_shuf);
zmin_shuf = min(allZ_shuf); zmax_shuf = max(allZ_shuf);
axmax_shuf = max(-xmin_shuf,xmax_shuf);
aymax_shuf = max(-ymin_shuf,ymax_shuf);
azmax_shuf = max(-zmin_shuf,zmax_shuf);

nSlices = floor((zmax-zmin)/60);
%nSlices = round(nSlices/(2*5.25945));

%% display slices in figure (optional)

% fig17 = figure(17);
% scsz = get(0,'ScreenSize'); % scsz = [left botton width height]
% fig17.Position = [scsz(3)/4 scsz(4)/4 scsz(3)/2 scsz(4)/2];
% 
% for i = 1:nSlices
%     % entire set, 3D
%     z0 = zmin+(i-1)/nSlices*(zmax-zmin);
%     z1 = zmin+i/nSlices*(zmax-zmin);
%     idx = find(PZ <= z0 | PZ >= z1);
%     subplot(1,3,1)
%     plot3(PX(idx),PY(idx),PZ(idx),'.r','MarkerSize',1), hold on
%     idx = find(PZ > z0 & PZ < z1);
%     plot3(PX(idx),PY(idx),PZ(idx),'.g','MarkerSize',1), hold off
%     axis([xMin xMax yMin yMax zMin zMax])
%     view(0,0)
%     
%     % subset, 3D
%     subplot(1,3,2)
%     
%     xyz = sqpairs{1,1};
%     xyz = [xyz; sqpairs{1,2}];
%     [x,y,z] = symplanecoord(xyz,sp,V);
%     idx = find(z <= z0 | z >= z1);
%     plot3(x(idx),y(idx),z(idx),'.','MarkerSize',1,'color',rgbs(1,:)), hold on
%     idx = find(z > z0 & z < z1);
%     plot3(x(idx),y(idx),z(idx),'.k','MarkerSize',1)
%     for pair = 2:npairs
%         xyz = sqpairs{pair,1};
%         xyz = [xyz; sqpairs{pair,2}];
%         [x,y,z] = symplanecoord(xyz,sp,V);
%         idx = find(z <= z0 | z >= z1);
%         plot3(x(idx),y(idx),z(idx),'.','MarkerSize',1,'color',rgbs(pair,:)), hold on
%         idx = find(z > z0 & z < z1);
%         plot3(x(idx),y(idx),z(idx),'.k','MarkerSize',1)
%     end
%     hold off
%     axis([xMin xMax yMin yMax zMin zMax])
%     view(0,0)
% 
%     s3 = subplot(1,3,3);
%     cla(s3)
%     hold on
%     for pair = 1:npairs
%         xyz = sqpairs{pair,1};
%         xyz = [xyz; sqpairs{pair,2}];
%         [x,y,z] = symplanecoord(xyz,sp,V);
%         plot(x,y,'.','MarkerSize',1,'color',hsv2rgb([1 0.25 1].*hsvs(pair,:)))
%     end
%     for pair = 1:npairs
%         sqA = sqpairs{pair,1};
%         sqB = sqpairs{pair,2};
%         if ~isempty(sqA)
%             [x,y,z] = symplanecoord(sqA,sp,V);
%             idx = find(z > z0 & z < z1);
%             if ~isempty(idx)
%                 plot(mean(x(idx)),mean(y(idx)),'.','MarkerSize',10,'color',rgbs(pair,:))
%             end
%         end
%         if ~isempty(sqB)
%             [x,y,z] = symplanecoord(sqB,sp,V);
%             idx = find(z > z0 & z < z1);
%             if ~isempty(idx)
%                 plot(mean(x(idx)),mean(y(idx)),'.','MarkerSize',10,'color',rgbs(pair,:))
%             end
%         end
%     end
%     axis([xmin xmax ymin ymax])
%     hold off
%     %saveas(gcf,sprintf('~/Desktop/Slices/Slice%05d.png',i));
%     pause(0.01)
% end

%% display slices on image

dodraw = 0;
if dodraw
    I = zeros(1000,1000,3);
    figure(18)
end
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
            idx = find(z >= z0 & z < z1);
            if ~isempty(idx)
                points(i,pair,1) = mean(x(idx));
                points(i,pair,2) = mean(y(idx));
                x = round((mean(x(idx))-xmin)/(xmax-xmin)*999+1);
                y = round((mean(y(idx))-ymin)/(ymax-ymin)*999+1);
                if dodraw
                    I = insertShape(I,'circle',[y x 10],'LineWidth',1,'Color',rgbs(pair,:));
                end
            end
        end
        if ~isempty(sqB)
            [x,y,z] = symplanecoord(sqB,sp,V);
            idx = find(z >= z0 & z < z1);
            if ~isempty(idx)
                points(i,pair,3) = mean(x(idx));
                points(i,pair,4) = mean(y(idx));
                x = round((mean(x(idx))-xmin)/(xmax-xmin)*999+1);
                y = round((mean(y(idx))-ymin)/(ymax-ymin)*999+1);
                if dodraw
                    I = insertShape(I,'circle',[y x 10],'LineWidth',1,'Color',rgbs(pair,:));
                end
            end
        end
    end
    if dodraw
        I(500:501,:,:) = 0.5;
        imshow(imrotate(I,90))
    end
    %pause(0.01)
end

% display slices on image for shuffle
dodraw = 0;
if dodraw
    I_shuf = zeros(1000,1000,3);
    figure(20)
end
points_shuf = nan(nSlices,npairs_shuf,4);
for i_shuf = 1:nSlices
    disp(i_shuf/nSlices)
    I_shuf = zeros(1000,1000,3);
    z0_shuf = zmin_shuf+(i_shuf-1)/nSlices*(zmax_shuf-zmin_shuf);
    z1_shuf = zmin_shuf+i_shuf/nSlices*(zmax_shuf-zmin_shuf);
    for pair_shuf = 1:npairs_shuf
        sqA_shuf = sqpairs_shuf{pair_shuf,1};
        sqB_shuf = sqpairs_shuf{pair_shuf,2};
        rr_shuf = 10;
        if ~isempty(sqA_shuf)
            [x_shuf,y_shuf,z_shuf] = symplanecoord(sqA_shuf,sp,V);
            idx_shuf = find(z_shuf >= z0_shuf & z_shuf < z1_shuf);
            if ~isempty(idx_shuf)
                points_shuf(i_shuf,pair_shuf,1) = mean(x_shuf(idx_shuf));
                points_shuf(i_shuf,pair_shuf,2) = mean(y_shuf(idx_shuf));
                x_shuf = round((mean(x_shuf(idx_shuf))-xmin_shuf)/(xmax_shuf-xmin_shuf)*999+1);
                y_shuf = round((mean(y_shuf(idx_shuf))-ymin_shuf)/(ymax_shuf-ymin_shuf)*999+1);
                if dodraw
                    I_shuf = insertShape(I_shuf,'circle',[y_shuf x_shuf 10],'LineWidth',1,'Color',rgbs_shuf(pair_shuf,:));
                end
            end
        end
        if ~isempty(sqB_shuf)
            [x_shuf,y_shuf,z_shuf] = symplanecoord(sqB_shuf,sp,V);
            idx_shuf = find(z_shuf >= z0_shuf & z_shuf < z1_shuf);
            if ~isempty(idx_shuf)
                points_shuf(i_shuf,pair_shuf,3) = mean(x_shuf(idx_shuf));
                points_shuf(i_shuf,pair_shuf,4) = mean(y_shuf(idx_shuf));
                x_shuf = round((mean(x_shuf(idx_shuf))-xmin_shuf)/(xmax_shuf-xmin_shuf)*999+1);
                y_shuf = round((mean(y_shuf(idx_shuf))-ymin_shuf)/(ymax_shuf-ymin_shuf)*999+1);
                if dodraw
                    I_shuf = insertShape(I_shuf,'circle',[y_shuf x_shuf 10],'LineWidth',1,'Color',rgbs_shuf(pair_shuf,:));
                end
            end
        end
    end
    if dodraw
        I_shuf(500:501,:,:) = 0.5;
        imshow(imrotate(I_shuf,90))
    end
    %pause(0.01)
end

%% analyze slices

% compute distance normalization
mx = 8000; % use 8um instead of computing
% mx = 0;
% for i = 1:nSlices
%     orgpoints = zeros(npairs,2);
%     refpoints = zeros(npairs,2);
%     for pair = 1:npairs
%         x1 = points(i,pair,1);
%         y1 = points(i,pair,2);
%         x2 = points(i,pair,3);
%         y2 = points(i,pair,4);
%         orgpoints(pair,:) = [x1 y1];
%         refpoints(pair,:) = [-x2 y2];
%     end
%     npoints = npairs;
%     for ii = 1:npoints-1
%         for jj = ii+1:npoints
%             p = orgpoints(ii,:);
%             q = orgpoints(jj,:);
%             r1 = q-p;
%             p = refpoints(ii,:);
%             q = refpoints(jj,:);
%             r2 = q-p;
%             n1 = norm(r1);
%             n2 = norm(r2);
%             nn = abs(n2-n1);
%             if nn > mx
%                 mx = nn;
%             end
%         end
%     end
% end

% compute distance normalization for shuffle
mx_shuf = 8000; % use 8um instead of computing
% mx_shuf = 0;
% for i_shuf = 1:nSlices
%     orgpoints_shuf = zeros(npairs_shuf,2);
%     refpoints_shuf = zeros(npairs_shuf,2);
%     for pair_shuf = 1:npairs_shuf
%         x1_shuf = points_shuf(i_shuf,pair_shuf,1);
%         y1_shuf = points_shuf(i_shuf,pair_shuf,2);
%         x2_shuf = points_shuf(i_shuf,pair_shuf,3);
%         y2_shuf = points_shuf(i_shuf,pair_shuf,4);
%         orgpoints_shuf(pair_shuf,:) = [x1_shuf y1_shuf];
%         refpoints_shuf(pair_shuf,:) = [-x2_shuf y2_shuf];
%     end
%     npoints_shuf = npairs_shuf;
%     for ii_shuf = 1:npoints_shuf-1
%         for jj_shuf = ii_shuf+1:npoints_shuf
%             p_shuf = orgpoints_shuf(ii_shuf,:);
%             q_shuf = orgpoints_shuf(jj_shuf,:);
%             r1_shuf = q_shuf-p_shuf;
%             p_shuf = refpoints_shuf(ii_shuf,:);
%             q_shuf = refpoints_shuf(jj_shuf,:);
%             r2_shuf = q_shuf-p_shuf;
%             n1_shuf = norm(r1_shuf);
%             n2_shuf = norm(r2_shuf);
%             nn_shuf = abs(n2_shuf-n1_shuf);
%             if nn_shuf > mx_shuf
%                 mx_shuf = nn_shuf;
%             end
%         end
%     end
% end

% % skeleton subset map (grayscale)
% xyz = sqpairs{1,1};
% xyz = [xyz; sqpairs{1,2}];
% [xs,~,zs] = symplanecoord(xyz,sp,V);
% for pair = 2:npairs
%     xyz = sqpairs{pair,1};
%     xyz = [xyz; sqpairs{pair,2}];
%     [x,~,z] = symplanecoord(xyz,sp,V);
%     xs = [xs; x];
%     zs = [zs; z];
% end
% nr = 2000; nc = 400;
% map = zeros(nr,nc,3);
% for j = 1:length(xs)
%     x = floor((xs(j)-xmin)/(xmax-xmin)*(nc-1))+1;
%     z = nr-floor((zs(j)-zmin)/(zmax-zmin)*(nr-1));
%     map(z,x,:) = 1;
% end
% map = imrotate(map,-90);
% map = imresize(map,[100 805]);
% %imshow(map)
% 
% % skeleton subset map (grayscale) for shuffle
% xyz_shuf = sqpairs_shuf{1,1};
% xyz_shuf = [xyz_shuf; sqpairs_shuf{1,2}];
% [xs_shuf,~,zs_shuf] = symplanecoord(xyz_shuf,sp,V);
% for pair_shuf = 2:npairs_shuf
%     xyz_shuf = sqpairs_shuf{pair_shuf,1};
%     xyz_shuf = [xyz_shuf; sqpairs_shuf{pair_shuf,2}];
%     [x_shuf,~,z_shuf] = symplanecoord(xyz_shuf,sp,V);
%     xs_shuf = [xs_shuf; x_shuf];
%     zs_shuf = [zs_shuf; z_shuf];
% end
% nr_shuf = 2000; nc_shuf = 400;
% map_shuf = zeros(nr_shuf,nc_shuf,3);
% for j_shuf = 1:length(xs_shuf)
%     x_shuf = floor((xs_shuf(j_shuf)-xmin_shuf)/(xmax_shuf-xmin_shuf)*(nc_shuf-1))+1;
%     z_shuf = nr_shuf-floor((zs_shuf(j_shuf)-zmin_shuf)/(zmax_shuf-zmin_shuf)*(nr_shuf-1));
%     map_shuf(z_shuf,x_shuf,:) = 1;
% end
% map_shuf = imrotate(map_shuf,-90);
% map_shuf = imresize(map_shuf,[100 805]);
% %imshow(map_shuf)

dodraw = 0;
% skeleton subset map (rgb)
% nr = 2000; nc = 400;
nr = 805; nc = 100;
map = zeros(nr,nc,3);
for pair = 1:npairs
    xyz = sqpairs{pair,1};
    [xs,~,zs] = symplanecoord(xyz,sp,V);
    x = (xs-xmin)./(xmax-xmin)*(nc-1)+1;
    z = nr-(zs-zmin)/(zmax-zmin)*(nr-1);
    x = smooth(x(1:10:end)); z = smooth(z(1:10:end));
    map = insertShape(map,'Line',reshape([x z]',1,[]),'Color',rgbs(pair,:),'Opacity',1);
    
    xyz = [sqpairs{pair,2}];
    [xs,~,zs] = symplanecoord(xyz,sp,V);
    x = (xs-xmin)./(xmax-xmin)*(nc-1)+1;
    z = nr-(zs-zmin)/(zmax-zmin)*(nr-1);
    x = smooth(x(1:10:end)); z = smooth(z(1:10:end));
    map = insertShape(map,'Line',reshape([x z]',1,[]),'Color',rgbs(pair,:),'Opacity',1);
end
map = imrotate(map,-90);
% map = imresize(map,[100 805]);
if dodraw
    fig21 = figure(21);
    scsz = get(0,'ScreenSize');
    fig21.Position = [scsz(3)/4 scsz(4)/4 scsz(3)/2 scsz(4)/2];
    imshow(map)
end

% skeleton subset map (rgb) for shuffle
% nr_shuf = 2000; nc_shuf = 400;
nr_shuf = 805; nc_shuf = 100;
map_shuf = zeros(nr_shuf,nc_shuf,3);
for pair_shuf = 1:npairs_shuf
    xyz_shuf = sqpairs_shuf{pair_shuf,1};
    [xs_shuf,~,zs_shuf] = symplanecoord(xyz_shuf,sp,V);
    x_shuf = (xs_shuf-xmin_shuf)./(xmax_shuf-xmin_shuf)*(nc_shuf-1)+1;
    z_shuf = nr_shuf-(zs_shuf-zmin_shuf)/(zmax_shuf-zmin_shuf)*(nr_shuf-1);
    x_shuf = smooth(x_shuf(1:10:end)); z_shuf = smooth(z_shuf(1:10:end));
    map_shuf = insertShape(map_shuf,'Line',reshape([x_shuf z_shuf]',1,[]),'Color',rgbs_shuf(pair_shuf,:),'Opacity',1);
    
    xyz_shuf = [sqpairs_shuf{pair_shuf,2}];
    [xs_shuf,~,zs_shuf] = symplanecoord(xyz_shuf,sp,V);
    x_shuf = (xs_shuf-xmin_shuf)./(xmax_shuf-xmin_shuf)*(nc_shuf-1)+1;
    z_shuf = nr_shuf-(zs_shuf-zmin_shuf)/(zmax_shuf-zmin_shuf)*(nr_shuf-1);
    x_shuf = smooth(x_shuf(1:10:end)); z_shuf = smooth(z_shuf(1:10:end));
    map_shuf = insertShape(map_shuf,'Line',reshape([x_shuf z_shuf]',1,[]),'Color',rgbs_shuf(pair_shuf,:),'Opacity',1);
end
map_shuf = imrotate(map_shuf,-90);
% map_shuf = imresize(map_shuf,[100 805]);
if dodraw
    fig22 = figure(22);
    scsz = get(0,'ScreenSize');
    fig22.Position = [scsz(3)/4 scsz(4)/4 scsz(3)/2 scsz(4)/2];
    imshow(map_shuf)
end

dodraw = 0;
doshow = 0;
dosave = 0;

% visualize
nStackRows = (npairs)*(npairs-1);
stack = zeros(nStackRows,nSlices,3);
stackNaN = nan(nStackRows,nSlices);
stackNaNum = nan(nStackRows,nSlices);
nStackRows_shuf = (npairs_shuf)*(npairs_shuf-1);
stack_shuf = zeros(nStackRows_shuf,nSlices,3);
stackNaN_shuf = nan(nStackRows_shuf,nSlices);
stackNaNum_shuf = nan(nStackRows_shuf,nSlices);
RDum_all = nan(npairs,npairs,nSlices);
RD_all = nan(npairs,npairs,nSlices);
RDum_shuf_all = nan(npairs,npairs,nSlices);
RD_shuf_all = nan(npairs,npairs,nSlices);
frameindex = 0;
a1122 = [];
a1122_shuf = [];
for ni = 1:nSlices
    i = nSlices-ni+1;
    disp(i/nSlices)
    
    z0 = zmin+(i-1)/nSlices*(zmax-zmin);
    z1 = zmin+i/nSlices*(zmax-zmin);
    z0_shuf = zmin_shuf+(i-1)/nSlices*(zmax_shuf-zmin_shuf);
    z1_shuf = zmin_shuf+i/nSlices*(zmax_shuf-zmin_shuf);
    
    if dodraw
        I = zeros(400,400,3);
        I(:,200:201,:) = 0.5;
        I_shuf = zeros(400,400,3);
        I_shuf(:,200:201,:) = 0.5;
    end
    
    orgpoints = zeros(npairs,2);
    refpoints = zeros(npairs,2);
    for pair = 1:npairs
        x1 = points(i,pair,1);
        y1 = points(i,pair,2);
        x2 = points(i,pair,3);
        y2 = points(i,pair,4);
        orgpoints(pair,:) = [x1 y1];
        refpoints(pair,:) = [-x2 y2];

        if dodraw
            if ~isnan(x1) && ~isnan(y1)
                %x1 = (x1-xmin)/(xmax-xmin)*400;
                x1 = 200+x1/axmax*199;
                y1 = 400-(y1-ymin)/(ymax-ymin)*400;
                I = insertShape(I,'circle',[x1 y1 5],'LineWidth',2,'Color',rgbs(pair,:));
            end
            if ~isnan(x2) && ~isnan(y2)
                %x2 = (x2-xmin)/(xmax-xmin)*400;
                x2 = 200+x2/axmax*199;
                y2 = 400-(y2-ymin)/(ymax-ymin)*400;
                I = insertShape(I,'circle',[x2 y2 5],'LineWidth',2,'Color',rgbs(pair,:));
            end
            I = insertText(I, [2 1], sprintf('z = %06.2f um',(mean([z0 z1])-zmin)/1000),'TextColor','white','BoxOpacity',0.0);
            I = insertText(I, [170 1], 'left','TextColor','white','BoxOpacity',0.0);
            I = insertText(I, [208 1], 'right','TextColor','white','BoxOpacity',0.0);
        end
    end
    
    orgpoints_shuf = zeros(npairs_shuf,2);
    refpoints_shuf = zeros(npairs_shuf,2);
    for pair_shuf = 1:npairs_shuf
        x1_shuf = points_shuf(i,pair_shuf,1);
        y1_shuf = points_shuf(i,pair_shuf,2);
        x2_shuf = points_shuf(i,pair_shuf,3);
        y2_shuf = points_shuf(i,pair_shuf,4);
        orgpoints_shuf(pair_shuf,:) = [x1_shuf y1_shuf];
        refpoints_shuf(pair_shuf,:) = [-x2_shuf y2_shuf];

        if dodraw
            if ~isnan(x1_shuf) && ~isnan(y1_shuf)
                %x1_shuf = (x1_shuf-xmin_shuf)/(xmax_shuf-xmin_shuf)*400;
                x1_shuf = 200+x1_shuf/axmax_shuf*199;
                y1_shuf = 400-(y1_shuf-ymin_shuf)/(ymax_shuf-ymin_shuf)*400;
                I_shuf = insertShape(I_shuf,'circle',[x1_shuf y1_shuf 5],'LineWidth',2,'Color',rgbs_shuf(pair_shuf,:));
            end
            if ~isnan(x2_shuf) && ~isnan(y2_shuf)
                %x2_shuf = (x2_shuf-xmin_shuf)/(xmax_shuf-xmin_shuf)*400;
                x2_shuf = 200+x2_shuf/axmax_shuf*199;
                y2_shuf = 400-(y2_shuf-ymin_shuf)/(ymax_shuf-ymin_shuf)*400;
                I_shuf = insertShape(I_shuf,'circle',[x2_shuf y2_shuf 5],'LineWidth',2,'Color',rgbs_shuf(pair_shuf,:));
            end
            I_shuf = insertText(I_shuf, [2 1], sprintf('z = %06.2f um',(mean([z0_shuf z1_shuf])-zmin_shuf)/1000),'TextColor','white','BoxOpacity',0.0);
            I_shuf = insertText(I_shuf, [170 1], 'left','TextColor','white','BoxOpacity',0.0);
            I_shuf = insertText(I_shuf, [208 1], 'right','TextColor','white','BoxOpacity',0.0);
        end
    end

    npoints = npairs;
    RDum = zeros(npoints,npoints); % real displacement in nm
    RD = zeros(npoints,npoints); % relative displacement
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
                if ii == 1 && jj == 2
                    v12 = r1/n1;
                    a1122 = [a1122; atan2(v12(2),v12(1))];
                end
            else
                RD(ii,jj) = NaN;
            end
            if jj > ii
                RD(jj,ii) = nn;
            end
        end
    end
    for ii = 1:npoints-1
        for jj = ii+1:npoints
            RDum(jj,ii) = RD(jj,ii);
            RD(jj,ii) = RD(jj,ii)/mx;
        end
    end

    npoints_shuf = npairs_shuf;
    RDum_shuf = zeros(npoints_shuf,npoints_shuf); % real displacement in nm
    RD_shuf = zeros(npoints_shuf,npoints_shuf); % relative displacement
    for ii_shuf = 1:npoints_shuf
        for jj_shuf = ii_shuf:npoints_shuf
            p_shuf = orgpoints_shuf(ii_shuf,:);
            q_shuf = orgpoints_shuf(jj_shuf,:);
            r1_shuf = q_shuf-p_shuf;
            p_shuf = refpoints_shuf(ii_shuf,:);
            q_shuf = refpoints_shuf(jj_shuf,:);
            r2_shuf = q_shuf-p_shuf;
            n1_shuf = norm(r1_shuf);
            n2_shuf = norm(r2_shuf);
            nn_shuf = abs(n2_shuf-n1_shuf);
            if n1_shuf > 0 && n2_shuf > 0
                RD_shuf(ii_shuf,jj_shuf) = (1-dot(r1_shuf/n1_shuf,r2_shuf/n2_shuf))/2;
                if ii_shuf == 1 && jj_shuf == 2
                    v12_shuf = r1_shuf/n1_shuf;
                    a1122_shuf = [a1122_shuf; atan2(v12_shuf(2),v12_shuf(1))];
                end
            else
                RD_shuf(ii_shuf,jj_shuf) = NaN;
            end
            if jj_shuf > ii_shuf
                RD_shuf(jj_shuf,ii_shuf) = nn_shuf;
            end
        end
    end
    for ii_shuf = 1:npoints_shuf-1
        for jj_shuf = ii_shuf+1:npoints_shuf
            RDum_shuf(jj_shuf,ii_shuf) = RD_shuf(jj_shuf,ii_shuf);
            RD_shuf(jj_shuf,ii_shuf) = RD_shuf(jj_shuf,ii_shuf)/mx_shuf;
        end
    end
    
    frameindex = frameindex+1;

    % colormap for heatmaps
    colormap('winter')
    WinterMap = colormap;
    
    rdum = [];
    rd = [];
    for ii = 1:size(RD,1)-1
        for jj = ii+1:size(RD,2)
            rdum = [rdum RDum(ii,jj)];
            rd = [rd RD(ii,jj)];
        end
    end
    % rd = nonzeros(triu(RD));
    for ii = 1:size(RD,1)-1
        for jj = ii+1:size(RD,2)
            rdum = [rdum RDum(jj,ii)];
            rd = [rd RD(jj,ii)];
        end
    end
    % rd = [rd; nonzeros(tril(RD))];
    for ii = 1:length(rd)
        if isnan(rd(ii))
            stack(ii,i,1) = 0.0;
        else
            %stack(ii,i,:) = rd(ii);
            graylevel = rd(ii);
            if graylevel > 1 % saturate if needed based on manually set mx
                graylevel = 1;
            end
            mapindex = round(graylevel*(size(WinterMap,1)-1))+1;
            mapcolor = WinterMap(mapindex,:);
            stack(ii,i,:) = reshape(mapcolor,[1 1 3]);
        end
    end
    stackNaNum(:,i) = rdum';
    RDum_all(:,:,i) = RDum;
    stackNaN(:,i) = rd';
    RD_all(:,:,i) = RD;

    rd_shuf = [];
    rdum_shuf = [];
    for ii_shuf = 1:size(RD_shuf,1)-1
        for jj_shuf = ii_shuf+1:size(RD_shuf,2)
            rdum_shuf = [rdum_shuf RDum_shuf(ii_shuf,jj_shuf)];
            rd_shuf = [rd_shuf RD_shuf(ii_shuf,jj_shuf)];
        end
    end
    % rd = nonzeros(triu(RD));
    for ii_shuf = 1:size(RD_shuf,1)-1
        for jj_shuf = ii_shuf+1:size(RD_shuf,2)
            rdum_shuf = [rdum_shuf RDum_shuf(jj_shuf,ii_shuf)];
            rd_shuf = [rd_shuf RD_shuf(jj_shuf,ii_shuf)];
        end
    end
    % rd = [rd; nonzeros(tril(RD))];
    for ii_shuf = 1:length(rd_shuf)
        if isnan(rd_shuf(ii_shuf))
            stack_shuf(ii_shuf,i,1) = 0.0;
        else
            %stack(ii,i,:) = rd(ii);
            graylevel_shuf = rd_shuf(ii_shuf);
            if graylevel_shuf > 1 % saturate if needed based on manually set mx
                graylevel_shuf = 1;
            end
            mapindex_shuf = round(graylevel_shuf*(size(WinterMap,1)-1))+1;
            mapcolor_shuf = WinterMap(mapindex_shuf,:);
            stack_shuf(ii_shuf,i,:) = reshape(mapcolor_shuf,[1 1 3]);
        end
    end
    stackNaNum_shuf(:,i) = rdum_shuf';
    RDum_shuf_all(:,:,i) = RDum_shuf;
    stackNaN_shuf(:,i) = rd_shuf';
    RD_shuf_all(:,:,i) = RD_shuf;
    
    if dodraw
        if doshow
            fig23 = figure(23);
            scsz = get(0,'ScreenSize');
            fig23.Position = [scsz(3)/4 scsz(4)/4 scsz(3)/2 scsz(4)/2];
        end
        RD2 = zeros(4*npoints,4*npoints,3);
        for ii = 1:npoints
            for jj = 1:npoints
                i0 = 4*(ii-1);
                j0 = 4*(jj-1);
                if ~isnan(RD(ii,jj))
                    %RD2(i0+1:i0+4,j0+1:j0+4,:) = RD(ii,jj);
                    graylevel = RD(ii,jj);
                    if graylevel > 1 % saturate if needed based on manually set mx
                        graylevel = 1;
                    end
                    mapindex = round(graylevel*(size(WinterMap,1)-1))+1;
                    mapcolor = WinterMap(mapindex,:);
                    RD2(i0+1:i0+4,j0+1:j0+4,:) = repmat(reshape(mapcolor,[1 1 3]),[4 4]);
                end
            end
        end
        J = imresize(RD2,[400 400],'nearest');
        resizeFactorR = size(J,1)/size(RD,1);
        resizeFactorC = size(J,2)/size(RD,2);
        for ii = 1:npoints
            for jj = 1:npoints
                i0 = round((ii-1)*resizeFactorR+0.5*size(J,1)/npoints);
                j00 = round((jj-1)*resizeFactorC+0.38*size(J,2)/npoints);
                j01 = round((jj-1)*resizeFactorC+0.62*size(J,2)/npoints);
                if ~isnan(RD(ii,jj))
                    RD2(i0+1:i0+4,j0+1:j0+4,:) =  RD(ii,jj);
                    if ii < jj
                        J = insertShape(J,'FilledCircle',[j00 i0 8],'LineWidth',2,'Color','white','Opacity',1);
                        J = insertShape(J,'FilledCircle',[j01 i0 8],'LineWidth',2,'Color','white','Opacity',1);
                        J = insertShape(J,'FilledCircle',[j00 i0 5],'LineWidth',2,'Color',rgbs(ii,:),'Opacity',1);
                        J = insertShape(J,'FilledCircle',[j01 i0 5],'LineWidth',2,'Color',rgbs(jj,:),'Opacity',1);
                    else
                        J = insertShape(J,'FilledCircle',[j00 i0 8],'LineWidth',2,'Color','white','Opacity',1);
                        J = insertShape(J,'FilledCircle',[j01 i0 8],'LineWidth',2,'Color','white','Opacity',1);
                        J = insertShape(J,'FilledCircle',[j00 i0 5],'LineWidth',2,'Color',rgbs(jj,:),'Opacity',1);
                        J = insertShape(J,'FilledCircle',[j01 i0 5],'LineWidth',2,'Color',rgbs(ii,:),'Opacity',1);
                    end
                end
            end
        end
        J = insertText(J, [4 380], 'distance difference','TextColor','white','BoxOpacity',0.0);
        J = insertText(J, [298 1], 'angle difference','TextColor','white','BoxOpacity',0.0);
        for ii = 1:400
            J(ii,ii,:) = 0.25;
        end
        W = 0.25*ones(size(I,1),5,3);
        RGB = cat(2,cat(2,I,W),J);
        stackTop = stack(1:size(stack,1)/2,:,:);
        stackBot = stack(size(stack,1)/2+1:end,:,:);
        stack2 = imresize(stackTop,[50 805],'nearest');
        stack2 = insertText(stack2, [2 1], sprintf('aggregate angle\ndifference'),'TextColor','white','BoxOpacity',0.0);
        W = 0.25*ones(5,size(RGB,2),3);
        RGB = cat(1,RGB,W);
        RGB = cat(1,RGB,stack2);
        stack2 = imresize(stackBot,[50 805],'nearest');
        stack2 = insertText(stack2, [2 1], sprintf('aggregate distance\ndifference'),'TextColor','white','BoxOpacity',0.0);
        W = 0.25*ones(1,size(RGB,2),3);
        RGB = cat(1,RGB,W);
        RGB = cat(1,RGB,stack2);
        W = 0.25*ones(5,size(RGB,2),3);
        RGB = cat(1,RGB,W);
        z0map = floor((z0-zmin)/(zmax-zmin)*(size(RGB,2)-1))+1;
        z1map = floor((z1-zmin)/(zmax-zmin)*(size(RGB,2)-1))+1;
        map2 = map;
        map0 = rgb2gray(map) > 0;
        map2(:,z0map:z1map,:) = 0;
        map2(:,z0map:z1map,:) = repmat(map0(:,z0map:z1map),[1 1 3]);
        map2 = insertText(map2, [2 3], 'location','TextColor','white','BoxOpacity',0.0);
        RGB = cat(1,RGB,map2);
        if doshow
            pause(0.01)
            imshow(RGB)
        end
        if dosave
            %mkdir(strcat(DataPath,filesep,Prefix,'video_900slices_rtol_xcent'))
            imwrite(RGB,strcat(DataPath,filesep,Prefix,'video_900slices_rtol_xcent',filesep,...
                sprintf('idx%05d.png',frameindex)))
        end
        
        if doshow
            fig24 = figure(24);
            scsz = get(0,'ScreenSize');
            fig24.Position = [scsz(3)/4 scsz(4)/4 scsz(3)/2 scsz(4)/2];
        end
        RD2_shuf = zeros(4*npoints_shuf,4*npoints_shuf,3);
        for ii_shuf = 1:npoints_shuf
            for jj_shuf = 1:npoints_shuf
                i0_shuf = 4*(ii_shuf-1);
                j0_shuf = 4*(jj_shuf-1);
                if ~isnan(RD_shuf(ii_shuf,jj_shuf))
                    %RD2_shuf(i0_shuf+1:i0_shuf+4,j0_shuf+1:j0_shuf+4,:) = RD_shuf(ii_shuf,jj_shuf);
                    graylevel_shuf = RD_shuf(ii_shuf,jj_shuf);
                    if graylevel_shuf > 1 % saturate if needed based on manually set mx
                        graylevel_shuf = 1;
                    end
                    mapindex_shuf = round(graylevel_shuf*(size(WinterMap,1)-1))+1;
                    mapcolor_shuf = WinterMap(mapindex_shuf,:);
                    RD2_shuf(i0_shuf+1:i0_shuf+4,j0_shuf+1:j0_shuf+4,:) = repmat(reshape(mapcolor_shuf,[1 1 3]),[4 4]);
                end
            end
        end
        J_shuf = imresize(RD2_shuf,[400 400],'nearest');
        resizeFactorR_shuf = size(J_shuf,1)/size(RD_shuf,1);
        resizeFactorC_shuf = size(J_shuf,2)/size(RD_shuf,2);
        for ii_shuf = 1:npoints_shuf
            for jj_shuf = 1:npoints_shuf
                i0_shuf = round((ii_shuf-1)*resizeFactorR_shuf+0.5*size(J_shuf,1)/npoints_shuf);
                j00_shuf = round((jj_shuf-1)*resizeFactorC_shuf+0.38*size(J_shuf,2)/npoints_shuf);
                j01_shuf = round((jj_shuf-1)*resizeFactorC_shuf+0.62*size(J_shuf,2)/npoints_shuf);
                if ~isnan(RD_shuf(ii_shuf,jj_shuf))
                    RD2_shuf(i0_shuf+1:i0_shuf+4,j0_shuf+1:j0_shuf+4,:) = RD(ii_shuf,jj_shuf);
                    if ii_shuf < jj_shuf
                        J_shuf = insertShape(J_shuf,'FilledCircle',[j00_shuf i0_shuf 8],'LineWidth',2,'Color','white','Opacity',1);
                        J_shuf = insertShape(J_shuf,'FilledCircle',[j01_shuf i0_shuf 8],'LineWidth',2,'Color','white','Opacity',1);
                        J_shuf = insertShape(J_shuf,'FilledCircle',[j00_shuf i0_shuf 5],'LineWidth',2,'Color',rgbs_shuf(ii_shuf,:),'Opacity',1);
                        J_shuf = insertShape(J_shuf,'FilledCircle',[j01_shuf i0_shuf 5],'LineWidth',2,'Color',rgbs_shuf(jj_shuf,:),'Opacity',1);
                    else
                        J_shuf = insertShape(J_shuf,'FilledCircle',[j00_shuf i0_shuf 8],'LineWidth',2,'Color','white','Opacity',1);
                        J_shuf = insertShape(J_shuf,'FilledCircle',[j01_shuf i0_shuf 8],'LineWidth',2,'Color','white','Opacity',1);
                        J_shuf = insertShape(J_shuf,'FilledCircle',[j00_shuf i0_shuf 5],'LineWidth',2,'Color',rgbs_shuf(jj_shuf,:),'Opacity',1);
                        J_shuf = insertShape(J_shuf,'FilledCircle',[j01_shuf i0_shuf 5],'LineWidth',2,'Color',rgbs_shuf(ii_shuf,:),'Opacity',1);
                    end
                end
            end
        end
        J_shuf = insertText(J_shuf, [4 380], 'distance difference','TextColor','white','BoxOpacity',0.0);
        J_shuf = insertText(J_shuf, [298 1], 'angle difference','TextColor','white','BoxOpacity',0.0);
        for ii_shuf = 1:400
            J_shuf(ii_shuf,ii_shuf,:) = 0.25;
        end
        W_shuf = 0.25*ones(size(I_shuf,1),5,3);
        RGB_shuf = cat(2,cat(2,I_shuf,W_shuf),J_shuf);
        stackTop_shuf = stack_shuf(1:size(stack_shuf,1)/2,:,:);
        stackBot_shuf = stack_shuf(size(stack_shuf,1)/2+1:end,:,:);
        stack2_shuf = imresize(stackTop_shuf,[50 805],'bilinear');
        stack2_shuf = insertText(stack2_shuf, [2 1], sprintf('aggregate angle\ndifference'),'TextColor','white','BoxOpacity',0.0);
        W_shuf = 0.25*ones(5,size(RGB_shuf,2),3);
        RGB_shuf = cat(1,RGB_shuf,W_shuf);
        RGB_shuf = cat(1,RGB_shuf,stack2_shuf);
        stack2_shuf = imresize(stackBot_shuf,[50 805],'nearest');
        stack2_shuf = insertText(stack2_shuf, [2 1], sprintf('aggregate distance\ndifference'),'TextColor','white','BoxOpacity',0.0);
        W_shuf = 0.25*ones(1,size(RGB_shuf,2),3);
        RGB_shuf = cat(1,RGB_shuf,W_shuf);
        RGB_shuf = cat(1,RGB_shuf,stack2_shuf);
        W_shuf = 0.25*ones(5,size(RGB_shuf,2),3);
        RGB_shuf = cat(1,RGB_shuf,W_shuf);
        z0map_shuf = floor((z0_shuf-zmin_shuf)/(zmax_shuf-zmin_shuf)*(size(RGB_shuf,2)-1))+1;
        z1map_shuf = floor((z1_shuf-zmin_shuf)/(zmax_shuf-zmin_shuf)*(size(RGB_shuf,2)-1))+1;
        map2_shuf = map_shuf;
        map0_shuf = rgb2gray(map_shuf) > 0;
        map2_shuf(:,z0map_shuf:z1map_shuf,:) = 0;
        map2_shuf(:,z0map_shuf:z1map_shuf,:) = repmat(map0_shuf(:,z0map_shuf:z1map_shuf),[1 1 3]);
        map2_shuf = insertText(map2_shuf, [2 3], 'location','TextColor','white','BoxOpacity',0.0);
        RGB_shuf = cat(1,RGB_shuf,map2_shuf);
        if doshow
            imshow(RGB_shuf)
            pause(0.01)
        end
        if dosave
            %mkdir(strcat(DataPath,filesep,Prefix,'video_900slices_rtol_xcent_shuf'))
            imwrite(RGB,strcat(DataPath,filesep,Prefix,'video_900slices_rtol_xcent_shuf',filesep,...
                sprintf('idx%05d.png',frameindex)))
        end
    end
end

stackang = stackNaN(1:nStackRows/2,:);
stackdst = stackNaNum(nStackRows/2+1:end,:);

stackang_shuf = stackNaN_shuf(1:nStackRows_shuf/2,:);
stackdst_shuf = stackNaNum_shuf(nStackRows_shuf/2+1:end,:);

%% consider absolute rather than relative relations within bundle
% angle analysis

scsz = get(0,'ScreenSize'); % scsz = [left botton width height]
figure('Position',scsz)

pair1 = 1;
pair2 = 5;

x1 = nan(1,nSlices);
y1 = nan(1,nSlices);
for slice = 1:nSlices
    x1(slice) = points(slice,pair1,1); % x of first point of pair
    y1(slice) = points(slice,pair1,2); % y of first point of pair
end
x2 = nan(1,nSlices);
y2 = nan(1,nSlices);
for slice = 1:nSlices
    x2(slice) = points(slice,pair2,1); % x of first point of pair
    y2(slice) = points(slice,pair2,2); % y of first point of pair
end

a = nan(1,nSlices);
for slice = 1:nSlices
    v1 = [x1(slice) y1(slice)];
    v2 = [x2(slice) y2(slice)];
    v = v2-v1;
    v = v/norm(v);
    a(slice) = atan2(v(2),v(1));
end

I = zeros(1000,1000,3);

for i = 1:nSlices
    disp(i/nSlices)
    I = zeros(1000,1000,3);
    z0 = zmin+(i-1)/nSlices*(zmax-zmin);
    z1 = zmin+i/nSlices*(zmax-zmin);
    for pair = [pair1 pair2]
        sqA = sqpairs{pair,1};
        sqB = sqpairs{pair,2};
        rr = 10;
        if ~isempty(sqA)
            [x,y,z] = symplanecoord(sqA,sp,V);
            idx = find(z >= z0 & z < z1);
            if ~isempty(idx)
                x = round((mean(x(idx))-xmin)/(xmax-xmin)*999+1);
                y = round((mean(y(idx))-ymin)/(ymax-ymin)*999+1);
                I = insertShape(I,'circle',[y x 20],'LineWidth',2,'Color',rgbs(pair,:));
            end
        end
    end
    
    subplot(1,2,1)
    plot(1:nSlices,a), hold on
    plot(i,a(i),'o'), hold off
    axis([1 nSlices -pi pi])
    title('angle')
    
    subplot(1,2,2)
    imshow(imrotate(I,90))
    
    pause(0.1)
end

%% generate summary heatmaps and plots

% fig70 = figure(70); clf;
% subplot(2,1,1), imshow(flipud(imresize(stackang,10,'nearest')))
% subplot(2,1,2), imshow(flipud(imresize(stackdst,10,'nearest')))

% heatmap angle difference
fig71 = figure(71); clf;
fig71.Position = [0 0 size(stackang,2)+200 size(stackang,1)+500];
heatmapper(stackang,[],[],[],'ColorBar',1,'UseLogColormap',false,...
    'Colormap',winter,'MaxColorValue',1,'MinColorValue',0,...
    'NaNColor',[0 0 0]);
ax = gca;
ax.TickLength = [0 0];
set(gcf,'PaperPositionMode','auto')
% saveas(gcf,strcat(DataPath,filesep,Prefix,'2D_SUBSETallMLF_AngleDiff'),'epsc');
% saveas(gcf,strcat(DataPath,filesep,Prefix,'2D_SUBSETallMLF_AngleDiff'),'svg');
% saveas(gcf,strcat(DataPath,filesep,Prefix,'2D_SUBSETallMLF_AngleDiff'),'png');

% heatmap shuffle angle difference
fig72 = figure(72); clf;
fig72.Position = [0 0 size(stackang_shuf,2)+200 size(stackang_shuf,1)+500];
heatmapper(stackang_shuf,[],[],[],'ColorBar',1,'UseLogColormap',false,...
    'Colormap',winter,'MaxColorValue',1,'MinColorValue',0,...
    'NaNColor',[0 0 0]);
ax_shuf = gca;
ax_shuf.TickLength = [0 0];
set(gcf,'PaperPositionMode','auto')
% saveas(gcf,strcat(DataPath,filesep,Prefix,'PairwiseCosts_Ang_Shuf'),'epsc');
% saveas(gcf,strcat(DataPath,filesep,Prefix,'PairwiseCosts_Ang_Shuf'),'svg');

% heatmap distance difference
fig73 = figure(73); clf;
fig73.Position = [0 0 size(stackdst,2)+200 size(stackdst,1)+500];
%stackdst2 = imresize(stackdst,[100 (7.20472441 * 300)],'nearest');
heatmapper(stackdst,[],[],{},'ColorBar',1,'UseLogColormap',false,...
    'Colormap',winter,'MaxColorValue',8000,'MinColorValue',0,...
    'NaNColor',[0 0 0]);
ax = gca;
ax.TickLength = [0 0];
set(gcf,'PaperPositionMode','auto')
% saveas(gcf,strcat(DataPath,filesep,Prefix,'2D_SUBSETallMLF_DistDiff'),'epsc');
% saveas(gcf,strcat(DataPath,filesep,Prefix,'2D_SUBSETallMLF_DistDiff'),'svg');
% saveas(gcf,strcat(DataPath,filesep,Prefix,'2D_SUBSETallMLF_DistDiff'),'png');

% heatmap shuffle distance difference
fig74 = figure(74); clf;
fig74.Position = [0 0 size(stackdst_shuf,2)+200 size(stackdst_shuf,1)+500];
%stackdst2 = imresize(stackdst,[100 (7.20472441 * 300)],'nearest');
heatmapper(stackdst_shuf,[],[],{},'ColorBar',1,'UseLogColormap',false,...
    'Colormap',winter,'MaxColorValue',8000,'MinColorValue',0,...
    'NaNColor',[0 0 0]);
ax_shuf = gca;
ax_shuf.TickLength = [0 0];
set(gcf,'PaperPositionMode','auto')
% saveas(gcf,strcat(DataPath,filesep,Prefix,'PairwiseCosts_Dist_Shuf'),'epsc');
% saveas(gcf,strcat(DataPath,filesep,Prefix,'PairwiseCosts_Dist_Shuf'),'svg');

% projection map
fig75 = figure(75); clf;
fig75.Position = [0 0 3*size(stackdst,2)+400 3*size(stackdst,1)+500];
xyz = sqpairs{1,1};
xyz = [xyz; sqpairs{1,2}];
[z,y,x] = symplanecoord(xyz,sp,V);
idx = find(z <= z0 | z >= z1);
plot3(x(idx),y(idx),z(idx),'.','MarkerSize',1,'color',rgbs(1,:)), hold on
idx = find(z > z0 & z < z1);
plot3(x(idx),y(idx),z(idx),'.k','MarkerSize',1)
for pair = 2:npairs
    xyz = sqpairs{pair,1};
    xyz = [xyz; sqpairs{pair,2}];
    [z,y,x] = symplanecoord(xyz,sp,V);
    idx = find(z <= z0 | z >= z1);
    plot3(x(idx),y(idx),z(idx),'.','MarkerSize',1,'color',rgbs(pair,:)), hold on
    idx = find(z > z0 & z < z1);
    plot3(x(idx),y(idx),z(idx),'.k','MarkerSize',1)
end
hold off
box off
axis on
axis equal
% set(gcf,'PaperPositionMode','auto')
view(0,180)
set(gca,'YTickLabel',{'-20000','0','20000','40000'})
set(gca,'YTickLabel',{'-250000','-200000','-150000','-100000','-50000','0','50000','100000','150000','200000','250000'})
set(gca,'XGrid','off','YGrid','off','TickDir','out','TickLength',[0.01,0.01])
set(gca,'Color',[240/255 240/255 240/255])
hold off
% saveas(gcf,strcat(DataPath,filesep,Prefix,'2D_SUBSETallMLF_SubsetProj'),'epsc');
% saveas(gcf,strcat(DataPath,filesep,Prefix,'2D_SUBSETallMLF_SubsetProj'),'svg');
% print(strcat(DataPath,filesep,Prefix,'2D_SUBSETidentMLF_SubsetProj'),'-dpng','-r1200');
% plot2svg(strcat(DataPath,filesep,Prefix,'2D_SUBSETidentMLF_SubsetProjB.svg'),fig75,'png');
% saveas(gcf,strcat(DataPath,filesep,Prefix,'2D_SUBSETallMLF_SubsetProj'),'png');

% shuffle projection map
fig76 = figure(76); clf;
fig76.Position = [0 0 size(stackdst,2)+400 size(stackdst,1)+500];
xyz_shuf = sqpairs_shuf{1,1};
xyz_shuf = [xyz_shuf; sqpairs_shuf{1,2}];
[z_shuf,y_shuf,x_shuf] = symplanecoord(xyz_shuf,sp,V);
idx_shuf = find(z_shuf <= z0_shuf | z_shuf >= z1_shuf);
plot3(x_shuf(idx_shuf),y_shuf(idx_shuf),z_shuf(idx_shuf),'.','MarkerSize',1,'color',rgbs_shuf(1,:)), hold on
idx_shuf = find(z_shuf > z0_shuf & z_shuf < z1_shuf);
plot3(x_shuf(idx_shuf),y_shuf(idx_shuf),z_shuf(idx_shuf),'.k','MarkerSize',1)
for pair_shuf = 2:npairs_shuf
    xyz_shuf = sqpairs_shuf{pair_shuf,1};
    xyz_shuf = [xyz_shuf; sqpairs_shuf{pair_shuf,2}];
    [z_shuf,y_shuf,x_shuf] = symplanecoord(xyz_shuf,sp,V);
    idx_shuf = find(z_shuf <= z0_shuf | z_shuf >= z1_shuf);
    plot3(x_shuf(idx_shuf),y_shuf(idx_shuf),z_shuf(idx_shuf),'.','MarkerSize',1,'color',rgbs_shuf(pair_shuf,:)), hold on
    idx_shuf = find(z_shuf > z0_shuf & z_shuf < z1_shuf);
    plot3(x_shuf(idx_shuf),y_shuf(idx_shuf),z_shuf(idx_shuf),'.k','MarkerSize',1)
end
hold off
box off
axis equal
set(gcf,'PaperPositionMode','auto')
view(0,180)
set(gca,'YTickLabel',{'-20000','0','20000','40000'})
set(gca,'XGrid','on','YGrid','on','TickDir','out')
% saveas(gcf,strcat(DataPath,filesep,Prefix,'PairwiseCosts_SubsetProj_Shuf'),'epsc');
% saveas(gcf,strcat(DataPath,filesep,Prefix,'PairwiseCosts_SubsetProj_Shuf'),'svg');

% plot of summed normalzied angle and normalized dist differences
stackNaN_nsum = nansum(stackNaN,1);
stackNaN_sum = smooth(stackNaN_nsum,0.005,'rloess'); %stackNaN_nsum;
stackNaN_shuf_nsum = nansum(stackNaN_shuf,1);
stackNaN_shuf_sum = smooth(stackNaN_shuf_nsum,0.005,'rloess'); %stackNaN_shuf_nsum;
Xax = 1:size(stackNaN,2);
fig77 = figure(77); clf;
fig77.Position = [0 0 size(stackdst,2)+400 size(stackdst,1)+500];
hold on
plot(Xax,stackNaN_sum,'b')
plot(Xax,stackNaN_shuf_sum,'r')
box off
% axis equal
set(gcf,'PaperPositionMode','auto')
ylim([0 (size(stackNaN,1)/16)*10]); xlim([1 size(stackNaN,2)]);
% saveas(gcf,strcat(DataPath,filesep,Prefix,'PairwiseComp_StackSumPlot'),'epsc');
% saveas(gcf,strcat(DataPath,filesep,Prefix,'PairwiseComp_StackSumPlot'),'svg');

% variance along z
% hmcolmap = flipud(double(parula));
fig81 = figure(81); clf;
fig81.Position = [0 0 500 500];
RD_var = nanvar(RD_all,1,3);
heatmapcust(RD_var,[],[],[],'ColorBar',1,'GridLines','-',...
    'TickAngle',270,'ShowAllTicks',1,'UseLogColormap',false,...
    'Colormap',winter,'MaxColorValue',0.1250,'MinColorValue',0);
hold on
ax = gca;
ax.TickLength = [0 0];
for c=1:size(RD_var,1)
    for r=1:size(RD_var,2)
        pos1 = [(c-0.25) (r-0.1) 0.25 0.25];
        pos2 = [(c) (r-0.1) 0.25 0.25];
        if c == r
            continue
        end
        if c < r
            rectangle('Position',pos1,'Curvature',[1 1],...
                'FaceColor',rgbs(c,:),'LineWidth',0.1);
            rectangle('Position',pos2,'Curvature',[1 1],...
                'FaceColor',rgbs(r,:),'LineWidth',0.1);
        else
            rectangle('Position',pos1,'Curvature',[1 1],...
                'FaceColor',rgbs(r,:),'LineWidth',0.1);
            rectangle('Position',pos2,'Curvature',[1 1],...
                'FaceColor',rgbs(c,:),'LineWidth',0.1);
        end
    end
end
hold off
axis square
% saveas(gcf,strcat(DataPath,filesep,Prefix,'2D_SUBSETidentMLF_VarianceOverZ'),'epsc');
% saveas(gcf,strcat(DataPath,filesep,Prefix,'2D_SUBSETidentMLF_VarianceOverZ'),'svg');
% saveas(gcf,strcat(DataPath,filesep,Prefix,'2D_SUBSETnucMLFMauth_VarianceOverZ'),'epsc');
% saveas(gcf,strcat(DataPath,filesep,Prefix,'2D_SUBSETnucMLFMauth_VarianceOverZ'),'svg');

% shuffle variance along z
fig82 = figure(82); clf;
fig82.Position = [0 0 500 500];
RD_shuf_var = nanvar(RD_shuf_all,1,3);
heatmapcust(RD_shuf_var,[],[],[],'ColorBar',1,'GridLines','-',...
    'TickAngle',270,'ShowAllTicks',1,'UseLogColormap',false,...
    'Colormap',winter,'MaxColorValue',0.1,'MinColorValue',0);
axis square
ax = gca;
ax.TickLength = [0 0];

% instantaneous heatmap at specific z
fig83 = figure(83); clf;
fig83.Position = [0 0 500 500];
i = 3547; disp(i+2013)
RD_oneZ = RD_all(:,:,i);
heatmapcust(RD_oneZ,[],[],[],'ColorBar',1,'GridLines','-',...
    'TickAngle',270,'ShowAllTicks',1,'UseLogColormap',false,...
    'Colormap',winter,'MaxColorValue',1,'MinColorValue',0);
hold on
ax = gca;
ax.TickLength = [0 0];
for c=1:size(RD_oneZ,1)
    for r=1:size(RD_oneZ,2)
        pos1 = [(c-0.25) (r-0.1) 0.25 0.25];
        pos2 = [(c) (r-0.1) 0.25 0.25];
        if c == r
            continue
        end
        if c < r
            rectangle('Position',pos1,'Curvature',[1 1],...
                'FaceColor',rgbs(c,:),'LineWidth',0.1);
            rectangle('Position',pos2,'Curvature',[1 1],...
                'FaceColor',rgbs(r,:),'LineWidth',0.1);
        else
            rectangle('Position',pos1,'Curvature',[1 1],...
                'FaceColor',rgbs(r,:),'LineWidth',0.1);
            rectangle('Position',pos2,'Curvature',[1 1],...
                'FaceColor',rgbs(c,:),'LineWidth',0.1);
        end
    end
end
hold off
axis square
% saveas(gcf,strcat(DataPath,filesep,Prefix,'2D_SUBSETnucMLFMauth_DiffMatSingleZ5560'),'epsc');
% saveas(gcf,strcat(DataPath,filesep,Prefix,'2D_SUBSETnucMLFMauth_DiffMatSingleZ5560'),'svg');
% fig84 = figure(84); clf;
% fig84.Position = [0 0 500 500];
% i = 3547; disp(i+2013)
% RD_oneZ = RDum_all(:,:,i);
% heatmapcust(RD_oneZ,[],[],[],'ColorBar',1,'GridLines','-',...
%     'TickAngle',270,'ShowAllTicks',1,'UseLogColormap',false,...
%     'Colormap',winter,'MaxColorValue',8000,'MinColorValue',0);
% hold on
% ax = gca;
% ax.TickLength = [0 0];
% hold off
% axis square

% 2D plot of points
fig86 = figure(86); clf;
fig86.Position = [0 0 500 500];
dodraw = 1;
% for ik = 1:400
i = 3500-13+60; disp(i+2013)
for pair = 1:npairs
    x1 = points(i,pair,1);
    y1 = points(i,pair,2);
    x2 = points(i,pair,3);
    y2 = points(i,pair,4);
    orgpoints(pair,:) = [x1 y1];
    refpoints(pair,:) = [-x2 y2];
    if dodraw
        if ~isnan(x1) && ~isnan(y1)
            x1 = (x1-xmin)/(xmax-xmin)*400;
            y1 = 400-(y1-ymin)/(ymax-ymin)*400;
            pos1 = [x1 y1 5 5];
            rectangle('Position',pos1,'Curvature',[1 1],...
                'FaceColor',rgbs(pair,:),'LineWidth',0.1);
        end
        if ~isnan(x2) && ~isnan(y2)
            x2 = (x2-xmin)/(xmax-xmin)*400;
            y2 = 400-(y2-ymin)/(ymax-ymin)*400;
            pos2 = [x2 y2 5 5];
            rectangle('Position',pos2,'Curvature',[1 1],...
                'FaceColor',rgbs(pair,:),'LineWidth',0.1);
        end
    end
end
axis equal
% pause(0.1)
% end
% saveas(gcf,strcat(DataPath,filesep,Prefix,'2D_SUBSETnucMLFMauth_PointsSingleZ5560'),'epsc');
% saveas(gcf,strcat(DataPath,filesep,Prefix,'2D_SUBSETnucMLFMauth_PointsSingleZ5560'),'svg');
