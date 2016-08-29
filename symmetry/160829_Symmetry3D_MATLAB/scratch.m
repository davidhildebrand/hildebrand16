clear
clc

% --------------------------------------------------
fprintf('loading data\n');
% --------------------------------------------------

% parameters
load('/home/mc457/files/CellBiology/IDAC/Marcelo/Hildebrand/sp.mat'); % symmetry plane points
load('/home/mc457/files/CellBiology/IDAC/Marcelo/Hildebrand/V.mat'); % symmetry plane perpendicular vector

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
fprintf('plot\n');
% --------------------------------------------------

plot3(D(1:100:end,4),D(1:100:end,5),D(1:100:end,6),'.'), hold on
fill3(sp(:,1),sp(:,2),sp(:,3),'k'), alpha(0.1), hold off, axis off