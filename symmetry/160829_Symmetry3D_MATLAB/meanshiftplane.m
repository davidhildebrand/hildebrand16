function [p1,p2,p3,p4,V] = meanshiftplane(pperps,displs,ms)

V = pperps;
md = mean(displs);
sd = std(displs);
D = 0.1*(displs-md)/sd;
R = [V D];

% old version
% n = length(displs);
% V0 = V(round(n/2),:);
% D0 = D(round(n/2));
% m = [V0 D0];
% [K,N] = meanshift_vd_3D_param(R);
% er = [];
% count = 0;
% while count < 30
%     count = count+1;
%     m = meanshift_vd_3D(R,m,K,N);
%     er = [er norm(m-[1 0 0 0])];
% end
% figure
% plot(er)
% title('mean shift error')
% meanshift_vd_3D_plot(R,m);
% V = m(1:3);
% D = m(4)*10*sd+md;

% new version
k = 4;
s = 0.1;
m = go_evmf(R,k,s);
meanshift_vd_3D_plot(R,m);
figure, subplot(1,1,1)
[x,y,z] = sphere;
surf(x,y,z), alpha(0.1), hold on
DV = V+repmat(D,[1 3]).*V;
plot3(DV(:,1),DV(:,2),DV(:,3),'.g')
dgm = m(1:3)+m(4)*m(1:3);
plot3(dgm(:,1),dgm(:,2),dgm(:,3),'sr')
hold off
axis equal, axis([-2 2 -2 2 -2 2])
view(0,0)
title('final vec/disp')
% zoom(2)
V = m(1:3); % plane's perpendicular vector
D = m(4)*10*sd+md; % distance to origin

mn2 = min(ms(:,2)); mx2 = max(ms(:,2));
mn3 = min(ms(:,3)); mx3 = max(ms(:,3));

% 4 points on plane
pp = @(y) y+(D-dot(y,V))*V;
p1 = pp([0 mn2 mn3]);
p2 = pp([0 mx2 mn3]);
p3 = pp([0 mx2 mx3]);
p4 = pp([0 mn2 mx3]);

end