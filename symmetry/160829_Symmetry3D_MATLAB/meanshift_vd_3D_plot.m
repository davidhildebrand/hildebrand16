function meanshift_vd_3D_plot(R,m)
    V = R(:,1:3);
    D = R(:,4);
    V0 = m(1:3);
    D0 = m(4);
    
    figure
    subplot(1,2,1)
    [x,y,z] = sphere;
    surf(x,y,z), alpha(0.8), hold on
    plot3(V(:,1),V(:,2),V(:,3),'.g','MarkerSize',0.5)
    plot3(V0(1),V0(2),V0(3),'or'), hold off
    axis equal, axis([-1 1 -1 1 -1 1])
    view(86,4)
    title(sprintf('final vec: (%.02f, %.02f, %.02f)', V0(1), V0(2), V0(3)))

    subplot(1,2,2)
    plot(D,zeros(size(D)),'.g','MarkerSize',0.5), hold on
    plot(D0,0,'or'), hold off
    axis equal, axis([-1 1 -1 1])
    title(sprintf('final disp: %f', D0))
end