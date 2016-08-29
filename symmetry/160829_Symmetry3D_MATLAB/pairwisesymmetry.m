function [m,s] = pairwisesymmetry(p,tp,q,tq)
    m = (p+q)/2;
    
    pq = ((q-p)/norm(q-p));

    R_tp = tp-2*dot(tp,pq)*pq; % reflection of tp w.r.t. 'plane through the origin orthogonal to pq'

    s = dot(tq,R_tp); % mirror symmetry (the closer to +1, the more mirror symmetric)
    
%     plot3(p(1),p(2),p(3),'or'), hold on
%     plot3([p(1) p(1)+tp(1)], [p(2) p(2)+tp(2)], [p(3) p(3)+tp(3)], '-r')
%     plot3(q(1),q(2),q(3),'og'), hold on
%     plot3([q(1) q(1)+tq(1)], [q(2) q(2)+tq(2)], [q(3) q(3)+tq(3)], '-g')
%     plot3([q(1) q(1)+R_tp(1)], [q(2) q(2)+R_tp(2)], [q(3) q(3)+R_tp(3)], '-r'), hold off
%     xlabel('x'), ylabel('y'), zlabel('z')
%     axis([-2 2 -2 2 -2 2]), grid on
%     view(0,0)
%     title(sprintf('%f',s))
%     pause
end