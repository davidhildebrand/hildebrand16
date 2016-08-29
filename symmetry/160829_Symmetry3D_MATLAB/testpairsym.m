p = [1 0 0];
q = [0 0 0];

for a = linspace(0,2*pi,16)
    tp = [cos(a) 0 sin(a)];
    tq = [cos(pi/3) 0 sin(pi/3)];

    [m,s] = pairwisesymmetry(p,tp,q,tq);
end
close all