clear all;
close all;
load('geometry_added_handle.mat');

load('dense_points.mat');
l = length(x);
for g=1:l
    F = @(q)mFun(x(g),y(g),1,q(1),q(2),q(3));
    options = optimset('Display','none','MaxFunEvals',300,'MaxIter',800);
    [val fval exit] = fsolve(F,[1,1,0],options);
    t(1,g) = 1;
    t(2,g) = val(1);
    t(3,g) = val(2);
    beta(g) = rad2deg(val(3));
    fprintf('%i / %i\nTensions = (%f,%f,%f)\nBeta=%f degrees\n\n',g,l,t(1,g),t(2,g),t(3,g),beta(g));
end
save('solution_blackbox.mat');