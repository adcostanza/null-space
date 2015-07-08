clear all;
close all;
load('matlabFunction_withGeometry.mat');
load('unequal_points.mat');
load('min_beta.mat');
for g=1:length(x)
    S = @(t)fun(x(g),y(g),min_beta(g))*[t(1);t(2);t(3)];
    t_val_temp = fsolve(S,[1,1,1]);
    t_val(:,g) = t_val_temp/min(abs(t_val_temp))';
end
t_val
save('tensions.mat','t_val');