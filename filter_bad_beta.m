clear all;
close all;
load('matlabFunction_withGeometry.mat');
load('dense_points.mat');
load('dense_beta.mat');
a = 1;
for g=1:length(x)
    S = @(t)fun(x(g),y(g),min_beta(g))*[t(1);t(2);t(3)];
    t_val_temp = fsolve(S,[1,1,1]);
    if((t_val_temp(1) < 0 && t_val_temp(2) < 0 && t_val_temp(3) < 0)&& isreal(t_val_temp))
        t_val_temp = -t_val_temp;
    elseif ((t_val_temp(1) > 0 && t_val_temp(2) > 0 && t_val_temp(3) > 0) && isreal(t_val_temp))
        t_val_temp = t_val_temp;
    else
        t_val(:,a) = [NaN;NaN;NaN];
        a=a+1;
        x_(a) = x(g);
        y_(a) = y(g);
        beta_(a) = NaN;
        continue;
    end
    t_val(:,a) = t_val_temp/min(abs(t_val_temp))';
    a=a+1;
    x_(a) = x(g);
    y_(a) = y(g);
    beta_(a) = min_beta(g);
end
t_val
save('beta_filtered.mat','t_val','x_','y_','beta_');