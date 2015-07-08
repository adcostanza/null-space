clear all;
close all;
load('new_s.mat');
load('lengths.mat');

figure;
s(:,29) = s2(:,29);
for d=1:3
    subplot(3,1,d);
    plot(s(d,:));
    hold on;
    plot(s2(d,:));
    
    legend('Old length','New length');
end
p_s = (s-s2)./max(max(s))*100;
figure;
plot(p_s');


load('old_tensions.mat');
load('new_tensions.mat');
old_t(:,29) = t_val(:,29);
figure;
for d=1:3
    subplot(3,1,d);
    o_t(d,:) = old_t(d,:)*250./max(old_t);
    plot(o_t(d,:));
    hold on;
    n_t(d,:) = t_val(d,:)*250./max(t_val);
    plot(n_t(d,:));
    legend('Old tension','New tension');
end
p_t = (o_t-n_t)./250*100;
figure;
plot(p_t');

figure;
load('beta.mat');
load('min_beta.mat');
min_beta = min_beta';
min_beta(29) = deg2rad(beta(29));
plot(beta);
hold on;
plot(rad2deg(min_beta));
