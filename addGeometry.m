clear all;
load('fun.mat');
syms r Al Ar dx0 dy0 dy1 alpha theta t1 t2 injx injy
r_ = 12/2; % radius of injector, 8mm
Al_ = 100; Ar_=125; %Arm left, Arm right lengths, 100mm
dx0_ = 5; dy0_=5; dy1_=5; %dimensions of injector, assuming symmetric left and right sides
alpha_ = 0.5236; %angle of arms for workspace
injx_ = 0;
injy_ = 13.4;

Sjac = subs(Sjac,r,r_);
Sjac = subs(Sjac,Al,Al_);
Sjac = subs(Sjac,Ar,Ar_);
Sjac = subs(Sjac,dx0,dx0_);
Sjac = subs(Sjac,dy0,dy0_);
Sjac = subs(Sjac,dy1,dy1_);
Sjac = subs(Sjac,alpha,alpha_);

save('geometry_added.mat','Sjac','injx','injy','theta');

