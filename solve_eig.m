clear all;
close all;
load('geometry_added.mat');
Sjac = subs(Sjac,injx,0);
Sjac = subs(Sjac,injy,20);
mFun = matlabFunction(Sjac);
syms lam;
z = det(mFun(10)-lam*eye(3));
eig(mFun(10))