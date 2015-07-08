clear all; 
close all;
load('beta_filtered.mat');
load('dense_points.mat');
Al = 100;
Ar = 125;
alpha = radtodeg(0.5236); %angle of arms for workspace
r = 12/2; % radius of injector, 8mm
dx0 = 7.17; dy0=4.14; dy1=8.28; %dimensions of injector, assuming symmetric left and right sides
%calculate coordinates of left and right base centers
bl = [-Al*sind(alpha);Al*cosd(alpha)]; %left base coordinate
br = [Ar*sind(alpha);Ar*cosd(alpha)]; %right base coordinate

xinj = x_;
yinj =y_;
F = TriScatteredInterp(xinj',yinj',beta_');
[qx,qy]=meshgrid(min(xinj):.05:max(xinj),0:.05:max(yinj));
qz = F(qx,qy);
figure;
colormap(jet);
mesh(qx,qy,qz);
hold on;
plot3(xinj,yinj,beta_,'o');
axis([bl(1)-r-5, br(1)+r+5, 0, br(2)+r+5 ]);
set(get(gca,'xlabel'),'FontSize', 13);
set(get(gca,'ylabel'),'FontSize', 13);
set(get(gca,'title'),'FontSize', 13);
set(gca,'FontSize',11);
set(gcf,'color','w');
%%
figure;

colormap(jet);
contourf(qx,qy,qz,10);
colorbar;
hold on;
plot(0,0,'.k','MarkerSize',50);
circle(bl(1),bl(2),r);
circle(br(1),br(2),r);
hold on;
plot([0,bl(1)],[0,bl(2)],'k','LineWidth',3);
hold on;
plot([0,br(1)],[0,br(2)],'k','LineWidth',3);
xlabel('X (mm)');
ylabel('Y (mm)');
axis([bl(1)-r-5, br(1)+r+5, 0, br(2)+r+5 ]);
set(get(gca,'xlabel'),'FontSize', 13);
set(get(gca,'ylabel'),'FontSize', 13);
set(get(gca,'title'),'FontSize', 13);
set(gca,'FontSize',11);
set(gcf,'color','w');
figure;

colormap(jet);
pcolor(qx,qy,qz)
colorbar;
set(get(gca,'xlabel'),'FontSize', 13);
set(get(gca,'ylabel'),'FontSize', 13);
set(get(gca,'title'),'FontSize', 13);
set(gca,'FontSize',11);
set(gcf,'color','w');

%plot what is already known and accept injector location
hold on;
plot(0,0,'.k','MarkerSize',50);
circle(bl(1),bl(2),r);
circle(br(1),br(2),r);
hold on;
plot([0,bl(1)],[0,bl(2)],'k','LineWidth',3);
hold on;
plot([0,br(1)],[0,br(2)],'k','LineWidth',3);
axis([bl(1)-r-5, br(1)+r+5, 0, br(2)+r+5 ]);
set(get(gca,'xlabel'),'FontSize', 13);
set(get(gca,'ylabel'),'FontSize', 13);
set(get(gca,'title'),'FontSize', 13);
set(gca,'FontSize',11);
set(gcf,'color','w');
xlabel('X (mm)');
ylabel('Y (mm)');