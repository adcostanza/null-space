inj = [10;90];

 r = 12/2; % radius of injector, 8mm
Al = 100; Ar=125; %Arm left, Arm right lengths, 100mm
dx0 = 7.17; dy0=4.14; dy1=8.28; %dimensions of injector, assuming symmetric left and right sides
alpha = radtodeg(0.5236); %angle of arms for workspace


%calculate coordinates of left and right base centers
bl = [-Al*sind(alpha);Al*cosd(alpha)]; %left base coordinate
br = [Ar*sind(alpha);Ar*cosd(alpha)]; %right base coordinate

%plot what is already known and accept injector location
figure(2);
hold off;
plot(0,0,'.k','MarkerSize',50);

circle(bl(1),bl(2),r);
circle(br(1),br(2),r);
hold on;
plot([0,bl(1)],[0,bl(2)],'k','LineWidth',3);
hold on;
plot([0,br(1)],[0,br(2)],'k','LineWidth',3);

axis equal;


%calculate h vectors to knot from injector center
h0 = [-dx0;dy0]; %left string h vector
h1 = [0;-dy1]; % middle
h2=[dx0;dy0]; % right

%create rotation matrix from injector coord frame to base
 %orientation of injector
 x = sym('x',[3,1]);
Rot = [cos(x(1)),-sin(x(1));sin(x(1)), cos(x(1))];

%calculate h vectors in base coordinate space
h0b = (Rot*h0);
h1b = (Rot*h1);
h2b = (Rot*h2);

%calculate knot locations in base coordinate frame
k0 = inj+h0b; %left
k1 = inj+h1b; %middle
k2 = inj+h2b; %right

%calculate l0, r0 vector from knot to center of bases
l0 = bl-k0;
l2 = br-k2;
l1=-k1; %no tangent point for middle, this is the actual S vector

%calculate tangent angle
kappal = acos(r/norm(l0));
kappar = acos(r/norm(l2));

rltan = [cos(kappal), -sin(kappal);sin(kappal), cos(kappal)]*-l0;
rrtan = [cos(-kappar), -sin(-kappar);sin(-kappar), cos(-kappar)]*-l2;

rltan = r*rltan/norm(rltan);
rrtan = r*rrtan/norm(rrtan);

rrtan = rrtan(1:2); %remove z component
rltan = rltan(1:2); %remove z component

%calculate rperp on each base, its the fixed location that the string will
%never reach and we add the angle between it and the rtan vectors. it is
%perpendicular to the arm
rlperp = cross([bl;0],[0;0;-1]); %arbitrary magnitude in correct direciton
rrperp = cross([br;0],[0;0;1]); %arbitrary magnitude in correct direction

rlperp = r*rlperp/sqrt(power(rlperp(1),2)+power(rlperp(2),2)); %make length of radius
rrperp = r*rrperp/sqrt(power(rrperp(1),2)+power(rrperp(2),2)); %make magnitude radius

rrperp = rrperp(1:2); %remove z component
rlperp = rlperp(1:2); %remove z component

%get angle of string around pulleys

rrperpsq = sqrt(power(rrperp(1),2)+power(rrperp(2),2));
rrtansq = sqrt(power(rrtan(1),2)+power(rrtan(2),2));
ralpha = acos(dot(rrperp,rrtan)/rrperpsq/rrtansq);


rlperpsq = sqrt(power(rlperp(1),2)+power(rlperp(2),2));
rltansq = sqrt(power(rltan(1),2)+power(rltan(2),2));
lalpha = acos(dot(rlperp,rltan)/rlperpsq/rltansq);

%calculate tangent points on bases
ltan = bl+rltan;
rtan = br+rrtan;

%calculate string vectors S
s0 = ltan-k0;
s2 = rtan-k2;
s1 = l1;

%create unit vectors
s0u = s0/sqrt(power(s0(1),2)+power(s0(2),2));
s1u = s1/sqrt(power(s1(1),2)+power(s1(2),2));
s2u = s2/sqrt(power(s2(1),2)+power(s2(2),2));

Sjac = [s0u,s1u,s2u];

%calculate bottom part of statics jacobian
q0 = s0u(1)*h0b(2)-s0u(2)*h0b(1);
q1 = s1u(1)*h1b(2)-s1u(2)*h1b(1);
q2 = s2u(1)*h2b(2)-s2u(2)*h2b(1);

Sjac = [Sjac; q0,q1,q2];
fun = Sjac * [1;x(2);x(3)];
%pretty(simplify(fun))

val = vpasolve(fun);
theta = eval(val.x1);
t1 = 1;
t2=eval(val.x2); 
t3 = eval(val.x3);

fprintf('Orientation: %f rad\nTension Left: %f N\nTension Middle: %f N\nTension Right: %f N\n',theta,t1,t2,t3);



hold on;
ltana = subs(ltan,x(1),eval(val.x1));
plot(ltana(1),ltana(2),'.r','MarkerSize',20);
hold on;
rtana = subs(rtan,x(1),eval(val.x1));
plot(rtana(1),rtana(2),'.r','MarkerSize',20);
hold on;
set(get(gca,'xlabel'),'FontSize', 13);
set(get(gca,'ylabel'),'FontSize', 13);
set(get(gca,'title'),'FontSize', 13);
set(gca,'FontSize',11);
set(gcf,'color','w');

k0a = eval(subs(k0,x(1),eval(val.x1)));
k1a = eval(subs(k1,x(1),eval(val.x1)));
k2a = eval(subs(k2,x(1),eval(val.x1)));
lalpha = eval(subs(lalpha,x(1),eval(val.x1)));
ralpha = eval(subs(ralpha,x(1),eval(val.x1)));
%full string length calculation
sl0 = eval(norm(subs(l0,x(1),eval(val.x1)))+lalpha*r);
sl2 = eval(norm(subs(l2,x(1),eval(val.x1)))+ralpha*r);
sl1 = eval(norm(subs(l1,x(1),eval(val.x1))));

plot([k0a(1),ltana(1)],[k0a(2),ltana(2)]);
plot([k1a(1),0],[k1a(2),0]);
plot([k2a(1),rtana(1)],[k2a(2),rtana(2)]);
plot(inj(1),inj(2),'.k','MarkerSize',10);
ylim([-5,120]);
str = sprintf('(x,y,th) = (%.2f,%.2f,%.2f)\n(t1, t2, t3) = (%.2f,%.2f,%.2f)\n(s1, s2, s3) = (%.2f,%.2f,%.2f)',inj(1),inj(2),theta,t1,t2,t3,sl0,sl1,sl2);
text(8,6,str);
xlabel('X (mm)');
ylabel('Y (mm)');
title('Numerically Solved Statics and Inverse Kinematics Example');

hold off;