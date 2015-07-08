load('min_beta.mat');
load('unequal_points.mat');
injx = x;
clear x;
for g=1:length(min_beta)
    inj = [injx(g); y(g)];
     r = 12/2; % radius of injector, 8mm
    Al = 100; Ar=125; %Arm left, Arm right lengths, 100mm
    dx0 = 7.17; dy0=4.14; dy1=8.28; %dimensions of injector, assuming symmetric left and right sides
    alpha = radtodeg(0.5236); %angle of arms for workspace


    %calculate coordinates of left and right base centers
    bl = [-Al*sind(alpha);Al*cosd(alpha)]; %left base coordinate
    br = [Ar*sind(alpha);Ar*cosd(alpha)]; %right base coordinate

  
    %calculate h vectors to knot from injector center
    h0 = [-dx0;dy0]; %left string h vector
    h1 = [0;-dy1]; % middle
    h2=[dx0;dy0]; % right

    %create rotation matrix from injector coord frame to base
     %orientation of injector
     x = sym('x',[3,1]);
    Rot = [cos(min_beta(g)),-sin(min_beta(g));sin(min_beta(g)), cos(min_beta(g))];

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
    %%
    Sjac = [s0u,s1u,s2u];

    %calculate bottom part of statics jacobian
    q0 = s0u(1)*h0b(2)-s0u(2)*h0b(1);
    q1 = s1u(1)*h1b(2)-s1u(2)*h1b(1);
    q2 = s2u(1)*h2b(2)-s2u(2)*h2b(1);

    Sjac = [Sjac; q0,q1,q2];
    fun = Sjac * [x(1);x(2);x(3)];
    %pretty(simplify(fun))

    val = vpasolve(fun);
    t1 = eval(val.x1); 
    t2=eval(val.x2); 
    t3 = eval(val.x3);

    fprintf('Orientation: %f rad\nTension Left: %f N\nTension Middle: %f N\nTension Right: %f N\n',min_beta(g),t1,t2,t3);



   
    ltana = subs(ltan,min_beta(g),eval(val.x1));
  
    rtana = subs(rtan,min_beta(g),eval(val.x1));


    k0a = eval(subs(k0,min_beta(g),eval(val.x1)));
    k1a = eval(subs(k1,min_beta(g),eval(val.x1)));
    k2a = eval(subs(k2,min_beta(g),eval(val.x1)));
    lalpha = eval(subs(lalpha,min_beta(g),eval(val.x1)));
    ralpha = eval(subs(ralpha,min_beta(g),eval(val.x1)));
    %full string length calculation
    sl0(g) = eval(norm(subs(l0,min_beta(g),eval(val.x1)))+lalpha*r);
    sl2(g) = eval(norm(subs(l2,min_beta(g),eval(val.x1)))+ralpha*r);
    sl1(g) = eval(norm(subs(l1,min_beta(g),eval(val.x1))));
    fprintf('\n(sl0,sl1,sl2) = (%f,%f,%f)\n',sl0(g),sl1(g),sl2(g));

    hold off;
end