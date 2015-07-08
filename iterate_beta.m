clear all;
close all;
load('matlabFunction_withGeometry.mat');

load('unequal_points.mat');
l = length(x);
for g=1:l
    nFun = @(theta)fun(x(g),y(g),theta);
    all_e = [];
    pspan = pi()/2;
    nspan = -pspan;
    num = 10000;
    b = [];
        min_val = 100;
        oldmin = 100;
        for beta = linspace(nspan,pspan,num)
            Sjaca = nFun(beta);
            e = eig(Sjaca);
            %if(isreal(e))
                all_e = [all_e min(abs(e))];

                min_val = min(abs([e(1),e(2),e(3),min_val]));
                if(min_val ~= oldmin) 
                    m_beta = beta;
                end
                b = [b beta];
                oldmin = min_val;
            %end
        end
    fprintf('(x,y) = (%f,%f)\nMinimum eigenvalue: %f, Beta: %f\n',x(g),y(g),min_val,m_beta);
        figure;
        scatter(b, all_e);
        xlabel('Beta (rad)');
        ylabel('min(abs(e))');
        t = sprintf('(x,y) = (%i,%i)',x(g),y(g));
        title(t);
        min_beta(g) = m_beta;
end