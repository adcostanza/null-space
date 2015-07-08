clear all;
close all;
load('matlabFunction_withGeometry.mat');

load('dense_points.mat');
l = length(x);
for g=1:l
    nFun = @(theta)fun(x(g),y(g),theta);
    all_e = [];
    pspan = pi()/4;
    nspan = -pspan;
    num = 1000;
    b = [];
    while(1) 
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
                if(min_val < 0.001), break; end
                oldmin = min_val;
            %end
        end
        fprintf('%i / %i\n',g,l);
    fprintf('(x,y) = (%f,%f)\nMinimum eigenvalue: %f, Beta: %f\n',x(g),y(g),min_val,m_beta);
        %figure;
        %scatter(b, all_e);
        %xlabel('Beta (rad)');
        %ylabel('min(abs(e))');
        %t = sprintf('(x,y) = (%i,%i)',x(g),y(g));
        %title(t);
        min_beta(g) = m_beta;
        if(min_val < 0.001), break; end
        if(min_val >= 0.001)
            nspan = m_beta-0.1;
            pspan = m_beta+0.1;
            num = num+1000;
        end
        if(num >= 6000), break; end;
    end
end
save('dense_beta.mat','min_beta');