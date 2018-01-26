% this version has nicer outputs

close all;
clc;
clear all;

% set parameters
t_max = 30 * 24;   % 30 days = 720 hours
dt = 0.1; 
T = 0:dt:t_max; 

intravasation_rate = 0.0000001; % rate at which cancer cells invade vasculature
seeding_rate = 0.0001;
vasc_rate = 0.0014;  %?????
birth_rate = 1/24;
K0 = 15000; 
Kmax = 2500^3; 
Kdelta = (Kmax - K0)/(1-0);
CTC_death_rate = 1/24;  % Circulating Tumor Cells
r_MN = 0.001;
r_NM_min = 0.0001;

% create a single tumor site
tumor_site.birth_rate = birth_rate;
tumor_site.intravasation_rate = intravasation_rate;
tumor_site.vasc_rate = vasc_rate;
tumor_site.N = zeros(size(T)); % tumor cells
tumor_site.M = zeros(size(T)); % M tumor cells
tumor_site.K = K0*ones(size(T)); % carrying capacity 
tumor_site.B = zeros(size(T)); % blood vessels: range [0,1]; start at 0
tumor_site.r_NM = r_NM_min*ones(size(T));

clear tumor_sites;
tumor_sites(1) = tumor_site; 

CTC = zeros(size(T));    % Circulating Tumor Cells

number_of_mets = zeros(size(T)); 

% set its initial values
tumor_sites(1).N(1) = 1; 
tumor_sites(1).M(1) = 0; 
tumor_sites(1).K(1) = K0; 
tumor_sites(1).B(1) = 0.0; 

CTC(1) = 0; 

% step through time until done
for i=2:length(T)
    t = T(i); 
    
    % update CTC count
    CTC(i) = CTC(i-1) - dt*CTC_death_rate*CTC(i-1); 

    for j=1:length(tumor_sites)
        % update N cell count 
        tumor_sites(j).N(i) = tumor_sites(j).N(i-1) + dt*tumor_sites(j).birth_rate* ...
            tumor_sites(j).N(i-1)*(1 - tumor_sites(j).N(i-1)/tumor_sites(j).K(i-1))- ...
            dt*tumor_site.r_NM(i)*tumor_sites(j).N(i-1) + dt*r_MN*tumor_sites(j).M(i-1);
        
        % update M cell count 
        tumor_sites(j).M(i) = tumor_sites(j).M(i-1) + dt*tumor_site.r_NM(i)*...
            tumor_sites(j).N(i-1)- dt*r_MN*tumor_sites(j).M(i-1);
        
        % intravasation
            % calculate number of cells lost from cell M
            loss = dt*tumor_sites(j).intravasation_rate * ...
                tumor_sites(j).B(i-1) * tumor_sites(j).M(i-1);
            
            % update M cell count 
            tumor_sites(j).M(i) = tumor_sites(j).M(i) - loss;  
            CTC(i) = CTC(i) + loss;  % circ tumor cells += loss

        % update blood vessels
        tumor_sites(j).B(i) = tumor_sites(j).B(i-1) + dt*tumor_sites(j).vasc_rate;
        tumor_sites(j).B(i) = min(tumor_sites(j).B(i), 1.0);
        
        % update carrying capacity
        % tumor_sites(j).K(i) = K0 + tumor_sites(j).K(i) * (Kmax - K0);  %%%%% maybe need a coefficient
        tumor_sites(j).K(i) = K0 + Kdelta*tumor_sites(j).B(i);  %%%%% maybe need a coefficient
        
        % rate of vascularization increases when the tumor cell population
        % is closer to the current car-rying capacity
        if (tumor_sites(j).N(i) > 0.9*tumor_sites(j).K(i))   %%%%%%%%% 0.9
            % tumor_sites(j).vasc_rate = tumor_sites(j).vasc_rate + vasc_rate; %%%%% maybe need a coefficient
            tumor_sites(j).vasc_rate = tumor_sites(j).vasc_rate + 0.0001*tumor_sites(j).vasc_rate;
        end
        
        %  r_NM increases when (M+N) is close to carrying capacity,
        if ((tumor_sites(j).N(i)+ tumor_sites(j).M(i)) > 0.9*tumor_sites(j).K(i))
            tumor_site.r_NM(i) = tumor_site.r_NM(i) + 0.001*(tumor_sites(j).N(i)+...
                tumor_sites(j).M(i))/tumor_sites(j).K(i)*tumor_site.r_NM(i);
        end
                      
    end
    
    % check if we seed a new met
    r = rand(); 
    prob = dt*seeding_rate*CTC(i);
    if( r < prob )
        tumor_sites(end+1) = tumor_site; 
        tumor_sites(end).N(i) = 1; 
        tumor_sites(end).M(i) = 0;
        CTC(i) = CTC(i) - 1; 
        disp( sprintf('new metastasis at time %3.1f days! %u tumor sites' , T(i)/24 , length(tumor_sites) ))
    end
    
    % calculate the number of metastases 
    number_of_mets(i) = length( tumor_sites ) - 1; 
    
end

figure(1)
plot(T/24,tumor_sites(1).N,'LineWidth',2);
hold on;
plot(T/24,tumor_sites(1).M,'LineWidth',2);

