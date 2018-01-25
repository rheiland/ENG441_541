% this version has nicer outputs

% set parameters
t_max = 30 * 24;
dt = 0.1; 
T = 0:dt:t_max; 

intravasation_rate = 0.001; 
seeding_rate = 0.001; 
birth_rate = 1/24; 
K0 = 15000; 
CTC_death_rate = 1/24;  % Circulating Tumor Cells

% create a single tumor site
tumor_site.birth_rate = birth_rate;
tumor_site.intravasation_rate = intravasation_rate; 
tumor_site.N = zeros(size(T)); % tumor cells
tumor_site.K = K0*ones(size(T)); % carrying capacity 
tumor_site.B = 0.1*ones(size(T)); % blood vessels 

clear tumor_sites;
tumor_sites(1) = tumor_site; 

CTC = zeros(size(T)); 

number_of_mets = zeros(size(T)); 

% set its initial values
tumor_sites(1).N(1) = 1; 
tumor_sites(1).K(1) = K0; 
tumor_sites(1).B(1) = 0.1; 

CTC(1) = 0; 

% step through time until done
for i=2:length(T)
    t = T(i); 
    
    % update CTC count
    CTC(i) = CTC(i-1) - dt*CTC_death_rate*CTC(i-1); 

    for j=1:length(tumor_sites)
        % update cell count 
        tumor_sites(j).N(i) = tumor_sites(j).N(i-1) + dt*tumor_sites(j).birth_rate* ...
            tumor_sites(j).N(i-1)*(1 - tumor_sites(j).N(i-1)/tumor_sites(j).K(i-1)); 

        % intravasation
            % calculate number of cells lost
            loss = dt*tumor_sites(j).intravasation_rate * ...
                tumor_sites(j).B(i-1) * tumor_sites(j).N(i-1);
            % update cell count 
            tumor_sites(j).N(i) = tumor_sites(j).N(i) - loss;
            CTC(i) = CTC(i) + loss; 

        % update blood vessels

        % update carrying capacity 
    end
    
    % check if we seed a new met
    r = rand(); 
    prob = dt*seeding_rate*CTC(i);
    if( r < prob )
        tumor_sites(end+1) = tumor_site; 
        tumor_sites(end).N(i) = 1; 
        CTC(i) = CTC(i) - 1; 
        disp( sprintf('new metastasis at time %3.1f days! %u tumor sites' , T(i)/24 , length(tumor_sites) ))
    end
    
    
    % calculate the nubmer of metastases 

    number_of_mets(i) = length( tumor_sites ) - 1; 
    
    % update time
end
