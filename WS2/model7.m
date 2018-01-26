%function retval = model7(a,b)
%display(a);

% fixed parameters
t_max = 30 * 24;   % 30 days = 720 hours
dt = 0.1; 
T = 0:dt:t_max; 

birth_rate = 1/24;
K0 = 15000; 
Kmax = 9.4*10^9;
Kdelta = Kmax - K0;
CTC_death_rate = 1/24;  % Circulating Tumor Cells

% varying params (from command line)
% to test this, execute 'clearvars'
set_params_from_cmdline = true;
if (set_params_from_cmdline == false)
    intravasation_rate = 1.e-7;
    seeding_rate = 0.0001;
    vasc_rate = 0.0014;
end

% create a single tumor site
tumor_site.birth_rate = birth_rate;
tumor_site.intravasation_rate = intravasation_rate;
tumor_site.vasc_rate = vasc_rate;
tumor_site.N = zeros(size(T)); % tumor cells
tumor_site.K = K0*ones(size(T)); % carrying capacity 
tumor_site.B = zeros(size(T)); % blood vessels: range [0,1]; start at 0

clear tumor_sites;
tumor_sites(1) = tumor_site; 

CTC = zeros(size(T));    % Circulating Tumor Cells

number_of_mets = zeros(size(T)); 

% set its initial values
tumor_sites(1).N(1) = 1; 
tumor_sites(1).K(1) = K0; 
tumor_sites(1).B(1) = 0.0; 

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
            CTC(i) = CTC(i) + loss;  % circ tumor cells += loss

        % update blood vessels
        tumor_sites(j).B(i) = tumor_sites(j).B(i-1) + dt*tumor_sites(j).vasc_rate;
        tumor_sites(j).B(i) = min(tumor_sites(j).B(i), 1.0);
        
        % update carrying capacity
        %tumor_sites(j).K(i) = K0 + tumor_sites(j).K(i) * (Kmax - K0);
        tumor_sites(j).K(i) = K0 + Kdelta * tumor_sites(j).B(i);
        
        % rate of vascularization increases when the tumor cell population
        % is closer to the current car-rying capacity
        if (tumor_sites(j).N(i) > 0.9*tumor_sites(j).K(i))
            %tumor_sites(j).vasc_rate = tumor_sites(j).vasc_rate + 0.1*vasc_rate;
            tumor_sites(j).vasc_rate = tumor_sites(j).vasc_rate + 0.0001*tumor_sites(j).vasc_rate;
        end
    end
    
    % check if we seed a new met
    r = rand(); 
    prob = dt*seeding_rate*CTC(i);
    if( r < prob )
        tumor_sites(end+1) = tumor_site; 
        tumor_sites(end).N(i) = 1; 
        CTC(i) = CTC(i) - 1; 
        %disp( sprintf('new metastasis at time %3.1f days! %u tumor sites' , T(i)/24 , length(tumor_sites) ))
    end
    
    % calculate the number of metastases 
    number_of_mets(i) = length( tumor_sites ) - 1;  % don't count primary site
    
    if (number_of_mets(i) > 5000)
        disp('breaking out due to exceeded max mets')
        break
    end
end

for idx=1:2 
    plot(T/24,tumor_sites(idx).N)
    xlabel('T(days)','FontSize',20)
    ystr = sprintf('sites(%d).N',idx)
    ylabel(ystr,'FontSize',20)
    png_fname = sprintf('T_N%d',idx)
    print(png_fname,'-dpng')
    
    plot(T/24,tumor_sites(idx).B)
    xlabel('T(days)','FontSize',20)
    ystr = sprintf('sites(%d).B',idx)
    ylabel(ystr,'FontSize',20)
    png_fname = sprintf('T_B%d',idx)
    print(png_fname,'-dpng')
    
    plot(T/24,tumor_sites(idx).K)
    xlabel('T(days)','FontSize',20)
    ystr = sprintf('sites(%d).K',idx)
    ylabel(ystr,'FontSize',20)
    png_fname = sprintf('T_K%d',idx)
    print(png_fname,'-dpng')
end
plot(T/24,number_of_mets)
xlabel('T(days)','FontSize',20)
ylabel('# mets','FontSize',20)
png_fname = sprintf('T_num_mets')
print(png_fname,'-dpng')

%retval = 0;
