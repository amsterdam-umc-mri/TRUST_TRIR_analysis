function [selected_voxels, T2_map, RMS_map, rsquare_map,IE_map, incl_voxels] = TRIR_Get_best_n_voxels(mask, data, nr_echos, TIs, use_nr_voxels, test_nr_of_voxels, IE_threshold, method, plots)
% Function to select n voxels with lowest RMS or highest rsquare in TRIR fit
% Copyright 2021, Koen Baas

% Pre-allocate some variables
T2_map = zeros(size(mask));
RMS_map = zeros(size(mask));
rsquare_map = zeros(size(mask));
IE_map = zeros(size(mask));

N   =   size(data,3);

masked_slice = data(:,:,nr_echos*4).*mask;

% Find highest voxels
highest_values = maxk(reshape(masked_slice,[numel(masked_slice),1]),test_nr_of_voxels);
incl_voxels = find(ismember(masked_slice, highest_values));

% Fit for each voxel individually
for k = 1:test_nr_of_voxels
    
    %% Calculate averages
    tmp_mask = incl_voxels(k);
    
    for i = 1:4
        for echo = 1:(N/4)
            
            slice = data(:,:,(i-1)*nr_echos + echo);
            Averages(echo,i) =  mean(slice(tmp_mask));
            
        end
    end
    
    % Do the fit and save the residual error
    x1 = TIs; x1(nr_echos+1:2*nr_echos) = TIs'; x1(2*nr_echos+1:3*nr_echos) = TIs'; x1(3*nr_echos+1:4*nr_echos) = TIs'; x1 = x1';
    z1 = Averages(:,1); z1(nr_echos+1:2*nr_echos) = Averages(:,2); z1(2*nr_echos+1:3*nr_echos) = Averages(:,3); z1(3*nr_echos+1:4*nr_echos) = Averages(:,4);
    y1(1:nr_echos) = 0; y1(nr_echos+1:2*nr_echos) = 40; y1(2*nr_echos+1:3*nr_echos) = 80; y1(3*nr_echos+1:4*nr_echos) = 160;
    Seedpoint = Averages(60,1);
    
    % Fit function
    H = @(M0,T2,IE,T1,x,y)(M0*abs(1-(1+exp(-(y/T2))*IE).*exp(-(x/T1))));
    
    % Perform the fit
    [thefit,fit_info] = fit([x1,y1'],z1,H,'StartPoint',[Seedpoint,60,0.8,1800],'Upper',[Inf,300,1,3000],'Lower',[Seedpoint*0.8,30,0.2,1400]);
    
    
    RMS(k,1) = fit_info.rmse;
    T2_values(k,1) = thefit.T2;
    rsquare(k,1) = fit_info.rsquare;
    IE_values(k,1) = thefit.IE;
    
    T2_map(incl_voxels(k)) = thefit.T2;
    RMS_map(incl_voxels(k)) = fit_info.rmse;
    rsquare_map(incl_voxels(k)) = fit_info.rsquare;
    IE_map(incl_voxels(k)) = thefit.IE;
end

% Find best n voxels
if strcmp(method,'rmse')
    tmp_RMS = RMS;
    tmp_RMS(rsquare<0.90) = NaN;
    tmp_RMS(T2_values<40) = NaN;
    tmp_RMS(IE_values<IE_threshold) = NaN;
    
    for i=1:use_nr_voxels
        [val(i),selected_voxels(i)] = min(tmp_RMS);
        % remove for the next iteration the last smallest value:
        tmp_RMS(selected_voxels(i)) = NaN;
    end
    
    
elseif strcmp(method,'rsquare')
    
    tmp_rsquare = rsquare;
    tmp_rsquare(rsquare<0.90) = NaN;
    tmp_rsquare(T2_values<40) = NaN;
    tmp_rsquare(IE_values<IE_threshold) = NaN;
    
    for i=1:use_nr_voxels
        [val(i),selected_voxels(i)] = max(tmp_rsquare);
        % remove for the next iteration the last smallest value:
        tmp_rsquare(selected_voxels(i)) = NaN;
    end
    
end

selected_voxels = incl_voxels(selected_voxels);
selected_voxels(find(isnan(val))) = [];

%Optional: plot the separate fits of best voxels
if plots
    
    figure,
    for k = 1: use_nr_voxels
        
        %% Calculate averages
        tmp_mask = selected_voxels(k);
        
        for i = 1:4
            for echo = 1:(N/4)
                
                slice = data(:,:,(i-1)*nr_echos + echo);
                Averages(echo,i) =  mean(slice(tmp_mask));
                
            end
        end
        
        % Now all the data is used to simultaneously fit M0, IE, T1 and T2.
        x1 = TIs; x1(nr_echos+1:2*nr_echos) = TIs'; x1(2*nr_echos+1:3*nr_echos) = TIs'; x1(3*nr_echos+1:4*nr_echos) = TIs'; x1 = x1';
        z1 = Averages(:,1); z1(nr_echos+1:2*nr_echos) = Averages(:,2); z1(2*nr_echos+1:3*nr_echos) = Averages(:,3); z1(3*nr_echos+1:4*nr_echos) = Averages(:,4);
        y1(1:nr_echos) = 0; y1(nr_echos+1:2*nr_echos) = 40; y1(2*nr_echos+1:3*nr_echos) = 80; y1(3*nr_echos+1:4*nr_echos) = 160;
        
        % Fit function
        H = @(M0,T2,IE,T1,x,y)(M0*abs(1-(1+exp(-(y/T2))*IE).*exp(-(x/T1))));
        Seedpoint = Averages(60,1);
        
        % Perform the fit
        [thefit,fit_info] = fit([x1,y1'],z1,H,'StartPoint',[Seedpoint,60,0.8,1800],'Upper',[Inf,300,1,3000],'Lower',[Seedpoint*0.8,30,0.2,1400]);
        
        subplot(5,2,k)
        plot(thefit,[x1,y1'],z1);
        yticks([0 40 80 120 160]);
        %set(gca,'FontSize',32)
        xlabel('TI [ms]')
        ylabel('eTE [ms]')
        zlabel('Magnetization [a.u.]')
        title(['Voxel ' num2str(k) ', T2 = ' num2str(thefit.T2)])
    end
    
end


end

