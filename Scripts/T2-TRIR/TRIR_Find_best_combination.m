function [Best_mask] = TRIR_Find_best_combination(Combinations, data, nr_echos, TIs, method, plots)
% This function finds the optimal combination from the previously derived
% combinations.
% Copyright 2021, Koen Baas

N   =   size(data,3);

% Pre-allocate some variables
T2_combinations = NaN(length(Combinations),1);
T1_combinations = NaN(length(Combinations),1);
IE_combinations = NaN(length(Combinations),1);
Errors_rmse = NaN(length(Combinations),1);
Errors_sse = NaN(length(Combinations),1);
Errors_rsquare  =NaN(length(Combinations),1);
Errors_dfe = NaN(length(Combinations),1);
Errors_adjrsquare = NaN(length(Combinations),1);


for k = 1: length(Combinations)
    
    %% Calculate averages
    mask = nonzeros(Combinations(k,:));
    
    for i = 1:4
        for echo = 1:(N/4)
            
            slice = data(:,:,(i-1)*nr_echos + echo);
            Averages(echo,i) =  mean(slice(mask));
            
        end
    end
    
    x1 = TIs; x1(nr_echos+1:2*nr_echos) = TIs'; x1(2*nr_echos+1:3*nr_echos) = TIs'; x1(3*nr_echos+1:4*nr_echos) = TIs'; x1 = x1';
    z1 = Averages(:,1); z1(nr_echos+1:2*nr_echos) = Averages(:,2); z1(2*nr_echos+1:3*nr_echos) = Averages(:,3); z1(3*nr_echos+1:4*nr_echos) = Averages(:,4);
    y1(1:nr_echos+1) = 0; y1(nr_echos+1:2*nr_echos+1) = 40; y1(2*nr_echos+1:3*nr_echos+1) = 80; y1(3*nr_echos+1:4*nr_echos) = 160;
    Seedpoint = z1(1);
   
    % Do the fit and save the residual error
    try
        % fit function
        H = @(M0,T2,IE,T1,x,y)(M0*abs(1-(1+exp(-(y/T2))*IE).*exp(-(x/T1))));

        % Perform the fit
        [thefit,fit_info] = fit([x1,y1'],z1,H,'StartPoint',[Seedpoint,60,0.8,1800],'Upper',[Inf,300,1,3000],'Lower',[Seedpoint*0.8,30,0.2,1400]);

        confidenceint = confint(thefit);
        T2_CI(k) = (confidenceint(4)-confidenceint(3));
        
        T2_combinations(k) = thefit.T2;
        T1_combinations(k) = thefit.T1;
        IE_combinations(k) = thefit.IE;
        
        Errors_rmse(k) = fit_info.rmse;
        Errors_sse(k) = fit_info.sse;
        Errors_rsquare(k) = fit_info.rsquare;
        Errors_dfe(k) = fit_info.dfe;
        Errors_adjrsquare(k) = fit_info.adjrsquare;
       
    catch
        T2_combinations(k) = NaN;
        Errors_rmse(k) = Inf;
        
    end
    
end


if strcmp(method, 'rmse')
    Best_combination = find(Errors_rmse == min(Errors_rmse));   
elseif strcmp(method, 'rsquare')
    Best_combination = find(Errors_rsquare == max(Errors_rsquare));
elseif stcmp(method, 'T2_CI')
     Best_combination = find(T2_CI == min(T2_CI));
end

Best_mask = nonzeros(Combinations(Best_combination,:));
Lowest_rmse = Errors_rmse(Best_combination);
Highest_rsquare = Errors_rsquare(Best_combination);

if strcmp(method, 'rmse')
display(['Optimal ROI, based on rmse, consists of voxel(s): ' num2str(Best_mask') ' and results in a rmse of ' num2str(Lowest_rmse)]);
elseif strcmp(method, 'rsquare')
display(['Optimal ROI, based on rsquare, consists of voxel(s): ' num2str(Best_mask') ' and results in a rsquare of ' num2str(Highest_rsquare)]);
end

%Optional: plot the separate fits of best voxels
if plots
selected_voxels = Best_mask;
figure,
for k = 1: length(selected_voxels)

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

subplot(4,4,k)
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



