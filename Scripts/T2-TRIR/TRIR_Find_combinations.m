function [Combinations] = TRIR_Find_combinations(mask, data, IE_threshold, nr_echos, TIs)
% This function finds all possible voxels combinations after testing if
% individual voxels result in reliable fitting of the data.
% Copyright 2021, Koen Baas

% Pre-allocate some variables
T2_map = zeros([size(data,1), size(data,2)]);
IE_map = zeros([size(data,1), size(data,2)]);
rsquare_map = zeros([size(data,1), size(data,2)]);
rmse_map = zeros([size(data,1), size(data,2)]);
T1_map = zeros([size(data,1), size(data,2)]);
T2_CI_map = zeros([size(data,1), size(data,2)]);
T1_CI_map = zeros([size(data,1), size(data,2)]);

% First, exclude voxels that are individually not reliable
voxels = find(mask == 1);
for voxel = 1:length(voxels)
    
     %% Calculate averages from best ROI
            N   =   size(data,3);
            
            for i = 1:4
                for echo = 1:(N/4)
                    slice = data(:,:,(i-1)*nr_echos + echo);
                    Averages(echo,i) =  mean(slice(voxels(voxel)));
                end
            end
            
            %% Do simultaneous fit
            x1 = TIs; x1(nr_echos+1:2*nr_echos) = TIs'; x1(2*nr_echos+1:3*nr_echos) = TIs'; x1(3*nr_echos+1:4*nr_echos) = TIs'; x1 = x1';
            z1 = Averages(:,1); z1(nr_echos+1:2*nr_echos) = Averages(:,2); z1(2*nr_echos+1:3*nr_echos) = Averages(:,3); z1(3*nr_echos+1:4*nr_echos) = Averages(:,4);
            y1(1:nr_echos) = 0; y1(nr_echos+1:2*nr_echos) = 40; y1(2*nr_echos+1:3*nr_echos) = 80; y1(3*nr_echos+1:4*nr_echos) = 160;
            
            % Fit function
            H = @(M0,T2,IE,T1,x,y)(M0*abs(1-(1+exp(-(y/T2))*IE).*exp(-(x/T1))));
            Seedpoint = Averages(60,1);
            
            % Perform the fit
            [thefit,fit_info] = fit([x1,y1'],z1,H,'StartPoint',[Seedpoint,60,0.8,1800],'Upper',[Inf,300,1,3000],'Lower',[Seedpoint*0.8,30,0.2,1400]);
            confidenceint = confint(thefit);
                        
            T2_CI = (confidenceint(4)-confidenceint(3));
            T1_CI = (confidenceint(8)-confidenceint(7));
            
            T2_map(voxels(voxel)) = thefit.T2;
            IE_map(voxels(voxel)) = thefit.IE;
            rsquare_map(voxels(voxel)) = fit_info.rsquare;
            rmse_map(voxels(voxel)) = fit_info.rmse;
            T1_map(voxels(voxel)) = thefit.T1;
            T2_CI_map(voxels(voxel)) =  T2_CI;
            T1_CI_map(voxels(voxel)) = T1_CI;
            
            %% Remove voxel from list if rsquare is too low
            
            if (fit_info.rsquare < 0.9 || thefit.T2 < 40 || thefit.IE < IE_threshold || isnan(T2_CI) || isnan(T1_CI))     
              voxels(voxel) = 0;
            end
            
end

% Make combinations with remaining voxels
masked_voxels = voxels(find(voxels>0));

Combinations = [];

for i = 1:length(masked_voxels)
Combinations(1+length(Combinations):length(Combinations) + size(combnk(masked_voxels,i),1),1:size(combnk(masked_voxels,i),2)) = combnk(masked_voxels,i);
end


figure,
subplot(3,3,1)
imagesc(T2_map)
title('T2 map')
subplot(3,3,2)
imagesc(IE_map)
title('IE map')
subplot(3,3,3)
imagesc(rsquare_map)
title('rsquare map')
subplot(3,3,4)
imagesc(rmse_map)
title('rmse map')
subplot(3,3,5)
imagesc(T1_map)
title('T1 map')
subplot(3,3,6)
imagesc(T2_CI_map)
title('T2 CI map')
subplot(3,3,7)
imagesc(T1_CI_map)
title('T1 CI map')

end






    




 