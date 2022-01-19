% Script to process TRUST data
% Author: Koen Baas
% Contact: k.p.baas@amsterdamumc.nl
% Last change: 9 Dec 2021
% Written in MATLAB R2019b
% Copyright 2021, Koen Baas
% Based on Lu et al 2008 (10.1002/mrm.21627)

clear all
close all
clc

%% Getting your folder structure
cd(fileparts(matlab.desktop.editor.getActiveFilename))
cd(['..' filesep]); Scripts_folder = pwd;
cd(['..' filesep]); Main_folder = pwd;
cd(fileparts(matlab.desktop.editor.getActiveFilename))

addpath(genpath(Scripts_folder));

%% Initialize run
use_nr_voxels = 4;       % Nr of voxels used for spatial averaging (default = 4)
Analysis_ID = 'Example'; % ID to separate output folders e.g. the current date

Main_dir = [Main_folder filesep 'Example_data' filesep 'Example_TRUST_TRIR_001' filesep 'PARREC' filesep ];   % Folder where scans are stored
TRUST_file = 'TRUST_TRIR_001_WIP_TRUST_ses1_scan1_3_1.PAR';      % File name (.PAR or .dcm)
Mask_file = '';                                                  % Mask file if available (.mat file)

% Constants
T1b = 1624;              % Blood T1 (ms)
eTEs = [0,40,80,160];    % Effective echo times (ms)


%% Load the data
[path, name, ext] = fileparts([Main_dir TRUST_file]);

% Load TRUST scan
if strcmp(ext, '.PAR')
    [data,INFO] = readrec_V4_2([Main_dir TRUST_file],'quiet');
elseif strcmp(ext,'.dcm')
    files = dir(fullfile(path, '*.dcm'));
    for k = 1:length(files)
        data(:,:,k) = dicomread([Main_dir files(k).name]);
    end
    info = dicominfo([Main_dir files(k).name]);
end

data = squeeze(data);

% Load large mask (if available)
if ~isempty(Mask_file)
    mask = load(Mask_file);
    mask = mask.mask;
else
    mean_eTE0 = mean(data(:,:,2:2:7),3) - mean(data(:,:,1:2:6),3);
    figure, imshow(mean_eTE0,[])
    set(gcf, 'Position', get(0, 'Screensize'));
    title('Roughly select the sagital sinus')
    mask = roipoly;
    close all
end

%% Subtract label images from control images
x_dim = size(data,1);
y_dim = size(data,2);
z_dim = size(data,3);

Av_cont_im(:,:,1) = mean(data(:,:,2:2:6),3);
Av_cont_im(:,:,2) = mean(data(:,:,8:2:12),3);
Av_cont_im(:,:,3) = mean(data(:,:,14:2:18),3);
Av_cont_im(:,:,4) = mean(data(:,:,20:2:24),3);

Av_lab_im(:,:,1) = mean(data(:,:,1:2:5),3);
Av_lab_im(:,:,2) = mean(data(:,:,7:2:11),3);
Av_lab_im(:,:,3) = mean(data(:,:,13:2:17),3);
Av_lab_im(:,:,4) = mean(data(:,:,19:2:23),3);

Av_diff_im = Av_cont_im - Av_lab_im;
masked_diff_im = Av_diff_im(:,:,1) .* mask(:,:);

%% Select voxels with largest difference signal at eTE = 0
highest_values = maxk(reshape(masked_diff_im(:,:,1),x_dim*y_dim,1),use_nr_voxels);

% Prevent crash if two voxels have the exact same value
for i = 1:use_nr_voxels
    if length(find(masked_diff_im == highest_values(i))) == 1
        incl_voxels(i) =  find(masked_diff_im == highest_values(i));
    else
        voxels_with_same_values = find(masked_diff_im == highest_values(i));
        
        for j = 1:length(voxels_with_same_values)
            
            if mask(voxels_with_same_values(j)) == 1
                incl_voxels(i) = masked_diff_im(voxels_with_same_values(j));
            end
        end
    end
end

% Create difference image with only the selected voxels
mask = zeros(size(mask));
mask(incl_voxels) = 1;
for i = 1:length(eTEs)
   masked_diff_im(:,:,i) =  Av_diff_im(:,:,i) .* mask;
end

%% Calculate average of masked voxels
S1 = nanmean(nonzeros(masked_diff_im(:,:,1)));
S2 = nanmean(nonzeros(masked_diff_im(:,:,2)));
S3 = nanmean(nonzeros(masked_diff_im(:,:,3)));
S4 = nanmean(nonzeros(masked_diff_im(:,:,4)));

%% Plot relaxation curve
figure, plot(eTEs,[S1 S2 S3 S4],'o')
title('Fit through T2 decay');
xlabel('TI [ms]')
ylabel('Magnetization [a.u.]');

%Copy values for fit
x  = eTEs';
y  = [S1 S2 S3 S4]';

%% Fit data to 3-parameter model (Li et al. 2017)
fo = fitoptions('Method','NonlinearLeastSquares', ...
    'Lower',[0,0], ...
    'Upper',[Inf,Inf], ...
    'StartPoint',[1000,80]);

% a = initial magnitude
% b = C -> T2b = 1 / ((1/T1b)-C)

TRUST_model                = fittype('a*exp(-x/b)','options',fo);


[fitted_parms,fit_info] = fit(x,y,TRUST_model);
fit_parms = coeffvalues(fitted_parms);
Confidence = confint(fitted_parms);


%% Calculate blood T2
C = -1/fit_parms(2);
T2b = 1/((1/T1b)-C);

CI_T2b = (1/((1/T1b)-(-1/Confidence(2,2)))) - (1/((1/T1b)-(-1/Confidence(1,2))));


%% Plot signal values and fitted curve
close all
figure, plot(fitted_parms,x,[S1 S2 S3 S4],'o');
xticks([0 40 80 120 160]);
title(['Signal and fitted curves -> T2b = ' num2str(T2b),' ï¿½ ', num2str(CI_T2b)]);
set(gca,'FontSize',20)
title('Fit through T2 decay');
display(['T2b = ' num2str(T2b),' ± ', num2str(CI_T2b)])
xlabel('eTE [ms]')
ylabel('Magnetization [a.u.]');
axis([0 160 0 1.1*S1])

%Show data with mask overlayed
tmp_mask = zeros(size(mask)); tmp_mask(incl_voxels) = 1;
Slice = Av_diff_im(:,:,1);
scale = 0.3*max(Slice(:));
check_mask_later = zeros(size(Slice));
check_mask_later(incl_voxels) = 1;
fig = figure;
imshow(Slice(:,:),[0 scale],'Border','tight')
title('ROI used for fitting (in green)')
pause(0.5)
green = cat(3, zeros(size(check_mask_later)), ones(size(Slice)), zeros(size(Slice)));
hold on
h = imshow(green);
hold off
set(h, 'AlphaData', 0.6.*check_mask_later)


%% Save the output
Output_dir = [Main_dir 'Data_analysis' filesep 'Output_' name '_' Analysis_ID];
mkdir(Output_dir)

fileID = fopen([Output_dir filesep 'TRUST_output.txt'],'w');
fprintf(fileID, 'Blood T2 = %3.1f ± %3.1f ms \r\n', T2b, CI_T2b);
fprintf(fileID, '%1.0f voxels were used for signal fitting \r\n', use_nr_voxels);
fprintf(fileID, 'T1b = %4.0f ms was assumed', T1b);
fclose(fileID);

save([Output_dir filesep 'Final_ROI.mat'], 'mask');
display(['Output was saved at: ' Output_dir]);