% Script to process T2-TRIR data
% Author: Koen Baas
% Contact: k.p.baas@amsterdamumc.nl
% Last change: 9 Dec 2021
% Written in MATLAB R2019b
% Copyright 2021, Koen Baas

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
test_nr_of_voxels = 10;   % Nr of voxels evaluated for ROI selection (default = 10)
use_nr_voxels = 4;        % Optional: use a fixed nr of voxels (default = 4)
find_best_comb = true;    % Find combination of voxels (default = true). This ignores 'use_nr_voxels'.
method = 'rsquare';       % Method to evaluate voxel combinations ('rsquare' or 'rmse')
IE_threshold = 0.6;       % Inversion efficiency threshold. Voxels with IE < threshold are not used. (default = 0.6)
fit_separate_TIs = false; % Fit T2 for each TI individually (takes relatively long, default = false)

Analysis_ID = 'Example';  % ID to separate output folders e.g. the current date

Main_dir  = [Main_folder filesep 'Example_data' filesep 'Example_TRUST_TRIR_001' filesep 'PARREC' filesep ''];       % Folder where scans are stored
TRIR_file = 'TRUST_TRIR_001_WIP_TRIR_adiab_ses1_scan1_7_1.PAR';  % File name (.PAR or .dcm)
Mask_file = '';                                                  % Mask file if available (.mat file)

% Other constant parameters
eTEs = [0,40,80,160];                   % Effective echo times (ms)
nr_echos = 60;                          % Nr of TIs in the TRIR sequence
TI1 = 250;                              % TI1 [ms]
TI = 150;                               % delta TI [ms]
TIs = (TI1:TI:(TI1+TI*(nr_echos-1)));   % calculate all TIs

%% Load data
[path, name, ext] = fileparts([Main_dir TRIR_file]);

% Load TRIR data
if strcmp(ext, '.PAR')
    [data,INFO] = readrec_V4_2([Main_dir TRIR_file],'quiet');
elseif strcmp(ext,'.dcm')
    files = dir(fullfile(path, '*.dcm'));
    for k = 1:length(files)
        data(:,:,k) = dicomread([Main_dir files(k).name]);
    end
    info = dicominfo([Main_dir files(k).name]);
    data = double(data);
    
    %Dicom need some specific restructuring
    data = reshape(data, [size(data,1),size(data,2), length(eTEs), nr_echos]);
end

data = squeeze(data);

%Restructure data
if length(size(data)) == 4
    tmp = zeros(size(data,1),size(data,2),size(data,3)*size(data,4));
    tmp(:,:,1:nr_echos)= data(:,:,1,:);
    tmp(:,:,nr_echos+1:nr_echos*2) = data(:,:,2,:);
    tmp(:,:,2*nr_echos+1:nr_echos*3) = data(:,:,3,:);
    tmp(:,:,3*nr_echos+1:nr_echos*4) = data(:,:,4,:);
    data = tmp;
end


% Load large mask (if available)
if ~isempty(Mask_file)
    mask = load(Mask_file);
    mask = mask.mask;
else
    figure, imshow(data(:,:,nr_echos),[])
    set(gcf, 'Position', get(0, 'Screensize'));
    title('Roughly select the sagital sinus')
    mask = roipoly;
    close all
end
save_mask = mask; % Mask for this specific run will be saved as .mat file

%% Find best combination of voxels

if find_best_comb
    
    % Isolate  masked voxels of last TI from eTE = 0 ms series
    masked_slice = data(:,:,nr_echos).*mask;
    
    % Find highest voxels
    highest_values = maxk(reshape(masked_slice,[numel(masked_slice),1]),test_nr_of_voxels);
    incl_voxels = find(ismember(masked_slice, highest_values));
    
    % Find all possible combinations of voxels for optimal fitting
    mask = zeros(80,80); mask(incl_voxels) = 1; % Update mask
    Combinations = TRIR_Find_combinations(mask, data, IE_threshold, nr_echos, TIs); % Find combinations with reliable voxels
    
    % Now find the best combination
    [Best_mask] = TRIR_Find_best_combination(Combinations, data, nr_echos, TIs, method, true);
    
elseif ~find_best_comb
    
    % Or simply use the n best voxels
    selected_voxels = TRIR_Get_best_n_voxels(mask, data, nr_echos, TIs, use_nr_voxels, test_nr_of_voxels, IE_threshold, method, true);
    Best_mask = selected_voxels;
    
end

%% Calculate averages from best ROI
N   =   size(data,3);

for i = 1:4
    for echo = 1:(N/4)
        slice = data(:,:,(i-1)*nr_echos + echo);
        Averages(echo,i) =  mean(slice(Best_mask));
    end
end

%% Fit the data with the selected voxels
x1 = TIs; x1(nr_echos+1:2*nr_echos) = TIs'; x1(2*nr_echos+1:3*nr_echos) = TIs'; x1(3*nr_echos+1:4*nr_echos) = TIs'; x1 = x1';
z1 = Averages(:,1); z1(nr_echos+1:2*nr_echos) = Averages(:,2); z1(2*nr_echos+1:3*nr_echos) = Averages(:,3); z1(3*nr_echos+1:4*nr_echos) = Averages(:,4);
y1(1:nr_echos) = 0; y1(nr_echos+1:2*nr_echos) = 40; y1(2*nr_echos+1:3*nr_echos) = 80; y1(3*nr_echos+1:4*nr_echos) = 160;
Seedpoint = Averages(60,1);

% Fit function
H = @(M0,T2,IE,T1,x,y)(M0*abs(1-(1+exp(-(y/T2))*IE).*exp(-(x/T1))));

% Perform the fit
[thefit,fit_info] = fit([x1,y1'],z1,H,'StartPoint',[Seedpoint,60,0.8,1800],'Upper',[Inf,300,1,3000],'Lower',[Seedpoint*0.8,30,0.2,1400]);

% Determine confidence intervals
confidenceint = confint(thefit);
T2_CI = (confidenceint(4)-confidenceint(3));
T1_CI = (confidenceint(8)-confidenceint(7));
IE_CI = (confidenceint(6)-confidenceint(5));

% Analyze TIs separately
if fit_separate_TIs
    Sep_T2_combinations = TRIR_Inspect_separate_TIs(data, Best_mask, nr_echos, TIs, thefit.T1, thefit.M0, thefit.IE);
end

% Show the selected voxels on top of a TRIR scan
scale = 600000;
tmp_mask = zeros(size(mask)); tmp_mask(incl_voxels) = 1;
Slice = data(:,:,240);
check_mask_later = zeros(size(Slice));
check_mask_later(Best_mask) = 1;
figure, imshow(Slice(:,:),[],'Border','tight')
title('ROI used for fitting (in green)')
pause(0.5)
green = cat(3, zeros(size(check_mask_later)), ones(size(Slice)), zeros(size(Slice)));
hold on
h = imshow(green);
hold off
set(h, 'AlphaData', 0.6.*check_mask_later)

%% Save the output
Output_dir = [Main_dir 'Data_analysis' filesep 'Output' name '_' Analysis_ID];
mkdir(Output_dir)
fileID = fopen([Output_dir filesep 'TRIR_output.txt'],'w');
fprintf(fileID, 'Blood T2 = %3.1f ± %3.1f ms \r\n', thefit.T2, T2_CI);
fprintf(fileID, 'Blood T1 = %3.1f ± %3.1f ms \r\n', thefit.T1, T1_CI);
fprintf(fileID, 'IE = %3.2f ± %3.2f ms \r\n', thefit.IE, IE_CI);
fprintf(fileID, '%1.0f voxel(s) were used for signal fitting \r\n', length(Best_mask));
fclose(fileID);

save([Output_dir filesep 'Initial_mask.mat'], 'save_mask');
save([Output_dir filesep 'Final_ROI.mat'], 'check_mask_later');

display(['Finished processing. T2b = ' num2str(thefit.T2) '. T1b = ' num2str(thefit.T1)]);
display(['Output was saved at: ' Output_dir]);


