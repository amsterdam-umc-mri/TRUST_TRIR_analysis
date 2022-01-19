% Script to process Phase Contrast data
% Author: Koen Baas
% Contact: k.p.baas@amsterdamumc.nl
% Last change: 9 Dec 2021
% Written in MATLAB R2019b
% Copyright 2021, Koen Baas

clear all
clc
close all

%% Getting your folder structure
cd(fileparts(matlab.desktop.editor.getActiveFilename))
cd(['..' filesep]); Scripts_folder = pwd;
cd(['..' filesep]); Main_folder = pwd;
cd(fileparts(matlab.desktop.editor.getActiveFilename))

addpath(genpath(Scripts_folder));

%% Initialize run
Main_dir  = [Main_folder filesep 'Example_data' filesep 'Example_TRUST_TRIR_001' filesep 'PARREC' filesep ''];       % Folder where scans are stored
Phase_contrast_file = 'TRUST_TRIR_001_WIP_QFlow_ses1_9_1.PAR'; % File name (.PAR or .dcm)
Segmentation = [Main_folder filesep 'Example_data' filesep 'Example_TRUST_TRIR_001' filesep 'Masks' filesep 'Seg_001_ses1.nii.gz']; % Optional: corresponding segmentation (if not supplied, a manual ROI can be drawn)
Venc = 80; % Velocity encoding [cm/s] -> needed for PARREC rescaling
enhanced = false; % Option for enhanced dicom

%% Unzip if necessary
[Path,Fileame,ext] = fileparts(Segmentation);

if strcmp(ext,'.gz')
    file_was_zipped = true;
    gunzip (Segmentation)
else
    file_was_zipped = false;
end

%% Load data
[Path,Filename,ext] = fileparts([Main_dir Phase_contrast_file]);

if strcmp(ext, '.dcm')
    
    Velocity_file = [Main_dir '00002.dcm'];
    Velocity_info = dicominfo(fullfile([Main_dir '00002.dcm']));
    Velocity_IM  = dicomread(Velocity_info);
    [tmp, tmp2, file_type] = fileparts(Velocity_file);
    
    Magnitude_info = dicominfo(fullfile([Main_dir '00001.dcm']));
    Magnitude_IM = dicomread(Magnitude_info);
    Magnitude_IM = double(Magnitude_IM);
    
    Velocity_info = dicominfo(fullfile([Main_dir '00003.dcm']));
    Velocity_IM  = dicomread(Velocity_info);
    Velocity_IM = double(Velocity_IM);
    
    file_type = '.dcm';
    
    Velocity_IM = Velocity_IM*Velocity_info.RescaleSlope + Velocity_info.RescaleIntercept;
    Magnitude_IM = Magnitude_IM*Magnitude_info.RescaleSlope + Magnitude_info.RescaleIntercept;
    
elseif strcmp(ext, '.PAR')
    [data, info] = readrec_V4_2([Main_dir Phase_contrast_file],'noscale'); % Phase contrast requires specific rescaling so 'noscale'
    data = squeeze(data);
    
    Velocity_IM = double(data(:,:,2));
    Magnitude_IM = data(:,:,1);
    file_type = '.PAR';
    
    % Rescaling
    Velocity_RS = info.tbl(3, info.tblcols.rescale_slope);
    Velocity_RI = info.tbl(3, info.tblcols.rescale_int);
    Velocity_SS = info.tbl(3, info.tblcols.scale_slope);
    
    Velocity_RS = Velocity_RS / (Venc/(1000*pi));
    Velocity_RI = Velocity_RI / (Venc/(1000*pi));
    
    Magnitude_RS = info.tbl(1, info.tblcols.rescale_slope);
    Magnitude_RI = info.tbl(1, info.tblcols.rescale_int);
    Magnitude_SS = info.tbl(1, info.tblcols.scale_slope);
    
    Velocity_IM = (Velocity_IM*Velocity_RS + Velocity_RI)/ (Velocity_RS * Velocity_SS);
    Magnitude_IM = (Magnitude_IM*Magnitude_RS + Magnitude_RI)/ (Magnitude_RS * Magnitude_SS);
end

%% Dimension in which segmentation is saved varies depending on which image
% was used to make the segmentation in ITK snap.
if isempty(Segmentation)
    manual = true
else
    manual = false;
end

if ~manual
seg_nii = spm_vol(Segmentation);
Seg = spm_read_vols(seg_nii);
if length(size(Seg)) > 2
    if sum(sum(Seg(:,:,1))) > 0
        Seg = Seg(:,:,1);
    elseif sum(sum(Seg(:,:,2))) > 0
        Seg = Seg(:,:,2);
    elseif sum(sum(Seg(:,:,3))) > 0
        Seg = Seg(:,:,3);
    end
    
end
Seg = permute(Seg,[2 1]);
else
   figure, imshow(Velocity_IM,[])
   pause(10) % to zoom in without roipoly bein activated
   Seg = roipoly;
end

% Show overlay for check
figure,imshow(Velocity_IM,[-70,70])
green = cat(3, zeros(size(Velocity_IM)), ones(size(Velocity_IM)), zeros(size(Velocity_IM)));
hold on
h = imshow(green);
hold off
set(h, 'AlphaData', Seg)

%% Calculate average velocity
CC = bwconncomp(Seg); % Only used if multiple vessels are segmented
Sum = 0;

for i = 1:CC.NumObjects
    Pixels = CC.PixelIdxList{i};
    
    for j = 1:length(Pixels)
        Velocity(j) = double(Velocity_IM(Pixels(j)));
    end
    Velocities{i} = Velocity';
    Sum = Sum + sum(Velocity);
    Seperate_sums(i) = sum(Velocity);
    Averages(i) = sum(Velocity)/length(Velocity);
    Sizes(i) = length(Pixels);
    
    clear Velocity
    clear Pixels
end

% Flow-weighted mean velocity:
% (1/sum(V)) * sum (V^2)/N

All_velocities = Velocity_IM(find(Seg == 1));

Flow_weighted_mean = (sum(All_velocities.^2)/sum(All_velocities));

%% Save output and cleanup

% Rename textfile if already exists
if exist([Main_dir 'Blood_velocity.txt'])
    copyfile([Main_dir 'Blood_velocity.txt'], [Main_dir 'Blood_velocity_old.txt'])
    delete([Main_dir 'Blood_velocity.txt'])
end

if (strcmp(file_type, '.dcm') && enhanced)
    fileID = fopen([Main_dir 'Data_analysis' filesep 'Blood_velocity.txt'],'w');
    fprintf(fileID, 'Flow weighted mean = %3.1f cm/s \r\n', Flow_weighted_mean);
    fprintf(fileID, ['Enhanced Dicom format was used in processing \r\n']);
    display(['Output was saved at: ' Main_dir 'Blood_velocity.txt']);
elseif (strcmp(file_type, '.dcm') && ~enhanced)
    fileID = fopen([Main_dir 'Data_analysis' filesep 'Blood_velocity.txt'],'w');
    fprintf(fileID, 'Flow weighted mean = %3.1f cm/s \r\n', Flow_weighted_mean);
    fprintf(fileID, ['Normal Dicom format was used in processing \r\n']);
    display(['Output was saved at: ' Main_dir 'Blood_velocity.txt']);
elseif strcmp(file_type, '.PAR')
    fileID = fopen([Main_dir 'Data_analysis' filesep Filename '_Blood_velocity.txt'],'w');
    fprintf(fileID, 'Flow weighted mean = %3.1f cm/s \r\n', Flow_weighted_mean);
    fprintf(fileID, ['PARREC format was used \r\n']);
    if manual
        fprintf(fileID, ['No segmentation was available so ROI was darwn manually in MATLAB']);
        save([Main_dir 'Data_analysis' filesep Filename '_Manually_drawn_ROI.mat'])
    end
    display(['Output was saved at: ' Main_dir 'Data_analysis' filesep Filename '_Blood_velocity.txt']);
end

fclose(fileID);

% Delete unzipped nii or zip if file was not zipped in the first place
if ~manual
if file_was_zipped
    delete(Segmentation(1:length(Segmentation)-3))
elseif ~file_was_zipped
    gzip(Segmentation, Seg_dir)
    delete(Segmentation)
end
end

display(['Average flow velocity is: ' num2str(Flow_weighted_mean) ' cm/s']);
