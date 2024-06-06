%% Distances between FISH spots - QS 
% Last update: 2022/07/22

% This script calculates 3D distances between mutual nearest neighbor FISH spots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Define MATLAB working directory
folder = cd;
% - Add the path of the folder containing functions
addpath([folder '/functions/']);

% Input data:
% - Make a folder named 'input' in your MATLAB working directory
% - Inside the input folder, deposit folders containing .tif files (channel separated images)

% Output:
% - Image (maximum projection)
% - Image (maximum projection) overlayed with DAPI and FISH segmentations, IDs of nuclei,
% FISH centroids and dashed lines between mutual nearest neighbor FISH spots
% - Summary file: image ID, DAPI ID, DAPI volume, number of FISH spots
% within each DAPI segmented object (nFC1, nFC2 and nFC3 for FISH channel 1, 2 and 3, respectively)
% distance between Mutual Nearest Neighbor FISH spots (MNN_D_FC1FC2, MNN_D_FC1FC3 and MNN_D_FC12FC3),
% IDs of FISH spots (per nucleus) and coordinates (x, y, z) for each FISH channel
% - Parameters file recording the dialog box inputs

%% Input & output paths
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
inputfolder = [folder '/input/'];
outputfolder = [folder '/output/'];
if ~exist(outputfolder, 'dir')
  mkdir(outputfolder);
end
list = dir(inputfolder);
list(strncmp({list.name}, '.', 1)) = [];

%% Dialog boxes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% General parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%
prompt = {'Enter xy pixel size (in µm)',...
          'Enter z-step size (in µm)',...
          'Enter DAPI channel name end',...
          'Enter FISH channel 1 name end',...
          'Enter FISH channel 2 name end',...
          'Enter FISH channel 3 name end', ...
          'Clear borders: type ''yes'' to discard nuclei touching image border'};
title = 'General parameters';
dims = [1 75];
definput = {'0.06', '0.3','C1', 'C2', 'C3', 'C4', 'NA'};
generalP = inputdlg(prompt,title,dims,definput);
if size(generalP,1) == 0
    return
end
xypixelsize = str2double(cell2mat(generalP(1)));
zsize = str2double(cell2mat(generalP(2)));
DAPIName = char(generalP(3));
FC1Name = char(generalP(4));
FC2Name = char(generalP(5));
FC3Name = char(generalP(6));
aclearborders = generalP{7};
clearborders = strcmp(aclearborders, 'yes');
nameGeneralP = {'xy_size'; 'z_size'; 'DAPI_channel'; 'FISH_channel_1'; 'FISH_channel_2'; 'FISH_channel_3'; 'clear_borders'};

% Gaussian filters
%%%%%%%%%%%%%%%%%%%%%%%%%%
    prompt = {'Enter  Gaussian filter sigma for DAPI'; ...
              'Enter Gaussian filter sigma for FISH C1'; ...
              'Enter Gaussian filter sigma for FISH C2'; ...
              'Enter Gaussian filter sigma for FISH C3';};
    title = 'Gaussian filters';
    dims = [1 75];
    definput = {'3'; '1'; '1'; '1'};
    filters = inputdlg(prompt,title,dims,definput);
    if size(filters,1) == 0
        return
    end
    filtDAPIV = str2double(cell2mat(filters(1)));
    filtFC1 = str2double(cell2mat(filters(2)));
    filtFC2 = str2double(cell2mat(filters(3)));
    filtFC3 = str2double(cell2mat(filters(4)));
    nameFilters = {'filter_DAPI'; 'filter_FISH_C1'; 'filter_FISH_C2'; 'filter_FISH_C3'};

% FISH segmentation thresholds
%%%%%%%%%%%%%%%%%%%%%%%%%%
prompt = {'Enter adaptative threshold value for FISH channel 1',...
          'Enter adaptative threshold value for FISH channel 2',...
          'Enter adaptative threshold value for FISH channel 3'};
title = 'Thresholds for FISH segmentation';
dims = [1 75];
definput = {'4', '4', '4'};
thresholds = inputdlg(prompt,title,dims,definput);
if size(thresholds,1) == 0
    return
end
threshVFC1 = str2double(cell2mat(thresholds(1)));
threshVFC2 = str2double(cell2mat(thresholds(2)));
threshVFC3 = str2double(cell2mat(thresholds(3)));
nameThresholds = {'threshold_FISH_C1'; 'threshold_FISH_C2'; 'threshold_FISH_C3'};


% Define min/max segmentation sizes
%%%%%%%%%%%%%%%%%%%%%%%%%%
    voxelvol = xypixelsize * xypixelsize * zsize;
    prompt = {'Enter min DAPI value (min DAPI volume = value/voxel volume)'; ...
              'Enter min FISH C1 value (min FISH C1 volume = value/voxel volume)'; ...
              'Enter max FISH C1 value (max FISH C1 volume = value/voxel volume)'; ...
              'Enter min FISH C2 value (min FISH C2 volume = value/voxel volume)'; ...
              'Enter max FISH C2 value (max FISH C2 volume = value/voxel volume)'; ...
              'Enter min FISH C3 value (min FISH C3 volume = value/voxel volume)'; ...
              'Enter max FISH C3 value (max FISH C3 volume = value/voxel volume)'};

    title = 'Size filtering for segmentation';
    dims = [1 75];
    definput = {'50'; '0.2'; '1.5'; '0.2'; '1.5'; '0.2'; '1.5'; };
    sizeP = inputdlg(prompt,title,dims,definput);
    if size(sizeP,1) == 0
        return
    end
    minvolumeDAPI = str2double(cell2mat(sizeP(1)))/voxelvol;
    minvolumeFC1 = str2double(cell2mat(sizeP(2)))/voxelvol;
    maxvolumeFC1 = str2double(cell2mat(sizeP(3)))/voxelvol;
    minvolumeFC2 = str2double(cell2mat(sizeP(4)))/voxelvol;
    maxvolumeFC2 = str2double(cell2mat(sizeP(5)))/voxelvol;
    minvolumeFC3 = str2double(cell2mat(sizeP(6)))/voxelvol;
    maxvolumeFC3 = str2double(cell2mat(sizeP(7)))/voxelvol;
    nameSizeP = {'min_DAPI_value';...
                 'min_FISH_C1_value'; 'max_FISH_C1_value'; ...
                 'min_FISH_C2_value'; 'max_FISH_C2_value'; ...
                 'min_FISH_C3_value'; 'max_FISH_C3_value'};

% Adjust contrast for output images
%%%%%%%%%%%%%%%%%%%%%%%%%%
    voxelvol = xypixelsize * xypixelsize * zsize;
    prompt = {'Enter low normalized intensity for FISH C1 contrast adjustment'; ...
              'Enter high normalized intensity for FISH C1 contrast adjustment'; ...
              'Enter low normalized intensity for FISH C2 contrast adjustment'; ...
              'Enter high normalized intensity for FISH C2 contrast adjustment'; ...
              'Enter low normalized intensity for FISH C3 contrast adjustment'; ...
              'Enter high normalized intensity for FISH C3 contrast adjustment'};

    title = 'Contrast adjustment fot output images';
    dims = [1 75];
    definput = {'0'; '0.75'; '0'; '0.75'; '0'; '0.75'};
    adjContrast = inputdlg(prompt,title,dims,definput);
    if size(adjContrast,1) == 0
        return
    end
    lowFC1 = str2double(cell2mat(adjContrast(1)));
    highFC1 = str2double(cell2mat(adjContrast(2)));
    lowFC2 = str2double(cell2mat(adjContrast(3)));
    highFC2 = str2double(cell2mat(adjContrast(4)));
    lowFC3 = str2double(cell2mat(adjContrast(5)));
    highFC3 = str2double(cell2mat(adjContrast(6)));
    nameAdjContrast = {'low_FISH_C1'; 'high_FISH_C1'; ...
                       'low_FISH_C2'; 'high_FISH_C2'; ...
                       'low_FISH_C3'; 'high_FISH_C3'};

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Summary parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%
parameter = vertcat(nameGeneralP, nameFilters, nameThresholds, nameSizeP, nameAdjContrast);
value = vertcat(generalP, filters, thresholds, sizeP, adjContrast);
parameters = table(parameter, value);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
for k = 1:length(list)
%% Define paths, files and create summary table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputsubfolder = [inputfolder list(k).name];
outputsubfolder = [outputfolder list(k).name];
if ~exist(outputsubfolder, 'dir')
  mkdir(outputsubfolder);
end

filepatternDAPI=fullfile(inputsubfolder, ['*' DAPIName '.tif']);
filesDAPI=dir(filepatternDAPI);
filepatternFC1=fullfile(inputsubfolder, ['*' FC1Name '.tif']);
filesFC1=dir(filepatternFC1);
filepatternFC2=fullfile(inputsubfolder, ['*' FC2Name '.tif']);
filesFC2=dir(filepatternFC2);
filepatternFC3=fullfile(inputsubfolder, ['*' FC3Name '.tif']);
filesFC3=dir(filepatternFC3);

nfiles=length(filesDAPI);
summary = table; summaryname = 'summary.csv';

%%
for i=1:nfiles
%% Import and read images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	inputnameDAPI=filesDAPI(i).name;
	inputsubfoldernameDAPI=fullfile(inputsubfolder, inputnameDAPI);
	inputnameFC1=filesFC1(i).name;
	inputsubfoldernameFC1=fullfile(inputsubfolder, inputnameFC1);
	inputnameFC2=filesFC2(i).name;
	inputsubfoldernameFC2=fullfile(inputsubfolder, inputnameFC2);
	inputnameFC3=filesFC3(i).name;
	inputsubfoldernameFC3=fullfile(inputsubfolder, inputnameFC3);

	imdata = imfinfo(inputsubfoldernameDAPI);
	width = [imdata.Width];
	nxpixels = width(1);
	height = [imdata.Height];
    nypixels = height(1);
	nzslices = numel(imdata);
     
    % 3D image
	im3dDAPI = zeros(nypixels, nxpixels, nzslices);
	im3dFC1 = zeros(nypixels, nxpixels, nzslices);
	im3dFC2 = zeros(nypixels, nxpixels, nzslices);
	im3dFC3 = zeros(nypixels, nxpixels, nzslices);

	for zslice = 1 : nzslices

        imslice = imread(inputsubfoldernameDAPI, zslice);
        im3dDAPI(:,:,zslice) = imslice;

        imslice = imread(inputsubfoldernameFC1, zslice);
        im3dFC1(:,:,zslice) = imslice;

        imslice = imread(inputsubfoldernameFC2, zslice);
        im3dFC2(:,:,zslice) = imslice;

        imslice = imread(inputsubfoldernameFC3, zslice);
        im3dFC3(:,:,zslice) = imslice;

	end

%% Image filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    imfDAPI = imgaussfilt3(im3dDAPI,filtDAPIV);
	imDAPI =uint16(imfDAPI); 
	projimDAPI = max(mat2gray(imfDAPI), [], 3);

    imfFC1 = imgaussfilt3(im3dFC1,filtFC1);
    imFC1 =uint16(imfFC1); 
    projimFC1 = max(mat2gray(imfFC1), [], 3);

    imfFC2 = imgaussfilt3(im3dFC2,filtFC2);
    imFC2 =uint16(imfFC2); 
    projimFC2 = max(mat2gray(imfFC2), [], 3);

    imfFC3 = imgaussfilt3(im3dFC3,filtFC3);
    imFC3 =uint16(imfFC3); 
    projimFC3 = max(mat2gray(imfFC3), [], 3);

%% DAPI segmentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	threshDAPI = graythresh(imDAPI);
	maskDAPI = imbinarize(imDAPI,threshDAPI); 

    connectedcDAPI = bwconncomp(maskDAPI);
    prestatsDAPI = regionprops3(connectedcDAPI, 'Volume');
    idxDAPI = find([prestatsDAPI.Volume] > minvolumeDAPI);    
    mask2DAPI = ismember(labelmatrix(connectedcDAPI), idxDAPI);
    mask2DAPI = imfill(mask2DAPI, 'holes');
    if clearborders == 1
         mask2DAPI = imclearborder(mask2DAPI);
    end
    projmask2DAPI = max(mask2DAPI, [], 3);
    statsDAPI = regionprops3(mask2DAPI, 'Volume', 'Centroid');
    
%% Single nucleus analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	nuclei = struct;
    connectedciDAPI = bwconncomp(mask2DAPI);
	for iDAPI = 1:size(statsDAPI,1)
        % DAPI
    	nuclei(iDAPI).nucleus_ID = iDAPI;
        nucleusMask = ismember(labelmatrix(connectedciDAPI), iDAPI);
        nuclei(iDAPI).nucleusProjMask = max(nucleusMask, [], 3);
        nuclei(iDAPI).statsDAPI = regionprops3(nucleusMask, 'Volume', 'Centroid');

        % FISH Channel 1
        imROIFC1 = imFC1; imROIFC1(~nucleusMask) = 0;
        intMaskFC1 = regionprops3(nucleusMask, imROIFC1, 'VoxelValues');
        allIntMaskFC1 = double(intMaskFC1.VoxelValues{1});
        threshFC1 = mean(allIntMaskFC1) + threshVFC1*std(allIntMaskFC1);
        mergeROIFC1 = imROIFC1 > threshFC1;        
        connectedcFC1 = bwconncomp(mergeROIFC1);
        prestatsFC1 = regionprops3(connectedcFC1, 'Volume');
        idxFC1 = find([prestatsFC1.Volume] > minvolumeFC1 & [prestatsFC1.Volume] < maxvolumeFC1);    
        nucleusMaskFC1 = ismember(labelmatrix(connectedcFC1), idxFC1);
        nuclei(iDAPI).nucleusProjMaskFC1 = max(nucleusMaskFC1, [], 3);
        nuclei(iDAPI).statsFC1 = regionprops3(nucleusMaskFC1, 'Volume', 'Centroid');
        nuclei(iDAPI).nFC1 = size(nuclei(iDAPI).statsFC1,1);

        % FISH Channel 2
        imROIFC2 = imFC2; imROIFC2(~nucleusMask) = 0;
        intMaskFC2 = regionprops3(nucleusMask, imROIFC2, 'VoxelValues');
        allIntMaskFC2 = double(intMaskFC2.VoxelValues{1});
        threshFC2 = mean(allIntMaskFC2) + threshVFC2*std(allIntMaskFC2);
        mergeROIFC2 = imROIFC2 > threshFC2;        
        connectedcFC2 = bwconncomp(mergeROIFC2);
        prestatsFC2 = regionprops3(connectedcFC2, 'Volume');
        idxFC2 = find([prestatsFC2.Volume] > minvolumeFC2 & [prestatsFC2.Volume] < maxvolumeFC2);    
        nucleusMaskFC2 = ismember(labelmatrix(connectedcFC2), idxFC2);
        nuclei(iDAPI).nucleusProjMaskFC2 = max(nucleusMaskFC2, [], 3);
        nuclei(iDAPI).statsFC2 = regionprops3(nucleusMaskFC2, 'Volume', 'Centroid');
        nuclei(iDAPI).nFC2 = size(nuclei(iDAPI).statsFC2,1);
        
        % FISH Channel 3
        imROIFC3 = imFC3; imROIFC3(~nucleusMask) = 0;
        intMaskFC3 = regionprops3(nucleusMask, imROIFC3, 'VoxelValues');
        allIntMaskFC3 = double(intMaskFC3.VoxelValues{1});
        threshFC3 = mean(allIntMaskFC3) + threshVFC3*std(allIntMaskFC3);
        mergeROIFC3 = imROIFC3 > threshFC3;        
        connectedcFC3 = bwconncomp(mergeROIFC3);
        prestatsFC3 = regionprops3(connectedcFC3, 'Volume');
        idxFC3 = find([prestatsFC3.Volume] > minvolumeFC3 & [prestatsFC3.Volume] < maxvolumeFC3);    
        nucleusMaskFC3 = ismember(labelmatrix(connectedcFC3), idxFC3);
        nuclei(iDAPI).nucleusProjMaskFC3 = max(nucleusMaskFC3, [], 3);
        nuclei(iDAPI).statsFC3 = regionprops3(nucleusMaskFC3, 'Volume', 'Centroid');
        nuclei(iDAPI).nFC3 = size(nuclei(iDAPI).statsFC3,1);
        
        % Distances between C2 and C3
        if nuclei(iDAPI).nFC1 == 0 || nuclei(iDAPI).nFC2 == 0 || nuclei(iDAPI).nFC3 == 0
           nuclei(iDAPI).distancesFC1FC2FC3 = NaN(1,3);
           nuclei(iDAPI).IDsFC1FC2FC3 = NaN(1,3);
           nuclei(iDAPI).coordFC1 = NaN(1,3);
           nuclei(iDAPI).coordFC2 = NaN(1,3);
           nuclei(iDAPI).coordFC3 = NaN(1,3);
        else
           [nuclei(iDAPI).distancesFC1FC2FC3, nuclei(iDAPI).IDsFC1FC2FC3] = f_3D_distances_3cols_QS_01(xypixelsize, zsize,...
            nuclei(iDAPI).statsFC1.Centroid, nuclei(iDAPI).statsFC2.Centroid, nuclei(iDAPI).statsFC3.Centroid);    
            nuclei(iDAPI).coordFC1 = nuclei(iDAPI).statsFC1.Centroid(nuclei(iDAPI).IDsFC1FC2FC3(:,1), :);
            nuclei(iDAPI).coordFC2 = nuclei(iDAPI).statsFC2.Centroid(nuclei(iDAPI).IDsFC1FC2FC3(:,2), :);
            nuclei(iDAPI).coordFC3 = nuclei(iDAPI).statsFC3.Centroid(nuclei(iDAPI).IDsFC1FC2FC3(:,3), :);
        end
    end

%% Summary file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    summaryimage = [];
    for iDAPI = 1:size(nuclei,2)
            summaryimage = [summaryimage; ...
                             repmat(nuclei(iDAPI).nucleus_ID, size(nuclei(iDAPI).distancesFC1FC2FC3,1),1), ...
                             repmat(nuclei(iDAPI).statsDAPI.Volume*voxelvol, size(nuclei(iDAPI).distancesFC1FC2FC3,1),1), ...
                             repmat(nuclei(iDAPI).nFC1, size(nuclei(iDAPI).distancesFC1FC2FC3,1),1), ...
                             repmat(nuclei(iDAPI).nFC2, size(nuclei(iDAPI).distancesFC1FC2FC3,1),1), ...
                             repmat(nuclei(iDAPI).nFC3, size(nuclei(iDAPI).distancesFC1FC2FC3,1),1), ...
                             nuclei(iDAPI).distancesFC1FC2FC3, nuclei(iDAPI).IDsFC1FC2FC3, nuclei(iDAPI).coordFC1, nuclei(iDAPI).coordFC2, nuclei(iDAPI).coordFC3];
    end
    summaryimage = array2table(summaryimage);
    imagename = string(repmat(strrep(inputnameDAPI, '_C1.tif', ''), size(summaryimage,1), 1));
    imagename = array2table(imagename);

   % summaryimage = table(imagename, summarynuclei(:,1), summarynuclei(:,2), summarynuclei(:,3),  summarynuclei(:,4), summarynuclei(:,5));
    summaryimage = [imagename, summaryimage];

    summaryimage.Properties.VariableNames= {'image_ID', 'DAPI_ID', 'DAPI_volume',...
        'nFC1', 'nFC2', 'nFC3', 'MNN_D_FC1FC2', 'MNN_D_FC1FC3', 'MNN_D_FC2FC3', 'ID_FC1', 'ID_FC2', 'ID_FC3', 'xFC1', 'yFC1', 'zFC1', 'xFC2', 'yFC2', 'zFC2', 'xFC3', 'yFC3', 'zFC3'};
    summary = [summary; summaryimage];
    

%% Images FISH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mergeimage = cat(3, projimFC2, projimFC1, projimFC3); % RGB
        I = imadjust(mergeimage,[lowFC1 lowFC2 lowFC3 ; highFC1 highFC2 highFC3]); % Adjust contrast
        lineWidth = .5;
        centroidSize = 6;
        FtSz = 4;
%% Image
        fig = figure('pos',[800 200 400 400]);
        set(gcf,'color','w'); set(gca,'xtick',[], 'ytick', [])
        ax = gca; outerpos = ax.OuterPosition; ti = ax.TightInset; 
        left = outerpos(1) + ti(1); bottom = outerpos(2) + ti(2);
        ax_width = outerpos(3) - ti(1) - ti(3); ax_height = outerpos(4) - ti(2) - ti(4);
        ax.Position = [left bottom ax_width ax_height];
        imshow(I)
        imname = [strrep(inputnameDAPI, ['_' DAPIName '.tif'],''), '.tif']; 
        print(fullfile(outputsubfolder, imname), '-dtiff', '-r300');
        close(fig);
      
%% Image + processing overlay
        fig = figure('pos',[800 200 400 400]);
        set(gcf,'color','w'); set(gca,'xtick',[], 'ytick', [])
        ax = gca; outerpos = ax.OuterPosition; ti = ax.TightInset; 
        left = outerpos(1) + ti(1); bottom = outerpos(2) + ti(2);
        ax_width = outerpos(3) - ti(1) - ti(3); ax_height = outerpos(4) - ti(2) - ti(4);
        ax.Position = [left bottom ax_width ax_height];
      
        imshow(I);
        hold on
        for iDAPI = 1:size(nuclei, 2)
            boundDAPI = bwboundaries(nuclei(iDAPI).nucleusProjMask);
            for lbDAPI = 1:length(boundDAPI)
                boundary = boundDAPI{lbDAPI};
                if ~isnan(nuclei(iDAPI).IDsFC1FC2FC3)
                plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', lineWidth/2)
                else
                plot(boundary(:,2), boundary(:,1), ':y', 'LineWidth', lineWidth/2)
                end
            end   
            text(nuclei(iDAPI).statsDAPI.Centroid(1), nuclei(iDAPI).statsDAPI.Centroid(2),...
            string(nuclei(iDAPI).nucleus_ID), 'Color', 'w', 'FontSize', FtSz + 2)
            
            boundFC1 = bwboundaries(nuclei(iDAPI).nucleusProjMaskFC1);
            for lb = 1:length(boundFC1)
                boundary = boundFC1{lb};
                plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', lineWidth/2)
            end   
            boundFC2 = bwboundaries(nuclei(iDAPI).nucleusProjMaskFC2);
            for lb = 1:length(boundFC2)
                boundary = boundFC2{lb};
                plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', lineWidth/2)
            end
            boundFC3 = bwboundaries(nuclei(iDAPI).nucleusProjMaskFC3);
            for lb = 1:length(boundFC3)
                boundary = boundFC3{lb};
                plot(boundary(:,2), boundary(:,1), 'b', 'LineWidth', lineWidth/2)
            end

            % Optional: plot all centroids 
            %{
            for nFC1 = 1:size(nuclei(iDAPI).statsFC1.Centroid,1)
                plot(nuclei(iDAPI).statsFC1.Centroid(nFC1,1), nuclei(iDAPI).statsFC1.Centroid(nFC1,2), 'g+', 'LineWidth',lineWidth/2, 'MarkerSize', centroidSize/2)
            end
            for nFC2 = 1:size(nuclei(iDAPI).statsFC2.Centroid,1)
                plot(nuclei(iDAPI).statsFC2.Centroid(nFC2,1), nuclei(iDAPI).statsFC2.Centroid(nFC2,2), 'r+', 'LineWidth', lineWidth/2, 'MarkerSize', centroidSize/2)
            end
            for nFC3 = 1:size(nuclei(iDAPI).statsFC3.Centroid,1)
                plot(nuclei(iDAPI).statsFC3.Centroid(nFC3,1), nuclei(iDAPI).statsFC3.Centroid(nFC3,2), 'b+', 'LineWidth', lineWidth/2, 'MarkerSize', centroidSize/2)
            end
            %}
            
        end

        % Plots MNN centroids
        plot(summaryimage.xFC1,summaryimage.yFC1, 'g+', 'LineWidth', lineWidth, 'MarkerSize', centroidSize)
        plot(summaryimage.xFC2,summaryimage.yFC2, 'r+', 'LineWidth', lineWidth, 'MarkerSize', centroidSize)
        plot(summaryimage.xFC3,summaryimage.yFC3, 'b+', 'LineWidth', lineWidth, 'MarkerSize', centroidSize)


        filtDAPI = summaryimage(~isnan(summaryimage.xFC1),:);
        for iFiltD = 1:size(filtDAPI,1)

                line([filtDAPI.xFC1(iFiltD) filtDAPI.xFC2(iFiltD)], [filtDAPI.yFC1(iFiltD) filtDAPI.yFC2(iFiltD)],...
                'Color', 'w', 'LineStyle', ':', 'LineWidth', lineWidth/2)  

                line([filtDAPI.xFC1(iFiltD) filtDAPI.xFC3(iFiltD)], [filtDAPI.yFC1(iFiltD) filtDAPI.yFC3(iFiltD)],...
                'Color', 'w', 'LineStyle', ':', 'LineWidth', lineWidth/2)  

                line([filtDAPI.xFC2(iFiltD) filtDAPI.xFC3(iFiltD)], [filtDAPI.yFC2(iFiltD) filtDAPI.yFC3(iFiltD)],...
                'Color', 'w', 'LineStyle', ':', 'LineWidth', lineWidth/2)  

                % Optional: write FISH IDs
                %{
                text(filtDAPI.xFC1(iFiltD)+2, filtDAPI.yFC1(iFiltD)-2, num2str(filtDAPI.ID_FC1(iFiltD)), 'FontSize', FtSz, 'Color', 'g')
                text(filtDAPI.xFC2(iFiltD)+2, filtDAPI.yFC2(iFiltD)-2, num2str(filtDAPI.ID_FC2(iFiltD)), 'FontSize', FtSz, 'Color', 'r')
                text(filtDAPI.xFC3(iFiltD)+2, filtDAPI.yFC3(iFiltD)-2, num2str(filtDAPI.ID_FC3(iFiltD)), 'FontSize', FtSz, 'Color', 'b')
                %}

        end

        imname = [strrep(inputnameDAPI, ['_' DAPIName '.tif'],'_MNN'), '.tif']; 
        print(fullfile(outputsubfolder, imname), '-dtiff', '-r300')
        close(fig);                    
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
end        

%% Save results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(summary,1) > 0
    outputsummaryname = fullfile(outputsubfolder, summaryname);
	writetable(summary, outputsummaryname);     
end
	outputparametersname = fullfile(outputsubfolder, 'parameters.txt');
	writetable(parameters, outputparametersname, 'WriteVariableNames', false,'Delimiter','\t');   


end

close all



