%% configure python use for cellpose
clear
SliceList = {'HB004-Slide50-Slice001','HB004-Slide50-Slice002','HB004-Slide50-Slice003','HB004-Slide50-Slice004'};%,'Slice002','Slice003','Slice004'};
MainFolder = repmat({'Z:\data_proc\HB\HB004\coppaFISH'},length(SliceList),1);
model_path = 'G:\Code\training_cellposeDAPI\models\CP_myDAPI_SB';  % where to find the trained model for DAPI segmentation
GenesetF = 'C:\Users\bugeon\Documents\MATLAB'; % where to find the scRNAseq reference dataset

%% manually pick tdtomato cells (if relevant) and adjust DAPI boundaries
GUI_findtdTomato(MainFolder{1})
%

for iSlice = 1:4
    SliceNb = SliceList{iSlice};
    OutputF = MainFolder{iSlice};
    F = fullfile(OutputF,SliceNb,'\output');
    cpt = cellpose(Model=model_path); % requires Medical Imaging Toolbox
    load(fullfile(F,'oCall_spots_OMP'));

    img = imadjust(imread(fullfile(F,'background_image_fixed.tif')));
    img = uint8(img./2^8);
    labels = segmentCells2D(cpt, img, ...
    'ImageCellDiameter', 80, ...
    'Tile', true, ...
    'GPUBatchSize', 1);   % Default is 8; lower this to 1 or 2 to drastically cut VRAM
    DapiBound = convert_cellpose_labels(labels);
    save(fullfile(F,'DAPI_Bound'), 'DapiBound', '-v7.3');
    
    % plot
    B = labeloverlay(img,labels);
    figure;imshow(B);
    
 
    %% rename output folder
    NewOut = fullfile(OutputF);
    o.TileFiles = cellfun(@(x) strrep(x,o.InputDirectory,NewOut), o.TileFiles, 'UniformOutput', false);
    o.OutputDirectory = strrep(o.OutputDirectory,o.InputDirectory,NewOut);
    %% call cells using reference geneSet
    % load('C:\Users\bugeon\Documents\MATLAB\gsetYao.mat') % whole hippocampus and cortex scRNAseq
    load(fullfile(GenesetF,'gsetTasic_123g.mat')) % Tasic cortex scRNAseq
    ExcludeGenes = {'Yjefn3'}; % gene to be excluded from cell calling
    Inefficiency = 0.005; % 0.1 for Yao?; 0.01 for Tasic
       
%     load(fullfile(F,'DAPI_Bound'));
    o.ExcludeGenes = ExcludeGenes;
    CellCalled = DoCellCalling(o,DapiBound,gSet,Inefficiency);
    save(fullfile(F,'Call_cells.mat'), 'CellCalled', '-v7.3');
    
end

%%
% I = img;
% 
% o.ompIntensityThresh = 0.5;
% o.ompIntensityThresh2 = 0.005;
% o.ompNeighbThresh = 18;
% o.ompNeighbThresh2 = 10;
% o.ompScoreThresh = 4.3;
% o.ompScoreThresh2 = 1.1;
% 
% o.MarkerSize = 5;
% o.PlotLineWidth = 1.2;
% o.MarkerType = 'GeneSymbols';
% 
% Roi = round([1, max(o.dpSpotGlobalYX(:,2)), ...
% 1, max(o.dpSpotGlobalYX(:,1))]);
% o.plot(I,Roi,'OMP');
% daspect([1 1 1])
% 
% hold on
% scatter(CellCalled.CellYX(:,2),CellCalled.CellYX(:,1),100,'or','LineWidth',2)
