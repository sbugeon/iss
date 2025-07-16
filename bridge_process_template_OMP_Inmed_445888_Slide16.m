% cellpose - requires Deep Learning Toolbox™, Computer Vision Toolbox™, and the Medical Imaging Toolbox™
% requires the bfmatlab package ! on a new matlab run: bfInitLogging('DEBUG');
addpath('C:\Users\bugeon\Desktop\bfmatlab')
addpath('C:\Users\bugeon\Documents\GitHub\iss') % main iss code: https://github.com/sbugeon/iss
addpath('C:\Users\bugeon\Documents\GitHub\Transcriptomics-master') % can be found at https://github.com/cortex-lab/Transcriptomics
addpath('C:\Users\bugeon\Documents\GitHub\RegSessPipeline\GUI_Slice_Curation') % https://github.com/sbugeon/RegSessPipeline

clear
% Folders and names
InputF = 'F:\';%'\\NETDATA\Apotome\Stephane\'; % where to find the raw data
OutputF = 'D:\ISS\'; % where to save the iss output
CodebookF = 'G:\Documents\C_doc\codebook_7rounds.txt'; % where to find the codebook for the genes
model_path = 'G:\Code\training_cellposeDAPI\models\CP_myDAPI_SB';  % where to find the trained model for DAPI segmentation
GenesetF = 'C:\Users\bugeon\Documents\MATLAB'; % where to find the scRNAseq reference dataset

An = '445888'; % animal name
SlideN = 'Slide16'; % slide name
SlicePrefix = 'Round-Anchor-'; % file prefix
DoISS = 1; % perform gene calling
DoCellCall = 1; % perform cell calling

MainFolder = fullfile(InputF,An,SlideN);
OutputF = fullfile(OutputF,An,[SlideN,'_test']);

%% Parameters that should be checked before each run
if DoISS
    o = iss_OMP;

    % Channel and round for Anchor, Gad, Gcamp and DAPI
    o.AnchorChannel = 3; % Channel where to find the anchor spots
    o.AnchorRound = 8; % Round where to find the anchor spots
    o.GadChannel = 2; % Channel where to find Gad spots
    o.GadRound = 8; % Round where to find Gad spots
    o.GcampChannel = 2; % Channel where to find Gcamp (or tdTomato) spots
    o.GcampRound = 8; % Round where to find Gcamp (or tdTomato) spots
    o.DapiChannel = 1; % Channel where to find DAPI
    o.DapiRound = 8; % Round where to find DAPI

    % File Names
    % FileBase{r} is the file name of the raw data of round r in o.InputDirectory
    o.FileBase = cell(1,1);
    o.FileBase{1} = strcat(SlicePrefix,'02');
    o.FileBase{2} = strcat(SlicePrefix,'03');
    o.FileBase{3} = strcat(SlicePrefix,'04');
    o.FileBase{4} = strcat(SlicePrefix,'05');
    o.FileBase{5} = strcat(SlicePrefix,'06');
    o.FileBase{6} = strcat(SlicePrefix,'07');
    o.FileBase{7} = strcat(SlicePrefix,'08');
    o.FileBase{8} = strcat(SlicePrefix,'09');

    % Other important parameters to be checked
    o.InitialShiftChannel = 5; % Channel to use to find initial shifts between rounds, needs to contain lots of spots in all rounds
    o.ReferenceRound = o.AnchorRound; % Global coordinate system is built upon o.ReferenceRound and
    o.ReferenceChannel = o.AnchorChannel; % o.ReferenceChannel. If RefRound = AnchorRound, this has to be AnchorChannel.
    o.RawFileExtension = '.czi'; % Format of raw data
    % order of bases in the color code, do not change unless you have misordered the channels on the microsocope
    o.bpLabels = {'0', '1', '2', '3', '4', '5', '6'};

    % these should be the same for all experiments
    o.LogToFile = 0;  % Set to 1 if you want to save command window to txt file, else set to 0.
    o.StripHack = true;
    o.InputDirectory = MainFolder;     % Folder path of raw data
    o.TileSz = 2048;        % Dimension of tile in pixels
    o.nBP = 7;              % Number of Channels
    o.nRounds = 7;          % Number of Imaging Rounds
    o.nExtraRounds = 1;     % Treat Anchor channel as extra round

    % indicate supplementary rounds done in case some z focus was wrong for
    % some rounds, allows to have less tiles than the main file
    % o.FileBase{9} = strcat(SliceNb,'06'); %supp r3
    % o.FileBase{10} = strcat(SliceNb,'11'); %supp r6
    % o.FileBase{11} = strcat(SliceNb,'12'); %supp r6
    % % give the round that was redone for each of the supplementary rounds
    % o.RefRounds = [4 7 7];

    %% Logging
    if o.LogToFile
        if isempty(o.LogFile)
            o.LogFile = fullfile(o.OutputDirectory,'Log.txt');
        end
    end

    o.TileDirectory = fullfile(OutputF,SlicePrefix,'\tiles');
    mkdir(o.TileDirectory);
    o.OutputDirectory = fullfile(OutputF,SlicePrefix,'\output');
    mkdir(o.OutputDirectory);

    %Codebook is a text file containing 2 columns - 1st is the gene name. 2nd is
    %the code, length o.nRounds and containing numbers in the range from 0 to o.nBP-1.
    o.CodeFile = CodebookF;
    %% extract and filter: z-project raw images and filter them
    %parameters
    o.FirstBaseChannel = 1;
    %These specify the dimensions of the filter. R1 should be approximately the
    %radius of the spot and R2 should be double this.
    o.ExtractR1 = 'auto';
    o.ExtractR2 = 'auto';
    o.ExtractScale = 'auto';
    % run code
    o = o.extract_and_filter;
    save(fullfile(o.OutputDirectory,'oExtract'), 'o', '-v7.3');

    % if there were some supplementary rounds, extract them and replace the former tiles with them
    if ~isempty(o.RefRounds)
        o.extract_and_filter_extraTiles(o.RefRounds);
    end
    oOut = o.OutputDirectory;

    % process the different regions (if there are some)
    for j = 1:length(o.TileConnectedID)
        load(fullfile(oOut,'oExtract'));
        o.AnchorChannel = 4; % Channel where to find the anchor spots
        o.ReferenceChannel = o.AnchorChannel;
        o.Graphics=1;
        o.EmptyTiles(:) = 1;
        gg = o.TilePosYX;
        for k = 1:length(o.TileConnectedID{j})
            Idx = gg(o.TileConnectedID{j}(k),:);
            o.EmptyTiles(Idx(1),Idx(2)) = 0;
        end

        o.OutputDirectory = fullfile(oOut,num2str(j));
        mkdir(fullfile(o.OutputDirectory))
        save(fullfile(o.OutputDirectory,'oExtract'), 'o', '-v7.3');
        %% register: roughly stitch overlapping tiles
        % parameters
        o.DetectionRadius = 2;
        o.SmoothSize = 0;
        o.IsolationRadius1 = 4;
        o.IsolationRadius2 = 14;
        o.DetectionThresh = 'auto';
        o.ThreshParam = 5;
        o.MinThresh = 10;
        o.minPeaks = 1;
        o.InitalShiftAutoMinScoreParam = 2;   %a lower value will make it quicker but more likely to fail
        o.DetectSpotsMaxSpots = 10000;
        o.RegMinScore = 'auto';
        o.RegStep = [5,5];
        o.RegSearch.South.Y = -1900:o.RegStep(1):-1700;
        o.RegSearch.South.X = -150:o.RegStep(2):150;
        o.RegSearch.East.Y = -150:o.RegStep(1):150;
        o.RegSearch.East.X = -1900:o.RegStep(2):-1700;
        o.RegWidenSearch = [50,50];

        % If a channel or round is faulty, you can ignore it by selecting only the
        % good ones in o.UseChannels and o.UseRounds.
        o.UseChannels = 1:o.nBP;
        o.UseRounds = 1:o.nRounds;

        % check channels
        % below will flag error if some channels are weak
        %     o = o.check_channels;
        % below will not flag error but remove weak channels automatically.
        % o = o.check_channels(true);

        % run code
        o = o.register2;
        save(fullfile(o.OutputDirectory,'oRegister'), 'o', '-v7.3');
        %% find spots: detect spots in each round/channel and align them to reference (anchor)

        % parameters
        o.FindSpotsMinScore = 'auto';
        o.FindSpotsStep = [5,5];
        % FindSpotsSearch can either be a 1x1 struct or a o.nRounds x 1 cell of
        % structs - have a different range for each round:
        o.FindSpotsSearch = struct();
        o.FindSpotsSearch.Y = -700:o.FindSpotsStep(1):700;
        o.FindSpotsSearch.X = -700:o.FindSpotsStep(2):700;
        % Make WidenSearch larger if you think you have a large shift between rounds
        o.FindSpotsWidenSearch = [250,250];
        o.DetectionThresh = 'auto';
        o.PcDist = 3;
        o.PointCloudMethod = 1;     %1 or 2, set to 2 if no anchor round.
        %2 assumes same scaling to each color channel across all rounds.

        % run code
        o = o.find_spots2;
        save(fullfile(o.OutputDirectory,'oFind_spots'), 'o', '-v7.3');
        %% call spots
        % run code
        o.CallSpotsCodeNorm = 'WholeCode';      % Other alternative is 'Round'
        o = o.call_spots;

        % OMP spot calling
        o.ompInitialNeighbThresh = 5;  % Increase to use less memory. Keep below 10.
        o = o.call_spots_omp;
        save(fullfile(o.OutputDirectory,'oCall_spots_OMP'), 'o', '-v7.3');
       

        %% save large stitched images for DAPI, Gad and Gcamp channels
        MakeDapi_images
         close all
    end
    clear o
end

%%% % to visualize results of gene calling and check diagnostics see
%%% Diagnostics_iss 

%% detect cells using DAPI (cellpose - requires Deep Learning Toolbox™, Computer Vision Toolbox™, and the Medical Imaging Toolbox™)
% F = fullfile(OutputF,SlicePrefix,'\output');
% cpt = cellpose(Model = model_path);
% load(fullfile(F,'oExtract'));
% 
% for j = 1:length(o.TileConnectedID)
%     Fj = fullfile(F,num2str(j));
%     img = imadjust(imread(fullfile(Fj,'background_image_fixed.tif')));
%     labels = segmentCells2D(cpt,img);
%     DapiBound = convert_cellpose_labels(labels);
%     save(fullfile(Fj,'DAPI_Bound'), 'DapiBound', '-v7.3');
%     
%     % plot 
%     B = labeloverlay(img,labels);
%     figure;imshow(B);title(num2str(j));
% end
% 
% %% manually pick tdtomato cells (if relevant) and adjust DAPI boundaries
% F = fullfile(OutputF,SlicePrefix,'\output');
% GUI_findtdTomato(F)
% 
% %% call cells using reference geneSet
% % load('C:\Users\bugeon\Documents\MATLAB\gsetYao.mat') % whole hippocampus and cortex scRNAseq
% load(fullfile(GenesetF,'gsetTasic_123g.mat')) % Tasic cortex scRNAseq
% ExcludeGenes = {'Akr1c18','Yjefn3'}; % gene to be excluded from cell calling
% Inefficiency = 0.01; % 0.1 for Yao?; 0.01 for Tasic
% 
% F = fullfile(OutputF,SlicePrefix,'\output');
% cpt = cellpose(Model=model_path);
% load(fullfile(F,'oExtract'));
% 
% for j = 1%:length(o.TileConnectedID)
%     Fj = fullfile(F,num2str(j));
%     load(fullfile(Fj,'oCall_spots_OMP'));
%     load(fullfile(Fj,'DAPI_Bound'));
%     o.ExcludeGenes = ExcludeGenes;
%     CellCalled = DoCellCalling(o,DapiBound,gSet,Inefficiency);
%     save(fullfile(Fj,'Call_cells.mat'), 'CellCalled', '-v7.3');
% end