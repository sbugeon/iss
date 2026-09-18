clear
SliceList = {'humansection-20260901','humansection-20260901-2'};
MainFolder = repmat({'D:\bugeon\data\Thomas\Thomas-20260901'},length(SliceList),1);
for iSlice = 2
    SliceNb = SliceList{iSlice};
    
    %     %% Parameters that should be checked before each run
    o = iss_OMP;
    o.AnchorChannel = 6;
    o.AnchorRound = 8;            %Channel that has most spots in o.AnchorRound
    o.GadChannel = 4;
    o.GadRound = 8;
    o.GcampChannel = 5;
    o.GcampRound = 8;
    o.DapiChannel = 1;
    o.DapiRound = 8;             %Channel in o.AnchorRound that contains Dapi images
    
    o.InitialShiftChannel = 5;      %Channel to use to find initial shifts between rounds
    o.ReferenceRound = o.AnchorRound;           %Global coordinate system is built upon o.ReferenceRound and
    o.ReferenceChannel = o.AnchorChannel;         %o.ReferenceChannel. If RefRound = AnchorRound, this has to be AnchorChannel.
    o.RawFileExtension = '.nd2';    %Format of raw data
    o.LogToFile = 0;                %Set to 1 if you want to save command window to txt file, else set to 0.
    o.StripHack = true;
    
    o.TileInitialPosYX = [3,3;3,2;3,1;2,1;2,2;2,3;1,3;1,2;1,1]; % hack for this experiment tile positions
    %% File Names
    %CHECK BEFORE EACH RUN
    o.InputDirectory = MainFolder{iSlice};     %Folder path of raw data
    
    o.TileSz = 3200;        %Dimension of tile in pixels
    o.nBP = 7;              %Number of Channels
    o.nRounds = 7;          %Number of Imaging Rounds
    o.nExtraRounds = 1;     %Treat Anchor channel as extra round
    
    %FileBase{r} is the file name of the raw data of round r in o.InputDirectory
    o.FileBase = cell(1,1);
    o.FileBase{1} = strcat(SliceNb,'R0');
    o.FileBase{2} = strcat(SliceNb,'R1');
    o.FileBase{3} = strcat(SliceNb,'R2');
    o.FileBase{4} = strcat(SliceNb,'R3');
    o.FileBase{5} = strcat(SliceNb,'R4');
    o.FileBase{6} = strcat(SliceNb,'R5');
    o.FileBase{7} = strcat(SliceNb,'R6');
    o.FileBase{8} = strcat(SliceNb,'anchordapi');

    o.TileDirectory = fullfile(o.InputDirectory,SliceNb,'\tiles');
    mkdir(o.TileDirectory);
    o.OutputDirectory = fullfile(o.InputDirectory,SliceNb,'\output');
    mkdir(o.OutputDirectory);
    
    o.RawFileExtension = '.nd2';
    %Codebook is a text file containing 2 columns - 1st is the gene name. 2nd is
    %the code, length o.nRounds and containing numbers in the range from 0 to o.nBP-1.
%     o.CodeFile = 'C:\Users\bugeon\Documents\data_coppaFISH\codebook_73g_ctx.txt';
     o.CodeFile = 'C:\Users\bugeon\Documents\data_coppaFISH\codebook_BaptisteHumanChimp.txt';%%%%%%%%%%%%
    %% Logging
    if o.LogToFile
        if isempty(o.LogFile)
            o.LogFile = fullfile(o.OutputDirectory,'Log.txt');
        end
    end
    %% extract and filter
    
    %parameters
    o.FirstBaseChannel = 1;
%     o.bpLabels = {'0', '2', '1','3', '4','5','6'}; %order of bases %%%%%%%%%%%%%%%
    o.bpLabels = {'0', '1', '2','3', '4','5','6'}; %order of bases %%%%%%%%%%%%%%%
    %These specify the dimensions of the filter. R1 should be approximately the
    %radius of the spot and R2 should be double this.
    o.ExtractR1 = 'auto';
    o.ExtractR2 = 'auto';
    
    %     o.ExtractScale = 2.8258;
    o.ExtractScale = 'auto';
    o.TilePixelValueShift = 15000;
    
    %Max time (seconds) to wait for raw .nd2 files to be obtained
    o.MaxWaitTime1 = 60;      %Less time for round 1 incase name is wrong
    o.MaxWaitTime = 21600;
    
    %     o.ExtractScale = 12.9230;
    %run code
    o = o.extract_and_filter;
    
    save(fullfile(o.OutputDirectory, 'oExtract'), 'o', '-v7.3');
    
    %% 
%     rename_tiles(o)
    %% register
    %o.AutoThresh(:,o.AnchorChannel,o.AnchorRound) = o.AutoThresh(:,o.AnchorChannel,o.AnchorRound)*0.25;     %As Anchor Threshold seemed too high
    %parameters
    %Anchor spots are detected in register2
    o.DetectionRadius = 2;
    o.SmoothSize = 0;
    o.IsolationRadius1 = 4;
    o.IsolationRadius2 = 14;
    
    o.DetectionThresh = 'auto';
    o.ThreshParam = 5;
    o.MinThresh = 10;
    o.minPeaks = 1;
    o.InitalShiftAutoMinScoreParam=2;   %a lower value will make it quicker but more likely to fail
    
    %paramaters to find shifts between overlapping tiles
    o.RegMinScore = 'auto';
    o.RegStep = [5,5];
 o.RegSearch.South.Y = -3000:o.RegStep(1):-2800;
    o.RegSearch.South.X = -150:o.RegStep(2):150;
    o.RegSearch.East.Y = -150:o.RegStep(1):150;
    o.RegSearch.East.X = -3000:o.RegStep(2):-2800;

    o.RegWidenSearch = [50,50];
    
    %If a channel or round is faulty, you can ignore it by selecting only the
    %good ones in o.UseChannels and o.UseRounds.
    o.UseChannels = 1:o.nBP;
    o.UseRounds = 1:5;
    
    % check channels
    % below will flag error if some channels are weak
%     o = o.check_channels;
    % below will not flag error but remove weak channels automatically.
    % o = o.check_channels(true);
    %run code
    
%     o.EmptyTiles(:) = 1;
% o.EmptyTiles(3) = 0;

%     o.RegMethod='Fft';
%        o.RegMethod= 'Fft'; 
       o.DetectionThresh=700;
       o.IsolationThresh = 300;
    o = o.register2;
    save(fullfile(o.OutputDirectory, 'oRegister'), 'o', '-v7.3');
    
%     o.EmptyTiles(:) = 0;
%     o.EmptyTiles(4:end)  = 1;
%     
    %% find spots

    %  o.PcImageMatchesThresh = 100;
    %     o.MinPCMatchFract =  0.07;
    
    %Search paramaters
    o.FindSpotsMinScore = 'auto';
    o.FindSpotsStep = [5,5];
    %FindSpotsSearch can either be a 1x1 struct or a o.nRounds x 1 cell of
    %structs - have a different range for each round:
    o.FindSpotsSearch = struct();
    o.FindSpotsSearch.Y = -500:o.FindSpotsStep(1):500;
    o.FindSpotsSearch.X = -150:o.FindSpotsStep(2):150;
    %Make WidenSearch larger if you think you have a large shift between rounds
    o.FindSpotsWidenSearch = [50,50];
    
    o.PcDist = 3;
    o.PointCloudMethod = 1;     %1 or 2, set to 2 if no anchor round.
    %2 assumes same scaling to each color channel across all rounds.
    
    o.DetectionThresh='auto';
    o.IsolationThresh ='auto';
    %run code
    
%     nMatches<o.PcMinSpots | o.AllBaseSpotNo<o.PcMinSpots;
   % o.PcMinSpotsScaling
    o = o.find_spots2;
    save(fullfile(o.OutputDirectory, 'oFind_spots'), 'o', '-v7.3');
    
    %% call spots
    %run code
    o.CallSpotsCodeNorm = 'WholeCode';      %Other alternative is 'Round'
    o = o.call_spots;
    
    %OMP
    o.ompInitialNeighbThresh = 5;  % Increase to use less memory. Keep below 10.
    o = o.call_spots_omp;
    save(fullfile(o.OutputDirectory, 'oCall_spots_OMP'), 'o', '-v7.3');
end
%% plot results
% o = o.call_spots; % to plot bleed matrix
% iss_color_diagnostics(o);

OldOut = 'D:\bugeon\data\Thomas\Thomas-20260901\';
NewOut = 'Z:\shared\Stephane\Thomas\';

o.TileFiles = cellfun(@(x) strrep(x,OldOut,NewOut), o.TileFiles, 'UniformOutput', false);
o.OutputDirectory = strrep(o.OutputDirectory,OldOut,NewOut);

I = imadjust(imread(fullfile(o.OutputDirectory,'anchor_image.tif'))); % background image = Anchor

% Parameters to select spot to plot

%  QualOK = NeighbNonZeros>o.ompNeighbThresh | o.([pf,'SpotIntensity2'])>o.ompIntensityThresh |...
%          o.([pf,'SpotScore'])>o.ompScoreThresh;
%  QualOK = QualOK & o.([pf,'SpotIntensity2'])>o.ompIntensityThresh2 & ...
%       NeighbNonZeros>o.ompNeighbThresh2 & o.([pf,'SpotScore'])>o.ompScoreThresh2;

o.ompIntensityThresh = 0.5;
o.ompIntensityThresh2 = 0.005;
o.ompNeighbThresh = 18;
o.ompNeighbThresh2 = 10;
o.ompScoreThresh = 4.3;
o.ompScoreThresh2 = 1.1;

%Spots assigned to gene that is not largest coefficient for pixel
%have a stronger thresholding as given by these:
o.ompIntensityThresh3 = 0.01;
o.ompIntensityThresh3_CoefDiffFactor = 0.27;
o.ompNeighbThresh3 = 28;
o.ompScoreThresh3 = 6.9;

o.MarkerSize = 5;
o.PlotLineWidth = 1.2;

Roi = round([1, max(o.dpSpotGlobalYX(:,2)), ...
1, max(o.dpSpotGlobalYX(:,1))]);
o.plot(I,Roi,'OMP');
daspect([1 1 1])

o.iss_change_plot('OMP',[],o.GeneNames); % show all genes

% o.iss_change_plot('OMP',[],{'FGF13','FGF14','FGF14_isoform1','LRRC37B'}); % show some genes
% o.iss_change_plot('OMP',[],{'ALDOC','CALB1'}); 
%% diagnostics per spot

iss_view_spot_omp3(o,234321) % diagnostic showing spot image for each round and color

iss_view_omp(o,234321) % diagnostic showing color code and gene probabilities

%% 
% o.ompSpotGlobalYX:    2D coordinates of the spots
% o.ompSpotCodeNo:      Gene code number for each spots (index on o.GeneNames)
% o.GeneNames:          Name of the genes

% QualOK = quality_threshold(o,'omp'); % boolean selecting good spots using parameters above (o.ompIntensityThresh, o.ompNeighbThresh, o.ompScoreThresh)
% change_gene_symbols function: assign a color and marker for each gene for
% plotting3

%%
