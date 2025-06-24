% requires the bfmatlab package !
%%
clear
SliceList = {'Round-'};
MainFolder = {'\\NETDATA\Apotome\Stephane\445951\Slide21'};
OutputF = 'D:\ISS\445951\Slide21';
for iSlice = 1:length(SliceList)
    SliceNb = SliceList{iSlice};
    
    %% Parameters that should be checked before each run
    o = iss_OMP;
    o.AnchorChannel = 3; % in which round 
    o.AnchorRound = 8;            %Channel that has most spots in o.AnchorRound
    o.GadChannel = 2;
    o.GadRound = 8;
    o.GcampChannel = 2;
    o.GcampRound = 8;
    o.DapiChannel = 1;
    o.DapiRound = 8;             %Channel in o.AnchorRound that contains Dapi images
    
    o.InitialShiftChannel = 5;      %Channel to use to find initial shifts between rounds
    o.ReferenceRound = o.AnchorRound;           %Global coordinate system is built upon o.ReferenceRound and
    o.ReferenceChannel = o.AnchorChannel;         %o.ReferenceChannel. If RefRound = AnchorRound, this has to be AnchorChannel.
    o.RawFileExtension = '.czi';   %Format of raw data 
    o.LogToFile = 0;                %Set to 1 if you want to save command window to txt file, else set to 0.
    o.StripHack = true;
    
    % File Names
    %CHECK BEFORE EACH RUN
    o.InputDirectory = MainFolder{iSlice};     %Folder path of raw data
    
    o.TileSz = 2048;        %Dimension of tile in pixels
    o.nBP = 7;              %Number of Channels
    o.nRounds = 7;          %Number of Imaging Rounds
    o.nExtraRounds = 1;     %Treat Anchor channel as extra round
    
    %FileBase{r} is the file name of the raw data of round r in o.InputDirectory
    o.FileBase = cell(1,1);
    o.FileBase{1} = strcat(SliceNb,'01');
    o.FileBase{2} = strcat(SliceNb,'02');
    o.FileBase{3} = strcat(SliceNb,'03');
    o.FileBase{4} = strcat(SliceNb,'04');
    o.FileBase{5} = strcat(SliceNb,'05');
    o.FileBase{6} = strcat(SliceNb,'06');
    o.FileBase{7} = strcat(SliceNb,'07');
    o.FileBase{8} = strcat(SliceNb,'08');
    
    % indicate supplementary rounds done in case some z focus was wrong for
    % some rounds, allows to have less tiles than the main file
    o.FileBase{9} = strcat(SliceNb,'06'); %supp r3
    o.FileBase{10} = strcat(SliceNb,'11'); %supp r6
    o.FileBase{11} = strcat(SliceNb,'12'); %supp r6
    % give the round that was redone for each of the supplementary rounds
    RefRounds = [4 7 7];
    
    o.TileDirectory = fullfile(OutputF,SliceNb,'\tiles');
    mkdir(o.TileDirectory);
    o.OutputDirectory = fullfile(OutputF,SliceNb,'\output');
    mkdir(o.OutputDirectory);
      
    %Codebook is a text file containing 2 columns - 1st is the gene name. 2nd is
    %the code, length o.nRounds and containing numbers in the range from 0 to o.nBP-1.
    o.CodeFile = 'G:\Documents\C_doc\codebook_7rounds.txt'; 
    
    o.bpLabels = {'0', '1', '2','3', '4','5','6'}; %order of bases 
    %% Logging
    if o.LogToFile
        if isempty(o.LogFile)
            o.LogFile = fullfile(o.OutputDirectory,'Log.txt');
        end
    end
    %% extract and filter
    
    %parameters
    o.FirstBaseChannel = 1;

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
    
    %run code
    o = o.extract_and_filter;
    save(fullfile(o.OutputDirectory,'oExtract'), 'o', '-v7.3');
    
    % if there were some supplementary rounds, extract them and replace the
    % former tiles with them
    o.extract_and_filter_extraTiles(RefRounds);
    
     oOut = o.OutputDirectory;
    
    for j = 1:length(o.TileConnectedID)
        load(fullfile(oOut,'oExtract'));
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
        %% register
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
        o.InitalShiftAutoMinScoreParam = 2;   %a lower value will make it quicker but more likely to fail
        
        %paramaters to find shifts between overlapping tiles
        o.RegMinScore = 'auto';
        o.RegStep = [5,5];
        o.RegSearch.South.Y = -1900:o.RegStep(1):-1700;
        o.RegSearch.South.X = -150:o.RegStep(2):150;
        o.RegSearch.East.Y = -150:o.RegStep(1):150;
        o.RegSearch.East.X = -1900:o.RegStep(2):-1700;
        o.RegWidenSearch = [50,50];
        
        %If a channel or round is faulty, you can ignore it by selecting only the
        %good ones in o.UseChannels and o.UseRounds.
        o.UseChannels = 1:o.nBP;
        o.UseRounds = 1:o.nRounds;
        
        % check channels
        % below will flag error if some channels are weak
        %     o = o.check_channels;
        % below will not flag error but remove weak channels automatically.
        % o = o.check_channels(true);
        %run code
        o.DetectSpotsMaxSpots = 10000;
        o = o.register2;
        save(fullfile(o.OutputDirectory,'oRegister'), 'o', '-v7.3');
        %% find spots
       
        %Search paramaters
        o.FindSpotsMinScore = 'auto';
        o.FindSpotsStep = [5,5];
        %FindSpotsSearch can either be a 1x1 struct or a o.nRounds x 1 cell of
        %structs - have a different range for each round:
        %o.FindSpotsSearch = cell(o.nRounds,1);
        o.FindSpotsSearch.Y = -700:o.FindSpotsStep(1):700;
        o.FindSpotsSearch.X = -700:o.FindSpotsStep(2):700;
        %Make WidenSearch larger if you think you have a large shift between rounds
        o.FindSpotsWidenSearch = [250,250];
        
        o.PcDist = 3;
        o.PointCloudMethod = 1;     %1 or 2, set to 2 if no anchor round.
        %2 assumes same scaling to each color channel across all rounds.
        
        %run code
        o.DetectionThresh = 'auto';
        o = o.find_spots2;
        save(fullfile(o.OutputDirectory,'oFind_spots'), 'o', '-v7.3');
        
        %% call spots
        %run code
        o.CallSpotsCodeNorm = 'WholeCode';      %Other alternative is 'Round'
        o = o.call_spots;
        
        %OMP
        o.ompInitialNeighbThresh = 5;  % Increase to use less memory. Keep below 10.
        o = o.call_spots_omp;
        save(fullfile(o.OutputDirectory,'oCall_spots_OMP'), 'o', '-v7.3');
        close all
    end
end

%% plot results
% o = o.call_spots; % to plot bleed matrix
% iss_color_diagnostics(o);
I = []; % background image

o.ompScoreThresh = 5; % more stringent threshold, need to be true for only one of the three
o.ompScoreThresh2 = 2; % less stringent threshold, need to be true for all three
o.ompIntensityThresh = 0.01; % more stringent threshold, need to be true for only one of the three
o.ompIntensityThresh2 = 0.005; % less stringent threshold, need to be true for all three
o.ompNeighbThresh = 12; % more stringent threshold, need to be true for only one of the three
o.ompNeighbThresh2 = 10; % less stringent threshold, need to be true for all three

o.MarkerSize = 5;
o.PlotLineWidth = 1.2;

Roi = round([1, max(o.dpSpotGlobalYX(:,2)), ...
    1, max(o.dpSpotGlobalYX(:,1))]);
o.plot(I,Roi,'OMP');
daspect([1 1 1])
camroll(90)


%% adjust thresholds
o.ompScoreThresh = 5;
o.ompScoreThresh2 = 2;
o.ompIntensityThresh = 0.01;
o.ompIntensityThresh2 = 0.005;
o.ompNeighbThresh = 12;
o.ompNeighbThresh2 = 10;

o.iss_change_plot('OMP',[],o.GeneNames)
% show only some genes
o.iss_change_plot('OMP',[], {'Npy','Pvalb','Penk','Pcp4'});
%% diagnostics per spot
iss_view_omp(o,234321)
iss_view_spot_omp3(o, 234321)


show7by7(o,150,OutputF )