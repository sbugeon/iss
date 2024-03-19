clear
SliceList = {'Round_-'};
MainFolder = {'\\NETDATA\Apotome\Stephane\SBI002\2024-03-06'};
OutputF = 'D:\ISS';
for iSlice = 1:length(SliceList)
    SliceNb = SliceList{iSlice};
    
    %     %% Parameters that should be checked before each run
    o = iss_OMP;
    o.AnchorChannel = 6;
    o.AnchorRound = 8;            %Channel that has most spots in o.AnchorRound
    o.GadChannel = 7;
    o.GadRound = 8;
    o.GcampChannel = 4;
    o.GcampRound = 8;
    o.DapiChannel = 2;
    o.DapiRound = 8;             %Channel in o.AnchorRound that contains Dapi images
    o.nRegions = 2;
    
    o.InitialShiftChannel = 4;      %Channel to use to find initial shifts between rounds
    o.ReferenceRound = o.AnchorRound;           %Global coordinate system is built upon o.ReferenceRound and
    o.ReferenceChannel = o.AnchorChannel;         %o.ReferenceChannel. If RefRound = AnchorRound, this has to be AnchorChannel.
    o.RawFileExtension = '.nd2';    %Format of raw data
    o.LogToFile = 0;                %Set to 1 if you want to save command window to txt file, else set to 0.
    o.StripHack = true;
    %% File Names
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
    o.FileBase{8} = strcat(SliceNb,'Anchor');
    
    o.TileDirectory = fullfile(OutputF,SliceNb,'\tilesB');
    mkdir(o.TileDirectory);
    o.OutputDirectory = fullfile(OutputF,SliceNb,'\outputB');
    mkdir(o.OutputDirectory);
    
    o.RawFileExtension = '.czi';
    %Codebook is a text file containing 2 columns - 1st is the gene name. 2nd is
    %the code, length o.nRounds and containing numbers in the range from 0 to o.nBP-1.
    %     o.CodeFile = 'G:\Documents\C_doc\codebook_73g_ctx.txt';
    o.CodeFile = 'G:\Documents\C_doc\codebook_7rounds.txt';%%%%%%%%%%%%
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
    
    o.ExtractScale = 10;
    
    %Max time (seconds) to wait for raw .nd2 files to be obtained
    o.MaxWaitTime1 = 60;      %Less time for round 1 incase name is wrong
    o.MaxWaitTime = 21600;
    
    %     o.ExtractScale = 12.9230;
    %run code
    o = o.extract_and_filter;
    save(fullfile(o.OutputDirectory,'oExtract'), 'o', '-v7.3');
    
    oOut = o.OutputDirectory;
    
    for j = 2%1:length(o.TileConnectedID)
        load(fullfile(oOut,'oExtract'));
        o.EmptyTiles(:) = 1;
        for k = 1:length(o.TileConnectedID{j})
            Idx = o.TilePosYX(o.TileConnectedID{j}(k),:);
            o.EmptyTiles(Idx(1),Idx(2)) = 0;
        end
        o.OutputDirectory = fullfile(oOut,num2str(j));
        
        mkdir(fullfile(o.OutputDirectory))
        save(fullfile(o.OutputDirectory,'oExtract'), 'o', '-v7.3');
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
        o = o.register2;
        save(fullfile(o.OutputDirectory,'oRegister'), 'o', '-v7.3');
        %% find spots
        
        %  o.PcImageMatchesThresh = 100;
        %     o.MinPCMatchFract =  0.07;
        
        %Search paramaters
        o.FindSpotsMinScore = 'auto';
        o.FindSpotsStep = [5,5];
        %FindSpotsSearch can either be a 1x1 struct or a o.nRounds x 1 cell of
        %structs - have a different range for each round:
        %o.FindSpotsSearch = cell(o.nRounds,1);
        o.FindSpotsSearch.Y = -100:o.FindSpotsStep(1):100;
        o.FindSpotsSearch.X = -100:o.FindSpotsStep(2):100;
        %Make WidenSearch larger if you think you have a large shift between rounds
        o.FindSpotsWidenSearch = [50,50];
        
        o.PcDist = 3;
        o.PointCloudMethod = 1;     %1 or 2, set to 2 if no anchor round.
        %2 assumes same scaling to each color channel across all rounds.
        
        %run code
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
    end
end
%% plot results
% o = o.call_spots; % to plot bleed matrix
% iss_color_diagnostics(o);
% I = imread(fullfile(o.OutputDirectory,'background_image.tif'))*0; % background image
I=[];
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
%% adjust thresholds
o.ompScoreThresh = 5;
o.ompScoreThresh2 = 4;
o.ompIntensityThresh = 0.01;
o.ompIntensityThresh2 = 0.005;
o.ompNeighbThresh = 12;
o.ompNeighbThresh2 = 10;

o.ompScoreThresh = 0;
o.ompScoreThresh2 = 0;
o.ompIntensityThresh = 0;
o.ompIntensityThresh2 = 0;%\0.005;
o.ompNeighbThresh = 0;
o.ompNeighbThresh2 = 0;

o.iss_change_plot('OMP',[],o.GeneNames)
% o.iss_change_plot('OMP',[],{'Ddit4l'})
o.iss_change_plot('OMP',[], {'Rab3c'});

%% diagnostics per spot
iss_view_spot_omp3(o,234321)
iss_view_omp(o,234321)



%% Check gene efficiency
gg = o.GeneEfficiency<0.1;badGenes = o.GeneNames(sum(gg,2)>2);
o.iss_change_plot('OMP',[],badGenes)
%
SpotC = o.ompSpotColors;
SpotCode = o.GeneNames(o.ompSpotCodeNo); % dim2 = colors; dim3 = rounds
QualOK = o.quality_threshold('OMP');

for i=1:length(o.GeneNames)
    gg = QualOK & strcmp(SpotCode,o.GeneNames{i});
    This = SpotC(gg,:,:);
    
    M = squeeze(mean(This,1));
    CodeShape = size(This);
    SpotCoefs = full(o.ompCoefs(gg,:));
    PredCode = SpotCoefs*o.ompBledCodes(:,:);
    PredCode = squeeze(mean(reshape(PredCode,CodeShape),1));
    CodeShape = size(M);
    BledCode = o.BledCodes(i,:);
    BledCode = reshape(BledCode,CodeShape);
    
    figure(1)
    clf
    subplot(4,1,1)
    imagesc(BledCode);title('Bled Code')
    subplot(4,1,2)
    imagesc(PredCode);title('Normalized Predicted Code')
    subplot(4,1,3)
    imagesc(M);title('Mean Code')
    sgtitle(o.GeneNames{i})
    subplot(4,1,4)
    plot(o.GeneEfficiency(i,:))
    pause
end
load('G:\Data-Analysis\SB031\ISS\PixelCall\oCall_spots_pixel_B5S4_Slice001_OMP.mat')
dd = o.GeneEfficiency;
load('G:\Data-Analysis\SB031\ISS\PixelCall\oCall_spots_pixel_B4S5_Slice003_OMP.mat')
dd2= o.GeneEfficiency;
load('D:\ISS\Round_-\outputB\1\oCall_spots_OMP.mat')
dd0 = o.GeneEfficiency;

figure
scatter(dd(:),dd2(:))
ylim([0 3])
xlim([0 3])
hold on
plot([0 3], [0 3])

figure
scatter(dd0(:),dd2(:))
ylim([0 3])
xlim([0 3])
hold on
plot([0 3], [0 3])
figure
scatter(dd0(:),dd(:))
ylim([0 3])
xlim([0 3])
hold on
plot([0 3], [0 3])

edges = [0:0.05:3];

figure
subplot(3,1,1)
histogram(dd,edges)
ylim([0 50])
subplot(3,1,2)
histogram(dd2,edges)
ylim([0 50])
subplot(3,1,3)
histogram(dd0,edges)
ylim([0 50])

%% get oligo bridge concentration
CC = NaN(size(o.GeneEfficiency));
cc = readtable('C:\Users\bugeon\Desktop\bridge_concentration.txt');
N = table2cell(cc(:,1));
for i=1:size(cc,1)
    d = cc{i,1}{1};
    t = strfind(d,'_'); 
    geneN = find(strcmp(o.GeneNames, d(1:t-1)));
    round = str2num(d(t+2:end));
    CC(geneN,round+1) = cc{i,2};
end

figure
scatter(CC(:),o.GeneEfficiency(:))


