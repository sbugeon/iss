clear
SliceList = {'Round-'};
MainFolder = {'\\NETDATA\Apotome\Stephane\445951\Slide22'};
OutputF = 'D:\ISS\445951\Slide22';
for iSlice = 1:length(SliceList)
    SliceNb = SliceList{iSlice};
    
    %% Parameters that should be checked before each run
    o = iss_OMP;
    o.AnchorChannel = 6; % in which round 
    o.AnchorRound = 8;            %Channel that has most spots in o.AnchorRound
    o.GadChannel = 2;
    o.GadRound = 8;
    o.GcampChannel = 2;
    o.GcampRound = 8;
    o.DapiChannel = 1;
    o.DapiRound = 8;             %Channel in o.AnchorRound that contains Dapi images
    o.nRegions = 7;
    
    o.InitialShiftChannel = 5;      %Channel to use to find initial shifts between rounds
    o.ReferenceRound = o.AnchorRound;           %Global coordinate system is built upon o.ReferenceRound and
    o.ReferenceChannel = o.AnchorChannel;         %o.ReferenceChannel. If RefRound = AnchorRound, this has to be AnchorChannel.
    o.RawFileExtension = '.nd2';    %Format of raw data
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
    o.FileBase{1} = strcat(SliceNb,'02');
    o.FileBase{2} = strcat(SliceNb,'03');
    o.FileBase{3} = strcat(SliceNb,'04');
    o.FileBase{4} = strcat(SliceNb,'05');
    o.FileBase{5} = strcat(SliceNb,'06');
    o.FileBase{6} = strcat(SliceNb,'08');
    o.FileBase{7} = strcat(SliceNb,'09');
    o.FileBase{8} = strcat(SliceNb,'10');
    
    o.TileDirectory = fullfile(OutputF,SliceNb,'\tiles');
    mkdir(o.TileDirectory);
    o.OutputDirectory = fullfile(OutputF,SliceNb,'\output');
    mkdir(o.OutputDirectory);
    
    o.RawFileExtension = '.czi';
    
    %Codebook is a text file containing 2 columns - 1st is the gene name. 2nd is
    %the code, length o.nRounds and containing numbers in the range from 0 to o.nBP-1.
%         o.CodeFile = 'G:\Documents\C_doc\codebook_73g_ctx.txt';
    o.CodeFile = 'G:\Documents\C_doc\codebook_7rounds.txt';%%%%%%%%%%%%   
    
    %     o.bpLabels = {'0', '2', '1','3', '4','5','6'}; %order of bases %%%%%%%%%%%%%%%
    o.bpLabels = {'0', '1', '2','3', '4','5','6'}; %order of bases %%%%%%%%%%%%%%%
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
%     
%     o.ExtractScale = 10;
%     fprintf('\n ExtractScale set to 10 for now, problem with saturated pixels otherwise')
    
    %Max time (seconds) to wait for raw .nd2 files to be obtained
    o.MaxWaitTime1 = 60;      %Less time for round 1 incase name is wrong
    o.MaxWaitTime = 21600;
    
    %     o.ExtractScale = 12.9230;
    %run code
    o = o.extract_and_filter;
    save(fullfile(o.OutputDirectory,'oExtract'), 'o', '-v7.3');
    
    oOut = o.OutputDirectory;
    
    for j = 4:length(o.TileConnectedID)
        load(fullfile(oOut,'oExtract'));
        o.Graphics=1;
        o.EmptyTiles(:) = 1;
        gg = o.TilePosYX;
        for k = 1:length(o.TileConnectedID{j})
            Idx = gg(o.TileConnectedID{j}(k),:);
            o.EmptyTiles(Idx(1),Idx(2)) = 0;
        end
%         o.EmptyTiles = flipud(fliplr(o.EmptyTiles));
%         o.EmptyTiles
%     end
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
        o.DetectionThresh = 700;
        o.DetectSpotsMaxSpots = 10000;
        o.Graphics = 2;
        o = o.register2;
%         o = o.register_FFt;
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
        o.FindSpotsSearch = struct();
        o.FindSpotsSearch.Y = -100:o.FindSpotsStep(1):100;
        o.FindSpotsSearch.X = -100:o.FindSpotsStep(2):100;
        %Make WidenSearch larger if you think you have a large shift between rounds
        o.FindSpotsWidenSearch = [200,200]; 
        o.InitalShiftAutoMinScoreParam = 5; % increased threshold for initial alignment score
        o.PcDist = 3;
        o.PointCloudMethod = 1;     %1 or 2, set to 2 if no anchor round.
        %2 assumes same scaling to each color channel across all rounds.
        
        %run code
        o.DetectionThresh = 'auto';
%         o.MinPCMatchFract = 0.05;
        o.Graphics = 1; 
%          o.ToPlot = [1 1:7 1:7];
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
%         close all
    MakeDapi_images
    end
end
%% plot results
% o = o.call_spots; % to plot bleed matrix
% iss_color_diagnostics(o);
I = imadjust(imread(fullfile(o.OutputDirectory,'background_image.tif'))*0,[0 0.1]); % background image

o.ompScoreThresh = 5; % more stringent threshold, need to be true for only one of the three
o.ompScoreThresh2 = 2; % less stringent threshold, need to be true for all three
o.ompIntensityThresh = 0.01; % more stringent threshold, need to be true for only one of the three
o.ompIntensityThresh2 = 0.005; % less stringent threshold, need to be true for all three
o.ompNeighbThresh = 12; % more stringent threshold, need to be true for only one of the three
o.ompNeighbThresh2 = 10; % less stringent threshold, need to be true for all three

% o.ompScoreThresh = 0; % more stringent threshold, need to be true for only one of the three
% o.ompScoreThresh2 = 0; % less stringent threshold, need to be true for all three
% o.ompIntensityThresh = 0; % more stringent threshold, need to be true for only one of the three
% o.ompIntensityThresh2 = 0; % less stringent threshold, need to be true for all three
% o.ompNeighbThresh = 0; % more stringent threshold, need to be true for only one of the three
% o.ompNeighbThresh2 = 0; % less stringent threshold, need to be true for all three
% o.ompScoreThresh = 7;
% o.ompScoreThresh2 = 7;
% o.ompIntensityThresh = 0.1;
% o.ompIntensityThresh2 = 0.1;
% o.ompNeighbThresh = 15;
% o.ompNeighbThresh2 = 15;
% o.ompScoreThresh = 0; % more stringent threshold, need to be true for only one of the three
% o.ompScoreThresh2 = 0; % less stringent threshold, need to be true for all three
% o.ompIntensityThresh =0; % more stringent threshold, need to be true for only one of the three
% o.ompIntensityThresh2 = 0; % less stringent threshold, need to be true for all three
% o.ompNeighbThresh = 0; % more stringent threshold, need to be true for only one of the three
% o.ompNeighbThresh2 = 0; % less stringent threshold, need to be true for all three
o.MarkerSize = 5;
o.PlotLineWidth = 1.2;
o.CombiQualThresh = 0.5;
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
o.iss_change_plot('OMP',[], {'Calb1'});
%% diagnostics per spot
iss_view_omp(o,234321)
iss_view_spot_omp3(o, 234321)

%%
FindColorCode

load('D:\ISS\445951\Slide22\Round-\output\oExtract.mat')
oOut = o.OutputDirectory;
for k = 2:length(o.TileConnectedID)
        load(fullfile(oOut,num2str(k),'oCall_spots_OMP'));
figure
for i=1:7
    subplot(1,7,i)
    title(Colors{i})
    for j = 1:7
        gg = AllCodes(:,j) == i-1;
        hold on
        Violin(o.GeneEfficiency(gg,j),j);
    end
    ylim([0 3])
    set(gca,'XTick',1:7,'XTickLabel',{'r0','r1','r2','r3','r4','r5','r6'})
    if i==1
        ylabel('Gene Efficiency')
    end
end
set(gcf,'Position',[244,926,2305,336])
sgtitle(['Gene efficiency per round and per channel for region ', num2str(k)])
end
dd = o.TileConnectedID;
%% copy dapi image to Cellpose folder
clear
Ani = '445951';
Slide = 'Slide22';
OutputF = fullfile('D:\ISS\',Ani,Slide,'Round-','output');
DapiF = 'D:\DapiSegmentationCellPose';

D = dir(OutputF);
Dn = {D.name};
Dn = Dn([D.isdir]);
Dn = Dn(3:end);
% 
% for i=1:length(Dn)
%     copyfile(fullfile(OutputF,Dn{i},'background_image_fixed.tif'),fullfile(DapiF,[Ani,'_',Slide,'_',Dn{i},'.tif']))
% end
% 
% ConvertCellPose2Bound_folder(DapiF)
% 
% for i=1:length(Dn)
%     copyfile(fullfile(DapiF,['DAPI_Bound_',Ani,'_',Slide,'_',Dn{i},'.mat']),fullfile(OutputF,Dn{i},['DAPI_Bound_',Ani,'_',Slide,'_',Dn{i},'.mat']))
% end

addpath('C:\Users\bugeon\Documents\GitHub\Transcriptomics-master')
addpath('C:\Users\bugeon\Documents\GitHub\iss')
addpath('C:\Users\bugeon\Documents\GitHub\RegSess')
load('C:\Users\bugeon\Documents\MATLAB\gsetYao.mat') % CA1 gene set
gSet.GeneName(strcmp(gSet.GeneName,'Ctgf')) = {'Ccn2'};
gSet.GeneName(strcmp(gSet.GeneName,'Fam19a1')) = {'Tafa1'};


% cell calling
for i=1:length(Dn)
    LocalFolder = fullfile(OutputF,Dn{i});
    CallCell_445951
end