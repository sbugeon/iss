clear
SliceList = {'Round-','Round-','Round-'};
MainFolder = {'\\NETDATA\Apotome\Stephane\445888\Slide14','\\NETDATA\Apotome\Stephane\445884\Slide 17','\\NETDATA\Apotome\Stephane\445879\Slide10'};
OutputF = {'D:\ISS\445888\Slide14bis','D:\ISS\445884\Slide17bis','D:\ISS\445879\Slide10bis'};
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
    o.nRegions = 7;
    
    o.InitialShiftChannel = 5;      %Channel to use to find initial shifts between rounds
    o.ReferenceRound = o.AnchorRound;           %Global coordinate system is built upon o.ReferenceRound and
    o.ReferenceChannel = o.AnchorChannel;         %o.ReferenceChannel. If RefRound = AnchorRound, this has to be AnchorChannel.
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
    
    o.TileDirectory = fullfile(OutputF{iSlice},SliceNb,'\tiles');
    mkdir(o.TileDirectory);
    o.OutputDirectory = fullfile(OutputF{iSlice},SliceNb,'\output');
    mkdir(o.OutputDirectory);
    
    o.RawFileExtension = '.czi'; %Format of raw data
    
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
    
    for j = 1:length(o.TileConnectedID)
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
I = imadjust(imread(fullfile(o.OutputDirectory,'background_image_fixed.tif')),[0 0.1]); % background image
I2 = imadjust(imread(fullfile(o.OutputDirectory,'tdTomato_image_fixed.tif')),[0 0.005]); % background image

I = I2;
o.ompScoreThresh = 5; % more stringent threshold, need to be true for only one of the three
o.ompScoreThresh2 = 2; % less stringent threshold, need to be true for all three
o.ompIntensityThresh = 0.01; % more stringent threshold, need to be true for only one of the three
o.ompIntensityThresh2 = 0.005; % less stringent threshold, need to be true for all three
o.ompNeighbThresh = 12; % more stringent threshold, need to be true for only one of the three
o.ompNeighbThresh2 = 10; % less stringent threshold, need to be true for all three


o.MarkerSize = 5;
o.PlotLineWidth = 1.2;
o.CombiQualThresh = 0.5;
Roi = round([1, max(o.dpSpotGlobalYX(:,2)), ...
    1, max(o.dpSpotGlobalYX(:,1))]);
o.plot(I,Roi,'OMP');
daspect([1 1 1])
camroll(90)
%% adjust thresholds
% o.ompScoreThresh = 5;
% o.ompScoreThresh2 = 2;
% o.ompIntensityThresh = 0.01;
% o.ompIntensityThresh2 = 0.005;
% o.ompNeighbThresh = 12;
% o.ompNeighbThresh2 = 10;
%
% o.iss_change_plot('OMP',[],o.GeneNames)
% o.iss_change_plot('OMP',[], {'Calb1'});
% %% diagnostics per spot
% iss_view_omp(o,234321)
% iss_view_spot_omp3(o, 234321)
%
% %%
% FindColorCode
%
% load('D:\ISS\445951\Slide22\Round-\output\oExtract.mat')
% oOut = o.OutputDirectory;
% for k = 2:length(o.TileConnectedID)
%         load(fullfile(oOut,num2str(k),'oCall_spots_OMP'));
% figure
% for i=1:7
%     subplot(1,7,i)
%     title(Colors{i})
%     for j = 1:7
%         gg = AllCodes(:,j) == i-1;
%         hold on
%         Violin(o.GeneEfficiency(gg,j),j);
%     end
%     ylim([0 3])
%     set(gca,'XTick',1:7,'XTickLabel',{'r0','r1','r2','r3','r4','r5','r6'})
%     if i==1
%         ylabel('Gene Efficiency')
%     end
% end
% set(gcf,'Position',[244,926,2305,336])
% sgtitle(['Gene efficiency per round and per channel for region ', num2str(k)])
% end
% dd = o.TileConnectedID;

%% run cellpose in matlab! 
model_path = 'G:\Code\training_cellposeDAPI\models\CP_myDAPI_SB';
img = imread("AT3_1m4_01.tif");
cp = cellpose(Model="cyto2");
labels = segmentCells2D(cp,img,ImageCellDiameter=56);
B = labeloverlay(img,labels);
imshow(B)
OutputModelFile = "diamondDetectionModel";
trainCellpose(trainingFolderName,OutputModelFile,...
    PreTrainedModel="cyto2", ...
    MaxEpoch=2, ... 
    ImageSuffix="_im", ...
    LabelSuffix="_mask");
cpt = cellpose(Model=OutputModelFile);
cpt.TrainingCellDiameter

%% copy dapi image to Cellpose folder
clear
MainFolder = {'\\NETDATA\Apotome\Stephane\445888\Slide14','\\NETDATA\Apotome\Stephane\445884\Slide 17','\\NETDATA\Apotome\Stephane\445879\Slide10'};

AniL = {'445888','445884','445879'};
SlideL = {'Slide14bis','Slide17bis','Slide10bis'};
for j=1:length(AniL)
    OutputF = fullfile('D:\ISS\',AniL{j},SlideL{j},'Round-','output');
    DapiF = 'D:\DapiSegmentationCellPose';
    
    D = dir(OutputF);
    Dn = {D.name};
    Dn = Dn([D.isdir]);
    Dn = Dn(3:end);
    
    for i=1:length(Dn)
        copyfile(fullfile(OutputF,Dn{i},'background_image_fixed.tif'),fullfile(DapiF,[Ani{j},'_',SlideL{j},'_',Dn{i},'.tif']))
    end
end


ConvertCellPose2Bound_folder(DapiF)


for j=1:length(AniL)
    OutputF = fullfile('D:\ISS\',AniL{j},SlideL{j},'Round-','output');
    DapiF = 'D:\DapiSegmentationCellPose';
    
    D = dir(OutputF);
    Dn = {D.name};
    Dn = Dn([D.isdir]);
    Dn = Dn(3:end);
    for i=1:length(Dn)
        copyfile(fullfile(DapiF,['DAPI_Bound_',AniL{j},'_',SlideL{j},'_',Dn{i},'.mat']),fullfile(OutputF,Dn{i},['DAPI_Bound_',AniL{j},'_',SlideL{j},'_',Dn{i},'.mat']))
    end
end

addpath('C:\Users\bugeon\Documents\GitHub\RegSessPipeline\GUI_Slice_Curation')
GUI_findtdTomato('D:\ISS\445879\Slide10bis\Round-\output')

%%
clear

AniL = {'445888','445884','445879'};
SlideL = {'Slide14bis','Slide17bis','Slide10bis'};
addpath('C:\Users\bugeon\Documents\GitHub\Transcriptomics-master')
addpath('C:\Users\bugeon\Documents\GitHub\iss')
addpath('C:\Users\bugeon\Documents\GitHub\RegSess')
% load('C:\Users\bugeon\Documents\MATLAB\gsetYao.mat') % CA1 gene set
load('C:\Users\bugeon\Documents\MATLAB\gsetTasic_123g.mat') % CA1 gene set
gSet.GeneName(strcmp(gSet.GeneName,'Ctgf')) = {'Ccn2'};
gSet.GeneName(strcmp(gSet.GeneName,'Fam19a1')) = {'Tafa1'};

for j=1:length(AniL)
    OutputF = fullfile('D:\ISS\',AniL{j},SlideL{j},'Round-','output');
    DapiF = 'D:\DapiSegmentationCellPose';
    D = dir(OutputF);
    Dn = {D.name};
    Dn = Dn([D.isdir]);
    Dn = Dn(3:end);
    Ani = AniL{j};
    Slide = SlideL{j};
    % cell calling
    for i=1:length(Dn)
        LocalFolder = fullfile(OutputF,Dn{i});
        CallCell_dlxGFP
    end
end
%%
clear

AniL = {'445888','445884','445879'};
SlideL = {'Slide14bis','Slide17bis','Slide10bis'};
AllSessionROI=struct();
AllSessionROI.pCellClass = [];
AllSessionROI.CellYX = [];
AllSessionROI.NSlice ={};
AllSessionROI.CellGeneCount=[];
AllSessionROI.tdTomato = [];
SliceCallInfo = struct();
AllSessionROI.tdTomatoIntensity = [];
AllSessionROI.tdTomato = [];
id=1;
for j=1:length(AniL)
    
    Ani = AniL{j};
    Slide = SlideL{j};
    
    OutputF = fullfile('D:\ISS\',Ani,Slide,'Round-','output');
    DapiF = 'D:\DapiSegmentationCellPose';
    
    D = dir(OutputF);
    Dn = {D.name};
    Dn = Dn([D.isdir]);
    Dn = Dn(3:end);
    
    Gn = readtable('G:\Documents\C_doc\codebook_73g_ctx.txt');
    Gn = Gn.Var1;Gn(strcmp(Gn,'Rab3c'))=  {'Rab3b'};
    
    % AllSessionROI.CellGeneCount2={};
    % recover cell calling results and tdTomato assignment
    for i=1:length(Dn)
        LocalFolder = fullfile(OutputF,Dn{i});
        load(fullfile(LocalFolder,'Call_cells.mat'))
        load(fullfile(LocalFolder,'oCall_spots_OMP.mat'))
        DapiImg = imread(fullfile(LocalFolder,'background_image_fixed.tif'));
        load(fullfile(LocalFolder,['DAPI_Bound_',Ani,'_',Slide,'_',Dn{i},'.mat']))
        
        AllSessionROI.pCellClass = [AllSessionROI.pCellClass; CellCalled.pCellClass(1:end-1,:)];
        AllSessionROI.CellYX = [AllSessionROI.CellYX ; CellCalled.CellYX];
        
        
        GeneCount = NaN(size(CellCalled.CellGeneCount,1),length(Gn));
        for k=1:length(Gn)
            DD = find(strcmp(Gn{k},CellCalled.GeneNames));
            if ~isempty(DD)
                GeneCount(:,k) = CellCalled.CellGeneCount(:,DD);
            else
                fprintf(['\n ' ,Gn{k}])
            end
        end
        
        AllSessionROI.CellGeneCount = [AllSessionROI.CellGeneCount; GeneCount(1:end-1,:)];
        for iSlice = 1:size(CellCalled.CellYX,1)
            AllSessionROI.NSlice = [AllSessionROI.NSlice; [Ani,'_',Slide,'_',Dn{i}]];
        end
        TomatoDetected = zeros(size(CellCalled.RabiesI,1),1);
        TomatoDetected(Tomato) = 1;
        
        AllSessionROI.tdTomatoIntensity = [AllSessionROI.tdTomatoIntensity;CellCalled.RabiesI];
        AllSessionROI.tdTomato = [AllSessionROI.tdTomato;TomatoDetected];
        %     AllSessionROI.CellGeneCount2 = [AllSessionROI.CellGeneCount2; CellCalled.CellGeneCount];
        if max(Tomato)>length(CellCalled.RabiesI)
            dd=0;
        end
        for k=1:length(Tomato)
            plotISScells(DapiImg,o,CellCalled,Tomato(k),'D:\ISS\plots_Bertrand',id)
            id=id+1;
        end
    end
    
    
end
SliceCallInfo.GeneNames = CellCalled.GeneNames;
SliceCallInfo.ClassNames = CellCalled.ClassNames;
SliceCallInfo.ClassNames = cellfun(@(x) x(2:end),SliceCallInfo.ClassNames,'UniformOutput',false);
SliceCallInfo.ClassNames{end} = 'Zero';

C = readtable('C:\Users\bugeon\Documents\Yao_classcolors.xlsx');
C.cluster_label = cellfun(@(x) strrep(x,'-','_'),C.cluster_label,'UniformOutput',false);
C.cluster_label = cellfun(@(x) strrep(x,' ','_'),C.cluster_label,'UniformOutput',false);
C.cluster_label = cellfun(@(x) strrep(x,'/','_'),C.cluster_label,'UniformOutput',false);
C.subclass_label(contains(C.cluster_label,'Mossy')) = {'Mossy'};

H = {'class','neighborhood','subclass','cluster'};
AllSessionROI.Subtype = {};
% assign subtype
for i=1:size(AllSessionROI.CellYX,1)
    p = AllSessionROI.pCellClass(i,:);
    if p(end)<0.1
        [~,m] = max(p);
        AllSessionROI.Subtype{i,1} = SliceCallInfo.ClassNames{m};
    else
        AllSessionROI.Subtype{i,1} = 'Zero';
    end
end

for k=1:length(H)
    AllSessionROI.(H{k}) = {};
    for i=1:size(AllSessionROI.CellYX,1)
        dd = find(strcmp(AllSessionROI.Subtype{i,1},C.cluster_label));
        if ~strcmp(AllSessionROI.Subtype{i,1},'Zero')
            AllSessionROI.(H{k}){i,1} = C.([H{k},'_label']){dd};
        else
            AllSessionROI.(H{k}){i,1} = 'Zero';
        end
    end
end
% C = C(1:364,:);
AllSessionROI.CellGeneCount(isnan(AllSessionROI.CellGeneCount))= 0; %%%%%%%%%%%%%% remove NaNs

dd = AllSessionROI.Subtype(contains(AllSessionROI.Subtype,'Sst'));
% dd = AllSessionROI.Subtype;
U = unique(dd);
gg=[];
for i=1:length(U)
   gg(i) = sum(strcmp(dd,U{i}));
    
end

close all
plot(gg,'-o')
set(gca,'XTick',1:length(U),'XTickLabel',U)
camroll(-90)

GE = AllSessionROI.CellGeneCount(logical(AllSessionROI.tdTomato),:);

figure
imagesc(log(1+GE(:,1:73)))
set(gca,'XTick',1:72,'XTickLabel',CellCalled.GeneNames)
set(gca,'YTick',1:10,'YTickLabel',strrep(AllSessionROI.Subtype(logical(AllSessionROI.tdTomato)),'_','-'))
camroll(-90)
ytickangle(45)