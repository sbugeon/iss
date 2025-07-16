clear
SliceList = {'Round-'};
MainFolder = {'F:\446003\Slide5'};
OutputF = 'D:\ISS\446003\Slide5_test';
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
    o.FileBase{1} = strcat(SliceNb,'01'); %r0
    o.FileBase{2} = strcat(SliceNb,'03'); %r1
    o.FileBase{3} = strcat(SliceNb,'02'); %r2
    o.FileBase{4} = strcat(SliceNb,'04'); %r3
    o.FileBase{5} = strcat(SliceNb,'05'); %r4
    o.FileBase{6} = strcat(SliceNb,'06'); %r5
    o.FileBase{7} = strcat(SliceNb,'09'); %r6
    o.FileBase{8} = strcat(SliceNb,'08'); %anchor/dapi % 14;
    
%     o.FileBase{9} = strcat(SliceNb,'06'); %supp r3
%     o.FileBase{10} = strcat(SliceNb,'11'); %supp r6
%     o.FileBase{11} = strcat(SliceNb,'12'); %supp r6
%     RefRounds = [4 7 7];
    
    o.TileDirectory = fullfile(OutputF,SliceNb,'\tiles');
    mkdir(o.TileDirectory);
    o.OutputDirectory = fullfile(OutputF,SliceNb,'\output');
    mkdir(o.OutputDirectory);
    
    o.RawFileExtension = '.czi';
    
    %Codebook is a text file containing 2 columns - 1st is the gene name. 2nd is
    %the code, length o.nRounds and containing numbers in the range from 0 to o.nBP-1.
    o.CodeFile = 'G:\Documents\C_doc\codebook_73g_ctx.txt';
    %     o.CodeFile = 'G:\Documents\C_doc\codebook_7rounds.txt';%%%%%%%%%%%%
    
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
%     o.extract_and_filter_extraTiles(RefRounds);
    oOut = o.OutputDirectory;
    
    for j = 1:13%length(o.TileConnectedID)
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
        o.DetectionThresh = 'auto';
        o.DetectSpotsMaxSpots = 10000;
%         o.Graphics = 2;
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
        o.FindSpotsSearch.Y = -700:o.FindSpotsStep(1):700;
        o.FindSpotsSearch.X = -700:o.FindSpotsStep(2):700;
        %Make WidenSearch larger if you think you have a large shift between rounds
        o.FindSpotsWidenSearch = [250,250];
        
        o.PcDist = 3;
        o.PointCloudMethod = 1;     %1 or 2, set to 2 if no anchor round.
        %2 assumes same scaling to each color channel across all rounds.
        
        %run code
        o.DetectionThresh = 'auto';
        %         o.MinPCMatchFract = 0.05;
%         o.Graphics = 2;
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
                close all
        MakeDapi_images
    end
end

%% plot results
% o = o.call_spots; % to plot bleed matrix
% iss_color_diagnostics(o);
I = imadjust(imread(fullfile(o.OutputDirectory,'background_image_fixed.tif'))); % background image

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
camroll(180)
% hold on
% scatter(CellCalled.CellYX(:,2),CellCalled.CellYX(:,1),'ob','LineWidth',2)
%% adjust thresholds
o.ompScoreThresh = 5;
o.ompScoreThresh2 = 2;
o.ompIntensityThresh = 0.01;
o.ompIntensityThresh2 = 0.005;
o.ompNeighbThresh = 12;
o.ompNeighbThresh2 = 10;

o.iss_change_plot('OMP',[],o.GeneNames)
o.iss_change_plot('OMP',[], {'Prkca','Pthlh','Rasgfr2','Tac1'});
%% diagnostics per spot
iss_view_omp(o,234321)
iss_view_spot_omp3(o, 234321)
%% make nice gene legend
GeneLegendPlot(4,o.GeneNames)
set(gcf,'Position',[844,321,539,764.6666666666667])
print(gcf, '-dsvg', fullfile('C:\Users\bugeon\Documents\BugeonLabWebsite',['Gene Legends 124g.svg']));
print(gcf, '-dpng', '-r800', fullfile('C:\Users\bugeon\Documents\BugeonLabWebsite',['Gene Legends 124g.png']));
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
Slide = 'Slide23';
OutputF = fullfile('D:\ISS\',Ani,Slide,'Round-','output');
DapiF = 'D:\DapiSegmentationCellPose';

D = dir(OutputF);
Dn = {D.name};
Dn = Dn([D.isdir]);
Dn = Dn(3:end);

for i=1:length(Dn)
%     copyfile(fullfile(OutputF,Dn{i},'background_image_fixed.tif'),fullfile(DapiF,[Ani,'_',Slide,'_',Dn{i},'.tif']))
    
     I = imadjust(imread(fullfile(OutputF,Dn{i},'background_image_fixed.tif')),[0 0.003]);
     imwrite(I,fullfile(DapiF,[Ani,'_',Slide,'_',Dn{i},'.tif']))
end

fprintf('\n now run DAPI cellpose')
%%
ConvertCellPose2Bound_folder(DapiF)

for i=1:length(Dn)
    copyfile(fullfile(DapiF,['DAPI_Bound_',Ani,'_',Slide,'_',Dn{i},'.mat']),fullfile(OutputF,Dn{i},['DAPI_Bound_',Ani,'_',Slide,'_',Dn{i},'.mat']))
end

%% do cell calling
addpath('C:\Users\bugeon\Documents\GitHub\Transcriptomics-master')
addpath('C:\Users\bugeon\Documents\GitHub\iss')
addpath('C:\Users\bugeon\Documents\GitHub\RegSess')
load('C:\Users\bugeon\Documents\MATLAB\gsetYao.mat') % CA1 gene set
gSet.GeneName(strcmp(gSet.GeneName,'Ctgf')) = {'Ccn2'};
gSet.GeneName(strcmp(gSet.GeneName,'Fam19a1')) = {'Tafa1'};

for i=1:length(Dn)
    LocalFolder = fullfile(OutputF,Dn{i});
    CallCell_445951
end

%% get cell classes
clear
Ani = '445951';
Slide = 'Slide23';
OutputF = fullfile('D:\ISS\',Ani,Slide,'Round-','output');
DapiF = 'D:\DapiSegmentationCellPose';

D = dir(OutputF);
Dn = {D.name};
Dn = Dn([D.isdir]);
Dn = Dn(3:end);

Gn = readtable('G:\Documents\C_doc\codebook_73g_ctx.txt');
Gn = Gn.Var1;Gn(strcmp(Gn,'Rab3c'))=  {'Rab3b'};
AllSessionROI=struct();
AllSessionROI.pCellClass = [];
AllSessionROI.CellYX = [];
AllSessionROI.NSlice ={};
AllSessionROI.CellGeneCount=[];
AllSessionROI.tdTomato = [];
SliceCallInfo = struct();
AllSessionROI.tdTomatoIntensity = [];
AllSessionROI.tdTomato = [];
% AllSessionROI.CellGeneCount2={};
% recover cell calling results and tdTomato assignment
for i=1:length(Dn)
    LocalFolder = fullfile(OutputF,Dn{i});
    load(fullfile(LocalFolder,'Call_cells.mat'))
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
end

SliceCallInfo.GeneNames = CellCalled.GeneNames;
SliceCallInfo.ClassNames = CellCalled.ClassNames;
SliceCallInfo.ClassNames = cellfun(@(x) x(2:end),SliceCallInfo.ClassNames,'UniformOutput',false);
SliceCallInfo.ClassNames{end} = 'Zero';

% assign subtype and type

Subi = find( contains(CellCalled.ClassNames,'SUB') | contains(CellCalled.ClassNames,'POST') | contains(CellCalled.ClassNames,'PRE') | contains(CellCalled.ClassNames,'PAR') ...
    | contains(CellCalled.ClassNames,'APr'));
PC = 124:317;
PC_Ctx = find(contains(CellCalled.ClassNames,'CTX') | contains(CellCalled.ClassNames,'RSPv') | contains(CellCalled.ClassNames,'Car3'));
PC_hipp = setdiff(PC,PC_Ctx); PC_hipp = setdiff(PC_hipp,Subi);

ClassCollapse = {CellCalled.ClassNames([1:4,365:387]), 'Not_neuron', [138, 134, 134]./255,[1:4,365:387] ; ...% gray
    CellCalled.ClassNames(PC_Ctx), 'PC.Cortex', [233, 237, 26]./255,PC_Ctx ; ...% yellow
    CellCalled.ClassNames(PC_hipp), 'PC.Hippocampus', [237, 26, 160]./255,PC_hipp ; ... %pink
    CellCalled.ClassNames(5:123), 'Interneurons', [240, 2, 2]./255,5:123; ...% red
    CellCalled.ClassNames(Subi), 'Subiculum', [237, 26, 160]./255,Subi; ... %pink
    CellCalled.ClassNames(329:348), 'PC.CA1', [26, 234, 237]./255,329:348; ...% cyan
    CellCalled.ClassNames(359:360), 'PC.CA2', [26, 110, 237]./255,359:360; ... % blue
    CellCalled.ClassNames(351:358), 'PC.CA3', [146, 26, 237]./255,351:358; ...% purple
    CellCalled.ClassNames(349:350), 'MossyCells', [237, 170, 26]./255,349:350; ... % orange
    CellCalled.ClassNames(361:364), 'DG', [26, 237, 75]./255,361:364; ... % green
    CellCalled.ClassNames(388), 'Zero', [0 0 0],388};%black

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
   C = C(1:364,:);
AllSessionROI.CellGeneCount(isnan(AllSessionROI.CellGeneCount))= 0; %%%%%%%%%%%%%% remove NaNs
%%
% MeanGC = nansum(AllSessionROI.CellGeneCount,2);
D = fieldnames(AllSessionROI);
AllSessionROI2 = struct();
for i=1:length(D)
    AllSessionROI2.(D{i}) = AllSessionROI.(D{i})(~strcmp( AllSessionROI.Type,'Zero') & ~strcmp( AllSessionROI.Type,'Not_neuron'),:);%& MeanGC>-Inf
end

addpath('G:\Code\Matlab_OneDrive\MATLAB\umap')
id=1;

yy = makeUMAP(log(1+AllSessionROI2.CellGeneCount),0.2,10);
%%
close all

figure
XY = [];
for k=1:length(H)
    figure
    [ClassCollapse,ia] = unique(C.([H{k},'_label']));
%     ia = ia(randperm(length(ia),length(ia)));
    for j=1:size(ClassCollapse,1)
        CC_color = hex2rgb(C.([H{k},'_color']){ia(j)});
        hold on
        dd = find(strcmp(AllSessionROI2.(H{k}),ClassCollapse{j}));
        if length(dd)>4
            
            nn = nanmedian(yy(dd,:));
            % find cell the closest to mean of point cloud 
            D = squareform(pdist( [nn;yy(dd,:)]));D = D(1,2:end);
            [~,m] = min(D);
            XY(j,:) = yy(dd(m),:);
            
        else
            XY(j,:) = [NaN NaN];
        end
        hold on
        scatter(yy(dd,1),yy(dd,2),'.','MarkerEdgeColor',CC_color)
        
    end
    scatter(yy(logical(AllSessionROI2.tdTomato),1),yy(logical(AllSessionROI2.tdTomato),2),'sk','LineWidth',1.5)
    title(['Min dist ',num2str(0.2),' Neighb ',num2str(10)])
    text(XY(:,1),XY(:,2),strrep(ClassCollapse,'_','-'),'FontSize',10,'LineWidth',2,'Color','k')
    daspect([1 1 1])
    set(gcf,'Position',[118,61,1589,1295])
    xlabel('UMAP axis 1')
    ylabel('UMAP axis 1')
    saveas(gcf,fullfile('D:\ISS\445951',['UMAP-',H{k},'.png']))
      saveas(gcf,fullfile('D:\ISS\445951',['UMAP-',H{k},'.fig']))
    % histogram showing percentage of Tomato per subtype and type
    U = unique(AllSessionROI2.(H{k}));
    
    [~,NewO] = ismember(U,C.([H{k},'_label']));
    [~,S] =sort(NewO);
    U = U(S);
    dd =[];BarCol = [];
    GoodU = {};
    for i = 1:length(U)
        if sum(strcmp(U{i}, AllSessionROI2.(H{k})))>10
            dd = [dd;sum(strcmp(U{i},  AllSessionROI2.(H{k})) & logical(AllSessionROI2.tdTomato)) / sum(strcmp(U{i}, AllSessionROI2.(H{k})) )];
            GoodU = [GoodU;[U{i},', n = ',num2str(sum(strcmp(U{i}, AllSessionROI2.(H{k}))))]];
            col = hex2rgb(C.([H{k},'_color'])(strcmp(U{i},C.([H{k},'_label']))));
            BarCol = [BarCol;col(1,:)];
        end
    end
    GoodU = strrep(GoodU,'_','-');
    figure
    for i=1:length(dd)
        hold on
        bar(i,100*dd(i),'FaceColor',BarCol(i,:))
    end
    set(gca,'XTick',1:length(dd),'XTickLabel',GoodU)
%     xtickangle(90)
    ylabel('Pourcentage de cellules tdTomato positives')
    camroll(90)
%     set(gca,'XAxisLocation','top')
    set(gca,'XDir','reverse')
    if k==4
        set(gcf,'Position',[1000,42,560,1296])
    end
    saveas(gcf,fullfile('D:\ISS\445951',['Barplot-',H{k},'.png']))
end



