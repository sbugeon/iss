%% First, add the load the iss code to matlab path, then load the 'oCall_spots_OMP' file in matlab

%% Adjust parameters
NewOut = 'Z:\shared\Stephane\Thomas\'; % new folder where the processed data is stored
BackIm = 'DAPI'; % which image to use for background( 'anchor', 'DAPI' or '')

% plotting parameters
o.MarkerSize = 5; % size of markers for plot
o.PlotLineWidth = 1.2; % Line width on the markers for plot
o.MarkerType = 'GeneSymbols';   % which type of markers to use 'GeneSymbols','Letters','Dots'

% which genes to plot
Gene2P = o.GeneNames; % if all genes are to be plotted
% Gene2P = {'ALDOC','CALB1','FGF13'}; % only plots these genes

% Quality thresholds
% to be kept, gene spots must pass at least one of these thresholds
o.ompIntensityThresh = 0.5;
o.ompNeighbThresh = 18;
o.ompScoreThresh = 4.3;

% to be kept, gene spots must pass all of these thresholds
o.ompScoreThresh2 = 1.1;
o.ompIntensityThresh2 = 0.005;
o.ompNeighbThresh2 = 10;

%% Make the gene spot plot from scratch
% adjust input folder name
o.TileFiles = cellfun(@(x) strrep(x,o.InputDirectory,NewOut), o.TileFiles, 'UniformOutput', false);
o.OutputDirectory = strrep(o.OutputDirectory,o.InputDirectory,NewOut);

% load background image
I=[];
switch BackIm
    case 'DAPI'
        I = imadjust(imread(fullfile(o.OutputDirectory,'background_image.tif')));
    case 'anchor'
        I = imadjust(imread(fullfile(o.OutputDirectory,'anchor_image.tif')));
end

% do the plotting
Roi = round([1, max(o.dpSpotGlobalYX(:,2)), ...
1, max(o.dpSpotGlobalYX(:,1))]);
o.plot(I,Roi,'OMP');
daspect([1 1 1])

o.iss_change_plot('OMP',[],Gene2P); % show all genes

%% Only update which genes to plot
Gene2P = {'ALDOC','CALB1','FGF13'}; % only plots these genes
o.iss_change_plot('OMP',[],Gene2P); % show all genes

%% Diagnostics per spot (requires the gene spot plot)

% diagnostic showing color code and gene probabilities
iss_view_omp(o,234321) 

% diagnostic showing spot image for each round and color
iss_view_spot_omp3(o,234321) 

%% Outputs:
% o.ompSpotGlobalYX:    2D coordinates of the spots
% o.ompSpotCodeNo:      Gene code number for each spots (index on o.GeneNames)
% o.GeneNames:          Name of the genes
% QualOK = quality_threshold(o,'omp'); % boolean selecting good spots using parameters above (o.ompIntensityThresh, o.ompNeighbThresh, o.ompScoreThresh)

% Example usages:

% find the coordinates of all gene spots for ALDOC and CALB1:
SpotYX = {};
GeneN = {'ALDOC','CALB1'};
QualOK = quality_threshold(o,'OMP');
figure
hold on
for i=1:length(GeneN)
    ThisG = find(strcmp(GeneN{i},o.GeneNames));
    KeepSpots = ismember(o.ompSpotCodeNo,ThisG) & QualOK;
    SpotYX{i} = o.ompSpotGlobalYX(KeepSpots,:);
    
    scatter(SpotYX{i}(:,1),SpotYX{i}(:,2))
end
legend(GeneN)


%% To change the gene symbols:
% Modify change_gene_symbols.m function: assign a color and marker for each
% gene for plotting
