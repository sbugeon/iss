% load('D:\ISS\445888\Slide14_test\Round-\output\1\oCall_spots_OMP.mat')
%% to plot bleed matrix
o = o.call_spots; 

%% color diagnostics
iss_color_diagnostics(o);

%% plot gene calling
I = []; % background image

o.ompScoreThresh = 5; % more stringent threshold, need to be true for only one of the three
o.ompScoreThresh2 = 2; % less stringent threshold, need to be true for all three
o.ompIntensityThresh = 0.01; % more stringent threshold, need to be true for only one of the three
o.ompIntensityThresh2 = 0.005; % less stringent threshold, need to be true for all three
o.ompNeighbThresh = 12; % more stringent threshold, need to be true for only one of the three
o.ompNeighbThresh2 = 10; % less stringent threshold, need to be true for all three

o.MarkerSize = 5; % size of markers in the plot
o.PlotLineWidth = 1.2; % line width of marker edges
o.MarkerType = 'GeneSymbols';   % which type of markers to use 'GeneSymbols','Dots'

Roi = round([1, max(o.dpSpotGlobalYX(:,2)), ...
    1, max(o.dpSpotGlobalYX(:,1))]);
o.plot(I,Roi,'OMP');
daspect([1 1 1])
camroll(90)

% show all genes
% o.iss_change_plot('OMP',[],o.GeneNames)

% show only some genes
% o.iss_change_plot('OMP',[], {'Npy','Pvalb','Penk','Pcp4'});

%% diagnostics per spot 
% light
iss_view_omp(o,234321)
% heavy
iss_view_spot_omp3(o, 234321)

%% this script will produce a 7 color image for each round after alignment, 
% and the final spot calling results around a given point of interest selected with the cross-hair
RegionSize = 150;
OutputF = fullfile(o.OutputDirectory,'ColorRounds');mkdir(OutputF)
show7by7(o,RegionSize,OutputF,I )

