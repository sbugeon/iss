function CellCalled = DoCellCalling(o,DapiBound,gSet,Ineff)

% cell calling parameters
o.CellCallMaxIter = 100;
o.nNeighbors = 3;
o.CellCallTolerance = 0.05;
o.ompScoreThresh = 5;
o.ompScoreThresh2 = 2;
o.ompIntensityThresh = 0.01;
o.ompIntensityThresh2 = 0.005;
o.ompNeighbThresh = 12;
o.ompNeighbThresh2 = 10;
o.rSpot = 2;
o.Inefficiency = Ineff;

%
LocalFolder = o.OutputDirectory;
BigDapi = imadjust(imread(fullfile(LocalFolder,'background_image_fixed.tif')));
GcampImg = imadjust(imread(fullfile(LocalFolder,'tdTomato_image_fixed.tif')));
ypoly = [size(BigDapi,1); size(BigDapi,1); 1; 1; size(BigDapi,1) ];
xpoly = [1; size(BigDapi,2); size(BigDapi,2); 1; 1];
o.CellCallRegionYX = round([ypoly xpoly]);
Boundaries = {};
TomatoIntensity = NaN( length(DapiBound),1);
CellMap = zeros(size(BigDapi,1),size(BigDapi,2));

% convert dapi boundaries to cellmap for cell calling
for iROI = 1: length(DapiBound)
    if ~isempty(DapiBound{iROI}) & ~iscell(DapiBound{iROI})
        BB = boundary(DapiBound{iROI}(:,1),DapiBound{iROI}(:,2),0.1);
        Boundaries{iROI} = DapiBound{iROI}(BB,:);
        X1R = min(Boundaries{iROI}(:,1));
        X2R = max(Boundaries{iROI}(:,1));
        Y1R = min(Boundaries{iROI}(:,2));
        Y2R = max(Boundaries{iROI}(:,2));
        Y1 = round(Y1R );
        Y2 = round(Y2R );
        X1 = round(X1R );
        X2 = round(X2R );
        X2 = min(X2,size(CellMap,2));
        Y2 = min(Y2,size(CellMap,1));
        X1 = max(X1,1);
        Y1 = max(Y1,1);
        CellMap(Y1:Y2,X1:X2) = iROI;
        TomatoIntensity(iROI) = nanmean(GcampImg(Y1:Y2,X1:X2),'all');
        % hold on
        % plot(Boundaries{iROI}(:,1),Boundaries{iROI}(:,2))
    end
end
% figure
% imshow(cat(3,imadjust(BigDapi),imadjust(uint16(CellMap)),imadjust(GcampImg)))

% perform cell calling
[o,CellCalled] = o.call_cells(gSet,Boundaries,CellMap,'OMP');
CellCalled.TomatoI = TomatoIntensity;

%% PLOTING
if length(gSet.Class) > 110
    o.ClassCollapse = {o.ClassNames([1:4,365:387]), 'Not neuron', [138, 134, 134]./255 ; ...% gray
        o.ClassNames(124:317), 'PC.Cortex', [233, 237, 26]./255 ; ...% yellow
        o.ClassNames(5:123), 'Interneurons', [240, 2, 2]./255 ; ...% red
        o.ClassNames(find(contains(CellCalled.ClassNames,'ProS') | contains(CellCalled.ClassNames,'SUB'))), 'Subiculum', [237, 26, 160]./255 ; ... %pink
        o.ClassNames(329:348), 'PC.CA1', [26, 234, 237]./255 ; ...% cyan
        o.ClassNames(359:360), 'PC.CA2', [26, 110, 237]./255 ; ... % blue
        o.ClassNames(351:358), 'PC.CA3', [146, 26, 237]./255 ; ...% purple
        o.ClassNames(349:350), 'MossyCells', [237, 170, 26]./255 ; ... % orange
        o.ClassNames(361:364), 'DG', [26, 237, 75]./255 ; ... % green
        o.ClassNames(388), 'Zero', [0 0 0]};%black
else

    o.ClassCollapse = {o.ClassNames(95:109), 'Not neuron', [0 0 0]./255 ; ...
        o.ClassNames(61:63), 'Layer2-3', [250, 169, 17]./255 ; ...
        o.ClassNames(64), 'Layer4', [5, 207, 230]./255 ; ...
        o.ClassNames([65:69,75:81]), 'Layer5', [166, 77, 235]./255 ; ...
        o.ClassNames([70:74,82:92]), 'Layer6', [250, 17, 130]./255 ; ...
        o.ClassNames(1:60), 'Interneurons', [0 0 0] ; ...
        o.ClassNames(93:94), 'others', [0 0 0] ; ...
        o.ClassNames(110), 'Zero', [0 0 0]};
end

o = pie_plot(o,Boundaries);
daspect([1 1 1])
set(gcf,'Position',[28,185,2362,1153])
set(gcf,'InvertHardCopy','off');

