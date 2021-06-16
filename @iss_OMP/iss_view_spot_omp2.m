function SpotNo = iss_view_spot_omp2(o, FigNo, ImSz, SpotLocation, ScoreMethod, Track, SpotNum)
%% iss_view_spot_omp2(o, FigNo, ImSz, SpotLocation,ScoreMethod, Track, SpotNum)
%
% Gives the final OMP coefficients as produced by call_spots_omp about a
% chosen location. 
%
% FigNo: figure number (default, current figure)
% ImSz: radius of image that is plotted for each round and channel.
% Default value is 7 pixels.
% SpotLocation: logical,  if true, will use location of spot closest to
% crosshair, otherwise will use actual position of crosshair. Default is false.
% Track: gives plots of residual and gene coefficients at each stage of
% iteration for central pixel. 
% SpotNum: spot to look at is o.pfSpotGlobalYX(SpotNum,:) where pf
% corresponds to ScoreMethod. Can also be yx location of interest.
% You can change o.ResidualThreshParam, o.ResidualThreshMin and
% o.ResidualThreshMax to produce different coefficients. 


%%
if nargin<3 || isempty(ImSz)
    ImSz = 7;
end
if ImSz>100
    warning('ImSz too large, setting to 7');
    ImSz = 7;
end

if nargin<4 || isempty(SpotLocation)
    SpotLocation = false;
end

if nargin<6 || isempty(Track)
    Track = false;
end

if nargin>=7
    if length(SpotNum)==2
        SpotLocation = false;
        xy = [SpotNum(2),SpotNum(1)];
        S.SpotYX = o.([o.CallMethodPrefix(ScoreMethod),'SpotGlobalYX']);
        [Dist,SpotNo] = min(sum(abs(S.SpotYX-[xy(2),xy(1)]),2));
        if round(Dist)==0
            SpotLocation=true;
        end
    else
        SpotLocation = true;
        SpotNo = SpotNum;
        xy = o.([o.CallMethodPrefix(ScoreMethod),'SpotGlobalYX'])(SpotNo,[2,1]);
    end
else
    if nargin>=2
        figure(FigNo);
    end
    CrossHairColor = [1,1,1];   %Make white as black background
    xy = ginput_modified(1,CrossHairColor);
    try
        S = evalin('base', 'issPlot2DObject');
        if nargin<5 || isempty(ScoreMethod)
            ScoreMethod = S.CallMethod;
        elseif ~strcmpi(S.CallMethod,ScoreMethod)
            S.QualOK = 1;
        end
    catch
        S.QualOK = 1;
        S.Roi = [1,inf,1,inf];
    end
    S.SpotYX = o.([o.CallMethodPrefix(ScoreMethod),'SpotGlobalYX']);
    if size(S.SpotYX,1)~=size(S.QualOK,1)
        S.QualOK = 1;
    end
    %Only consider spots that can be seen in current plot
    InRoi = all(int64(round(S.SpotYX))>=S.Roi([3 1]) & round(S.SpotYX)<=S.Roi([4 2]),2);
    PlotSpots = find(InRoi & S.QualOK);        
    [Dist,SpotIdx] = min(sum(abs(S.SpotYX(PlotSpots,:)-[xy(2),xy(1)]),2));
    SpotNo = PlotSpots(SpotIdx);
    if SpotLocation || round(Dist)==0
        SpotLocation = true;
        xy = S.SpotYX(SpotNo,[2,1]);
    end   
end

if ~ismember({ScoreMethod},o.CallMethods)
    error('Method invalid, must be member of o.CallMethods.');
end
pf = o.CallMethodPrefix(ScoreMethod);


fprintf('loading channel/round images...');
if SpotLocation == false
    %Find tile that the point is on and local centered coordinates in reference round
    t = o.get_local_tile([xy(2),xy(1)]);
    SpotNo = [];
else
    t = o.([pf,'LocalTile'])(SpotNo);
end
LocalYX = [xy(2),xy(1)]-o.TileOrigin(t,:,o.ReferenceRound);

o.UseRounds = o.UseRounds(o.UseRounds<=o.nRounds);
[RoundTile,~] = get_SpotTileEachRound(o,[xy(2),xy(1)],t);
load(fullfile(o.OutputDirectory, 'FindSpotsWorkspace.mat'), 'AllBaseLocalYX');
[SpotColor,PointCorrectedLocalYX] = get_spot_colors(o,LocalYX,t,...
    RoundTile,AllBaseLocalYX);
if SpotLocation==true
    if max(max(abs(double(o.([pf,'SpotColors'])(SpotNo,:,o.UseRounds(o.UseRounds<=o.nRounds)))...
            -SpotColor(:,:,o.UseRounds(o.UseRounds<=o.nRounds)))))>0.3
        warning('Spot Color found is different from than in o object');
    end
end
 

%Get spot colors for whole grid
SpotColors = zeros((ImSz*2+1)^2,o.nBP,o.nRounds);
for r=1:o.nRounds  
    for b=1:o.nBP        
        rbYX = round(PointCorrectedLocalYX(1,:,r,b));
        y0 = rbYX(1);
        x0 = rbYX(2);
        if y0>o.TileSz || y0<1 || x0>o.TileSz || x0<1
            continue;
        end
        y1 = max(1,y0 - ImSz);
        y2 = min(o.TileSz,y0 + ImSz);
        x1 = max(1,x0 - ImSz);
        x2 = min(o.TileSz,x0 + ImSz);
        BaseIm = int32(imread(o.TileFiles{r,t}, b, 'PixelRegion',...
            {[y1 y2], [x1 x2]}))-o.TilePixelValueShift;
        if o.SmoothSize
            SE = fspecial3('ellipsoid',o.SmoothSize);
            BaseImSm = imfilter(BaseIm, SE);
        else
            BaseImSm = BaseIm;
        end
        SpotColors(:,b,r) = BaseImSm(:);               
    end
end

%Get dot product with the gene bled codes for each pixel in image
S.ImSz = ImSz;
S.ImShape = [2*S.ImSz+1,2*S.ImSz+1];
S.NormSpotColors = (double(SpotColors)-o.z_scoreSHIFT)./o.z_scoreSCALE;
S.ResidualThresh = o.ResidualThreshParam*prctile(abs(S.NormSpotColors(:,:))',47.5*100/49.0)';
S.ResidualThresh(S.ResidualThresh<o.ResidualThreshMin) = o.ResidualThreshMin;
S.ResidualThresh(S.ResidualThresh>o.ResidualThreshMax) = o.ResidualThreshMax;
nSpots = size(S.NormSpotColors,1);

S.NormBledCodes = o.ompBledCodes(:,:)./vecnorm(o.ompBledCodes(:,:),2,2);
S.NormSpotColors = S.NormSpotColors(:,:);
S.nCodes = length(o.CharCodes);
S.nBackground = o.nBackground;
S.GeneNames = o.GeneNames;
for g=S.nCodes+1:S.nCodes+S.nBackground
    S.GeneNames{g}=['Bckgrnd',num2str(g-S.nCodes)];
end
S.x0 = ImSz+1;
S.y0 = ImSz+1;

coefs = zeros(nSpots,S.nCodes+S.nBackground);
for s=1:nSpots
    coefs(s,:) = omp_free_background(S.NormBledCodes',S.NormSpotColors(s,:)',...
        o.ompMaxGenes,S.ResidualThresh(s),S.nCodes+1:S.nCodes+S.nBackground,1:S.nCodes)';
end
fprintf('\n');

PlotGenes = find(sum(coefs~=0,1)>0);
climits = [min(coefs(:)),max(coefs(:))];

%Get Spots Found by Algorithm
y_range = xy(2)-ImSz:xy(2)+ImSz;
x_range = xy(1)-ImSz:xy(1)+ImSz;
[A,B] = meshgrid(y_range,x_range);
c=cat(2,A',B');
GlobalYX=reshape(c,[],2);
InRoi = all(int64(round(o.ompSpotGlobalYX))>=min(GlobalYX,[],1) &...
    round(o.ompSpotGlobalYX)<=max(GlobalYX,[],1),2);
AlgFoundGenes = unique(o.ompSpotCodeNo(InRoi))';
PlotGenes = unique([PlotGenes,AlgFoundGenes]);
ompQualOK = o.quality_threshold('OMP');

%Get Ground Truth Spots
gtInRoi = cell(max(o.gtRounds),o.nBP);
for r=o.gtRounds
    for b=o.UseChannels
        if o.gtGeneNo(r,b)~=0
            gInRoi = all(int64(round(o.gtSpotGlobalYX{r,b}))>=min(GlobalYX,[],1) &...
                round(o.gtSpotGlobalYX{r,b})<=max(GlobalYX,[],1),2);
            gtInRoi{r,b} = gInRoi;
            if sum(gInRoi)>0
                PlotGenes = unique([PlotGenes,o.gtGeneNo(r,b)]);
            end
        end
    end
end

Small = 1e-6;
se1 = strel('disk', o.PixelDetectRadius);     %Needs to be bigger than in detect_spots
ImYX = GlobalYX-min(GlobalYX,[],1)+1;
ImInd = sub2ind(S.ImShape,ImYX(:,1),ImYX(:,2));

%Indicate if local maxima is duplicate
[AllLocalTile, ~] = which_tile(GlobalYX,...
    o.TileOrigin(:,:,o.ReferenceRound), o.TileSz);
NotDuplicate = (AllLocalTile==t);

Fig = figure(38102);
set(Fig,'Position',[336,123,820,677]);
nRows = floor((length(PlotGenes)-0.000001)/7)+1;
nCols = min(length(PlotGenes),7);
for g=1:length(PlotGenes)
    subplot(nRows, nCols, g);
    CoefIm = reshape(coefs(:,PlotGenes(g)),S.ImShape);
    imagesc(x_range,y_range,CoefIm);
    caxis(climits);
    colormap(gca,bluewhitered);
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
    hold on
    l1=scatter(o.ompSpotGlobalYX(InRoi&o.ompSpotCodeNo==PlotGenes(g)&ompQualOK,2),...
        o.ompSpotGlobalYX(InRoi&o.ompSpotCodeNo==PlotGenes(g)&ompQualOK,1),100,'+',...
        'MarkerEdgeColor',[0.5,0.5,0.5],'LineWidth',2);
    l2=scatter(o.ompSpotGlobalYX(InRoi&o.ompSpotCodeNo==PlotGenes(g)&~ompQualOK,2),...
        o.ompSpotGlobalYX(InRoi&o.ompSpotCodeNo==PlotGenes(g)&~ompQualOK,1),100,'x',...
        'MarkerEdgeColor',[0.5,0.5,0.5],'LineWidth',2);
    plot(xlim, [S.y0 S.y0], 'k'); plot([S.x0 S.x0], ylim, 'k');
    if ismember(PlotGenes(g),o.gtGeneNo)
        [r,b] = ind2sub(size(o.gtGeneNo),find(o.gtGeneNo==PlotGenes(g)));
        l3=scatter(o.gtSpotGlobalYX{r,b}(gtInRoi{r,b}&o.omp_gtFound{r,b}==2,2),...
            o.gtSpotGlobalYX{r,b}(gtInRoi{r,b}&o.omp_gtFound{r,b}==2,1),...
            30,'kv','LineWidth',2);
        l4=scatter(o.gtSpotGlobalYX{r,b}(gtInRoi{r,b}&o.omp_gtFound{r,b}==1,2),...
            o.gtSpotGlobalYX{r,b}(gtInRoi{r,b}&o.omp_gtFound{r,b}==1,1),...
            30,'k^','LineWidth',2);
        l5=scatter(o.gtSpotGlobalYX{r,b}(gtInRoi{r,b}&o.omp_gtFound{r,b}==0,2),...
            o.gtSpotGlobalYX{r,b}(gtInRoi{r,b}&o.omp_gtFound{r,b}==0,1),...
            30,'kd','LineWidth',2);
    end
    if PlotGenes(g)<=S.nCodes
        Dilate = imdilate(CoefIm, se1);
        MaxPixels = find(CoefIm + Small >= Dilate & CoefIm > 0);
        PeakInd = find(ismember(ImInd,MaxPixels));
        ndPeak = NotDuplicate(PeakInd);
        l6=scatter(GlobalYX(PeakInd(ndPeak),2),GlobalYX(PeakInd(ndPeak),1),...
            30,'gs','LineWidth',2);
        l7=scatter(GlobalYX(PeakInd(~ndPeak),2),GlobalYX(PeakInd(~ndPeak),1),...
            30,'gx','LineWidth',2);
    end
    hold off
    set(gca, 'YDir', 'normal');
    title([num2str(PlotGenes(g)),': ',S.GeneNames{PlotGenes(g)}]);
end
try
    if sum(~NotDuplicate)>0
        hL = legend([l1,l2,l3,l4,l5,l6,l7],{'ompAlg: Pass QualOK',...
            'ompAlg: Fail QualOK','GT: Missed',...
            'GT: Found','GT: Not in Set','Local Maxima with Current Thresholds',...
            'Duplicate Local Maxima'});
    else
        hL = legend([l1,l2,l3,l4,l5,l6],{'ompAlg: Pass QualOK',...
            'ompAlg: Fail QualOK','GT: Missed',...
            'GT: Found','GT: Not in Set','Local Maxima with Current Thresholds'});
    end
catch
    if sum(~NotDuplicate)>0
        hL = legend([l1,l2,l6,l7],{'ompAlg: Pass QualOK',...
            'ompAlg: Fail QualOK','Local Maxima with Current Thresholds',...
            'Duplicate Local Maxima'});
    else
        hL = legend([l1,l2,l6],{'ompAlg: Pass QualOK',...
            'ompAlg: Fail QualOK','Local Maxima with Current Thresholds'});
    end
end
newPosition = [0.4 0.94 0.2 0.05];
newUnits = 'normalized';
set(hL,'Position', newPosition,'Units', newUnits,...
    'Orientation','horizontal','color','none');
set( hL, 'Box', 'off' ) ;

%% Give plots of omp parameters at each iteration. 
if Track
    NormSpotColor = (double(SpotColor)-o.z_scoreSHIFT)./o.z_scoreSCALE;
    ResidualThresh = o.ResidualThreshParam*prctile(abs(NormSpotColor(:,:))',47.5*100/49.0)';
    ResidualThresh(ResidualThresh<o.ResidualThreshMin) = o.ResidualThreshMin;
    ResidualThresh(ResidualThresh>o.ResidualThreshMax) = o.ResidualThreshMax;
    o.omp_free_background_track(S.NormBledCodes',NormSpotColor(:,:)',...
        o.ompMaxGenes,ResidualThresh,S.nCodes+1:S.nCodes+S.nBackground);
end
end

