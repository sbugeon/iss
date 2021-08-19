function [SpotNo,ViewOMP3_Info] = iss_view_spot_omp2(o, FigNo, ImSz, SpotLocation, ScoreMethod, Track, SpotNum)
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

if nargin<5
    ScoreMethod = [];
end

if nargin<7
    SpotNum = [];
end

[xy, SpotLocation, ScoreMethod, SpotNo, Dist]  = ...
    get_crosshair_location(o, FigNo, SpotLocation, ScoreMethod, SpotNum);
pf = o.CallMethodPrefix(ScoreMethod);
[SpotColors, PointCorrectedLocalYX] = get_spot_colors_grid(o, pf, xy, ImSz, SpotNo,...
    SpotLocation);


%Get dot product with the gene bled codes for each pixel in image
S.ImSz = ImSz;
S.ImShape = [2*S.ImSz+1,2*S.ImSz+1];
S.NormSpotColors = (double(SpotColors)-o.z_scoreSHIFT)./o.z_scoreSCALE;
S.NormSpotColors = S.NormSpotColors(:,:);
S.nCodes = length(o.CharCodes);
S.nBackground = o.nBackground;
S.GeneNames = o.GeneNames;
for g=S.nCodes+1:S.nCodes+S.nBackground
    S.GeneNames{g}=['Bckgrnd',num2str(g-S.nCodes)];
end
S.x0 = ImSz+1;
S.y0 = ImSz+1;


%PlotGenes = find(sum(coefs~=0,1)>0);
%climits = [min(coefs(:)),max(coefs(:))];

%Get Spots Found by Algorithm
y_range = xy(2)-ImSz:xy(2)+ImSz;
x_range = xy(1)-ImSz:xy(1)+ImSz;
[A,B] = meshgrid(y_range,x_range);
c=cat(2,A',B');
GlobalYX=reshape(c,[],2);

if isprop(o,'ompGetCoefMethod')
    if o.ompGetCoefMethod == 1
        coefs = o.get_omp_coefs(S.NormSpotColors);
    elseif o.ompGetCoefMethod == 2
        coefs = o.get_omp_coefs2(S.NormSpotColors, GlobalYX);
    end
else
    coefs = o.get_omp_coefs(S.NormSpotColors);
end
    
PlotGenes = find(sum(coefs~=0,1)>ImSz/2);
climits = [min(coefs(:)),max(coefs(:))];

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
if SpotLocation == false
    %Find tile that the point is on and local centered coordinates in reference round
    t = o.get_local_tile([xy(2),xy(1)]);
else
    t = o.([pf,'LocalTile'])(SpotNo);
end
[AllLocalTile, ~] = which_tile(GlobalYX,...
    o.TileOrigin(:,:,o.ReferenceRound), o.TileSz);
NotDuplicate = (AllLocalTile==t);

if ishandle(38102)
    Fig = figure(38103);
else 
    Fig = figure(38102);
end
set(Fig,'Position',[336,123,820,677]);
nRows = floor((length(PlotGenes)-0.000001)/7)+1;
nCols = min(length(PlotGenes),7);
for g=1:length(PlotGenes)
    subplot(nRows, nCols, g);
    CoefIm = reshape(coefs(:,PlotGenes(g)),S.ImShape);
    if nargout>1
        % called from iss_view_spot_omp3
        imagesc(x_range,y_range,CoefIm, 'ButtonDownFcn', {@update_view_omp3, PlotGenes(g), o});
    else
        imagesc(x_range,y_range,CoefIm);
    end
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
    try
        if sum(~NotDuplicate)>0
            hL = legend([l1,l2,l6,l7],{'ompAlg: Pass QualOK',...
                'ompAlg: Fail QualOK','Local Maxima with Current Thresholds',...
                'Duplicate Local Maxima'});
        else
            hL = legend([l1,l2,l6],{'ompAlg: Pass QualOK',...
                'ompAlg: Fail QualOK','Local Maxima with Current Thresholds'});
        end
    catch
    end
end
newPosition = [0.4 0.04 0.2 0.05];
newUnits = 'normalized';
try
    set(hL,'Position', newPosition,'Units', newUnits,...
        'Orientation','horizontal','color','none');
    set( hL, 'Box', 'off' ) ;
catch
end
iss_view_plot_title(o, ScoreMethod, SpotLocation, SpotNo);

%% Give plots of omp parameters at each iteration. 
if Track
    NormSpotColor = (double(SpotColor)-o.z_scoreSHIFT)./o.z_scoreSCALE;
    ResidualThresh = o.ResidualThreshParam*prctile(abs(NormSpotColor(:,:))',47.5*100/49.0)';
    ResidualThresh(ResidualThresh<o.ResidualThreshMin) = o.ResidualThreshMin;
    ResidualThresh(ResidualThresh>o.ResidualThreshMax) = o.ResidualThreshMax;
    o.omp_free_background_track(o.ompBledCodes',NormSpotColor(:,:)',...
        o.ompMaxGenes,ResidualThresh,S.nCodes+1:S.nCodes+S.nBackground);
end

ViewOMP3_Info.SpotCodeNo = o.([pf,'SpotCodeNo'])(SpotNo);
ViewOMP3_Info.SpotColors = SpotColors;
ViewOMP3_Info.PointCorrectedLocalYX = PointCorrectedLocalYX;
ViewOMP3_Info.ImSz = ImSz;
ViewOMP3_Info.Dist = Dist;
ViewOMP3_Info.coefs = coefs;
ViewOMP3_Info.nCodes = S.nCodes;
end

function S = update_view_omp3(aH,~,idx,o)
fig = ancestor(aH,'figure');
click_type = get(fig,'SelectionType');
S = evalin('base', 'issViewSpotOMP3Object');
if idx>S.nCodes
    % Add all background or remove all background all at once
    idx = S.nCodes+1:S.nCodes+o.nBackground;
end
if strcmp(click_type,'normal')
    S.SelectGenes = unique([S.SelectGenes,idx]);
elseif strcmp(click_type,'alt')
    S.SelectGenes = setdiff(S.SelectGenes, idx);
end
NormSpotColors = (double(S.SpotColors)-o.z_scoreSHIFT)./o.z_scoreSCALE;
for g=S.SelectGenes
    NormSpotColors = NormSpotColors - ...
        permute(repmat(reshape(o.ompBledCodes(g,:),[o.nBP,o.nRounds]),...
        1,1,size(S.coefs,1)),[3,1,2]).*S.coefs(:,g);
end
S.SpotColorsCurrent = NormSpotColors.*o.z_scoreSCALE + o.z_scoreSHIFT;
figure(S.FigNo);
if S.Norm
    plot_spot_colors_grid(o, NormSpotColors, S.PointCorrectedLocalYX,...
        S.ImSz, S.Dist, S.SpotCodeNo, S.Clim);
else
    plot_spot_colors_grid(o, S.SpotColorsCurrent, S.PointCorrectedLocalYX,...
        S.ImSz, S.Dist, S.SpotCodeNo, S.Clim);
end
GeneNames = o.GeneNames;
GeneNames{S.nCodes+1} = 'Background';
if isempty(S.SelectGenes)
    sgtitle('Original SpotColors');
else
    nGenes = length(S.SelectGenes(S.SelectGenes<=S.nCodes+1));
    title = ['SpotColors after removing ', GeneNames{S.SelectGenes(1)}];
    if nGenes>1
        for g=2:nGenes-1
            title = [title, ', ', GeneNames{S.SelectGenes(g)}];
        end
        title = [title, ' and ', GeneNames{S.SelectGenes(nGenes)}];
    end
    sgtitle(title);
end
assignin('base','issViewSpotOMP3Object',S);
end

