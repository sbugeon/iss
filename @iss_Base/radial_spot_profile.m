function [RadialProfile,PlotLogical,OutSpotImages] = radial_spot_profile(o,GeneNo,ImRad,Method,SetInput,...
    OutImageRound,OutImageChannel,AllBaseLocalYX)
%% [RadialProfile,PlotLogical,OutSpotImages] = o.radial_spot_profile(GeneNo,ImRad,Method,SetInput,...
%    OutImageRound,OutImageChannel,AllBaseLocalYX);
% This gives the radial intensity profile about a spot of a given gene in
% each round/channel.
% Inputs
% o: iss object
% GeneNo: list of genes want profile of (If have multiple Sets, only
% allowed one gene).
% ImRad: Max radius of radial profile plot.
% Method: 'DotProduct', 'Prob' or 'Pixel' i.e. which spot data the SetInput
% logical is referring to.
% SetInput: SetInput(:,i) is the ith logical set of spots to build plot
% from. Will produce plot for each i. If multiple genes, only allowed one
% Set.
% OutImageRound, OutImageChannel: Can specify a round/channel you would
% like to retrieve all the spot images from.
% AllBaseLocalYX: Saved in find_spots_workspace, need for registration.
% Outputs
% RadialProfile(s,b,r,:) is the radial profile for spot s.
% PlotLogical(s,i) indicates whether spot s was included in the ith plot. i
% refers either to different genes or different SetInput.
% OutSpotImages: the output spots if OutImageRound and OutImageChannel
% specified.


if nargin<3 || isempty(ImRad)
    ImRad = 10;
end
if nargin<4 || isempty(Method)
    Method = 'Pixel';
end
pf = o.CallMethodPrefix(Method);
if nargin<5 || isempty(SetInput)
    SetInput = o.quality_threshold(Method);
end
if nargin<8 || isempty(AllBaseLocalYX)
    load(fullfile(o.OutputDirectory, 'FindSpotsWorkspace.mat'), 'AllBaseLocalYX');
end
Set = max(SetInput,[],2)&ismember(o.([pf,'SpotCodeNo']),GeneNo);
AllBaseSpotNo = cell2mat(cellfun(@size,AllBaseLocalYX,'uni',false));
AllBaseSpotNo = AllBaseSpotNo(:,1:2:o.nBP*2,:);
nSpots = sum(Set);
nSets = size(SetInput,2);
nGenes = length(GeneNo);
if nSets>1 && nGenes>1
    error('Require either number of genes or number of thresholding sets to be 1');
end
nPlots = max(nGenes,nSets);

if nargin<6 || isempty(OutImageRound) || isempty(OutImageChannel)
    OutImage = false;
else
    OutImage = true;
    OutSpotImages = zeros(nSpots,2*ImRad+1,2*ImRad+1);
end

LocalTile = o.([pf,'LocalTile'])(Set);
GlobalYX = o.([pf,'SpotGlobalYX'])(Set,:);
ImageTiles = unique(LocalTile)';
RoundTile = o.get_SpotTileEachRound(GlobalYX,LocalTile);
LocalYX = [GlobalYX-o.TileOrigin(LocalTile,:,o.ReferenceRound)-o.TileCentre,ones(nSpots,1)];
RadialBins = 0:ImRad;
RadialProfile = nan(nSpots,o.nBP,o.nRounds,ImRad+1);
numCharCode = zeros(nGenes,o.nRounds);
PlotLogical = false(nSpots,nGenes);
LegendNames = cell(nPlots,1);
for i=1:nPlots
    if nGenes>1 || nSets==1
        numCharCode(i,:) = str2double(regexp(cell2mat(o.CharCodes(GeneNo(i))),'\d','match'))+1;
        PlotLogical(:,i) = o.([pf,'SpotCodeNo'])(Set)==GeneNo(i);
        LegendNames{i} = [o.GeneNames{GeneNo(i)}, ': ',...
            num2str(sum(PlotLogical(:,i))), ' Spots'];
    else
        numCharCode = str2double(regexp(cell2mat(o.CharCodes(GeneNo)),'\d','match'))+1;
        PlotLogical(:,i) = SetInput(Set,i);
        LegendNames{i} = [o.GeneNames{GeneNo},' ',num2str(i), ': ',...
            num2str(sum(PlotLogical(:,i))), ' Spots'];
    end
end

SpotImYX = zeros((ImRad*2+1)^2,2);
SpotImYX(:,1) = repelem(-ImRad:ImRad,1,ImRad*2+1);
SpotImYX(:,2) = repmat(-ImRad:ImRad,1,ImRad*2+1);

for t=ImageTiles
    for r=o.UseRounds
        % open file for this tile/round
        FileName = o.TileFiles{r,t};
        TifObj = Tiff(FileName);
        % find the home tile for all current spots in the ref round
        RefRoundHomeTiles = LocalTile(RoundTile(:,r)==t);
        MyRefTiles = unique(RefRoundHomeTiles);
        
        fprintf('Tile %d, Round %d: %d spots\n', t, r,length(RefRoundHomeTiles));
        
        for b=o.UseChannels
            TifObj.setDirectory(o.FirstBaseChannel + b - 1);
            BaseIm = int32(TifObj.read())-o.TilePixelValueShift;
            if o.SmoothSize
                BaseImSm = double(round(imfilter(double(BaseIm), fspecial('disk', o.SmoothSize))));
            else
                BaseImSm = double(BaseIm);
            end
            for t2 = MyRefTiles(:)'
                MyBaseSpots = (RoundTile(:,r)==t & LocalTile==t2);
                MyLocalYX = LocalYX(MyBaseSpots,:);
                if t == t2
                    MyPointCorrectedYX = MyLocalYX*o.D(:,:,t,r,b)+o.TileCentre;
                    MyPointCorrectedYX = round(MyPointCorrectedYX);
                else
                    [MyPointCorrectedYX, ~, nMatches_diff_tile] = ...
                        o.different_tile_transform(AllBaseLocalYX,o.RawLocalYX, MyLocalYX,t,t2,r,b);
                    if isempty(nMatches_diff_tile) || nMatches_diff_tile<o.MinPCMatchFract*AllBaseSpotNo(t,b,r)
                        MyPointCorrectedYX = nan(sum(MyBaseSpots),2);
                    end
                end
                PointImageYX = repelem(MyPointCorrectedYX,(ImRad*2+1)^2,1)+repmat(SpotImYX,sum(MyBaseSpots),1);
                SpotImages = reshape(IndexArrayNan(BaseImSm, PointImageYX'),ImRad*2+1,ImRad*2+1,sum(MyBaseSpots));
                RadialProfile(MyBaseSpots,b,r,:) = radialAverage(SpotImages,ImRad+1,ImRad+1,RadialBins);
                if OutImage && r == OutImageRound && b ==OutImageChannel
                    OutSpotImages(MyBaseSpots,:,:) = permute(SpotImages,[3,1,2]);
                end
            end
        end
        TifObj.close();
    end
    fprintf('\nDone %d/%d spots\n',sum(~isnan(RadialProfile(:,numCharCode(1,1),1,1))),nSpots);
end

MeanProfile = zeros(o.nBP,o.nRounds,ImRad+1,nPlots);
StdProfile = zeros(o.nBP,o.nRounds,ImRad+1,nPlots);
Colors = distinguishable_colors(nPlots,[1,1,1]);
for i=1:nPlots
    MeanProfile(:,:,:,i) = squeeze(mean(RadialProfile(PlotLogical(:,i),:,:,:),1));
    StdProfile(:,:,:,i) = squeeze(std(RadialProfile(PlotLogical(:,i),:,:,:),[],1));
end
YLim = [min(min(min(MeanProfile-StdProfile/2,[],4),[],3),[],2),...
    max(max(max(MeanProfile+StdProfile/2,[],4),[],3),[],2)];
LinePlots = cell(nPlots,1);

try
    clf(545686)
    figure(545686)
catch
    figure(545686)
end
set(gcf,'Position',[164,108,1621,805]);
Ylegends = {o.bpLabels{:}};
Xlegends = string(1:o.nRounds);
for r=o.UseRounds
    for b=o.UseChannels
        h = subplot(o.nBP, o.nRounds, (b-1)*o.nRounds + r);
        if r == 1 && b == 1
            Pos1 = get(h,'position');
        elseif r == 1 && b == o.nBP
            Pos2 = get(h,'position');
        elseif r == o.nRounds && b == o.nBP
            Pos3 = get(h,'position');
        end
        hold on       
        for i=1:nPlots
            brMean = squeeze(MeanProfile(b,r,:,i));
            brStd = squeeze(StdProfile(b,r,:,i))/2;
            fill([RadialBins';flipud(RadialBins')],...
                [brMean-brStd;flipud(brMean+brStd)],...
                Colors(i,:),'linestyle','none','facealpha',.03);
            LinePlots{i} = plot(RadialBins,brMean,'Color',Colors(i,:),...
                'LineWidth',1, 'DisplayName', LegendNames{i});
        end       
        hold off
        ylim(YLim(b,:));
        if sum(numCharCode(:,r)==b)>0
            ax = gca;
            ax.Box = true;
            ax.XColor = mean(Colors(find(numCharCode(:,r)==b),:),1);
            ax.YColor = mean(Colors(find(numCharCode(:,r)==b),:),1);
            ax.LineWidth=1;
        end
        if r==1
            ylabel(Ylegends{b},'Color',[0.15 0.15 0.15]); 
        else
            set(gca,'yticklabel',[]);
        end
        if b==o.nBP
            xlabel(Xlegends(r),'Color',[0.15 0.15 0.15]); 
        else
            set(gca,'xticklabel',[]);
        end
    end
end
PosDev = 0.02;
SuperAxisPos = [Pos2(1:2)-PosDev,Pos3(1)+Pos3(2)-Pos2(1)+PosDev*2,Pos1(2)+Pos1(4)-Pos3(2)+PosDev*2];
hSuper=axes('position',SuperAxisPos,'visible','off');
hSuper.XLabel.Visible='on';
hSuper.YLabel.Visible='on';
axes(hSuper);
ylabel('Channel');
xlabel('Round');
LegendData = [];
for i=1:nPlots
    LegendData = [LegendData,LinePlots{i}];
end
lg = legend(hSuper,LegendData,'Orientation','horizontal','Location','northoutside');
legend('boxoff');


% %Eigenvector analysis
% NanSum=sum(isnan(SpotImages(:,:)),2);
% EigSpotImages = nan(21,21,9043);
% for g=1:size(GeneLogical,2)
%     [TopEvec, s2] = eigs(double(SpotImages(GeneLogical(:,g)&NanSum==0,:)'*...
%         SpotImages(GeneLogical(:,g)&NanSum==0,:))/sum(GeneLogical(:,g)&NanSum==0));
%     if g<=2
%         Coef = TopEvec(:,2)'*SpotImages(GeneLogical(:,g),:)';
%         EigIm = reshape(TopEvec(:,2),[21,21]);
%     else
%         Coef = TopEvec(:,1)'*SpotImages(GeneLogical(:,g),:)';
%         EigIm = reshape(TopEvec(:,1),[21,21]);
%     end
%     EigSpotImages(:,:,GeneLogical(:,g))=permute(Coef'.*...
%         permute(repmat(EigIm,1,1,length(Coef)),[3,1,2]),[2,3,1]);
% end
% EigRadialProfile = radialAverage(EigSpotImages,11,11,0:10);
% EigRadialProfile(EigRadialProfile(:,1)<0,:) = nan;
%% Below replicates radial plot in OutSpotImages round/channel using 
% only eigenvector component
% SpotImages = OutSpotImages;
% NanSum=sum(isnan(SpotImages(:,:)),2);
% EigSpotImages = nan(size(SpotImages,[2,3,1]));
% GeneLogical = PlotLogical;
% EigUse = 1;
% for g=1:size(GeneLogical,2)
%     [TopEvec, ~] = eigs(double(SpotImages(NanSum==0,:)'*...
%         SpotImages(NanSum==0,:))/sum(NanSum==0));
%     Coef = TopEvec(:,EigUse)'*SpotImages(GeneLogical(:,g),:)';
%     EigIm = reshape(TopEvec(:,EigUse),[2*ImRad+1,2*ImRad+1]);
% %     if gEig<=0
% %         Coef = TopEvec(:,2)'*SpotImages(GeneLogical(:,g),:)';
% %         EigIm = reshape(TopEvec(:,2),[21,21]);
% %     else
% %         Coef = TopEvec(:,1)'*SpotImages(GeneLogical(:,g),:)';
% %         EigIm = reshape(TopEvec(:,1),[21,21]);
% %     end
%     EigSpotImages(:,:,GeneLogical(:,g))=permute(Coef'.*...
%         permute(repmat(EigIm,1,1,length(Coef)),[3,1,2]),[2,3,1]);
% end
% EigRadialProfile = radialAverage(EigSpotImages,ImRad+1,ImRad+1,0:ImRad);
% EigRadialProfile(EigRadialProfile(:,1)<0,:) = nan;
% MeanProfile = zeros(ImRad+1,4);
% StdProfile = zeros(ImRad+1,4);
% Colors = distinguishable_colors(4,[1,1,1]);
% for i=1:4
%     MeanProfile(:,i) = squeeze(nanmean(EigRadialProfile(GeneLogical(:,i),:,:,:),1));
%     StdProfile(:,i) = squeeze(nanstd(EigRadialProfile(GeneLogical(:,i),:,:,:),[],1));
% end
% figure;
% hold on
% RadialBins = 0:ImRad;
% for i=1:4
%     brMean = MeanProfile(:,i);
%     brStd = StdProfile(:,i)/2;
%     fill([RadialBins';flipud(RadialBins')],...
%         [brMean-brStd;flipud(brMean+brStd)],...
%         Colors(i,:),'linestyle','none','facealpha',.03);
%     LinePlots = plot(RadialBins,brMean,'Color',Colors(i,:),...
%         'LineWidth',1);
% end
end

