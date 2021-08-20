function [RadialProfile,PlotLogical] = radial_spot_profile_Z(o,GeneNo,ImRad,Method,SetInput,...
    r,b,ImRadZ,AllBaseLocalYX)
%% [RadialProfile,PlotLogical,OutSpotImages] = o.radial_spot_profile(GeneNo,ImRad,Method,SetInput,...
%    OutImageRound,OutImageChannel,AllBaseLocalYX);
% This gives the radial intensity profile about a spot of a given gene in
% each round/channel.
% Inputs
%   o: iss object
%   GeneNo: list of genes want profile of (If have multiple Sets, only
%   allowed one gene).
%   ImRad: Max radius of radial profile plot.
%   Method: 'DotProduct', 'Prob' or 'Pixel' i.e. which spot data the SetInput
%   logical is referring to.
%   SetInput: SetInput(:,i) is the ith logical set of spots to build plot
%       from. Will produce plot for each i. If multiple genes, only allowed one
%       Set. Set can also be list of spot index of interest. 
%   OutImageRound, OutImageChannel: Can specify a round/channel you would
%       like to retrieve all the spot images from.
%   AllBaseLocalYX: Saved in find_spots_workspace, need for registration.
% Outputs
%   RadialProfile(s,b,r,:) is the radial profile for spot s.
%   PlotLogical(s,i) indicates whether spot s was included in the ith plot. i
%       refers either to different genes or different SetInput.
%   OutSpotImages: the output spots if OutImageRound and OutImageChannel specified.


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
if nargin<8 || isempty(ImRadZ)
    ImRadZ=3;
end
if nargin<9 || isempty(AllBaseLocalYX)
    load(fullfile(o.OutputDirectory, 'FindSpotsWorkspace.mat'), 'AllBaseLocalYX');
end

if length(SetInput)~=length(o.([pf,'SpotCodeNo']))
    %This is if Set given as list of spots. 
    GeneNo = o.([pf,'SpotCodeNo'])(SetInput(1));
    SpotIndices = SetInput;
    SetInput = false(length(o.([pf,'SpotCodeNo'])),length(SpotIndices));
    for i=1:length(SpotIndices)
       SetInput(SpotIndices(i),i)=true;
    end   
    Set = max(SetInput,[],2);
else
    Set = max(SetInput,[],2)&ismember(o.([pf,'SpotCodeNo']),GeneNo);
end
    
AllBaseSpotNo = cell2mat(cellfun(@size,AllBaseLocalYX,'uni',false));
AllBaseSpotNo = AllBaseSpotNo(:,1:2:o.nBP*2,:);
nSpots = sum(Set);
nSets = size(SetInput,2);
nGenes = length(GeneNo);
if nSets>1 && nGenes>1
    error('Require either number of genes or number of thresholding sets to be 1');
end
nPlots = max(nGenes,nSets);

LocalTile = o.([pf,'LocalTile'])(Set);
GlobalYX = o.([pf,'SpotGlobalYX'])(Set,:);
ImageTiles = unique(LocalTile)';
RoundTile = o.get_SpotTileEachRound(GlobalYX,LocalTile);
LocalYX = [GlobalYX-o.TileOrigin(LocalTile,:,o.ReferenceRound)-o.TileCentre,ones(nSpots,1)];
PlotLogical = false(nSpots,nGenes);
LegendNames = cell(nPlots,1);
for i=1:nPlots
    if nGenes>1 || nSets==1
        PlotLogical(:,i) = o.([pf,'SpotCodeNo'])(Set)==GeneNo(i);
        LegendNames{i} = [o.GeneNames{GeneNo(i)}, ': ',...
            num2str(sum(PlotLogical(:,i))), ' Spots'];
    else
        PlotLogical(:,i) = SetInput(Set,i);
        LegendNames{i} = [o.GeneNames{GeneNo},' ',num2str(i), ': ',...
            num2str(sum(PlotLogical(:,i))), ' Spots'];
    end
    if exist('SpotIndices','var')
        LegendNames{i} = ['SpotNo: ', num2str(SpotIndices(i)),', Gene: ',...
            o.GeneNames{o.([pf,'SpotCodeNo'])(SpotIndices(i))}];
    end
end
%Spot coordinates about central spot
SpotImYX = zeros((ImRad*2+1)^2,2);
SpotImYX(:,1) = repelem(-ImRad:ImRad,1,ImRad*2+1);
SpotImYX(:,2) = repmat(-ImRad:ImRad,1,ImRad*2+1);

% construct a Bio-Formats reader with the Memoizer wrapper
imfile = fullfile(o.InputDirectory, [o.FileBase{r}, o.RawFileExtension]);  %raw data file name for round r 
bfreader = loci.formats.Memoizer(bfGetReader(), 0);
bfreader.setId(imfile);
% get some basic image metadata
[nSeries, nSerieswPos, nChannels, nZstacks, xypos, pixelsize] = ...
    get_ome_tilepos(bfreader);
scene = nSeries/nSerieswPos;
bfreader.close();

SpotImagesCell = cell(ImRadZ*2+1,1);
for z2=1:ImRadZ*2+1
    SpotImagesCell{z2} = zeros(ImRad*2+1,ImRad*2+1,nSpots);
end
SpotZPlane = zeros(nSpots,1);

for t=ImageTiles
        % find the home tile for all current spots in the ref round
        RefRoundHomeTiles = LocalTile(RoundTile(:,r)==t);
        MyRefTiles = unique(RefRoundHomeTiles);
        
        fprintf('Tile %d, Round %d: %d spots\n', t, r,length(RefRoundHomeTiles));
        
        for t2 = MyRefTiles(:)'
            t_rawdata = str2double(cell2mat(extractBetween(o.TileFiles{r,t2},'_t','.tif')));
            bfreader = javaObject('loci.formats.Memoizer', bfGetReader(), 0);
            % use the memo file cached before
            bfreader.setId(imfile);
            
            bfreader.setSeries(scene*t_rawdata-1);
            BaseImSm = zeros(o.TileSz,o.TileSz,nZstacks);
            for z = 1:nZstacks
                iPlane = bfreader.getIndex(z-1, b-1, 0)+1;
                BaseImSm(:,:,z) = bfGetPlane(bfreader, iPlane);
            end
            MyBaseSpots = (RoundTile(:,r)==t & LocalTile==t2);
            tSpotIndex = find(MyBaseSpots);
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
            %Find z-plane each spot is on
            YIndex = repelem(MyPointCorrectedYX(:,1),nZstacks,1);
            XIndex = repelem(MyPointCorrectedYX(:,2),nZstacks,1);
            ZIndex = repmat((1:nZstacks)',sum(MyBaseSpots),1);
            ZProfileInd = sub2ind(size(BaseImSm),YIndex,XIndex,ZIndex);
            ZProfile = reshape(BaseImSm(ZProfileInd),nZstacks,sum(MyBaseSpots))';
            [~,MaxZ] = max(ZProfile(:,1+ImRadZ:nZstacks-ImRadZ),[],2);
            MaxZ = MaxZ+ImRadZ;
            SpotZPlane(tSpotIndex) = MaxZ;
            
            %For each spot, save spot image between z=MaxZ-ImRadZ:MaxZ+ImRadZ 
            %z is actual z-plane. z2 is relative zplane of each spot.
            Z_Spots = min(MaxZ-ImRadZ):max(MaxZ+ImRadZ);            
            for z=Z_Spots                
                MyZSpotsZ2Index = z-MaxZ+ImRadZ+1;
                MyZSpots = MyZSpotsZ2Index>=1 & MyZSpotsZ2Index<=ImRadZ*2+1;
                PointImageYX = repelem(MyPointCorrectedYX(MyZSpots,:),(ImRad*2+1)^2,1)+repmat(SpotImYX,sum(MyZSpots),1);
                zSpotImages = reshape(IndexArrayNan(BaseImSm(:,:,z), PointImageYX'),ImRad*2+1,ImRad*2+1,sum(MyZSpots));
                for z2 = unique(MyZSpotsZ2Index(MyZSpots))'
                    MyZ2Spots = MyZSpotsZ2Index(MyZSpots)==z2;
                    SpotImagesCell{z2}(:,:,tSpotIndex(MyZSpotsZ2Index==z2)) = zSpotImages(:,:,MyZ2Spots);
                end
            end
        end
    fprintf('\nDone %d/%d spots\n',sum(SpotImagesCell{ImRadZ+1}(ImRad+1,ImRad+1,:)~=0),nSpots);
end
SpotImages = zeros(nSpots,ImRad*2+1,ImRad*2+1,ImRadZ*2+1);
for z2=1:ImRadZ*2+1
    SpotImages(:,:,:,z2) = permute(SpotImagesCell{z2},[3,1,2]);
end

%% Plot eigenvector

NanSum=sum(isnan(SpotImages(:,:)),2);
EigNo = 1;

if ~exist('SpotIndices','var')
    [TopEvec, s2] = eigs(double(SpotImages(NanSum==0,:)'*...
        SpotImages(NanSum==0,:))/sum(NanSum==0));
    EigImage = reshape(TopEvec(:,EigNo),ImRad*2+1,ImRad*2+1,ImRadZ*2+1)*sqrt(s2(EigNo,EigNo));
    EigImage = EigImage*sign(EigImage(ImRad+1,ImRad+1,ImRadZ+1));   %Make positive at centre.
    show_image_eigenvector(o,EigImage,'Best Eigenvector found using all genes');
end

if nPlots>1 || exist('SpotIndices','var')
    for i=1:nPlots
        [TopEvec, s2] = eigs(double(SpotImages(PlotLogical(:,i)&NanSum==0,:)'*...
            SpotImages(PlotLogical(:,i)&NanSum==0,:))/sum(PlotLogical(:,i)&NanSum==0));
        EigImage = reshape(TopEvec(:,EigNo),ImRad*2+1,ImRad*2+1,ImRadZ*2+1)*sqrt(s2(EigNo,EigNo));
        EigImage = EigImage*sign(EigImage(ImRad+1,ImRad+1,ImRadZ+1));   %Make positive at centre.
        show_image_eigenvector(o,EigImage,LegendNames{i});
    end
end

%% Plot radial profile in mid z-plane
% Also plot dependence on central pixel with z.
RadialBins = 0:ImRad;
ZBins = -ImRadZ:ImRadZ;
RadialProfile = radialAverage(permute(SpotImages(:,:,:,ImRadZ+1),[2,3,1]),ImRad+1,ImRad+1,RadialBins);
MeanProfile = zeros(ImRad+1,nPlots);
StdProfile = zeros(ImRad+1,nPlots);
ZProfile = squeeze(SpotImages(:,ImRad+1,ImRad+1,:));
MeanZProfile = zeros(ImRadZ*2+1,nPlots);
StdZProfile = zeros(ImRadZ*2+1,nPlots);
Colors = distinguishable_colors(nPlots,[1,1,1]);
for i=1:nPlots
    MeanProfile(:,i) = mean(RadialProfile(PlotLogical(:,i),:),1);
    StdProfile(:,i) = std(RadialProfile(PlotLogical(:,i),:),[],1);
    MeanZProfile(:,i) = mean(ZProfile(PlotLogical(:,i),:),1);
    StdZProfile(:,i) = std(ZProfile(PlotLogical(:,i),:),[],1);
end


figure;
subplot(2,1,1);
hold on
for i=1:nPlots
    fill([RadialBins';flipud(RadialBins')],...
        [MeanProfile(:,i)-StdProfile(:,i)/2;flipud(MeanProfile(:,i)+StdProfile(:,i)/2)],...
        Colors(i,:),'linestyle','none','facealpha',.03,'HandleVisibility','off');
    plot(RadialBins,MeanProfile(:,i),'Color',Colors(i,:),...
        'LineWidth',1, 'DisplayName', LegendNames{i});
end
hold off
xlabel('Distance From Spot');
ylabel(sprintf('Round %.0f, Channel %.0f Intensity',r,b));
legend;
title('Z=0 Radial Intensity Profile');

subplot(2,1,2);
hold on
for i=1:nPlots
    fill([ZBins';flipud(ZBins')],...
        [MeanZProfile(:,i)-StdZProfile(:,i)/2;flipud(MeanZProfile(:,i)+StdZProfile(:,i)/2)],...
        Colors(i,:),'linestyle','none','facealpha',.05,'HandleVisibility','off');
    plot(ZBins,MeanZProfile(:,i),'Color',Colors(i,:),...
        'LineWidth',1, 'DisplayName', LegendNames{i});
end
hold off
xlabel('Z Plane');
ylabel(sprintf('Round %.0f, Channel %.0f Intensity',r,b));
title('Mid Pixel Z Intensity Profile');

%% Plot dependence on Z-coordinate Each Spot Found
figure;
for i=1:nPlots
    subplot(nPlots,1,i);
    boxplot_x = [];
    boxplot_g = [];
    j=1;
    for z=unique(SpotZPlane)'
        bxplt = log10(ZProfile(PlotLogical(:,i)&SpotZPlane==z,ImRadZ+1));
        if length(bxplt)==0
            bxplt=nan;
        end
        boxplot_x = [boxplot_x;bxplt];
        boxplot_g = [boxplot_g;j*ones(size(bxplt))];
        j=j+1;
    end
    boxplot(boxplot_x,boxplot_g, 'Colors',Colors(i,:), 'plotstyle', 'compact');
    title(LegendNames{i});
    xticklabels(unique(SpotZPlane));
    xlabel('Z Plane');
    ylabel('Log10(Intensity)');
end
sgtitle(sprintf('Round %.0f, Channel %.0f',r,b));


end

function [Image2D,IFS] = show_image_eigenvector(o,Image3D,PlotTitle)
%% This plots the eigenvector Image_3D as well as the 2D projected and 3D versions of it.
%Filter properties
h = -hanning(o.ExtractR2*2+1);
h = -h/sum(h);
h(o.ExtractR2+1-o.ExtractR1:o.ExtractR2+1+o.ExtractR1) = ...
    h(o.ExtractR2+1-o.ExtractR1:o.ExtractR2+1+o.ExtractR1)+hanning(o.ExtractR1*2+1)/sum(hanning(o.ExtractR1*2+1));
SE = ftrans2(h');
SE = single(SE);

ImRadZ = (size(Image3D,3)-1)/2;
ImRad = (size(Image3D,1)-1)/2;
Image2D = cell(ImRadZ*2+1,1);
figure;
CaxLim = [min(Image3D(:)),max(Image3D(:))];
for z=1:ImRadZ*2+1
    subplot(ceil((ImRadZ*2+3)/ImRadZ),ImRadZ,z);
    imagesc(Image3D(:,:,z));
    title(sprintf('Z = %.0f',z-(ImRadZ+1)));
    caxis(CaxLim);
    Image2D{z} = Image3D(:,:,z);
    xticks([1,ImRad+1,ImRad*2+1]);
    xticklabels([-ImRad,0,ImRad]);
    yticks([1,ImRad+1,ImRad*2+1]);
    yticklabels([-ImRad,0,ImRad]);
end
%Add 2D projection and filtered image
Image2D = o.fstack_modified(Image2D);
subplot(ceil((ImRadZ*2+3)/ImRadZ),ImRadZ,z+1);
imagesc(Image2D);
caxis(CaxLim);
colorbar;
xticks([1,ImRad+1,ImRad*2+1]);
xticklabels([-ImRad,0,ImRad]);
yticks([1,ImRad+1,ImRad*2+1]);
yticklabels([-ImRad,0,ImRad]);
title('2D Projection');
IFS = single(padarray(Image2D,(size(SE)-1)/2,'replicate','both'));
IFS = convn(IFS,SE,'valid');
%IFS = IFS*o.ExtractScale;
subplot(ceil((ImRadZ*2+3)/ImRadZ),ImRadZ,z+2);
imagesc(IFS);
colormap(gca,bluewhitered);
colorbar;
xticks([1,ImRad+1,ImRad*2+1]);
xticklabels([-ImRad,0,ImRad]);
yticks([1,ImRad+1,ImRad*2+1]);
yticklabels([-ImRad,0,ImRad]);
title('Filtered Image');
if nargin>=3
    sgtitle(PlotTitle);
end

end

