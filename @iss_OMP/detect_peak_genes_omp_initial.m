function [PeakLocalYX,PeakSpotColors,PeakResOverBackground,...
    Peak2ndBestRes,PeakCoef,OriginalTile,PeakBestGene,BackgroundSpotColors] = ...
    detect_peak_genes_omp_initial(o,GoodSpotColors,GoodLocalYX,t)
%% [PeakLocalYX,PeakSpotColors,PeakResOverBackground,...
%    Peak2ndBestRes,PeakCoef,OriginalTile,PeakBestGene,BackgroundSpotColors] = ...
%    detect_peak_genes_omp_initial(o,LookupTable,GoodSpotColors,GoodLocalYX,t)
%
% This finds the local maxima in residual reduction compared to background
% for each gene.
% 
% Input
% o: iss object
% GoodSpotColors(S,b,r) is the intensity for spot S in channel b, round r.
% S should cover all pixel values that don't go off edge of tile in any b,r.
% GoodLocalYX(S,:) is the corresponding pixel location.
% t is the current tile of interest
%
% Output
% PeakLocalYX{G} contains the YX position of local maxima of gene G.
% PeakSpotColors{G} contains the corresponding spot colors.
% PeakResOverBackground{G} contains the corresponding 
% Residual Reduction relative to the background.
% Peak2ndBestRes{G} contains the Residual Reduction relative to the
% background for the second best match at that location.
% PeakCoef{G} is the omp coefficient of gene G at that location.
% OriginalTile{G} = t
% PeakBestGene{G} contains the best gene at the location of local maxima of
% gene G. I.e. most will be G but few will be overlapping. 
% BackgroundSpotColors contains the spot colors of all pixels far from any
% gene score maxima spots.

%% Use Channel strips as background for initial search
nCodes = length(o.CharCodes);
nChannels = length(o.UseChannels);
nRounds = length(o.UseRounds);
z_scoredSpotColors = (double(GoodSpotColors)-o.z_scoreSHIFT)./o.z_scoreSCALE;
z_scoredSpotColors = z_scoredSpotColors(:,o.UseChannels,o.UseRounds);
BledGeneCodes = reshape(o.iompBledCodes,[nCodes,o.nBP,o.nRounds]);
BledGeneCodes = BledGeneCodes(:,o.UseChannels,o.UseRounds);
BledGeneCodes(isnan(BledGeneCodes)) = 0;

if length(o.UseRounds)>2
    %Can only distinguish genes from background if have more than 2 rounds
    BackgroundStripCodes = zeros(o.nBP,o.nBP,o.nRounds);
    for b=1:o.nBP
        BackgroundStripCodes(b,b,:) = 1;
    end
    nBackground = length(o.UseChannels);
    BackgroundStripCodes = BackgroundStripCodes(o.UseChannels,o.UseChannels,o.UseRounds);
    BledCodes = zeros(nCodes+nBackground,length(o.UseRounds)*nBackground);
    BledCodes(1:nCodes,:) = BledGeneCodes(:,:);
    BledCodes(nCodes+1:nCodes+nBackground,:) = BackgroundStripCodes(:,:);
else
    BledCodes = BledGeneCodes(:,:);
    nBackground = 0;
end
BledCodes = BledCodes./vecnorm(BledCodes(:,:),2,2);
BledCodes(isnan(BledCodes)) = 0;
%% PROPER METHOD: Get residual reduction obtained by removing each gene in addition to 
%that achieved by removing background
nSpots = size(z_scoredSpotColors,1);
z_scoredSpotColors = z_scoredSpotColors(:,:)';
%AllSpotIntensity(i) is the 4th largest value of AllSpotColors(i,:) if
%o.nRounds=7, i.e. this is saying we need high intensity in at least
%one channel in half the number of rounds to be a gene.
AllSpotIntensity = prctile(z_scoredSpotColors,...
    (nRounds*nChannels-nRounds/2)*100/(nRounds*nChannels))';
Use = AllSpotIntensity>o.iompInitialIntensityThresh;
z_scoredSpotColors = z_scoredSpotColors(:,Use);
BackgroundResNorm = zeros(nSpots,1);
BackgroundResNorm(Use) = o.get_spot_residual(BledCodes',z_scoredSpotColors,...
    nCodes+1:nCodes+nBackground);
AllResOverBackground = zeros(nSpots,nCodes);
AllCoef = zeros(nSpots,nCodes);
fprintf('Finding residual after removing gene     ');
for GeneNo = 1:nCodes
    g_num = sprintf('%.6f', GeneNo);
    fprintf('\b\b\b\b%s',g_num(1:4));
    [gResNorm,gCoef] = o.get_spot_residual(BledCodes',z_scoredSpotColors,...
        [nCodes+1:nCodes+nBackground,GeneNo]);
    AllResOverBackground(Use,GeneNo) = BackgroundResNorm(Use) - gResNorm;
    AllCoef(Use,GeneNo) = gCoef(:,end);   
end
fprintf('\n');
clear BackgroundResNorm gCoef gResNorm z_scoredSpotColors

%% CHEAT/QUICKER METHOD: Update coefficient independently.
% Really all background eigenvectors should be orthogonal for this to work.

% DotProduct = z_scoredSpotColors(:,:)*BledCodes(nCodes+1:nCodes+o.nBP,:)';
% ToRemove = DotProduct*BledCodes(nCodes+1:nCodes+o.nBP,:);
% bSpotColors = z_scoredSpotColors(:,:) - ToRemove;
% clearvars z_scoredSpotColors ToRemove
% bNorm = vecnorm(bSpotColors,2,2);
% DotProduct = bSpotColors*BledCodes(1:nCodes,:)';
% fprintf('Finding residual after removing gene     ');
% for GeneNo = 1:nCodes
%     g_num = sprintf('%.6f', GeneNo);
%     fprintf('\b\b\b\b%s',g_num(1:4));
%     gSpotColors = bSpotColors - DotProduct(:,GeneNo)*BledCodes(GeneNo,:);
%     AllResOverBackground(:,GeneNo) = bNorm - vecnorm(gSpotColors,2,2);
%     AllCoef(:,GeneNo) = DotProduct(:,GeneNo);   
% end
% fprintf('\n');
% clearvars DotProduct bSpotColors bNorm
%% For each gene, find peaks in residual images. Keep these as spots going forward
PeakSpotColors = cell(nCodes,1);
PeakLocalYX = cell(nCodes,1);
PeakResOverBackground = cell(nCodes,1);
Peak2ndBestRes = cell(nCodes,1);
PeakCoef = cell(nCodes,1);
OriginalTile = cell(nCodes,1);
PeakBestGene = cell(nCodes,1);

GeneIm = zeros(max(GoodLocalYX));     %Y index is first in zeros
Ind = sub2ind(size(GeneIm),GoodLocalYX(:,1),GoodLocalYX(:,2));

fprintf('Finding peaks for gene     ');
for GeneNo = 1:nCodes    
    g_num = sprintf('%.6f', GeneNo);
    fprintf('\b\b\b\b%s',g_num(1:4));
    
    %Find local maxima in gene image
    GeneIm(Ind) = AllResOverBackground(:,GeneNo); 
    Small = 1e-6;
    se1 = strel('disk', o.PixelDetectRadius);     %Needs to be bigger than in detect_spots
    Dilate = imdilate(GeneIm, se1);
    MaxPixels = find(GeneIm + Small >= Dilate & GeneIm>o.ompInitialResThresh);
    
    %Get Indices of Good Global Spot Colors / YX
    PeakInd = find(ismember(Ind,MaxPixels));        %As position in Ind = ResOverBackground Index = Good Index
    nPeaks = length(PeakInd);
    %Save information for that gene
    PeakSpotColors{GeneNo} = GoodSpotColors(PeakInd,:,:);
    PeakLocalYX{GeneNo} = GoodLocalYX(PeakInd,:);
    peakPoverB = AllResOverBackground(PeakInd,:);
    PeakResOverBackground{GeneNo} = peakPoverB(:,GeneNo); 
    PeakCoef{GeneNo} = AllCoef(PeakInd,GeneNo);
    
    %Find 2nd best gene so can give score relative to it
    [~,PeakBestGene{GeneNo}] = max(peakPoverB,[],2);
    peakPoverB(sub2ind(size(peakPoverB),(1:nPeaks)',PeakBestGene{GeneNo}))=-inf;
    Peak2ndBestRes{GeneNo} = max(peakPoverB,[],2);
    %SortProb = sort(peakPoverB,2,'descend');
    %Peak2ndBestLogProb{GeneNo} = SortProb(PeakInd,2);        
    OriginalTile{GeneNo} = int16(ones(nPeaks,1)*t);
    PeakBestGene{GeneNo} = int16(PeakBestGene{GeneNo});
    
    if o.Graphics==2
        figure(50965467); clf;
        imagesc(GeneIm); hold on;
        colormap(gca,bluewhitered);
        plot(PeakLocalYX{GeneNo}(:,2), PeakLocalYX{GeneNo}(:,1), 'gx');
        hold off
        drawnow;
    end   
end
fprintf('\n');
clearvars AllResOverBackground AllCoef peakPoverB PeakInd MaxPixels Dilate
%Find background pixels as those far from peaks for all genes
tree = KDTreeSearcher(double(cell2mat(PeakLocalYX)));
[~,Dist] = tree.knnsearch(double(GoodLocalYX),'K',1);
BackgroundPixels = Dist>o.ompBackgroundDist;
BackgroundSpotColors = GoodSpotColors(BackgroundPixels,:,:);
if o.Graphics==2
    figure(50965468); clf;
    imagesc(GeneIm); hold on;
    colormap(gca,bluewhitered);
    plot(GoodLocalYX(BackgroundPixels,2), GoodLocalYX(BackgroundPixels,1), 'gx');
    hold off
end
end

