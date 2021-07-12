function o = call_spots_omp_initial(o)
%% o = o.call_spots_omp_initial(LookupTable)
% THIS ACCEPTS ALL GROUND TRUTH LOCAL MAXIMA EVEN DUPLICATES
% This is another probability method for gene calling.
% The difference here, is that an image is built up for each gene and then
% the spots are the local maxima on each gene image. This allows for
% multiple gene matches at each location i.e. overlapping spots. 
%
% o: iss object
% LookupTable: should be returned from call_spots_prob. It
% just gives the the probabilities that each spot score is explained by each
% gene. It saves calculating the probabilities explicitly each time.
%
% produces 
% pxSpotColors(Spot,b,r): intensity of Spot in channel b, round r
% pxSpotGlobalYX(Spot,:): global yx coordinate for each spot
% pxSpotCodeNo(Spot): gene index for each spot
% pxLogProbOverBackground(Spot): log of probability spot can be explained
% by gene relative to probability it can be explained by background.
% pxSpotScore(Spot): pxLogProbOverBackground of best gene match relative to
% second best gene match at that location.
% pxSpotScoreDev(Spot): standard deviation in spot scores across all genes
% at that location.
% pxSpotIntensity(Spot): intensity of the spot. Takes into account
% pxSpotCodeNo. Calculated by get_spot_intensity.
% pxSpotBestGene(Spot): gene index of gene with highest probability at this
% location.
% 
%% Logging
if o.LogToFile
    diary(o.LogFile);
    cleanup = onCleanup(@()diary('off'));
end

if o.ProbMethod==2 && o.ScoreScale~=1
    warning('o.ProbMethod=2 so changing o.ScoreScale to 1');
    o.ScoreScale = 1;
end

%% Load in images, make images for each gene, find local maxima in images
nCodes = length(o.CharCodes);
rr = o.ReferenceRound;      %Round to base coordinate system on
[nY, nX] = size(o.EmptyTiles);
nTiles = nY*nX;
NonemptyTiles = find(~o.EmptyTiles)';
if size(NonemptyTiles,2)==1
    NonemptyTiles = NonemptyTiles';
end

PeakSpotColors = cell(nCodes,1);
PeakLocalYX = cell(nCodes,1);
PeakResOverBackground = cell(nCodes,1);
Peak2ndBestRes = cell(nCodes,1);
PeakCoef = cell(nCodes,1);
OriginalTile = cell(nCodes,1);
PeakBestGene = cell(nCodes,1);
BackgroundSpotColors = cell(1,1);

%Get output file names so don't have data from all tiles in matlab at once
OutputTileNo = cumsum(equal_split(int32(length(NonemptyTiles)),...
    round(length(NonemptyTiles)/o.OMPFileMaxTiles)));
nFiles = length(OutputTileNo);
o.OMP_FileNames = cell(nFiles,1);  

for FileIdx=1:nFiles
    if isequal(NonemptyTiles,1:nTiles)
        o.OMP_FileNames(FileIdx) =...
            {fullfile(o.OutputDirectory,...
            strcat('OMP_Initial_Workspace',num2str(OutputTileNo(FileIdx)),'.mat'))};
    else
         o.OMP_FileNames(FileIdx) =...
            {fullfile(o.OutputDirectory,...
            strcat('Abridged_OMP_Initial_Workspace',num2str(OutputTileNo(FileIdx)),'.mat'))};
    end
end

FileIdx = 1;
for t=1:length(NonemptyTiles)  
    tile_no = NonemptyTiles(t);
    if exist(o.OMP_FileNames{FileIdx}, 'file')
        if ismember(t,OutputTileNo)
            FileIdx=FileIdx+1;
        end
        fprintf('Tile %d already done.\n', tile_no);
        continue;
    end
        
    %Get pixel colors
    [GoodAnchorLocalYX,GoodSpotColors] = o.get_spot_colors_all_pixels(tile_no);
    
    %Get local maxima log probabilities for each gene
    [tPeakLocalYX,tPeakSpotColors,tPeakResOverBackground,...
    tPeak2ndBestRes,tPeakCoef,tOriginalTile,tPeakBestGene,tBackgroundSpotColors] = ...
    o.detect_peak_genes_omp_initial(GoodSpotColors,GoodAnchorLocalYX,tile_no);
    clearvars GoodSpotColors GoodAnchorLocalYX;
    
    %Keep data for all tiles together
    PeakSpotColors = cellfun( @(x,y) [x;y], PeakSpotColors, tPeakSpotColors,...
        'UniformOutput', false );
    PeakLocalYX = cellfun( @(x,y) [x;y], PeakLocalYX, tPeakLocalYX, 'UniformOutput', false );
    clearvars tPeakSpotColors tPeakLocalYX;
    PeakResOverBackground = cellfun( @(x,y) [x;y], PeakResOverBackground, tPeakResOverBackground,...
        'UniformOutput', false );
    Peak2ndBestRes = cellfun( @(x,y) [x;y], Peak2ndBestRes, tPeak2ndBestRes,...
        'UniformOutput', false );
    PeakCoef = cellfun( @(x,y) [x;y], PeakCoef, tPeakCoef, 'UniformOutput', false );
    OriginalTile = cellfun( @(x,y) [x;y], OriginalTile, tOriginalTile, 'UniformOutput', false );
    PeakBestGene = cellfun( @(x,y) [x;y], PeakBestGene, tPeakBestGene, 'UniformOutput', false );
    BackgroundSpotColors{1} = [BackgroundSpotColors{1};tBackgroundSpotColors];
    clearvars tPeakLogProbOverBackground tPeak2ndBestRes tPeakCoef...
        tOriginalTile tPeakBestGene tBackgroundSpotColors;
    
    if ismember(t,OutputTileNo)
        save(o.OMP_FileNames{FileIdx}, 'PeakSpotColors','PeakLocalYX',...
            'PeakResOverBackground','Peak2ndBestRes','PeakCoef','OriginalTile',...
            'PeakBestGene','BackgroundSpotColors', '-v7.3');
        PeakSpotColors = cell(nCodes,1);
        PeakLocalYX = cell(nCodes,1);
        PeakResOverBackground = cell(nCodes,1);
        Peak2ndBestRes = cell(nCodes,1);
        PeakCoef = cell(nCodes,1);
        OriginalTile = cell(nCodes,1);
        PeakBestGene = cell(nCodes,1);
        BackgroundSpotColors = cell(1,1);
        FileIdx=FileIdx+1;
    end
end

%% Deal with each file one by one
SpotCodeNo = cell(1,1);
SpotColors = cell(1,1);
GlobalYX = cell(1,1);
ResOverBackground = cell(1,1);
SecondBestRes = cell(1,1);
Coef = cell(1,1);
SpotLocalTile = cell(1,1);
SpotBestGene = cell(1,1);
AllBackgroundSpotColors = cell(1,1);
    
nFiles = length(o.OMP_FileNames);
fprintf('\nGetting results from file      ');
for f = 1:nFiles
    if nFiles<10
        fprintf('\b\b\b%d/%d',f,nFiles);
    else
        if f<10; fprintf('\b\b\b\b%d/%d',f,nFiles);
        else; fprintf('\b\b\b\b\b%d/%d',f,nFiles); end
    end

    %Get global coordinates of peaks
    load(cell2mat(o.OMP_FileNames(f)));
    PeakGlobalYX = cell(nCodes,1);
    for GeneNo = 1:nCodes
        PeakGlobalYX{GeneNo} = bsxfun(@plus,double(PeakLocalYX{GeneNo}),...
            o.TileOrigin(OriginalTile{GeneNo},:,rr));
    end 
    % Clear memory by removing bad matches, pSpotScore<-5
    clearvars PeakLocalYX
    %Have to filter results initially so don't save too many and have
    %memory issues. Only take gene which are 1st or second best at thier location
    %unless probability relative to background is good.
    QualOK = cellfun(@(x1,x2,x3) x1-x2>=o.ompInitialScoreThresh & x3>o.ompInitialCoefThresh,...
        PeakResOverBackground,Peak2ndBestRes,PeakCoef,'UniformOutput',false);
    PeakSpotColors = cellfun(@(x1,x2) x1(x2,:,:),PeakSpotColors,QualOK,'UniformOutput',false);
    PeakGlobalYX = cellfun(@(x1,x2) x1(x2,:),PeakGlobalYX,QualOK,'UniformOutput',false);
    PeakResOverBackground = ...
        cellfun(@(x1,x2) x1(x2),PeakResOverBackground,QualOK,'UniformOutput',false);
    Peak2ndBestRes = cellfun(@(x1,x2) x1(x2),Peak2ndBestRes,QualOK,'UniformOutput',false);
    PeakCoef = cellfun(@(x1,x2) x1(x2),PeakCoef,QualOK,'UniformOutput',false);
    OriginalTile = cellfun(@(x1,x2) x1(x2),OriginalTile,QualOK,'UniformOutput',false);
    PeakBestGene = cellfun(@(x1,x2) x1(x2),PeakBestGene,QualOK,'UniformOutput',false);
    clearvars QualOK
        
    % Remove duplicates by keeping only spots detected on their home tile
    ndSpotColors = cell(nCodes,1);
    ndGlobalYX = cell(nCodes,1);
    ndResOverBackground = cell(nCodes,1);
    nd2ndBestRes = cell(nCodes,1);
    ndCoef = cell(nCodes,1);
    ndOriginalTile = cell(nCodes,1);
    ndBestGene = cell(nCodes,1);
    
    for GeneNo = 1:nCodes
        if o.Graphics==2
            figure(1001)
            plot(PeakGlobalYX{GeneNo}(:,2), PeakGlobalYX{GeneNo}(:,1), '.', 'markersize', 1);
            title('All global coords including duplicates');
            %set(gca, 'YDir', 'reverse');
        end
        
        [AllLocalTile, ~] = which_tile(PeakGlobalYX{GeneNo}, o.TileOrigin(:,:,rr), o.TileSz);
        NotDuplicate = (AllLocalTile==OriginalTile{GeneNo});
        
        ndSpotColors{GeneNo} = PeakSpotColors{GeneNo}(NotDuplicate,:,:);
        ndGlobalYX{GeneNo} = PeakGlobalYX{GeneNo}(NotDuplicate,:);
        ndResOverBackground{GeneNo} = PeakResOverBackground{GeneNo}(NotDuplicate);
        nd2ndBestRes{GeneNo} = Peak2ndBestRes{GeneNo}(NotDuplicate);
        ndCoef{GeneNo} = PeakCoef{GeneNo}(NotDuplicate);
        ndOriginalTile{GeneNo} = OriginalTile{GeneNo}(NotDuplicate); 
        ndBestGene{GeneNo} = PeakBestGene{GeneNo}(NotDuplicate);
        if o.Graphics==2
            figure(1002); clf
            plot(ndGlobalYX{GeneNo}(:,2), ndGlobalYX{GeneNo}(:,1), '.', 'markersize', 1);
            title('Global coords without duplicates');
            drawnow;
            %set(gca, 'YDir', 'reverse');
        end
    end
    
    %Free up memory
    clearvars PeakSpotColors PeakGlobalYX PeakResOverBackground Peak2ndBestRes PeakCoef...
        OriginalTile AllLocalTile NotDuplicate PeakBestGene;
    
    % Get final results
    nGeneSpots = cell2mat(cellfun(@length,ndCoef,'uni',false));    
    for GeneNo = 1:nCodes
        SpotCodeNo{1} = [SpotCodeNo{1};ones(nGeneSpots(GeneNo),1)*GeneNo];
        SpotColors{1} = [SpotColors{1};ndSpotColors{GeneNo}];
        GlobalYX{1} = [GlobalYX{1};ndGlobalYX{GeneNo}];
        ResOverBackground{1} = [ResOverBackground{1};ndResOverBackground{GeneNo}];
        SecondBestRes{1} = [SecondBestRes{1};nd2ndBestRes{GeneNo}];
        Coef{1} = [Coef{1};ndCoef{GeneNo}];
        SpotLocalTile{1} = [SpotLocalTile{1};ndOriginalTile{GeneNo}];
        SpotBestGene{1} = [SpotBestGene{1};ndBestGene{GeneNo}];
    end
    clearvars ndSpotColors ndGlobalYX ndResOverBackground nd2ndBestRes ndCoef
    
    AllBackgroundSpotColors{1} = [AllBackgroundSpotColors{1};BackgroundSpotColors{1}];
    clearvars BackgroundSpotColors
end
fprintf('\n');

%% Get background eigenvectors
% o.UseBackgroundEigenvectors = o.BackgroundEigenvalues>=o.BackgroundEigenvalues(...
%     min(o.Max_nBackground,length(o.BackgroundEigenvalues))) & ...
%     o.BackgroundMaxGeneDotProduct<o.BackgroundMaxGeneDotProductThresh;
[o.BackgroundEigenvectors,o.BackgroundEigenvalues,o.BackgroundMaxGeneDotProduct,...
    o.BackgroundMaxGeneDotProductGene] = ...
    o.get_background_codes(cell2mat(AllBackgroundSpotColors));
if o.ompBackgroundChannelStrips
    o.UseBackgroundEigenvectors = true(sum(o.UseChannels>0),1);
else
    o.UseBackgroundEigenvectors = o.BackgroundEigenvalues>o.ompBackgroundEigvalueThresh;
end
o.nBackground = sum(o.UseBackgroundEigenvectors);

%% Add results to iss object

o.iompSpotColors = cell2mat(SpotColors);
o.iompSpotCodeNo = int16(cell2mat(SpotCodeNo));
o.iompSpotGlobalYX = cell2mat(GlobalYX);
o.iompResOverBackground = cell2mat(ResOverBackground);
o.iompSpotScore = o.iompResOverBackground-cell2mat(SecondBestRes);
o.iompCoef = cell2mat(Coef);        
[o.iompSpotIntensity,o.iompSpotIntensity2] = ...
    o.get_spot_intensity(o.iompSpotCodeNo,o.iompSpotColors,o.z_scoreSCALE);
o.iompLocalTile = int16(cell2mat(SpotLocalTile));
o.iompSpotBestGene = int16(cell2mat(SpotBestGene));
end


