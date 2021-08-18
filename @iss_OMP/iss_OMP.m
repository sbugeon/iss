%% iss: code for processing of in situ sequencing
% Kenneth D. Harris and Xiaoyan Qian
% GPL 3.0 https://www.gnu.org/licenses/gpl-3.0.en.html
%
% to use:
% o = iss_OMP; % create structure, default parameters
% % change any parameters you want, and set o.FileN
% o = o.extract_and_filter; % create top-hat filtered tiffs for each tile
% o = o.find_spots; % find spot positions in global coordinates
% o = o.call_spots; % allocate each spot to a gene
% o = o.call_cells; % identify cells
%
% this current file iss.m just contains default parameter values
classdef iss_OMP < iss_GroundTruth
    
    properties
        %% OMP method variables
        %Spot colors are z scored by subtracting z_scoreSHIFT and dividing by
        %z_scoreSCALE
        z_scoreSHIFT;
        z_scoreSCALE;
        
        % ompBleedMatrixEigMethod == 'Mean': ScaledKMeans used to compute each
        % column of bleed matrix.
        % ompBleedMatrixEigMethod == 'Median': ScaledKMedians used to compute
        % each column of bleed matrix - just take median of all spots
        % assigned to each cluster.
        ompBleedMatrixEigMethod = 'Mean';
        
        % BleedMatrix used to estimate BledCodes in call_spots_omp.
        % Each round and channel is z-scored.
        z_scoreBleedMatrix
        
        %BackgroundEigenvectors are found from spots on tiles in
        %BackgroundEigenvectorTiles;
        BackgroundEigenvectorTiles;
                
        %BackgroundEigenvectors(i,b,r) is the intensity in channel b, round
        %r for the ith eigenvector of the covariance matrix comprised of
        %pixels a distance > ompBackgroundDist from any spot in iomp search
        %BackgroundEigenvalue(i) is the corresponding eigenvalue.
        ompBackgroundDist = 15;
        BackgroundEigenvectors;
        BackgroundEigenvalues;
        
        %BackgroundMaxGeneDotProduct(i) is the absolute dot product of
        %BackgroundEigenvectors(i,:) with
        %z_scoreBledCodes(BackgroundMaxGeneDotProductGene(i),:).
        %This is the largest dot product of any gene with BackgroundEigenvectors(i,:)
        %Idea is that large dot product means it is similar to a gene code
        %so reject any larger than BackgroundMaxGeneDotProductThresh.
        BackgroundMaxGeneDotProduct;
        BackgroundMaxGeneDotProductGene;
        BackgroundMaxGeneDotProductThresh = 0.5;
        
        %BackgroundEigenvectors(UseBackgroundEigenvectors,:) are appended
        %to z_scoreBledCodes to form ompBledCodes. Default is to choose
        %all with eigenvalue > ompBackgroundEigvalueThresh.
        ompBackgroundEigvalueThresh = 0.15;
        UseBackgroundEigenvectors;
        nBackground;    %size(UseBackgroundEigenvectors).
        Max_nBackground = 7; %there will be at most this number of background codes.
        ompBledCodes;
        
        %ompBackgroundChannelStrips = true means there will be 7 background
        %vectors, each being a single color channel in every round i.e.
        %CharCode = '0000000' for the first one.
        ompBackgroundChannelStrips = true;
        
        %OMP stops when reduction in residual drops below
        %ResidualThresh = ResidualThreshParam*SpotSecondLargestIntensity.
        ResidualThreshParam = 0.0612;
        %ResidualThresh is clamped between the two values below.
        ResidualThreshMin = 0.0100;
        ResidualThreshMax = 3.0;
        
        %ompMaxGenes is the maximum number of genes that can be assigned to
        %each pixel.
        ompMaxGenes = 6;
        %ompMaxGenes = 30; %THIS IS BEST FOR CONSTBACKGROUND_WEIGHTED VERSION
        
        % Only save to oObject spots with more than ompInitialNeighbThresh
        % non zero pixels near it.
        ompInitialNeighbThresh = 5;
        
        % To save time in call_spots_omp:
        % If BasePixels are pixels with 4th (ceil(o.nRounds/2))
        % largest z_scoredSpotColor greater than ompInitialIntensityThresh
        % Then OMP only carried out on pixels closer than
        % o.PixelDetectRadius to these.
        % Will probably have ompInitialIntensityThresh = ompIntensityThresh2
        ompInitialIntensityThresh = 0.0150;
        
        %For ompSpotScore:
        %GeneEfficiency(g,r)<ompScore_GeneEfficiencyThresh means don't use
        %round r for gene g when getting ompSpotScore.
        ompScore_GeneEfficiencyThresh = 0.1;
        %Score Preferentially weights rounds/channels with large error:
        %Rounds/Channels where error larger than ompScore_LargeErrorPrcntileThresh
        %of all rounds/channels for that spot are treated as large.
        ompScore_LargeErrorPrcntileThresh = 80;
        %Rounds/Channels where error is larger than
        %ompScore_LargeErrorMax*ompScore_LargeErrorPrcntileThreshValue is set to 
        %ompScore_LargeErrorMax*ompScore_LargeErrorPrcntileThreshValue
        ompScore_LargeErrorMax = 3;
        
        %For quality_threshold:
        %Spots must have  ompSpotIntensity2>ompIntensityThresh2 &
        %ompNeighbNonZeros > ompNeighbThresh2 &...
        %(ompSpotIntensity2>ompIntensityThresh | ompNeighbNonZeros >
        %ompNeighbThresh | ompSpotScore | ompScoreThresh);
        %ompNeighbNonZeros>ompNeighbThresh or ompSpotScore>ompScoreThresh.
        ompIntensityThresh = 0.5;       
        ompIntensityThresh2 = 0.005;
        ompNeighbThresh = 18;       
        ompNeighbThresh2 = 10;
        ompScoreThresh = 4.3;    
        ompScoreThresh2 = 1.1;
        
        %Spots assigned to gene that is not largest coefficient for pixel
        %have a stronger thresholding as given by these:
        ompIntensityThresh3 = 0.01;
        ompIntensityThresh3_CoefDiffFactor = 0.27;
        ompNeighbThresh3 = 28;
        ompScoreThresh3 = 6.9;
        
        %% OMP method outputs
        
        % ompSpotGlobalYX(Spot,1:2) contains y,x coordinates of every spot in
        % global coordiate system. Both combinatorial and extra spots
        % pixel spots are at different locations.
        ompSpotGlobalYX;
        
        % ompSpotColors(Spot, Base, Round) contains spot color on each base
        % and round. only for combinatorial splots
        ompSpotColors;
        
        %ompSpotIntensity is the mean intensity in the unbled code of gene
        %ompSpotCodeNo(s) in the z_scored ompSpotColors(s,:).
        ompSpotIntensity;
        
        %ompSpotIntensity2 is the median intensity in the unbled code of gene
        %ompSpotCodeNo(s) in the z_scored ompSpotColors(s,:).
        ompSpotIntensity2;
                      
        %ompScore is the based onreduction in error caused by introduction 
        %of gene in rounds/channels where there is not already a gene. 
        %Given by get_omp_score
        ompSpotScore;
                
        %ompSpotCodeNo is the gene found for each spot
        ompSpotCodeNo;
        
        %ompLocalTile(s) is the tile spot s was found on
        ompLocalTile;
        
        %ompCoef(s,g) is the weighting of spot s for gene g. Most are zero.
        %The last nBackground are background codes and are always non zero.
        ompCoefs;
        
        %ompNeighbNonZeros(s) is the number of pixels in the surrounding
        %region around spot s that have non zero coefficient for the gene
        %specified by ompSpotCodeNo(s). Min value is 1. Region size
        %specified by o.PixelDetectRadius. For default,
        %o.PixelDetectRadius=4, max value of this is 37.
        ompNeighbNonZeros;
        %% Initial OMP to get gene efficiencies
        
        %OMPFileNames contains the names of files in which initial omp
        %method data is stored
        OMP_FileNames;
        
        %OMPFileMaxTiles is approximately the maximum number of tiles
        %that can be stored in a single file. Output data to files so don't
        %get memory problems.
        OMPFileMaxTiles = 6;
        
        % To save time in call_spots_omp_initial:
        % iOMP only carried out on pixels with 4th (ceil(o.nRounds/2))
        % largest z_scoredSpotColor greater than iompInitialIntensityThresh
        % Will probably have iompInitialIntensityThresh =
        % GeneEfficiencyIntensityThresh.
        iompInitialIntensityThresh = 0.1;
        
        % iompSpotColors(Spot, Base, Round) contains spot color on each base
        % and round.
        iompSpotColors;
        
        % iompSpotCodeNo is the gene found for each spot
        iompSpotCodeNo;
        
        % iompSpotGlobalYX(Spot,1:2) contains y,x coordinates of every spot in
        % global coordiate system.
        iompSpotGlobalYX;
        
        % iompResOverBackground(s) is the norm of iompSpotColors(s,:) after
        % Background vectors have been fit to it minus the norm after
        % Background vectors and gene iompSpotCodeNo(s) have. 
        iompResOverBackground;
        
        % iompSpotScore is the difference between iompResOverBackground(s)
        % when gene iompSpotCodeNo(s) is used and when the 2nd best gene is
        % used. 
        iompSpotScore;
        
        % iompCoef(s) is the quantity of BledCodes(iompSpotCodeNo(s),:) fit
        % to iompSpotColors(s,:)
        iompCoef;
        
        %iompSpotIntensity is the mean intensity in the unbled code of gene
        %iompSpotCodeNo(s) in the z_scored iompSpotColors(s,:).
        iompSpotIntensity;
        
        %iompSpotIntensity2 is the median intensity in the unbled code of gene
        %iompSpotCodeNo(s) in the z_scored iompSpotColors(s,:).
        iompSpotIntensity2;
        
        %iompLocalTile(s) is the tile spot s was found on
        iompLocalTile;
        
        %iompBestGene(s) is gene number of best gene at location of spot s.
        iompSpotBestGene;
        
        % Only save to oObject spots with more than iompSpotScore,
        % iompResOverBackground, iompCoef exceeding these thresholds.
        ompInitialScoreThresh = 0;
        ompInitialResThresh = 0.005;
        ompInitialCoefThresh = 0;       
        
        %GroundTruth Data for initial OMP
        iomp_gtColor;
        iomp_gtIdentity;
        iomp_gtFound;
        
        %iompBledCodes are the bled codes used in the initial omp search.
        iompBledCodes;
        
        %The mean of iompSpots with Score,Coef,Intensity and ResOverBackground
        %all exceeding these thresholds will be used to find GeneEfficiency
        GeneEfficiencyScoreThresh = 0.01;
        GeneEfficiencyCoefThresh = 0.5;
        GeneEfficiencyIntensityThresh = 0.1;
        GeneEfficiencyResThresh = 0.1;
        
        %If less than GeneEfficiencyMinSpots satisfy the above thresholds
        %for gene g, GeneEfficiency(g,:) = 1. 
        GeneEfficiencyMinSpots = 10;
        
        %GeneEfficiency(g,r) is the strength of gene g in round r compared
        %to o.z_scoreBleedMatrix(:,CharCode{g}(r)+1) or
        GeneEfficiency;
        
        % z_scoreBledCodes(nCodes, nBP*nRounds): code vectors after modeling
        % crosstalk and z-scoring each round/channel. These bled codes
        % include the effect GeneEfficiency(g,r) on gene g, round r.
        % These are the gene bled codes used in the final omp search.
        z_scoreBledCodes;
        
    end
end

