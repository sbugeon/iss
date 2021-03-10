function [o,LookupTable] = call_spots_prob(o,LookupTableInput)
%% [o, LookupTable] = o.call_spots_prob
% calls spots to codes for in-situ sequencing. Run this after find_spots
%
% o: iss object
% LookupTableInput: Can skip making LookupTable if provide one e.g. if
% running but with different o.ScoreScale
% LookupTable: gives the the probabilities that each spot score is explained 
% by each gene. It saves calculating the probabilities explicitly each time.
%
% produces 
% pSpotCodeNo(Spot): gene index for each spot
% pLogProbOverBackground(Spot): log of probability spot can be explained
% by gene relative to probability it can be explained by background.
% pSpotScore(Spot): pLogProbOverBackground of best gene match relative to
% second best gene match at that location.
% pSpotScoreDev(Spot): standard deviation in spot scores across all genes
% at that location.
% pSpotIntensity(Spot): intensity of the spot. Takes into account
% pSpotCodeNo. Calculated by get_spot_intensity.

%% Logging
if o.LogToFile
    diary(o.LogFile);
    cleanup = onCleanup(@()diary('off'));
end

if o.ProbMethod==2 && o.ScoreScale~=1
    warning('o.ProbMethod=2 so changing o.ScoreScale to 1');
    o.ScoreScale = 1;
end

%% Assign prob variables that are equal to dot product method
o.pSpotGlobalYX = o.dpSpotGlobalYX;
o.pSpotColors = o.dpSpotColors;
o.pSpotIsolated = o.dpSpotIsolated;
o.pLocalTile = o.dpLocalTile;

%% Make Bleed Matrices
%Only using channels and rounds given by o.UseChannels and o.UseRounds
if isempty(o.UseChannels)
    o.UseChannels = 1:o.nBP;
end
    
if isempty(o.UseRounds)
    o.UseRounds = 1:o.nRounds;
end

AllChannels = 1:o.nBP;
IgnoreChannels = setdiff(AllChannels,o.UseChannels);
nChans = size(o.UseChannels,2);
nRounds = size(o.UseRounds,2);

%Normalise each colour channel by a percentile as to correct for weaker
%colour channels
p = o.BledCodesPercentile;
NormBleedMatrix = o.BleedMatrix;
if o.Graphics
    figure(98043715); clf
    for i=1:nRounds
        subplot(ceil(nRounds/3),3,i);
        imagesc(NormBleedMatrix(:,:,i));
        caxis([0 1]);
        title(sprintf('Round %d', o.UseRounds(i)));
        set(gca, 'xtick', 1:o.nBP);
        set(gca, 'XTickLabel', o.bpLabels(AllChannels));
        set(gca, 'ytick', 1:o.nBP);
        set(gca, 'yTickLabel', o.bpLabels(AllChannels));
        if i==4
            xlabel('Actual')
            ylabel('Measured');
        end
    end
    drawnow;
end

%Unnormalise Bleed matrices by multiplying rows by percentiles
BleedMatrix = zeros(o.nBP,o.nBP,nRounds);
for r=1:o.nRounds
    for b=1:o.nBP
        BleedMatrix(b,:,r) = p(:,b,r)*NormBleedMatrix(b,:,r);
    end
end

if o.Graphics
    figure(98043764); clf
    for i=1:nRounds
        subplot(ceil(nRounds/3),3,i); 
        imagesc(BleedMatrix(:,:,i)); 
        %caxis([0 1]); 
        title(sprintf('Cycle %d', o.UseRounds(i))); 
        set(gca, 'xtick', 1:o.nBP);
        set(gca, 'XTickLabel', o.bpLabels(AllChannels));
        set(gca, 'ytick', 1:o.nBP);
        set(gca, 'yTickLabel', o.bpLabels(AllChannels));
        if i==4
            xlabel('Actual')
            ylabel('Measured');
        end
    end
    drawnow;
%     subplot(2,3,6);
%     caxis([0 1]); 
%     axis off
%     colormap hot
% %     colorbar
end

%Set bad colour channels to 0 for all rounds
BleedMatrix(IgnoreChannels,:,:) = 0;
%o.BleedMatrix = NormBleedMatrix;
o.pBleedMatrix = BleedMatrix;
%BleedMatrix = o.pBleedMatrix;
%% Get BledCodes for each gene
% now load in the code book and apply bleeds to it
%codebook_raw = importdata(o.CodeFile);
%CharCode = codebook_raw.textdata(2:end,5);
%GeneName = codebook_raw.textdata(2:end,3);
GeneName = {};
CharCode = {};
fp = fopen(o.CodeFile, 'r');
tmp = textscan(fp, '%s %s', inf);
GeneName=tmp{1};
CharCode=tmp{2};
fclose(fp);

% bit of a hack to get rid of Sst and Npy (assume always in the end)
nCodes = size(CharCode,1) - nnz(cellfun(@(v) strcmp(v(1:2),'SW'), CharCode));

% put them into object o but without the extras
o.CharCodes=CharCode(1:nCodes);
o.GeneNames=GeneName(1:nCodes);

% create numerical code (e.g. 33244 for CCGAA)
NumericalCode = zeros(nCodes, o.nRounds);
for r=1:o.nRounds
    if r<=o.nRounds-o.nRedundantRounds
        for c=1:nCodes
            try
                [~, NumericalCode(c,r)] = ismember(CharCode{c}(r), o.bpLabels);
            catch
                warning('Code %s has no channel for round %.0f',GeneName{c},r);
                NumericalCode(c,r)=0;
            end
        end
    else
        % redundant round - compute codes automatically
        % find pseudobases for this code
        for c=1:nCodes
            PseudoCode = repmat('0',1,o.nRounds-o.nRedundantRounds);
            for p = 1:length(o.RedundantPseudobases)
                PseudoCode(1,ismember(CharCode{c}, o.RedundantPseudobases{p}))=('0'+p);
            end
            % now match them to the redundant codes
            for cc=1:o.nBP
                rrn = r-o.nRounds+o.nRedundantRounds;
                if ~isempty(regexp(PseudoCode, o.RedundantCodes{rrn,cc}, 'once'))
                    NumericalCode(c,r)=cc;
                end
            end
        end
    end
end

BledCodes = zeros(nCodes, o.nBP*o.nRounds);
UnbledCodes = zeros(nCodes, o.nBP*o.nRounds);
% make starting point using bleed vectors (means for each base on each day)
for i=1:nCodes
    for r=1:nRounds
        if any(AllChannels == NumericalCode(i,o.UseRounds(r))) == 0 continue; end
        BledCodes(i,AllChannels+o.nBP*(r-1)) = BleedMatrix(:, find(AllChannels == NumericalCode(i,o.UseRounds(r))), r);
        UnbledCodes(i,AllChannels(find(AllChannels == NumericalCode(i,o.UseRounds(r))))+o.nBP*(r-1)) = 1;
    end
end

o.pBledCodes = BledCodes;

%% Find probability of each spot to each gene
%Comes from P(f) = P_lambda(lambda)P_hist(f-lambda*g) as f = lambda*g+background

%Load histogram data - background prob distribution
ModHistCounts = o.SmoothHistCounts;     %Smooth extreme values of HistCounts
nBins = length(o.HistValues);
HistProbs = ModHistCounts./sum(ModHistCounts);
o.HistProbs = (HistProbs+o.alpha)./(1+nBins*o.alpha);
o.HistZeroIndex = find(o.HistValues==0);

%Get Lambda probability distribution for all genes
%Find length of array by ensuring that gamma distribution for max predicted gene code
%drops to 0 within it.
x_All = 1:10^6;
gMax = max(o.pBledCodes(:));
MaxProbDist = gampdf(x_All,o.GammaShape,gMax/(o.GammaShape-1));
Max_x = min(find(x_All>gMax & MaxProbDist<1e-7));
if o.ProbMethod == 2
    x = min(o.pSpotColors(:))-1:Max_x;    %subsitiution x=lambda*g, -1 due to matlab indexing
else
    x = 0:Max_x;
end
o.ZeroIndex = find(x==0);     %need when looking up conv values
xSz = length(x);
LambdaDistAll = zeros(xSz,nCodes,o.nBP,o.nRounds);


fprintf('\nGetting probability distributions for gene   ');
for GeneNo = 1:nCodes
    if GeneNo<10
        fprintf('\b%d',GeneNo);
    else
        fprintf('\b\b%d',GeneNo);
    end
    BledCode = reshape(o.pBledCodes(GeneNo,:),[o.nBP,o.nRounds]);
    numCharCode = str2double(regexp(cell2mat(o.CharCodes(GeneNo)),'\d','match'))+1;
    
    for b=1:o.nBP
        for r=1:o.nRounds
            g = BledCode(b,r);
            if g<1
                g=1;
            end
            if o.ProbMethod == 1
                %GammaScale is g/(o.GammaShape-1) so distribution peaks at g.
                LambdaDist = gampdf(x,o.GammaShape,g/(o.GammaShape-1));

            elseif o.ProbMethod == 2
                if numCharCode(r)==b
                    %for b/r in CharCodes, expect non zero lambda.
                    %g always >0 in this case
                    LambdaDist = raylpdf(x/g,o.RaylConst);
                else
                    %for b/r not in CharCodes, expect approx zero lambda.
                    LambdaDist = (o.ExpConst/2)*exp(-o.ExpConst*abs(x/g));
                end
            end
            %Normalise as for small g might not be normalised due to
            %discrete distribution.
            LambdaDistAll(:,GeneNo,b,r) = LambdaDist/sum(LambdaDist);
            %Ensure no zero values as log(0)=-inf.
            LambdaDistAll(:,GeneNo,b,r) = (LambdaDistAll(:,GeneNo,b,r)+o.alpha)./(1+xSz*o.alpha);
        end
    end
end
o.LambdaDist = LambdaDistAll;

if o.ProbMethod == 1
    %Get background Prob using gamma with g=1 i.e. very small prediction of
    %intensity.
    gBackground = 1.0;
    BackgroundGamma = gampdf(x,o.GammaShape,gBackground/(o.GammaShape-1));
    BackgroundGamma = BackgroundGamma/sum(BackgroundGamma);
elseif o.ProbMethod == 2
    %This just shifts o.HistProbs when convolved with it to where x=0.
    BackgroundGamma = dirac(x);
    BackgroundGamma(BackgroundGamma==inf)=1;
end
%Ensure no zero values as log(0)=-inf.
BackgroundGamma = (BackgroundGamma+o.alpha)./(1+xSz*o.alpha);
o.BackgroundLambdaDist = BackgroundGamma;
o.BackgroundProb = zeros(xSz+length(o.HistCounts)-1,o.nBP,o.nRounds);
for b=1:o.nBP
    for r=1:o.nRounds
        o.BackgroundProb(:,b,r) = conv(BackgroundGamma,o.HistProbs(:,b,r));
    end
end

%Store convolution results as look up table
BledCodes = reshape(o.BledCodes,[nCodes,o.nBP,o.nRounds]);
if nargin<2 || isempty(LookupTableInput)    
    LookupTable = nan(xSz+length(o.HistCounts)-1,nCodes,o.nBP,o.nRounds);
    fprintf('\nDoing convolutions for Channel  ');
    for b=1:o.nBP
        fprintf('\b%d',b);
        if ismember(b,o.UseChannels)
            for r=1:o.nRounds
                NonZeroGenes = BledCodes(:,b,r)>0;
                LookupTable(:,NonZeroGenes,b,r)=log(conv2(LambdaDistAll(:,NonZeroGenes,b,r),o.HistProbs(:,b,r)));
                %If bled code is 0, just set to background distribution.
                LookupTable(:,~NonZeroGenes,b,r)=log(repmat(o.BackgroundProb(:,b,r),[1,sum(~NonZeroGenes)]));
            end
        end
    end
else
    LookupTable = LookupTableInput;
end
fprintf('\n');
%Only keep LookupTable values within spot range.
LookupTable = LookupTable(1:o.HistZeroIndex+o.ZeroIndex-1+max(o.HistValues),:,:,:);

%Get log probs for each spot 
LogProbOverBackground = o.get_LogProbOverBackground(o.pSpotColors,LookupTable);
[LogProbOverBackground,SpotCodeNo] = sort(LogProbOverBackground,2,'descend');

o.pLogProbOverBackground = LogProbOverBackground(:,1);
o.pSpotCodeNo = SpotCodeNo(:,1);
o.pSpotScore = LogProbOverBackground(:,1)-LogProbOverBackground(:,2);
%Store deviation in spot scores - can rule out matches based on a low
%deviation.
o.pSpotScoreDev = std(LogProbOverBackground,[],2);
[o.pSpotIntensity,o.pSpotIntensity2] = o.get_spot_intensity(o.pSpotCodeNo,o.pSpotColors);
if nargin<2 || isempty(LookupTableInput)
    save(fullfile(o.OutputDirectory, ['LookupTable',num2str(o.ProbMethod),'.mat']),'LookupTable','-v7.3');
end

end
