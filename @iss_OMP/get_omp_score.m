function NormExpScore = get_omp_score(o,z_scoredSpotColors,OriginalCoefs,SpotCodeNo)
%% ModScore = o.get_omp_score(z_scoredSpotColors,OriginalCoefs,SpotCodeNo)
% omp score is based on the reduction in error caused by introduction of gene in
% rounds/channels where there is not already a gene.
% It is normalised so independent of spot intensity
% A perfect spot would have a score of about 1 in each round so a total score of 7. 
% (Can get a better than this but just to give a order of magnitude
% for good spots. Actual best score is o.ompScore_LargeErrorMax*o.nRounds). 
%
% z_scoredSpotColors: spot colors that have been z-scored by channel and
% round.
% OriginalCoefs(s,g) is the coefficient of spot s for gene g found by the
% OMP algorithm.
% SpotCodeNo(s): the gene g that OriginalCoefs(s,g) was a peak for in the
% whole tile image. 

if nargin==1
    pf = o.CallMethodPrefix('OMP'); 
    SpotColors = o.([pf,'SpotColors']);
    z_scoredSpotColors = (double(SpotColors)-o.z_scoreSHIFT)./o.z_scoreSCALE;
    OriginalCoefs = o.([pf,'Coefs']);
    SpotCodeNo = o.([pf,'SpotCodeNo']);
    clearvars SpotColors
end

%% Overall score is difference between the residual after all genes 
% added to when all but gene specified by SpotCodeNo(s) added.
nSpots = length(SpotCodeNo);
FinalPredCodes = OriginalCoefs*o.ompBledCodes(:,:);
BestGeneIndex = sub2ind(size(OriginalCoefs),(1:nSpots)',SpotCodeNo);
RemoveGeneCoefs = OriginalCoefs;
RemoveGeneCoefs(BestGeneIndex) = 0;
RemoveGenePredCodes = RemoveGeneCoefs*o.ompBledCodes(:,:);
clear RemoveGeneCoefs BestGeneIndex
FinalError = abs(z_scoredSpotColors(:,:)-FinalPredCodes);
clear FinalPredCodes 
RemoveGeneError = abs(z_scoredSpotColors(:,:)-RemoveGenePredCodes); 
clear z_scoredSpotColors RemoveGenePredCodes
OverallScore = RemoveGeneError - FinalError;

% Save memory
clearvars FinalError 
%% Get multiplier to modify OverallScore so only rounds/channels in 
% the unbled code for specfied gene contribute. Also, if round in
% unbled code but gene efficiency very low (i.e. RCPs failed in that
% round) then don't include. 

nCodes = length(o.CharCodes);
%Need AllChannels for nRoundsUsedNorm calculation - exclude rounds with
%channel not in UseChannels.
GeneMultiplierKeepAllChannels = zeros(o.nRounds*o.nBP,nCodes);
GeneMultiplier = zeros(o.nRounds*o.nBP,nCodes);
%Neglect channels/rounds not in UseChannels/UseRounds.
UnbledCodesReshape = reshape(o.UnbledCodes,nCodes,o.nBP,o.nRounds);
FailedChannels = ~ismember(1:o.nBP,o.UseChannels);
UnbledCodesReshape(:,FailedChannels,:) = 0;
FailedRounds = ~ismember(1:o.nRounds,o.UseRounds);
UnbledCodesReshape(:,:,FailedRounds) = 0;
for g=1:nCodes
    GeneMultiplier(:,g) = UnbledCodesReshape(g,:);
    GeneMultiplierKeepAllChannels(:,g) = o.UnbledCodes(g,:);
end
GeneMultiplier(:,nCodes+1) = 0;     %Use when no overlapping spots
GeneMultiplierKeepAllChannels(:,nCodes+1) = 0;
%Set Low Gene Efficiency Rounds to all not in unbled code 
%i.e. forget these rounds
GeneEfficiency = repelem(o.GeneEfficiency,1,o.nBP)';
GeneEfficiency(:,nCodes+1) = 0;
FailedGeneRound = GeneEfficiency<o.ompScore_GeneEfficiencyThresh;
%Bad gene efficiency or not in UseChannels:
GeneUnbledFailedRounds = FailedGeneRound+GeneMultiplier==2 | ...
    (GeneMultiplierKeepAllChannels==1 & GeneMultiplier == 0);
GeneMultiplier(FailedGeneRound)=0;
ScoreMultiplier = GeneMultiplier(:,SpotCodeNo);
ModScore = OverallScore.*ScoreMultiplier';
%Normalise score because want score metric independent of intensity to go
%with other metrics (SpotIntensity and NeighbNonZeros) that do depend on
%intensity).
NormModScore = ModScore./RemoveGeneError;
NormModScore(ScoreMultiplier'==0)=0;  %Set any Nans due to o.UseChannels to 0.
clear ModScore OverallScore
%% When looking at codes, the gene looks more likely if it explains rounds
% that have large error so give larger weightings to these through
% AbsErrorFactor. Precise form is through trial and error - 
% By larger error we mean more than 80th (o.ompScore_LargeErrorPrcntileThresh)
% percentile and don't want to
% overly weight one anomalously large round/channel hence upper bound.
AbsErrorFactor = RemoveGeneError./prctile(RemoveGeneError',...
    o.ompScore_LargeErrorPrcntileThresh)';
AbsErrorFactor(AbsErrorFactor<1)=1;
AbsErrorFactor(AbsErrorFactor>o.ompScore_LargeErrorMax)=o.ompScore_LargeErrorMax;
clear RemoveGeneError
%% Now want to modify this so neglect rounds which already have a
%positive gene in. Every gene needs to explain something on their own independent 
%of what other genes present
[CoefValues,SortedCoefs]=sort(OriginalCoefs(:,1:nCodes),2,'descend');
clear OriginalCoefs
SortedCoefs(CoefValues<=0) = nCodes+1;
clear CoefValues
%Deal with each overlapping genes in turn until no spot has any more genes
MinCodeNo = 1;
nOverlaps = 0;
ScoreMultiplierOverlaps = zeros(size(ScoreMultiplier));
while MinCodeNo<nCodes+1
    nOverlaps=nOverlaps+1;
    ScoreMultiplierOverlaps = ScoreMultiplierOverlaps+...
        GeneMultiplier(:,SortedCoefs(:,nOverlaps));
    MinCodeNo=min(SortedCoefs(:,nOverlaps));
end
clear SortedCoefs
% %Set genes that are by far the best to have no overlap i.e. use all rounds
% %for clear matches DOESN'T MAKE TOO MUCH DIFFERENCE
% SpotInd = sub2ind(size(OriginalCoefs),(1:length(SpotCodeNo))',SpotCodeNo);
% SpotCoef = OriginalCoefs(SpotInd);
% SortCoef = sort(abs(OriginalCoefs(:,1:73)),2,'descend');
% GoodMatch = SpotCoef-SortCoef(:,2)>0;
% ScoreMultiplierOverlaps(:,GoodMatch) = 0;

%Set rounds/channels overlapping with other gene = 0.
AbsErrorFactor(ScoreMultiplierOverlaps'>1&ScoreMultiplier'==1)=0;
%Set rounds where Gene RCPS failed to 0.
AbsErrorFactor(GeneUnbledFailedRounds(:,SpotCodeNo)')=0;
%For spots with overlapping genes or that match to genes with failed rounds,
%give more weight to other rounds to compensate for less rounds used in score.
nRoundsUsed = max(o.nRounds-sum(AbsErrorFactor==0,2),1);
nRoundsUsedNorm = o.nRounds./nRoundsUsed;
clear ScoreMultiplier ScoreMultiplierOverlaps SpotCodeNo

%% Combine factors
% Again, through trial and error. Use exponential as to put a lower bound
% on how bad single round can be i.e. min(ExpScore) = -AbsErrorFactor.
% Use log(2) factor so limit on how good single round can be is same as
% limit on how bad round can be. 
% I.e. max(NormModScore)=1 so max(ExpScore) = AbsErrorFactor
ExpScore = (exp(NormModScore*log(2))-1).*AbsErrorFactor;
ExpScore = ExpScore.*nRoundsUsedNorm;
ExpScore(ExpScore<-o.ompScore_LargeErrorMax) = -o.ompScore_LargeErrorMax;
ExpScore(ExpScore>o.ompScore_LargeErrorMax) = o.ompScore_LargeErrorMax;
if nSpots>1
    NormExpScore = sum(ExpScore,2);
    %NormExpScore = ExpScore.*nRoundsUsedNorm;
else
%     NormExpScore = reshape(ExpScore.*nRoundsUsedNorm,...
%         o.nRounds,o.nBP);
    NormExpScore = reshape(ExpScore,...
        o.nRounds,o.nBP);
end
end