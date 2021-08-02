function QualOK = quality_threshold(o,Method,MaxCoef,SpotCoef)
% QualOK = o.quality_threshold
% quick function that returns a binary saying which spots are above quality
% threshold
%Method = 'DotProduct','Prob','Pixel' or 'OMP' to consider gene assignments given
%by o.dpSpotCodeNo, o.pSpotCodeNo, o.pxSpotCodeNo, o.ompSpotCodeNo respectively.
%MaxCoef and SpotCoef can be provided to compute the QualOK('OMP') more
%quickly if running in a for loop e.g. gtAnalysis.m

if nargin<2 || isempty(Method)
    Method = 'DotProduct';
end

if ~ismember({Method},o.CallMethods)
    error('Method invalid, must be member of o.CallMethods.');
end
pf = o.CallMethodPrefix(Method);

if strcmpi('Prob',Method) || strcmpi('Pixel',Method) || strcmpi('GroundTruth',Method)
%     QualOK = (o.([pf,'SpotScore'])>o.pScoreThresh & o.([pf,'SpotIntensity'])>o.pIntensityThresh2 | ...
%     o.([pf,'SpotScore'])>o.pScoreThresh2 & o.([pf,'SpotScore'])+o.([pf,'LogProbOverBackground'])>o.pLogProbThresh2 &...
%     o.([pf,'SpotIntensity'])>o.pIntensityThresh); 
%     QualOK = QualOK & o.([pf,'SpotScore'])>=0;       %3rd best spots do more harm than good.
    %QualOK = QualOK & o.pSpotIntensity2 > o.pIntensity2Thresh;
    QualOK1 = o.([pf,'SpotScore'])>0 & o.([pf,'LogProbOverBackground'])+o.pQualParam1*o.([pf,'SpotScore'])>o.pQualThresh1;
    if strcmpi('Prob',Method) || strcmpi('GroundTruth',Method) || (strcmpi('Pixel',Method)&&isempty(o.([pf,'SpotScore2'])))
        QualOK2 = o.([pf,'SpotScore'])==0 & o.([pf,'LogProbOverBackground'])+o.pQualParam2*o.([pf,'SpotScore'])>o.pQualThresh2;
        QualOK3 = o.([pf,'SpotScore'])<0 & o.([pf,'LogProbOverBackground'])>o.pQualThresh3 & o.([pf,'SpotScore'])>o.pQualThresh4;
    elseif strcmpi('Pixel',Method) && ~isempty(o.([pf,'SpotScore2']))
        QualOK2 = o.([pf,'SpotScore2'])>0 & o.([pf,'LogProbOverBackground2'])+o.pQualParam2*o.([pf,'SpotScore2'])>o.pQualThresh2;
        QualOK3 = o.([pf,'SpotScore2'])==0 & o.([pf,'LogProbOverBackground2'])>o.pQualThresh3;
    end
    QualOK = (QualOK1 | QualOK2 | QualOK3) & o.([pf,'SpotIntensity'])>o.pIntensityThresh;
elseif strcmpi('DotProduct',Method)
    QualOK = o.([pf,'SpotCombi']) & o.([pf,'SpotScore'])>o.CombiQualThresh &...
        o.([pf,'SpotIntensity'])>o.CombiIntensityThresh & o.([pf,'SpotScoreDev'])>o.CombiDevThresh;
elseif strcmpi('OMP',Method)
    if nargin<3 || isempty(MaxCoef) || isempty(SpotCoef)
        nCodes = length(o.CharCodes);
        MaxCoef = max(abs(o.ompCoefs(:,1:nCodes )),[],2);
        SpotInd = sub2ind(size(o.ompCoefs),(1:length(o.ompSpotCodeNo))',o.ompSpotCodeNo);
        SpotCoef = o.ompCoefs(SpotInd);
    end
     %Old method below
     %QualOK = o.ompNeighbNonZeros>o.ompNeighbThresh | ...
     %    (o.ompSpotIntensity>o.ompIntensityThresh & o.ompNeighbNonZeros>o.ompNeighbThresh2);
     %QualOK = QualOK & o.ompSpotIntensity2 > o.ompIntensity2Thresh;
     %New method, found using PyTorch
     %Use median SpotIntensity2 not mean SpotIntensity as more robust to
     %cases where OMP is bad with one good round. 
     NeighbNonZeros = o.ompNeighbNonZeros(:,1)*o.ompNeighbNearPosNeighbMultiplier+o.ompNeighbNonZeros(:,2);
     QualOK = NeighbNonZeros>o.ompNeighbThresh | o.([pf,'SpotIntensity2'])>o.ompIntensityThresh |...
         o.([pf,'SpotScore'])>o.ompScoreThresh;
     QualOK = QualOK & o.([pf,'SpotIntensity2'])>o.ompIntensityThresh2 & ...
         NeighbNonZeros>o.ompNeighbThresh2 & o.([pf,'SpotScore'])>o.ompScoreThresh2;
     %Spots assigned to genes that are not max coefficient have stronger
     %thresholds:
     NotBestGene = MaxCoef>SpotCoef;
     QualOK(NotBestGene) = QualOK(NotBestGene) & (o.ompSpotIntensity2(NotBestGene)>...
         o.ompIntensityThresh3.*(1+round((MaxCoef(NotBestGene)-SpotCoef(NotBestGene))./o.ompIntensityThresh3_CoefDiffFactor)) | ...
         o.ompSpotScore(NotBestGene)>o.ompScoreThresh3 | NeighbNonZeros(NotBestGene)>o.ompNeighbThresh3);
%      CoefInd = sub2ind(size(o.ompCoefs),[1:length(o.([pf,'SpotCodeNo']))]',o.([pf,'SpotCodeNo']));
%      QualOK = o.([pf,'NeighbNonZeros'])>o.ompNeighbThresh | o.([pf,'Coefs'])(CoefInd)>o.ompIntensityThresh |...
%          o.([pf,'SpotScore'])>o.ompScoreThresh;
%      QualOK = QualOK & o.([pf,'Coefs'])(CoefInd)>o.pIntensityThresh & o.([pf,'NeighbNonZeros'])>o.ompNeighbThresh2;
    
     
elseif strcmpi('Spatial',Method)
    QualOK = o.([pf,'SpotScore'])>0;        %All spots at the moment.
else
    %If new method, just accept everything
    QualOK = o.([pf,'SpotCodeNo'])>0;   
end


% % HACK ALERT
% QualOK = QualOK & o.cSpotIsolated;

% nCombiCodes = sum(~strcmp(o.CharCodes, 'EXTRA'));
% 
% % now extras - they have their own thresholds, set manually for each type
% for i=1:size(o.ExtraCodes,1)
%     MySpots = (o.SpotCodeNo == nCombiCodes+i);
%     QualOK(MySpots) = o.SpotIntensity(MySpots)>o.ExtraCodes{i,4};
% end