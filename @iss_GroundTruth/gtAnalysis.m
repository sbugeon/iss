%% GadGroundTruthNewTransform 

%%
%Method = 'Pixel'; %OMP or Pixel
Method = 'OMP';
pf = o.CallMethodPrefix(Method);
        
if strcmpi('OMP',Method)
    %Get primary and secondary sets
    [~,SortedCoefs]=sort(o.ompCoefs(:,1:73)','descend');
    SortedCoefs = SortedCoefs';
    PrimarySet = o.ompSpotCodeNo==SortedCoefs(:,1);
    SecondarySet = o.ompSpotCodeNo==SortedCoefs(:,2);
    %o.ompNeighbThresh = 13;
    %o.ompIntensityThresh = 700;
elseif strcmpi('Pixel',Method)
    %Get primary and secondary sets
    PrimarySet = o.pxSpotScore>0;
    SecondarySet = o.pxSpotScore==0;
    %o.pScoreThresh = 50;
    %o.pScoreThresh2 = 10;
    %o.pIntensityThresh = 200;
    %o.pLogProbThresh2 = 0;
    %o.pIntensityThresh2 = 50;
elseif strcmpi('Spatial',Method)
    PrimarySet = o.spSpotCodeNo>0;
    SecondarySet = o.spSpotCodeNo==0;
end
QualOK = quality_threshold(o,Method);
%QualOK = o.pxNotDuplicate  & (o.pxSpotScore>=0|o.pxLogProbOverBackground>-5);
%QualOK = o.pxNotDuplicate & QualOK;
fprintf('Total Primary Spots: \t\t\t\t%d\n',sum(QualOK&PrimarySet));
fprintf('Total Primary or Secondary Spots: \t\t%d\n',sum(QualOK&(PrimarySet|SecondarySet)));
fprintf('Total Spots: \t\t\t\t\t%d\n',sum(QualOK));

for r=o.gtRounds
    for b=o.UseChannels
        if o.gtGeneNo(r,b)==0; continue; end
        pfTruePosSet = o.([pf,'_gtIdentity']){r,b}==1;
        pfFalsePosSet = o.([pf,'_gtIdentity']){r,b}==2;
        fprintf('There are %d %s peak spots\n', sum(o.gtTruePositiveSet{r,b}),...
            o.GeneNames{o.gtGeneNo(r,b)});
        fprintf('of which, we can achieve %d\n',sum(pfTruePosSet));
        fprintf('False positive set has %d spots.\n',sum(pfFalsePosSet));

        %print excel data
        fprintf('Total Primary True Positives: \t\t\t%d\n',sum(QualOK&PrimarySet&pfTruePosSet));
        fprintf('Total Primary or Secondary True Positives: \t%d\n',...
            sum(QualOK&(PrimarySet|SecondarySet)&pfTruePosSet));
        fprintf('Total True Positives: \t\t\t\t%d\n',sum(QualOK&pfTruePosSet));
        fprintf('Total Primary False Positives: \t\t\t%d\n',sum(QualOK&PrimarySet&pfFalsePosSet));
        fprintf('Total Primary or Secondary False Positives: \t%d\n',...
            sum(QualOK&(PrimarySet|SecondarySet)&pfFalsePosSet));
        fprintf('Total False Positives: \t\t\t\t%d\n',sum(QualOK&pfFalsePosSet));
        fprintf('\n');
    end
end

%% Save data
%Make Tables
MakeTables = false;
if MakeTables
    FileLocation = {''};
    Method = {''};
    Intensity_Method = {''};
    ScoreScale = 0;
    pScoreThresh = 0;
    pIntensityThresh = 0;
    pScoreThresh2 = 0;
    pIntensityThresh2 = 0;
    pLogProbThresh2 = 0;
    pQualThresh1 = 0;
    pQualParam1 = 0;
    pQualThresh2 = 0;
    pQualParam2 = 0;
    ompNeighbThresh = 0;
    ompIntensityThresh = 0;
    ompScoreThresh = 0;
    nPrimarySpots = 0;
    nPrimarySecondarySpots = 0;
    nTotalSpots = 0;
    Combined_Score = 0;
    Combined_nTP = 0;
    Combined_nTP_Missed = 0;
    Combined_nFP = 0;
    TruePosData.Summary = table(FileLocation,Method,Intensity_Method,ScoreScale,pScoreThresh,pIntensityThresh,...
        pScoreThresh2,pIntensityThresh2,pLogProbThresh2,pQualThresh1,pQualParam1,pQualThresh2,...
        pQualParam2,ompNeighbThresh,ompIntensityThresh,ompScoreThresh,nPrimarySpots,nPrimarySecondarySpots,...
        nTotalSpots,Combined_Score,Combined_nTP,Combined_nTP_Missed,Combined_nFP);
    nTP = 0;
    TP_Max = 0;
    TP_Primary = 0;
    TP_PrimarySecondary = 0;
    TP_Total = 0;
    FP_Max = 0;
    FP_Primary = 0;
    FP_PrimarySecondary = 0;
    FP_Total = 0;
    for r=o.gtRounds
        for b=o.UseChannels
            if o.gtGeneNo(r,b)==0; continue; end
            TruePosData.(cell2mat(o.GeneNames(o.gtGeneNo(r,b)))) = ...
                table(nTP,TP_Max,TP_Primary,TP_PrimarySecondary,TP_Total,...
                FP_Max,FP_Primary,FP_PrimarySecondary,FP_Total);
        end
    end
end
%% Add data to tables
% Run after the above so have QualOK and Primary.

i = size(TruePosData.Summary,1)+1;      %INDEX OF DATA TO BE ADDED
TruePosData.Summary(i,'FileLocation') = ...
    {fullfile(o.OutputDirectory, 'oOMP_NewPCR_NewSpotScore')};
%Method = 'Pixel: pLogThresh, ProbMethod = 1, GammaShape=3, NoFilter, Smooth, pQualThresh3=100, pQualThresh=-25';
%Method = 'Pixel: pQualThresh, ProbMethod = 1, pQualThresh3=51, Median BleedMatrix, Used GeneEfficiencies, Used get_secondary_gene_prob with remove_full_code';
%Method = 'MP_Weights: Use DotProduct>2 as only thresh in OMP. Stop MP iteration when any artificial gene found, SpotIntensity is median(Z_scored in Unbled code) and no subtraction, weight by round, FinalThresholding was QualOK = quality_threshold(o,Method), o.ompNeighbThresh2=5. Proper initial OMP, 7 Background Channel Strips, Gene Efficiencies found from Mean Code';
%Method = 'Spatial With ompBledCodes(ompBledCodes<0.1)=0';
Method = 'OMP: Use PCR with max 1500 spots on each imaging tile, o.ompNeighbThresh2=10, o.ompScoreThresh2=1.1. Include Higher Thresh for Non Max Genes';
%Method = 'Spatial';
Intensity_Method = 'Median Unbled Z_scored';
%Intensity_Method = 'Mean Unbled';
%Intensity_Method = 'Median Unbled';
%Intensity_Method = 'Mean Unbled - Mean Not Unbled, Z-Scored';

TruePosData.Summary(i,'Method') = {Method};
TruePosData.Summary(i,'Intensity_Method') = {Intensity_Method};
for k=1:size(TruePosData.Summary,2)
    try
        TruePosData.Summary(i,k) = {nan};
    end
end

if contains(Method,'Pixel: pLogThresh')
    if o.ProbMethod==1
        TruePosData.Summary(i,'ScoreScale') = {o.ScoreScale};
    end
    TruePosData.Summary(i,'pScoreThresh') = {o.pScoreThresh};
    TruePosData.Summary(i,'pIntensityThresh') = {o.pIntensityThresh};
    TruePosData.Summary(i,'pScoreThresh2') = {o.pScoreThresh2};
    TruePosData.Summary(i,'pIntensityThresh2') = {o.pIntensityThresh2};
    TruePosData.Summary(i,'pLogProbThresh2') = {o.pLogProbThresh2};
elseif contains(Method,'Pixel: pQualThresh')
    TruePosData.Summary(i,'Intensity_Method') = {'N/A'};
    if o.ProbMethod==1
        TruePosData.Summary(i,'ScoreScale') = {o.ScoreScale};
    end
    TruePosData.Summary(i,'pIntensityThresh') = {o.pIntensityThresh};
    TruePosData.Summary(i,'pQualThresh1') = {o.pQualThresh1};
    TruePosData.Summary(i,'pQualParam1') = {o.pQualParam1};
    TruePosData.Summary(i,'pQualThresh2') = {o.pQualThresh2};
    TruePosData.Summary(i,'pQualParam2') = {o.pQualParam2};
elseif contains(Method,'OMP')
    TruePosData.Summary(i,'ompNeighbThresh') = {o.ompNeighbThresh};
    TruePosData.Summary(i,'ompIntensityThresh') = {o.ompIntensityThresh};
    TruePosData.Summary(i,'ompScoreThresh') = {o.ompScoreThresh};
    TruePosData.Summary(i,'pIntensityThresh') = {o.ompIntensityThresh2};
end
TruePosData.Summary(i,'nPrimarySpots') = {sum(QualOK&PrimarySet)};
TruePosData.Summary(i,'nPrimarySecondarySpots') = {sum(QualOK&(PrimarySet|SecondarySet))};
TruePosData.Summary(i,'nTotalSpots') = {sum(QualOK)};
Combined_Score = 0;
Combined_nTP = 0;
Combined_nTP_Missed = 0;
Combined_nFP = 0;
for r=o.gtRounds
    for b=o.UseChannels
        if o.gtGeneNo(r,b)==0; continue; end
        pfTruePosSet = o.([pf,'_gtIdentity']){r,b}==1;
        TruePosData.(cell2mat(o.GeneNames(o.gtGeneNo(r,b))))(i,'nTP') = ...
            {sum(o.gtTruePositiveSet{r,b})};
        TruePosData.(cell2mat(o.GeneNames(o.gtGeneNo(r,b))))(i,'TP_Max') = ...
            {sum(pfTruePosSet)};
        TruePosData.(cell2mat(o.GeneNames(o.gtGeneNo(r,b))))(i,'TP_Primary') = ...
            {sum(QualOK&PrimarySet&pfTruePosSet)}; 
        TruePosData.(cell2mat(o.GeneNames(o.gtGeneNo(r,b))))(i,'TP_PrimarySecondary') = ...
            {sum(QualOK&(PrimarySet|SecondarySet)&pfTruePosSet)};
        TruePosData.(cell2mat(o.GeneNames(o.gtGeneNo(r,b))))(i,'TP_Total') = ...
            {sum(QualOK&pfTruePosSet)};
        
        pfFalsePosSet = o.([pf,'_gtIdentity']){r,b}==2;
        TruePosData.(cell2mat(o.GeneNames(o.gtGeneNo(r,b))))(i,'FP_Max') = ...
            {sum(pfFalsePosSet)};
        TruePosData.(cell2mat(o.GeneNames(o.gtGeneNo(r,b))))(i,'FP_Primary') = ...
            {sum(QualOK&PrimarySet&pfFalsePosSet)};
        TruePosData.(cell2mat(o.GeneNames(o.gtGeneNo(r,b))))(i,'FP_PrimarySecondary') = ...
            {sum(QualOK&(PrimarySet|SecondarySet)&pfFalsePosSet)};
        TruePosData.(cell2mat(o.GeneNames(o.gtGeneNo(r,b))))(i,'FP_Total') = ...
            {sum(QualOK&pfFalsePosSet)};
        
        Combined_nTP = Combined_nTP+sum(QualOK&pfTruePosSet);
        Combined_nTP_Missed = Combined_nTP_Missed+sum(~QualOK&pfTruePosSet)+...
            sum(o.gtTruePositiveSet{r,b})-sum(pfTruePosSet);
        Combined_nFP = Combined_nFP+sum(QualOK&pfFalsePosSet);
        %Combined_Score = Combined_Score+sum(QualOK&pfTruePosSet)/sum(pfTruePosSet)+...
        %    5*sum(~QualOK&pfFalsePosSet)/sum(pfFalsePosSet);    
        %Min Score is best
        Combined_Score = Combined_Score+sum(~QualOK&pfTruePosSet)+...
            sum(o.gtTruePositiveSet{r,b})-sum(pfTruePosSet)+2*sum(QualOK&pfFalsePosSet);
    end
end
TruePosData.Summary(i,'Combined_Score') = {Combined_Score};
TruePosData.Summary(i,'Combined_nTP') = {Combined_nTP};
TruePosData.Summary(i,'Combined_nTP_Missed') = {Combined_nTP_Missed};
TruePosData.Summary(i,'Combined_nFP') = {Combined_nFP};

%save(fullfile(o.OutputDirectory, 'GroundTruth_Data'), 'TruePosData', '-v7.3');
%% Find Best Params - PixelBased
%QualThresh method
% QualThresh2 = o.pQualThresh2;
% QualParam2 = o.pQualParam2;
% QualThresh1 = -300:100:100;
% QualParam1 =  0:0.5:3;
% QualThresh1 = -140:5:-60;
% QualParam1 =  0.8:0.05:1.3;

% QualThresh1 = o.pQualThresh1;
% QualParam1 = o.pQualParam1;
% % QualThresh2 = -250:50:150;
% % QualParam2 = 0:0.5:10;
% QualThresh2 = -40:5:40;
% QualParam2 = 1.6:0.05:2.4;
% %QualParam2 = 0;   %Doesn't affect it at all
% ScoreImage = zeros(length(QualThresh1),length(QualParam1),length(QualThresh2),length(QualParam2));
% for i=1:length(QualThresh1)
%     o.pQualThresh1 = QualThresh1(i);
%     for j=1:length(QualParam1)
%         o.pQualParam1 = QualParam1(j);
%         for k=1:length(QualThresh2)
%             o.pQualThresh2 = QualThresh2(k);
%             for k2=1:length(QualParam2)
%                 o.pQualParam2 = QualParam2(k2);
%                 QualOK = quality_threshold(o,'Pixel');
%                 for r=o.gtRounds
%                     for b=o.UseChannels
%                         if o.gtGeneNo(r,b)==0; continue; end
%                         pfTruePosSet = o.([pf,'_gtIdentity']){r,b}==1;
%                         pfFalsePosSet = o.([pf,'_gtIdentity']){r,b}==2;
%                         %ScoreImage(i,j,k,k2) = ScoreImage(i,j,k,k2)+sum(QualOK&pfTruePosSet)/sum(pfTruePosSet)+...
%                         %     5*sum(~QualOK&pfFalsePosSet)/sum(pfFalsePosSet);
%                         ScoreImage(i,j,k,k2) = ScoreImage(i,j,k,k2)+sum(~QualOK&pfTruePosSet)+...
%                             sum(o.gtTruePositiveSet{r,b})-sum(pfTruePosSet)+2*sum(QualOK&pfFalsePosSet);
%                     end
%                 end
%             end
%         end
%     end
% end
% [a,b] = min(ScoreImage(:));
% a
% [a,b,c,d]=ind2sub(size(ScoreImage),b);
% o.pQualThresh1=QualThresh1(a);
% o.pQualParam1=QualParam1(b);
% o.pQualThresh2=QualThresh2(c);
% o.pQualParam2=QualParam2(d);
% QualThresh1(a)
% QualParam1(b)
% QualThresh2(c)
% QualParam2(d)

%QualThresh3 = 0:25:300;
%QualThresh4 = -300:25:0;
QualThresh3 = 50:1:100;
%QualThresh4 = -50:5:0;
QualThresh4 = inf;
ScoreImage = zeros(length(QualThresh3),length(QualThresh4));
for i=1:length(QualThresh3)
    o.pQualThresh3 = QualThresh3(i);
    for j=1:length(QualThresh4)
        o.pQualThresh4 = QualThresh4(j);
        QualOK = quality_threshold(o,'Pixel');
        for r=o.gtRounds
            for b=o.UseChannels
                if o.gtGeneNo(r,b)==0; continue; end
                pfTruePosSet = o.([pf,'_gtIdentity']){r,b}==1;
                pfFalsePosSet = o.([pf,'_gtIdentity']){r,b}==2;
                %ScoreImage(i,j,k,k2) = ScoreImage(i,j,k,k2)+sum(QualOK&pfTruePosSet)/sum(pfTruePosSet)+...
                %     5*sum(~QualOK&pfFalsePosSet)/sum(pfFalsePosSet);
                ScoreImage(i,j) = ScoreImage(i,j)+sum(~QualOK&pfTruePosSet)+...
                    sum(o.gtTruePositiveSet{r,b})-sum(pfTruePosSet)+2*sum(QualOK&pfFalsePosSet);
            end
        end
    end
end
[a,b] = min(ScoreImage(:));
a
[a,b]=ind2sub(size(ScoreImage),b);
o.pQualThresh3=QualThresh3(a);
o.pQualThresh4=QualThresh4(b);
QualThresh3(a)
QualThresh4(b)

%figure; imagesc(ScoreImage);

%pLogThresh method
% % ScoreThresh = -100:100:100;
% % IntensityThresh = -10000:2000:2000;
% % o.pScoreThresh2 = -0.000001;
% % IntensityThresh2 = -10000:2000:2000;
% % LogProbThresh2 = -300:100:100;
% ScoreThresh = 50:5:100;
% IntensityThresh = -250:10:-200;
% o.pScoreThresh2 = -0.000001;
% IntensityThresh2 = -inf;
% LogProbThresh2 = -40:5:-20;
% ScoreImage = zeros(length(ScoreThresh),length(IntensityThresh),length(IntensityThresh2),length(LogProbThresh2));
% for i=1:length(ScoreThresh)
%     o.pScoreThresh = ScoreThresh(i);
%     for j=1:length(IntensityThresh)
%         o.pIntensityThresh = IntensityThresh(j);
%         for k=1:length(IntensityThresh2)
%             o.pIntensityThresh2 = IntensityThresh2(k);
%             for k2=1:length(LogProbThresh2)
%                 o.pLogProbThresh2 = LogProbThresh2(k2);
%                 QualOK = quality_threshold(o,'Pixel');
%                 for r=o.gtRounds
%                     for b=o.UseChannels
%                         if o.gtGeneNo(r,b)==0; continue; end
%                         pfTruePosSet = o.([pf,'_gtIdentity']){r,b}==1;
%                         pfFalsePosSet = o.([pf,'_gtIdentity']){r,b}==2;
%                         %ScoreImage(i,j,k,k2) = ScoreImage(i,j,k,k2)+sum(QualOK&pfTruePosSet)/sum(pfTruePosSet)+...
%                         %     5*sum(~QualOK&pfFalsePosSet)/sum(pfFalsePosSet);
%                         ScoreImage(i,j,k,k2) = ScoreImage(i,j,k,k2)+sum(~QualOK&pfTruePosSet)+...
%                             sum(o.gtTruePositiveSet{r,b})-sum(pfTruePosSet)+2*sum(QualOK&pfFalsePosSet);
%                     end
%                 end
%             end
%         end
%     end
% end
% %figure; imagesc(ScoreImage);
% [a,b] = min(ScoreImage(:));
% a
% [a,b,c,d]=ind2sub(size(ScoreImage),b);
% o.pScoreThresh=ScoreThresh(a);
% o.pIntensityThresh=IntensityThresh(b);
% o.pIntensityThresh2=IntensityThresh2(c);
% o.pLogProbThresh2=LogProbThresh2(d);
% ScoreThresh(a)
% IntensityThresh(b)
% IntensityThresh2(c)
% LogProbThresh2(d)

%% Testing OMP thresholds
%%
nCodes = length(o.CharCodes);
MaxCoef = max(abs(o.ompCoefs(:,1:nCodes )),[],2);
SpotInd = sub2ind(size(o.ompCoefs),(1:length(o.ompSpotCodeNo))',o.ompSpotCodeNo);
SpotCoef = o.ompCoefs(SpotInd);
%%        
o.ompNeighbThresh2 = 10;
o.ompIntensityThresh2 = 0.005;
o.ompScoreThresh2 = 1.1;
o.ompIntensityThresh3 = 0.01;
o.ompIntensityThresh3_CoefDiffFactor = 0.27;
o.ompNeighbThresh3 = 28;
o.ompScoreThresh3 = 6.9;
% NeighbThresh = 1:4:33;
% QualParam2 = -2:0.5:3;
%ScoreThresh = -4:4:24;
NeighbThresh = 16:1:20;
QualParam2 = 0.3:0.1:0.7;
ScoreThresh = 4.1:0.1:4.5;
% ScoreThresh = 3.2:0.1:3.8;
% ScoreThresh = 0:0.001:0.9;
%ScoreThresh = 4:0.1:6;
%NeighbThresh = o.ompNeighbThresh;
%QualParam2 = o.ompIntensityThresh;
% ScoreThresh = o.ompScoreThresh;
ScoreImage = zeros(length(NeighbThresh),length(QualParam2),length(ScoreThresh));
        
for n=1:length(NeighbThresh)
    o.ompNeighbThresh=NeighbThresh(n);
    for i=1:length(QualParam2)
        o.ompIntensityThresh=QualParam2(i);
        for s=1:length(ScoreThresh)
            o.ompScoreThresh=ScoreThresh(s);
            QualOK = quality_threshold(o,'OMP',MaxCoef,SpotCoef);
            %Use = MaxCoef-SpotCoef>0.26;
            %QualOK(Use) = QualOK(Use) & (o.ompSpotIntensity2(Use)>0.03 | o.ompSpotScore(Use)>o.ompScoreThresh*2);
       
            %Use = MaxCoef-SpotCoef>0;
            %QualOK(Use) = QualOK(Use) & (o.ompSpotIntensity2(Use)>0.01.*(1+round((MaxCoef(Use)-SpotCoef(Use))./0.27)) | o.ompSpotScore(Use)>o.ompScoreThresh*1.6 | o.ompNeighbNonZeros(Use)>o.ompNeighbThresh*1.3);
            for r=o.gtRounds
                for b=o.UseChannels
                    if o.gtGeneNo(r,b)==0; continue; end
                    pfTruePosSet = o.([pf,'_gtIdentity']){r,b}==1;
                    pfFalsePosSet = o.([pf,'_gtIdentity']){r,b}==2;
                    %ScoreImage(i,j,k,k2) = ScoreImage(i,j,k,k2)+sum(QualOK&pfTruePosSet)/sum(pfTruePosSet)+...
                    %     5*sum(~QualOK&pfFalsePosSet)/sum(pfFalsePosSet);
                    ScoreImage(n,i,s) = ScoreImage(n,i,s)+sum(~QualOK&pfTruePosSet)+...
                        sum(o.gtTruePositiveSet{r,b})-sum(pfTruePosSet)+2*sum(QualOK&pfFalsePosSet);
                end
            end
        end
    end
end
[a,b] = min(ScoreImage(:));
a
[a,b,c]=ind2sub(size(ScoreImage),b);
o.ompNeighbThresh=NeighbThresh(a);
o.ompIntensityThresh=QualParam2(b);
o.ompScoreThresh = ScoreThresh(c);
NeighbThresh(a)
QualParam2(b)
ScoreThresh(c)

%% Initial OMP thresholding
% Store initial ground truth data first
% o.CallMethods = {'DotProduct','Prob','Pixel','OMP','Spatial','GroundTruth','iOMP'};
% o.CallMethodPrefix = ...
%     containers.Map({'DotProduct','Prob','Pixel','OMP','Spatial','GroundTruth','iOMP'},...
%     {'dp','p','px','omp','sp','gt_px','iomp'});
% o = o.get_gtRoundColor('iOMP');
% o = get_pf_gtTruePositiveSets(o,'iOMP');
pf = o.CallMethodPrefix('iOMP');

%QualThresh method
% QualThresh2 = o.pQualThresh2;
% % QualThresh1 = 0:0.01:0.4;
% % QualParam1 =  0:0.5:3;
% QualThresh1 = 0:0.001:0.1;
% QualParam1 =  0:0.05:1.3;

% QualThresh1 = o.pQualThresh1;
% QualParam1 = o.pQualParam1;
% % QualThresh2 = -250:50:150;
% QualThresh2 = 0:0.01:1;
% ScoreImage = zeros(length(QualThresh1),length(QualParam1),length(QualThresh2));
% for i=1:length(QualThresh1)
%     o.pQualThresh1 = QualThresh1(i);
%     for j=1:length(QualParam1)
%         o.pQualParam1 = QualParam1(j);
%         for k=1:length(QualThresh2)
%             o.pQualThresh2 = QualThresh2(k);
%             QualOK1 = o.iompSpotScore>0 & o.iompResOverBackground+o.pQualParam1*o.iompSpotScore>o.pQualThresh1;
%             QualOK2 = o.iompSpotScore==0 & o.iompResOverBackground>o.pQualThresh2;
%             QualOK3 = o.iompSpotScore<0 & o.iompResOverBackground>o.pQualThresh3 & o.iompSpotScore>o.pQualThresh4;
%             QualOK = (QualOK1 | QualOK2 | QualOK3) & o.iompSpotIntensity>o.pIntensityThresh & o.iompCoef>0;
%             for r=o.gtRounds
%                 for b=o.UseChannels
%                     if o.gtGeneNo(r,b)==0; continue; end
%                     pfTruePosSet = o.([pf,'_gtIdentity']){r,b}==1;
%                     pfFalsePosSet = o.([pf,'_gtIdentity']){r,b}==2;
%                     %ScoreImage(i,j,k) = ScoreImage(i,j,k,k2)+sum(QualOK&pfTruePosSet)/sum(pfTruePosSet)+...
%                     %     5*sum(~QualOK&pfFalsePosSet)/sum(pfFalsePosSet);
%                     ScoreImage(i,j,k) = ScoreImage(i,j,k)+sum(~QualOK&pfTruePosSet)+...
%                         sum(o.gtTruePositiveSet{r,b})-sum(pfTruePosSet)+2*sum(QualOK&pfFalsePosSet);
%                 end
%             end
%         end
%     end
% end
% [a,b] = min(ScoreImage(:));
% a
% [a,b,c]=ind2sub(size(ScoreImage),b);
% o.pQualThresh1=QualThresh1(a);
% o.pQualParam1=QualParam1(b);
% o.pQualThresh2=QualThresh2(c);
% QualThresh1(a)
% QualParam1(b)
% QualThresh2(c)

QualThresh3 = 0:0.01:0.7;
QualThresh4 = -0.7:0.01:0;
% QualThresh3 = 50:1:100;
% QualThresh4 = -50:5:0;
% QualThresh4 = inf;
ScoreImage = zeros(length(QualThresh3),length(QualThresh4));
for i=1:length(QualThresh3)
    o.pQualThresh3 = QualThresh3(i);
    for j=1:length(QualThresh4)
        o.pQualThresh4 = QualThresh4(j);
        QualOK1 = o.iompSpotScore>0 & o.iompResOverBackground+o.pQualParam1*o.iompSpotScore>o.pQualThresh1;
        QualOK2 = o.iompSpotScore==0 & o.iompResOverBackground>o.pQualThresh2;
        QualOK3 = o.iompSpotScore<0 & o.iompResOverBackground>o.pQualThresh3 & o.iompSpotScore>o.pQualThresh4;
        QualOK = (QualOK1 | QualOK2 | QualOK3) & o.iompSpotIntensity>o.pIntensityThresh & o.iompCoef>0;
        for r=o.gtRounds
            for b=o.UseChannels
                if o.gtGeneNo(r,b)==0; continue; end
                pfTruePosSet = o.([pf,'_gtIdentity']){r,b}==1;
                pfFalsePosSet = o.([pf,'_gtIdentity']){r,b}==2;
                %ScoreImage(i,j,k,k2) = ScoreImage(i,j,k,k2)+sum(QualOK&pfTruePosSet)/sum(pfTruePosSet)+...
                %     5*sum(~QualOK&pfFalsePosSet)/sum(pfFalsePosSet);
                ScoreImage(i,j) = ScoreImage(i,j)+sum(~QualOK&pfTruePosSet)+...
                    sum(o.gtTruePositiveSet{r,b})-sum(pfTruePosSet)+2*sum(QualOK&pfFalsePosSet);
            end
        end
    end
end
[a,b] = min(ScoreImage(:));
a
[a,b]=ind2sub(size(ScoreImage),b);
o.pQualThresh3=QualThresh3(a);
o.pQualThresh4=QualThresh4(b);
QualThresh3(a)
QualThresh4(b)


%% Spatial
ScoreThresh = 0:0.01:0.5;
ScoreImage = zeros(size(ScoreThresh));
for s=1:length(ScoreThresh)
    QualOK = A>ScoreThresh(s);
    for r=o.gtRounds
        for b=o.UseChannels
            if o.gtGeneNo(r,b)==0; continue; end
            pfTruePosSet = o.([pf,'_gtIdentity']){r,b}==1;
            pfFalsePosSet = o.([pf,'_gtIdentity']){r,b}==2;
            %ScoreImage(i,j,k,k2) = ScoreImage(i,j,k,k2)+sum(QualOK&pfTruePosSet)/sum(pfTruePosSet)+...
            %     5*sum(~QualOK&pfFalsePosSet)/sum(pfFalsePosSet);
            ScoreImage(s) = ScoreImage(s)+sum(~QualOK&pfTruePosSet)+...
                sum(o.gtTruePositiveSet{r,b})-sum(pfTruePosSet)+2*sum(QualOK&pfFalsePosSet);
        end
    end
end
[a,b] = min(ScoreImage(:));
a
ScoreThresh(b)

%% What genes are near missed ground truth?
r = 8;
b = 4;
g = o.gtGeneNo(r,b);
MaxDist = 5;
MaxGenes = 12;
GT_MissedYX = o.gtSpotGlobalYX{r,b}(o.omp_gtFound{r,b}==2,:);
nMissedGT = size(GT_MissedYX,1);
GT_FoundYX = o.gtSpotGlobalYX{r,b}(o.omp_gtFound{r,b}==1,:);
nFoundGT = size(GT_FoundYX,1);
QualOK = o.quality_threshold('OMP');
QualOK_CodeNo = o.ompSpotCodeNo(QualOK);
tree = KDTreeSearcher(o.ompSpotGlobalYX(QualOK,:));
[Ind_Missed,Dist_Missed] = tree.knnsearch(GT_MissedYX,'K',MaxGenes);
[Ind_Found,Dist_Found] = tree.knnsearch(GT_FoundYX,'K',MaxGenes);
Dist_Found(QualOK_CodeNo(Ind_Found)==g) = inf;      %Don't show g for found ground truth
nNearCoefsMissed = sum(Dist_Missed<=MaxDist,2);
nNearCoefsFound = sum(Dist_Found<=MaxDist,2);
figure; histogram(nNearCoefsMissed,-0.5:MaxGenes+0.5,'Normalization','probability');
hold on;
histogram(nNearCoefsFound,-0.5:MaxGenes+0.5,'Normalization','probability');
legend(['Missed Ground Truth Spots: ',num2str(nMissedGT)],...
    ['Found Ground Truth Spots: ',num2str(nFoundGT)]);
title(sprintf('Number of Genes found near %s Ground Truth',o.GeneNames{g}));
figure; 
Hist = histogram(QualOK_CodeNo(Ind_Missed(Dist_Missed<MaxDist)),0.5:73.5);
GeneNoCounts = Hist.BinCounts;
hold on;
histogram(QualOK_CodeNo(Ind_Found(Dist_Found<MaxDist)),0.5:73.5);
legend(['Missed Ground Truth Spots: ',num2str(nMissedGT)],...
    ['Found Ground Truth Spots: ',num2str(nFoundGT)]);
title(sprintf('Genes found near %s Ground Truth',o.GeneNames{g}));
xticks(1:73);
xticklabels(o.GeneNames);
xtickangle(90);

%% Plot Gene Bled Codes of most common genes to missed ground truth
nModeNearGenes = 2;
[~,ModeGeneNumbers] = sort(GeneNoCounts,2,'descend');
%find which rounds/channels overlap
OverlapBledCodeThresh = 0.04;
BledCodeOverlap_br = cell(nModeNearGenes,1);
AllGlyphs = ['x','d','*','o','>','^'];
GlyphLineWidth = 2;
rShift = -0.25:0.5/(nModeNearGenes-1):0.25; % So glyphs don't overlap
for i=1:nModeNearGenes
    g_i = ModeGeneNumbers(i);
    IndOverlap = find(o.ompBledCodes(g,:)>OverlapBledCodeThresh &...
        o.ompBledCodes(g_i,:)>OverlapBledCodeThresh);
    [Overlap_b,Overlap_r] = ind2sub([o.nBP,o.nRounds],IndOverlap);
    Overlap_r = Overlap_r + rShift(i);
    BledCodeOverlap_br{i} = [Overlap_b',Overlap_r'];
end
CLims = [min(o.ompBledCodes(:)), max(o.ompBledCodes(:))];
figure; 
subplot(nModeNearGenes+1,1,1);
imagesc(reshape(o.ompBledCodes(g,:),[o.nRounds,o.nBP]));
title(sprintf('Ground Truth Bled Code, gene %.0f: %s',...
    g,o.GeneNames{g}));
xticks(1:o.nRounds);
xticklabels([]);
yticks(1:o.nBP);
yticklabels(o.bpLabels);
ylabel('Channel');
caxis(CLims);
hold on;
for i=1:nModeNearGenes
    scatter(BledCodeOverlap_br{i}(:,2),BledCodeOverlap_br{i}(:,1),...
        'r',AllGlyphs(i),'LineWidth',GlyphLineWidth);
end
for i=1:nModeNearGenes
    g_i = ModeGeneNumbers(i);
    subplot(nModeNearGenes+1,1,i+1);
    imagesc(reshape(o.ompBledCodes(g_i,:),[o.nRounds,o.nBP]));
    nGT_spots_near_gi = ...
        sum(sum(QualOK_CodeNo(Ind_Missed)==g_i & Dist_Missed<MaxDist,2)>0);
    title(sprintf('Gene %.0f: %s. %.0f/%.0f missed ground truth near this',...
        g_i,o.GeneNames{g_i},nGT_spots_near_gi,nMissedGT));
    xticks(1:o.nRounds);
    if i==nModeNearGenes
        xticklabels(1:7);
        xlabel('Round');
    else
        xticklabels([]);
    end
    yticks(1:o.nBP);
    yticklabels(o.bpLabels);
    ylabel('Channel');
    caxis(CLims);
    hold on
    scatter(BledCodeOverlap_br{i}(:,2),BledCodeOverlap_br{i}(:,1),...
        'r',AllGlyphs(i),'LineWidth',GlyphLineWidth);
   
end

