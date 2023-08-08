function o = get_gene_efficiencies(o)
%% o = o.get_gene_efficiencies;
% Finds Gene Efficiencies from average code of good genes as found by
% call_spots_omp_initial. Then uses these to update iompBledCodes to
% something more representative of genes encountered.
% Input
%   o: iss object
% Output
%   o.GeneEfficiency(g,r): mean of o.iompSpotColors corresponding to gene g
%       in round r equals o.GeneEfficiency(g,r) x o.iompBledCodes(g,r).
%   o.z_scoreBledCodes(g,r) = o.GeneEfficiency(g,r) x o.iompBledCodes(g,r)
%       and then normalised so norm of o.z_scoreBledCodes(g,:) equals 1. 
%%
nCodes = length(o.CharCodes);
o.GeneEfficiency = ones(nCodes,o.nRounds);
%Need median Intensity too to avoid spots with just one massive square.
GoodQualOK = ...
    o.iompSpotIntensity>o.GeneEfficiencyIntensityThresh &...
    o.iompSpotIntensity2>o.GeneEfficiencyIntensityThresh &...
    o.iompCoef>o.GeneEfficiencyCoefThresh & ...
    o.iompSpotScore>o.GeneEfficiencyScoreThresh & ...
    o.iompResOverBackground>o.GeneEfficiencyResThresh;
% GoodQualOK = o.iompCoef>0.7 &...
%     o.iompSpotScore>0.1 & o.iompResOverBackground>0.5 &...
%     o.iompSpotCodeNo==o.iompSpotBestGene;
BledCodes = zeros(nCodes,o.nBP,o.nRounds);
for g=1:nCodes
    Use = find(GoodQualOK&o.iompSpotCodeNo==g);
%     %Only keep genes not near other genes
%     tree = KDTreeSearcher(o.iompSpotGlobalYX(GoodQualOK&o.iompSpotCodeNo~=g,:));
%     [~,Dist2NearestOtherGene] = ...
%         tree.knnsearch(o.iompSpotGlobalYX(GoodQualOK&o.iompSpotCodeNo==g,:),'K',1);
%     Use = Use(Dist2NearestOtherGene>=5);
    MeanCode = mean(double(o.iompSpotColors(Use,:,:)));
    %MeanCode = median(double(o.iompSpotColors(Use,:,:))); %Median gave slightly better results
    NormMeanCode = (MeanCode-o.z_scoreSHIFT)./o.z_scoreSCALE;
    %Weight average so channels where no intensity don't
    %contribute much
    for r=1:o.nRounds
        %W = abs(o.z_scoreBleedMatrix(:,ChannelNo)');
        W = 1;
        % ChannelNo = str2num(o.CharCodes{g}(r))+1; %% wrong if the channel order is changed!!!!!
        [~, ChannelNo] = ismember(o.CharCodes{g}(r), o.bpLabels); % changed by SB 2023-07-10
        
        if sum(Use)>=o.GeneEfficiencyMinSpots && ismember(ChannelNo,o.UseChannels)
            o.GeneEfficiency(g,r) = (W.*NormMeanCode(:,o.UseChannels,r))/...
                (W.*o.z_scoreBleedMatrix(o.UseChannels,ChannelNo)');
        end
        BledCodes(g,:,r) = o.GeneEfficiency(g,r)*...
            o.z_scoreBleedMatrix(:,ChannelNo)';
    end
%     %Method to find eigenvectors of unbled code for each gene so get
%     %relation of intensity between rounds - Failed. 
%     %MAYBE TRY DIVIDING ALL INTENSITIES BY ROUND R INTENSITY FIRST??
%     NormSpotColors = (double(o.iompSpotColors(Use,:,:))-o.z_scoreSHIFT)./o.z_scoreSCALE;
%     UnbledSpotColors = zeros(sum(Use),o.nRounds);
%     for r=1:o.nRounds
%         ChannelNo = str2num(o.CharCodes{g}(r))+1;
%         UnbledSpotColors(:,r) = NormSpotColors(:,ChannelNo,r);
%     end
%     [NormTopEvec,s2] = eigs(double(UnbledSpotColors'*UnbledSpotColors)/sum(Use), 1);
%     NormTopEvec = NormTopEvec*sign(mean(NormTopEvec)); % make positive
%     TopEvec = NormTopEvec * sqrt(s2);
%     
%     
%     for r=1:o.nRounds
%         ChannelNo = str2num(o.CharCodes{g}(r))+1;
%         o.GeneEfficiency(g,r) = TopEvec(r)/...
%             o.z_scoreBleedMatrix(ChannelNo,ChannelNo)';
%         BledCodes(g,:,r) = o.GeneEfficiency(g,r)*...
%             o.z_scoreBleedMatrix(:,ChannelNo)';
%     end

    if o.Graphics==2 || sum(Use)<o.GeneEfficiencyMinSpots
        if sum(Use)<o.GeneEfficiencyMinSpots
            FigNo = 47280+g;
            warning(['GeneEfficiency set to 1.0 for all rounds of gene %.0f: %s ',...
                'because only %.0f spots pass thresholding. ',...
                'Need o.GeneEfficiencyMinSpots=%.0f spots to use ',...
                'mean code shown in figure %.0f'],...
                g,o.GeneNames{g},sum(Use),o.GeneEfficiencyMinSpots,FigNo);
        else
            FigNo = 47280;
        end
        CLims = [min(min(NormMeanCode(:)),min(BledCodes(g,:))),...
            max(max(NormMeanCode(:)),max(BledCodes(g,:)))];
        fig = figure(FigNo);
        fig.Position = [163,312,1150,414];
        subplot(1,3,1);
        imagesc(squeeze(NormMeanCode));
        caxis(CLims);
        set(gca, 'ytick', 1:o.nBP);
        yticklabels(o.bpLabels);
        ylabel('Channel');
        set(gca, 'xtick', 1:o.nRounds);
        xlabel('Round');
        title('Mean Bled Code');
        subplot(1,3,2);
        imagesc(squeeze(BledCodes(g,:,:)));
        caxis(CLims);
        colorbar;
        set(gca, 'ytick', 1:o.nBP);
        yticklabels(o.bpLabels);
        ylabel('Channel');
        set(gca, 'xtick', 1:o.nRounds);
        xticklabels(round(o.GeneEfficiency(g,:),2,'Significant'));
        xlabel('Gene Efficiency');
        title('BledCode From GeneEfficiency');
        subplot(1,3,3);
        imagesc(squeeze(NormMeanCode)-squeeze(BledCodes(g,:,:)));
        caxis([-0.1,0.1]);
        colormap(gca,bluewhitered);
        colorbar;
        set(gca, 'ytick', 1:o.nBP);
        yticklabels(o.bpLabels);
        ylabel('Channel');
        set(gca, 'xtick', 1:o.nRounds);
        xlabel('Round');
        title('MeanCode - GeneEfficiencyCode');
        sgtitle(sprintf('Gene %.0f: %s, %.0f Good Spots',g,o.GeneNames{g},...
            sum(Use)));
    end
end
o.z_scoreBledCodes = BledCodes(:,:);
% Set non o.UseChannels, o.UseRounds to 0 not nan.
% Both SpotColors and BledCodes are 0 in these rounds/channels so OMP won't
% try and fit them. 
o.z_scoreBledCodes(isnan(o.z_scoreBledCodes)) = 0;  
o.z_scoreBledCodes = o.z_scoreBledCodes./vecnorm(o.z_scoreBledCodes(:,:),2,2);
end

