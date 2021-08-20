function [BackgroundEigenvectors,BackgroundEigenvalues,BackgroundMaxGeneDotProduct,...
    BackgroundMaxGeneDotProductGene] = get_background_codes(o,SpotColors)
%% [BackgroundEigenvectors,BackgroundEigenvalues,BackgroundMaxGeneDotProduct,...
%   BackgroundMaxGeneDotProductGene] = o.get_background_codes(SpotColors);
% Find vectors to best represent background
% Finds covariance matrix of background pixels. Then takes eigenvectors of
% this. 
% If o.ompBackgroundChannelStrips is true, then will just return o.nBP background
% eigenvectors which are just strips in each color channel. This is the
% default.
%
% Input: 
% o: iss_OMP object.
% SpotColors: Raw SpotColors of pixels representative of background
%
% Output:
% BackgroundEigenvectors(i,b,r): intensity in channel b, round
%   r for the ith eigenvector of the covariance matrix comprised of
%   background pixels
% BackgroundEigenvalues(i): corresponding eigenvalue.
% BackgroundMaxGeneDotProduct(i): the absolute dot product of
%   BackgroundEigenvectors(i,:) with
%   o.iompBledCodes(BackgroundMaxGeneDotProductGene(i),:).
%   This is the largest dot product of any gene with BackgroundEigenvectors(i,:)
%%
if o.ompBackgroundChannelStrips 
    %Eigenvector b only has intensity in channel b for all rounds. 
    BackgroundEigenvectors = zeros(o.nBP,o.nBP,o.nRounds);
    for b=1:o.nBP
        BackgroundEigenvectors(b,b,:) = 1;
    end
    BackgroundEigenvectors = BackgroundEigenvectors(o.UseChannels,:,:);
    BackgroundEigenvectors = BackgroundEigenvectors./vecnorm(BackgroundEigenvectors(:,:),2,2);
    BackgroundEigenvalues = ones(length(o.UseChannels),1);
else
    z_scoredSpotColors = (double(SpotColors)-o.z_scoreSHIFT)./o.z_scoreSCALE;
    z_scoredSpotColors = z_scoredSpotColors(:,o.UseChannels,o.UseRounds);
    %CovMatrix = z_scoredSpotColors(:,:)'*z_scoredSpotColors(:,:);
    CovMatrix = cov(z_scoredSpotColors(:,:));
    [Eig,BackgroundEigenvalues] = get_eig(CovMatrix);
    BackgroundEigenvectors = zeros(length(BackgroundEigenvalues),o.nBP,o.nRounds);
    BackgroundEigenvectors(:,o.UseChannels,o.UseRounds) = ...
        reshape(Eig',[length(BackgroundEigenvalues),length(o.UseChannels),length(o.UseRounds)]);
    BackgroundEigenvalues = BackgroundEigenvalues/vecnorm(BackgroundEigenvalues);
end
%Find out how much background eigenvectors resemble genes
nCodes = length(o.CharCodes);
BledGeneCodes = reshape(o.iompBledCodes,[nCodes,o.nBP,o.nRounds]);
BledGeneCodes = BledGeneCodes(:,o.UseChannels,o.UseRounds);
BledGeneCodes(isnan(BledGeneCodes)) = 0;
NormBledGeneCodes = BledGeneCodes(:,:)./vecnorm(BledGeneCodes(:,:),2,2);
NormBledGeneCodes(isnan(NormBledGeneCodes)) = 0;
BackgroundEigenvectorsCrop = ...
    reshape(BackgroundEigenvectors,[size(BackgroundEigenvalues,1),o.nBP,o.nRounds]);
BackgroundEigenvectorsCrop = BackgroundEigenvectors(:,o.UseChannels,o.UseRounds);
BackgroundEigenvectorsCrop(isnan(BackgroundEigenvectorsCrop))=0;
BackgroundEigenvectorsCrop = ...
    BackgroundEigenvectorsCrop(:,:)./vecnorm(BackgroundEigenvectorsCrop(:,:),2,2);
BackgroundEigenvectorsCrop(isnan(BackgroundEigenvectorsCrop))=0;
[BackgroundMaxGeneDotProduct,BackgroundMaxGeneDotProductGene] = ...
    get_bled_code_max_dot_product(BackgroundEigenvectorsCrop,NormBledGeneCodes,nCodes);

if o.Graphics
    figure(180432); clf
    for i=1:min(9,length(BackgroundEigenvalues))
        subplot(3,3,i);
        if o.nRounds>1
            imagesc(squeeze(BackgroundEigenvectors(i,:,:)));
        else
            imagesc(squeeze(BackgroundEigenvectors(i,:,:))');
        end
        title(sprintf('Eigenvalue = %.3f\n%s DotProduct = %.3f',...
            BackgroundEigenvalues(i),o.GeneNames{BackgroundMaxGeneDotProductGene(i)},...
            BackgroundMaxGeneDotProduct(i)));
        set(gca, 'ytick', 1:o.nBP);
        set(gca, 'yTickLabel', o.bpLabels(1:o.nBP));
        set(gca, 'xtick', 1:o.nRounds);
        set(gca, 'XTickLabel', 1:o.nRounds);
        caxis([-1,1]);
        colormap(gca,bluewhitered);
        if i==4
            xlabel('Round')
            ylabel('Channel');
        end
    end
    drawnow;
    sgtitle('Background Eigenvectors');
end
     
end

function [MaxDotProduct,MaxGeneNumber,DotProductStd] = get_bled_code_max_dot_product(SpotColors,GeneBledCodes,nCodes)
%For each spot, s, finds dot product of SpotColors(s,:) with each gene.
%Then returns max absolute dot product and corresponding gene. 
%Both SpotColors and GeneBledCodes should have been normalised. 
GeneDotProduct = zeros(size(SpotColors,1),nCodes);
for GeneNo=1:nCodes
    GeneDotProduct(:,GeneNo) = SpotColors(:,:)*GeneBledCodes(GeneNo,:)';
end
[MaxDotProduct,MaxGeneNumber] = max(abs(GeneDotProduct),[],2);
DotProductStd = std(abs(GeneDotProduct),[],2);
end   
