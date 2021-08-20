function [SpotIntensity, MedianIntensity] = ...
    get_spot_intensity(o,SpotCodeNo,SpotColors,NormFactor)
%% [SpotIntensity, MedianIntensity] = o.get_spot_intensity(SpotCodeNo,SpotColors,NormFactor);
% This gives a modified spot intensity taking account of the gene it is
% assigned to.
% If don't give NormFactor:
%   For spot s, matched to gene g = SpotCodeNo(s),
%   SpotIntensity(s) =
%   mean(SpotColors(s,InCharCode_gene_g))-mean(SpotColors(s,NotInCharCode_gene_g))
%   MedianIntensity(s) = median(SpotColors(s,InCharCode_gene_g))
% If do give NormFactor:
%   SpotIntensity(s) = mean(NormSpotColors(s,InCharCode_gene_g))
%   MedianIntensity(s) = median(NormSpotColors(s,InCharCode_gene_g))
% Hence a high SpotIntensity indicates a high intensity and a good match.
% Inputs
%   o: iss object.
%   SpotCodeNo(s): gene that spot s was matched to.
%   SpotColors(s,b,r): gives intensity for spot s in channel b, round r.
%   NormFactor(1,b,r): the value that SpotColors(:,b,r) is divided by to
%       give Normalised SpotColors (Optional).

Norm = false;
if nargin>=4
    SpotColors = double(SpotColors)./NormFactor; %Normalise
    Norm = true;
end

nCodes = length(o.GeneNames);
CodeIndex = zeros(nCodes,o.nRounds);
NonCodeIndex = cell(nCodes,1);
for g=1:nCodes
    GeneChannels = str2double(regexp(cell2mat(o.CharCodes(g)),'\d','match'))+1; 
    GeneChannels(~ismember(GeneChannels,o.UseChannels)) = nan;
    CodeIndex(g,1:length(GeneChannels)) = sub2ind([o.nBP,o.nRounds],GeneChannels,1:length(GeneChannels));
    UnusedChannels = setdiff(o.UseChannels,GeneChannels);
    UnusedChannelIndex = sub2ind([o.nBP,o.nRounds],repelem(UnusedChannels,1,o.nRounds),...
        repmat(1:o.nRounds,1,length(UnusedChannels)));
    NonCodeIndex{g} = setdiff(1:o.nRounds*o.nBP,[CodeIndex(g,:),UnusedChannelIndex]);
end

nSpots = length(SpotCodeNo);
SpotIntensity = zeros(nSpots,1);
MedianIntensity = zeros(nSpots,1);

fprintf('Percentage of spot intensities found:       ');
for s=1:nSpots
    SpotCode = SpotColors(s,:);
    % NEED TO INVESTIGATE SPOT INTENSITY: Z-SCORE and NO SUBTRACTION BEST I THINK
    sCodeIndex = CodeIndex(SpotCodeNo(s),:);
    sCodeIndex = sCodeIndex(sCodeIndex>0);
    if Norm
        SpotIntensity(s) = mean(SpotCode(sCodeIndex));
    else
        SpotIntensity(s) = mean(SpotCode(sCodeIndex))-mean(SpotCode(NonCodeIndex{SpotCodeNo(s)}));
    end
    %SpotIntensity(s) = mean(SpotCode(CodeIndex(SpotCodeNo(s),:)));
    %SpotIntensity(s) = mean(SpotCode(NonCodeIndex{SpotCodeNo(s)}));
    MedianIntensity(s) = median(SpotCode(sCodeIndex));
    if mod(s,round(nSpots/100))==0
        Percent = sprintf('%.6f', round(s*100/nSpots));
        fprintf('\b\b\b\b\b%s%%',Percent(1:4));
    end
end
fprintf('\n');
end


% %This was original method used, differs slightly as when colour channel
% %appeared in more than one round, took mean of all rounds it appeared first. 
% %This method was much slower though and pretty much same results
% nSpots = length(SpotCodeNo);
% SpotIntensity = zeros(nSpots,1);
% RoundCode = 1:o.nRounds;
% for s=1:nSpots
%     SpotCode = o.cSpotColors(s,:,:);
%     numCharCode = str2double(regexp(cell2mat(o.CharCodes(SpotCodeNo(s))),'\d','match'))+1;    
%     sIntensity = zeros(7,1);
%     for b=1:o.nBP
%         UseRounds = find(numCharCode==b);
%         if ~isempty(UseRounds)
%             sIntensity(b)= mean(SpotCode(:,b,UseRounds))-mean(SpotCode(:,b,setdiff(RoundCode,UseRounds)));
%         end
%     end
%     SpotIntensity(s) = mean(sIntensity(sIntensity~=0));
% end
% 
