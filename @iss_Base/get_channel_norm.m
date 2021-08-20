function p = get_channel_norm(o)
%% p = o.get_channel_norm;
% p(1,b,r) is the normalisation for channel b, round r. 
% This finds the normalisations for each color channel such that if
% bNormSpotColor = SpotColor(:,b,:)./p(:,b,:), then the probability that
% bNormSpotColor>o.ChannelNormValue1 is approximately o.ChannelNormThresh1
% for each channel, b.
% Probabilities are based on all pixels in image so can penalise a channel
% if there are few spots.
% To ensure not sharp increase in probability at smaller intensity, also
% require that probability bNormSpotColor>o.ChannelNormValue2 is less than
% o.ChannelNormThresh2 where o.ChannelNormValue2<o.ChannelNormValue1.
% To ensure not lots of high intensity outliers, also
% require that probability bNormSpotColor>o.ChannelNormValue3 is less than
% o.ChannelNormThresh3 where o.ChannelNormValue3>o.ChannelNormValue1.

p = zeros(1,o.nBP,o.nRounds);
for b = 1:o.nBP
    bHistCounts = sum(o.HistCounts(:,b,:),3);
    nPixels = sum(bHistCounts);
    bCumSum = cumsum(bHistCounts);
    pb = o.HistValues(find(bCumSum>(1-o.ChannelNormThresh1)*nPixels, 1 )+1)/o.ChannelNormValue1;
    Prob2 = sum(bHistCounts(o.HistValues >= o.ChannelNormValue2*pb))/nPixels;
    Prob3 = sum(bHistCounts(o.HistValues >= o.ChannelNormValue3*pb))/nPixels;
    if Prob2>o.ChannelNormThresh2
        pbOld = pb;
        pb = o.HistValues(find(bCumSum>(1-o.ChannelNormThresh2)*nPixels, 1 )+1)/o.ChannelNormValue2;
        warning(['Channel %.0f: Probability of %.4f that intensity is ',...
            'larger than o.ChannelNormValue2*o.BledCodesPercentile(:,%.0f,:) ',...
            'exceeds threshold of %.4f. Norm Changed from %.0f to %.0f.'],...
            b,Prob2,b,o.ChannelNormThresh2,pbOld,pb);
    end
    if Prob3>o.ChannelNormThresh3
        pbOld = pb;
        pb = o.HistValues(find(bCumSum>(1-o.ChannelNormThresh3)*nPixels, 1 )+1)/o.ChannelNormValue3;
        warning(['Channel %.0f: Probability of %.6f that intensity is ',...
            'larger than o.ChannelNormValue3*o.BledCodesPercentile(:,%.0f,:) ',...
            'exceeds threshold of %.6f Norm Changed from %.0f to %.0f.'],...
            b,Prob3,b,o.ChannelNormThresh3,pbOld,pb);
    end
    p(:,b,:) = pb;
end

% %Old method
% p = zeros(1,o.nBP,o.nRounds);
% for b = 1:o.nBP
%     bSpotColors = o.dpSpotColors(:,b,:);
%     p(:,b,:) = prctile(bSpotColors(:), o.SpotNormPrctile);
% end

end

