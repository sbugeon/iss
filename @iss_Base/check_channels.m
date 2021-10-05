function o = check_channels(o, auto_adjust)
%% o = o.check_channels(auto_adjust);
% Checks to see if any of the color channels are weak.
% Weak color channels are those with norm factor less than 
% o.ChannelNormQualThresh
% I.e. normalisation factor is so small that background pixels can be
% boosted to seem like spots when compared to other channels with larger
% norm factors.
% If auto_adjust is true, o.UseChannels will change automatically.
% If false, an error will be raised if any channels are weak.
% default: false

if nargin < 2 || isempty(auto_adjust)
    auto_adjust = false;
end

if isempty(o.UseChannels)
    o.UseChannels = 1:o.nBP;
end

p = o.get_channel_norm;
channels_to_change = find(min(p,[],3)<o.ChannelNormQualThresh);
UseChannelsAdjust = setdiff(1:o.nBP, channels_to_change);
    
if size(channels_to_change,2) > 0  && ~isempty(setdiff(o.UseChannels,...
        UseChannelsAdjust))
    str_bad_channels = [sprintf('%d,', channels_to_change(1:end-1)), ...
        sprintf('%d', channels_to_change(end))];
    if size(UseChannelsAdjust,2) == 0
        str_UseChannelsAdjust = ('[]');
    else
        str_UseChannelsAdjust = [sprintf('%d,', UseChannelsAdjust(1:end-1)), ...
            sprintf('%d', UseChannelsAdjust(end))];
    end
    if auto_adjust
        warning(['Channels ', str_bad_channels ,...
            ' are weak. Setting o.UseChannels = ',...
            str_UseChannelsAdjust, '.']);
        o.UseChannels = UseChannelsAdjust;
        if size(o.UseChannels,2) < 4
            error('Only '+ string(size(o.UseChannels,2))+ ' Channels Remain')
        end
    else
        error(['Channels ', str_bad_channels ,...
            ' are weak. Recommend setting o.UseChannels = ',...
            str_UseChannelsAdjust, '.']);
    end
end

