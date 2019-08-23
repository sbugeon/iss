function [] = plotShiftSearch(shifts,score,name)
% Plots the Score used to find the best shift
shiftsNew = shifts-min(shifts)+1;

ScoreImage = zeros(max(shiftsNew));
for i = 1:size(shiftsNew,1)
    ScoreImage(shiftsNew(i,1),shiftsNew(i,2),shiftsNew(i,3)) = score(i);
end

S.BestShift = shifts(score == max(score),:);


S.fh = figure('units','pixels','position',[500 200 800 600],'name',name,'numbertitle','off');  %Left, Bottom, Width, Height
S.Cmap = ScoreImage;
clearvars Correlation logCorrel;
S.x = min(shifts(:,2)):max(shifts(:,2));
S.y = min(shifts(:,1)):max(shifts(:,1));
S.z = min(shifts(:,3)):max(shifts(:,3));
S.Background = imagesc('XData',S.x,'YData',S.y,'CData',S.Cmap(:,:,1));

%S.Data = LocalYXZ;
%ToUse = round(S.Data(:,3)) == 1;
%S.x = S.Data(ToUse,2);  % For plotting.  
%S.y = S.Data(ToUse,1);
%S.LN = scatter(S.x,S.y,'wx');    %Format dots at this point
title(['Z Shift = ' num2str(S.z(1))]);
xlim([min(shifts(:,2)) max(shifts(:,2))]);
ylim([min(shifts(:,1)) max(shifts(:,1))]);
xlabel('X');
ylabel('Y');
%S.caxis = [-5, max(S.Cmap(:))];
%caxis(S.caxis);
if S.BestShift(3) == S.z(1)
    hold on;
    scatter(S.BestShift(2),S.BestShift(1),80,'kx');
    hold off;
end

S.sl = uicontrol('style','slide',...
                 'unit','pix',...
                 'position',[60 8 693 18],...
                 'min',1,'max',size(S.z,2),'val',1,...
                 'sliderstep',[1/(size(S.z,2)-1) 1/(size(S.z,2)-1)],...
                 'callback',{@sl_call,S});  
set( findall( S.fh, '-property', 'Units' ), 'Units', 'Normalized' )


function [] = sl_call(varargin)
% Callback for the slider.
[h,S] = varargin{[1,3]};  % calling handle and data structure.
ZPlane = round(get(h,'value'));
S.Background = imagesc('XData',S.x,'YData',S.y,'CData',S.Cmap(:,:,ZPlane));
%caxis(S.caxis);
%ToUse = round(S.Data(:,3)) == ZPlane;
%S.x = S.Data(ToUse,2);  
%S.y = S.Data(ToUse,1);
%S.LN = scatter(S.x,S.y,'wx');
if S.BestShift(3) == S.z(ZPlane)
    hold on;
    scatter(S.BestShift(2),S.BestShift(1),80,'kx');
    hold off;
end
title(['Z Shift = ' num2str(S.z(ZPlane))]);

