function [] = plotCorrelation(Correlation,shift,name)
% Plots the correlation found from ImRegFft3D
% shift is the YXZ optimal shift found.

LogCorrel = zeros(size(Correlation));
TileSzY = size(Correlation,1)/2;
TileSzX = size(Correlation,2)/2;
TileSzZ = size(Correlation,3)/2;
%Recentre plot so 0 shift is in the middle
LogCorrel(1:TileSzY,1:TileSzX,1:TileSzZ) = Correlation(TileSzY+1:TileSzY*2,TileSzX+1:TileSzX*2,TileSzZ+1:TileSzZ*2);
LogCorrel(1:TileSzY,1:TileSzX,TileSzZ+1:TileSzZ*2) = Correlation(TileSzY+1:TileSzY*2,TileSzX+1:TileSzX*2,1:TileSzZ);
LogCorrel(1:TileSzY,TileSzX+1:TileSzX*2,1:TileSzZ) = Correlation(TileSzY+1:TileSzY*2,1:TileSzX,TileSzZ+1:TileSzZ*2);
LogCorrel(TileSzY+1:TileSzY*2,1:TileSzX,1:TileSzZ) = Correlation(1:TileSzY,TileSzX+1:TileSzX*2,TileSzZ+1:TileSzZ*2);
LogCorrel(1:TileSzY,TileSzX+1:TileSzX*2,TileSzZ+1:TileSzZ*2) = Correlation(TileSzY+1:TileSzY*2,1:TileSzX,1:TileSzZ);
LogCorrel(TileSzY+1:TileSzY*2,1:TileSzX,TileSzZ+1:TileSzZ*2) = Correlation(1:TileSzY,TileSzX+1:TileSzX*2,1:TileSzZ);
LogCorrel(TileSzY+1:TileSzY*2,TileSzX+1:TileSzX*2,1:TileSzZ) = Correlation(1:TileSzY,1:TileSzX,TileSzZ+1:TileSzZ*2);
LogCorrel(TileSzY+1:TileSzY*2,TileSzX+1:TileSzX*2,TileSzZ+1:TileSzZ*2) = Correlation(1:TileSzY,1:TileSzX,1:TileSzZ);

LogCorrel(LogCorrel<0) = 10^-200; %remove negative numbers before taking log
LogCorrel = log10(LogCorrel);

S.fh = figure('units','pixels','position',[500 200 800 600],'name',name,'numbertitle','off');  %Left, Bottom, Width, Height
S.Cmap = LogCorrel;
clearvars Correlation logCorrel;
S.x = -TileSzX:TileSzX-1;
S.y = -TileSzY:TileSzY-1;
S.z = -TileSzZ:TileSzZ-1;
S.Background = imagesc('XData',S.x,'YData',S.y,'CData',S.Cmap(:,:,1));

%S.Data = LocalYXZ;
%ToUse = round(S.Data(:,3)) == 1;
%S.x = S.Data(ToUse,2);  % For plotting.  
%S.y = S.Data(ToUse,1);
%S.LN = scatter(S.x,S.y,'wx');    %Format dots at this point
title(['Z Shift = ' num2str(S.z(1))]);
xlim([-(TileSzX-1) TileSzX]);
ylim([-(TileSzY-1) TileSzY]);
xlabel('X');
ylabel('Y');
S.caxis = [-5, max(S.Cmap(:))];
caxis(S.caxis);
S.shift = shift;
if S.shift(3) == S.z(1)
    hold on;
    scatter(S.shift(2),S.shift(1),80,'kx');
    hold off;
end
S.sl = uicontrol('style','slide',...
                 'unit','pix',...
                 'position',[60 8 693 18],...
                 'min',1,'max',TileSzZ*2,'val',1,...
                 'sliderstep',[1/(TileSzZ*2-1) 1/(TileSzZ*2-1)],...
                 'callback',{@sl_call,S});  
set( findall( S.fh, '-property', 'Units' ), 'Units', 'Normalized' )


function [] = sl_call(varargin)
% Callback for the slider.
[h,S] = varargin{[1,3]};  % calling handle and data structure.
ZPlane = round(get(h,'value'));
S.Background = imagesc('XData',S.y,'YData',S.x,'CData',S.Cmap(:,:,ZPlane));
caxis(S.caxis);
%ToUse = round(S.Data(:,3)) == ZPlane;
%S.x = S.Data(ToUse,2);  
%S.y = S.Data(ToUse,1);
%S.LN = scatter(S.x,S.y,'wx');
if S.shift(3) == S.z(ZPlane)
    hold on;
    scatter(S.shift(2),S.shift(1),80,'kx');
    hold off;
end
title(['Z Shift = ' num2str(S.z(ZPlane))]);

