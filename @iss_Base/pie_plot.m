function o = pie_plot(o,Boundaries)
% plot pie chart for each cell, showing probability of it belonging to all
% classes


nC = size(o.CellYX,1);
nK = size(o.pCellClass,2);

% find classes to collapse
CollapseMe = zeros(nK,1);
Colors = zeros(nK,3);
% DisplayName = o.ClassNames;
DisplayName = {};
for i=1:size(o.ClassCollapse,1)
    ClassList = o.ClassCollapse{i,1};
    for j=1:length(ClassList)
        MyClasses = strmatch(ClassList{j}, o.ClassNames);
        if length(MyClasses)==0; continue; end
        CollapseMe(MyClasses)=i;
        Colors(MyClasses,:) = repmat(o.ClassCollapse{i,3},length(MyClasses),1);
        DisplayName(MyClasses) = o.ClassCollapse(i,2);
    end
end

nColorWheel = sum(CollapseMe==0);

Colors0 = hsv(ceil(nColorWheel*1.2));
Colors(~CollapseMe,:) = Colors0(1:nColorWheel,:); % last is zero

figure(43908765)
clf; 
set(gcf, 'Color', 'k');
set(gca, 'color', 'k');
hold on

% Idx = CellMap>0;
% B = bwconncomp(Idx);
for c=1:nC
    hold on
    fprintf(['\n Cell Nb', num2str(c)]);

    pMy = o.pCellClass(c,:);
    [~ , WorthShowing] = max(pMy);
    if WorthShowing>0.9 && ~isempty(Boundaries{c})
        if sum(Colors(WorthShowing(1),:))>0
            patch(Boundaries{c}(:,1),Boundaries{c}(:,2),Colors(WorthShowing(1),:),'EdgeColor', 'none');
        end
    end
    
end
yMax = max(o.CellYX(:,1));
xMax = max(o.CellYX(:,2));
yMin = min(o.CellYX(:,1));
xMin = min(o.CellYX(:,2));

ClassShown = true(length(DisplayName),1);
ClassDisplayNameShown = DisplayName(ClassShown);
[uDisplayNames, idx] = unique(ClassDisplayNameShown, 'stable');
nShown = length(uDisplayNames);
for k=1:nShown
    h = text(xMax*1.1 - xMin*.5, yMin + k*(yMax-yMin)/nShown, uDisplayNames{k});
    set(h, 'color', Colors((idx(k)),:));
    h.FontSize = 35;
end
daspect([1 1 1])

