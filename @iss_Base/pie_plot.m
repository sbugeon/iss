function o = pie_plot(o,Boundaries,Alpha)
% plot pie chart for each cell, showing probability of it belonging to all
% classes

if nargin<3
    Alpha = 1;
end
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


set(gcf, 'Color', 'k');
set(gca, 'color', 'k');
hold on

% Idx = CellMap>0;
% B = bwconncomp(Idx);
for c=1:nC
    hold on
    % fprintf(['\n Cell Nb', num2str(c)]);
    
    pMy = o.pCellClass(c,:);
    
    % max prob for subclasses
    Group1 = o.ClassCollapse(:,2);
    Psubclass=[];
    for i = 1 : length(Group1)
        ClassList = o.ClassCollapse{i,1};
        MyClasses = [];
        for j=1:length(ClassList)
            MyClasses(j) = find(strcmp(ClassList{j}, o.ClassNames));
        end
        Psubclass(i) = sum(pMy(MyClasses));
    end
    
    [~ , WorthShowing] = max(pMy);
    if ~isempty(Boundaries{c}) %&& max(Psubclass)>0.8
        if sum(Colors(WorthShowing(1),:))>0
            if  pMy(WorthShowing(1))>0.6 || max(Psubclass)>0.8
                patch(Boundaries{c}(:,1),Boundaries{c}(:,2),Colors(WorthShowing(1),:),'EdgeColor', 'none','FaceAlpha',Alpha);
            else % cell not above threshold are in grey
                patch(Boundaries{c}(:,1),Boundaries{c}(:,2),[0.2 0.2 0.2],'EdgeColor', 'none','FaceAlpha',Alpha);
            end
        end
    end
    
end
yMax = max(o.CellYX(:,1));
xMax = max(o.CellYX(:,2));
yMin = min(o.CellYX(:,1));
xMin = min(o.CellYX(:,2));

ClassShown = true(length(DisplayName),1);
ClassDisplayNameShown = DisplayName(ClassShown);

ClassDisplayNameShown(end+1) = {'Not-assigned'};
Colors= [Colors;[0.2 0.2 0.2]];
[uDisplayNames, idx] = unique(ClassDisplayNameShown, 'stable');
nShown = length(uDisplayNames);
if Alpha ==1
    for k=1:nShown
        h = text(xMax*1.1 - xMin*.5, yMin + k*(yMax-yMin)/nShown, uDisplayNames{k});
        set(h, 'color', Colors((idx(k)),:));
        h.FontSize = 35;
    end
end
daspect([1 1 1])

