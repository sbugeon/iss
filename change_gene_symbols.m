function New_symbols = change_gene_symbols(MarkerSize, FontSize, MultiCol,LineWidth,MarkerType)
% ChangeGeneSymbols(MarkerSize, FontSize, nPerCol);
%
% changes gene symbols so in situ plots look nice.
% MarkerSize defaults to 6 - if 0, won't change existing sizes
%
% FontSize is font size for legend
%
% nPerCol says how many legend entries per column (0 for one-column)
%
% Kenneth D. Harris, 29/3/17
% GPL 3.0 https://www.gnu.org/licenses/gpl-3.0.en.html

if nargin<1 || isempty(MarkerSize)
    MarkerSize = 6;
    LineWidth=1;
    MarkerType = 'GeneSpots';
end
% MarkerSize = 6;

if nargin<2 || isempty(FontSize)
    FontSize = 10;
end

if nargin<3 || isempty(MultiCol)  || MultiCol ==0
    MultiCol = 69;
end

% colors
non_neuron = hsv2rgb([0 0 1]);
pc_or_in =   hsv2rgb([.4 .5 .5]);
less_active =   hsv2rgb([.3 .2 .7]);
pc =        hsv2rgb([1/3 1 1]);
pc2 =       hsv2rgb([.27 1 .7]);
in_general = hsv2rgb([2/3 1 1]);

sst =   hsv2rgb([.55 1 1]);
pvalb = hsv2rgb([.7 .8 1]);
ngf =   hsv2rgb([.85 1 1]);
cnr1 =  hsv2rgb([ 1 1 1]);
vip =   hsv2rgb([ .13 1 1]);
cxcl14= hsv2rgb([.1 1 .6]);
sstpv = [50 135 230]./255;

BadGenes = [0 0 0];

% All symbols: +o*.xsd^v<>ph - have them in that order to not get confused
New_symbols = {...
    
% 'Ddit4l',    BadGenes, '.'; ...

'Snca',     in_general, '+'; ...
'Kcnip1',    in_general, 'o'; ...
'Nxph1',    in_general, '*'; ...
'Cplx2',    in_general, '.'; ...
'Qpct',    in_general, 'x'; ...
'Lhx6',     in_general, 's'; ...
'Zcchc12',    in_general, 'd'; ...
'Col25a1',  in_general, '^'; ...
'Pnoc',     in_general, '>'; ...
'Rab3b',    in_general, '<'; ...
'Gad1',     in_general, 'p'; ...
'Slc6a1',   in_general, 'h'; ...
%      '',    in_general, 'v'; ...

% sst complet
'Th',       sst, '+'; ...
'Crhbp',    sst, 'o'; ...
'Sst',      sst, '*'; ...
'Npy',      sst, '.'; ...
'Synpr',    sst, 'x'; ...
'Chodl',    sst, 's';...
'Cort',     sst, 'd'; ...
'Reln',     sst, '^'; ...
'Serpini1', sst, '<'; ...
'Satb1',    sst, '>'; ...
'Grin3a',   sst, 'p'; ...
'Cdh9',   sst, 'h'; ...
'Spp1',   sst, 'v'; ...

%
'Lgals1',   sstpv, 'o'; ...

% pv complet
'Tac1',     pvalb, 'o'; ...
'Pvalb',    pvalb, '*'; ...
'Thsd7a',   pvalb, 'd'; ...
'Cox6a2',   pvalb, 'v'; ...
'Chrm2',    pvalb, 'p'; ...
'Akr1c18',    pvalb, '<'; ...
'Gabrg1',    pvalb, '>'; ...
'Gpr83',    pvalb, '.'; ...
'Lypd6',    pvalb, 'x'; ...
'Nek7',    pvalb, '+'; ...
'Pdlim3',    pvalb, '^'; ...
'Prss23',    pvalb, 's'; ...
'Tpbg',    pvalb, 'h'; ...


'Id2',      ngf, '+'; ...
'Hapln1',   ngf, 'o'; ...
'Gabrd',    ngf, '*'; ...
'Cryab',    ngf, 'x'; ...
'Kit',      ngf, 's'; ...
'Ndnf',     ngf, 'd'; ...
'Nos1',     ngf, '^'; ...
'Lamp5',    ngf, '>'; ...
'Fam19a1',    ngf, '.'; ...
'Tafa1',    ngf, '.'; ...
'Mia',    ngf, 'h'; ...
'Fos',    ngf, '<'; ...
'Arc',    ngf, 'v'; ...


'Cadps2',   cxcl14, 'o'; ...
'Cxcl14',   cxcl14, '*'; ...
'Ntng1',    cxcl14, 's'; ...
'Cpne5',    cxcl14, 'd'; ...
'Rgs12',    cxcl14, 'h'; ...

% cnr1 complet
'Sncg',     cnr1, 'o'; ...
'Cnr1',     cnr1, '*'; ...
'Cck',      cnr1, '.'; ...
'Trp53i11', cnr1, 'x'; ...
'Sema3c',   cnr1, 's'; ...
'Yjefn3',   cnr1, 'v'; ...
'Rgs10',    cnr1, '>'; ...
'Nov',      cnr1, '<'; ...
'Kctd12',   cnr1, 'p'; ...
'Slc17a8',  cnr1, 'h'; ...
'Cpne2',  cnr1, '+'; ...
'Tnfaip8l3',  cnr1, '^'; ...


% vip 1 spot left
'Tac2',     vip, '+'; ...
'Npy2r',    vip, 'o'; ...
'Calb2',    vip, '*'; ...
'Htr3a',    vip, '.'; ...
'Penk',     vip, 's';...
'Pthlh',    vip, '^'; ...
'Vip',      vip, 'v'; ...
'Crh',      vip, '>'; ...
'Car4',      vip, '<'; ...
'Mybpc1',      vip, 'h'; ...
'Zbtb20',      vip, 'd'; ...
'Crispld2',      vip, 'x'; ...
%'othergene',      vip, 'p'; ...
% All symbols: +o*.xsd^v<>ph - have them in that order to not get confused

'Gda',      pc_or_in, '+'; ...
'Bcl11b',   pc_or_in, 'o'; ...
'Rgs4',     pc_or_in, '*'; ...
'Prkca',    pc_or_in, 'd'; ...
'Cdh13',    pc_or_in, '^'; ...
'Pde1a',    pc_or_in, '<'; ...
'Amigo2',    pc_or_in, '>'; ...
'Ccn2',    pc_or_in, '.'; 
'Ctgf',    pc_or_in, '.'; ...
'Chrna6',    pc_or_in, 'v'; ...
'Cpne7',    pc_or_in, 'x'; ...
'Crym',    pc_or_in, 's'; ...
 'Ddit4l',    pc_or_in, 'h'; ...
'Dkkl1',    pc_or_in, 'h'; ...

'Fezf2',    less_active, '+';...
'Gnb4',    less_active, 'o';...
'Calb1',    less_active, '*';...
'Gpr88',    less_active, '.';...
'Grik1',    less_active, 'x';...
'Hs3st2',    less_active, 's';...
'Hsd11b1',    less_active, 'd';...
'Lypd1',    less_active, '^';...
'Myl4',    less_active, 'v';...
'Nptx2',    less_active, '<';...
'Parm1',    less_active, '>';...
'Plpp4',    less_active, 'p';'Ppapdc1a',    less_active, 'p';...
'Qrfpr',    less_active, 'h';...

'Plcxd2',   pc, 'v'; ...
'Nrn1',     pc, '*'; ...
'Pcp4',     pc, '.'; ...
'Enpp2',    pc, 'x';...
'Rasgrf2',  pc, 's'; ...
'Wfs1',     pc, 'd'; ...

'Neurod6',  pc2, '+'; ...
'Nr4a2',    pc2, 'o'; ...
'Rprm',    pc2, '*'; ...
'Kcnk2',    pc2, '.'; ...
'Scgn',    pc2, 'x'; ...
'Sema5a',    pc2, 's'; ...
'Tnfaip6',    pc2, 'd'; ...
'Tnmd',    pc2, '<'; ...
'Tshz2',    pc2, '>'; ...

'Plp1',     non_neuron, 'x'; ...
'Aldoc',    non_neuron, 'o'; ...
'Serpine2',  non_neuron, '.'; ...
};
New_symbols = flipud(New_symbols);
if nargout<1
    % delete any existing legend
    fc = get(gcf, 'Children');
    for i=1:length(fc)
        if strcmp(get(fc(i), 'UserData'), 'key')
            delete(fc(i));
        end
    end
    
    MainAxes = gca;
    
    n =  size(New_symbols,1);
    
    gc = get(MainAxes, 'children');
    MyChildren = [];
    for i=1:length(gc)
        if (strcmp(gc(i).Type, 'line') || strcmp(gc(i).Type, 'scatter')) ...
                && ~isempty(gc(i).DisplayName)
            MyChildren = [MyChildren; i];
        end
    end
    DisplayNames = {gc(MyChildren).DisplayName};
    % get first word of display name as gene
    GeneNames = cell(size(DisplayNames));
    for i=1:length(DisplayNames)
        GeneNames{i} = strtok(DisplayNames{i});
    end
    
    clear h s;
    j=1;
    Present = [];
    for i=1:n
        MyGeneName = New_symbols{i,1};
        l = find(strcmp(MyGeneName, GeneNames));
        if ~isempty(l)
            h(j) = gc(MyChildren(l));
            if strcmp(h(j).Type, 'line')
                set(h(j), 'Color', New_symbols{i,2});
            elseif strcmp(h(j).Type, 'scatter')
                set(h(j), 'CData', New_symbols{i,2});
            end
            if strcmp(MarkerType,'GeneSymbols')
                set(h(j), 'Marker', New_symbols{i,3},'LineWidth',LineWidth);
            elseif strcmp(MarkerType,'Letters')
                 set(h(j), 'Marker', New_symbols{i,3},'LineWidth',LineWidth);
            end
            
            if MarkerSize>0
                if strcmp(gc(l).Type, 'line')
                    set(h(j), 'MarkerSize', MarkerSize);
                elseif strcmp(gc(l).Type, 'scatter')
                    set(h(j), 'SizeData', MarkerSize);
                end
            end
            Present(j) = i;
            j=j+1;
        end
    end
    
    other_h = setdiff(gc(MyChildren), h);
    other_symbols = {other_h.DisplayName};
    
    all_h = [h(:); other_h(:)];
    all_sym = {New_symbols{Present,1}, other_symbols{:}};
    
    % lh = legend([h(:); other_h(:)], ...
    %     {New_symbols{s,1}, other_symbols{:}}, ...
    %     'color', 'k', 'textcolor', 'w', 'fontsize', FontSize);
    % set(lh, 'color', 'k');
    %
    % return;
    
    if MultiCol==0
        lh = legend(all_h, all_sym, 'color', 'k', 'textcolor', 'w', 'fontsize', FontSize);
        set(lh, 'color', 'k');
    else
        ah = axes('Position', [.75 .13 .15 .8]);
        set(ah, 'color', 'k'); cla; hold on; box off
        set(ah, 'UserData', 'key');
        for j=1:length(Present)
            i = Present(j);
            plot(ceil(j/MultiCol)+.1, mod(j-1,MultiCol), New_symbols{i,3}, 'Color', New_symbols{i,2});
            text(ceil(j/MultiCol)+.3, mod(j-1,MultiCol), New_symbols{i,1}, 'color', 'w', 'fontsize', FontSize);
        end
        if ~isempty(setdiff(GeneNames,New_symbols(:,1)))
            j=j+1;
            plot(ceil(j/MultiCol)+.1, mod(j-1,MultiCol), '.', 'Color', hsv2rgb([0,0,0.5]));
            text(ceil(j/MultiCol)+.3, mod(j-1,MultiCol), 'Non Neuron', 'color', 'w', 'fontsize', FontSize);
        end
        ylim([-1 MultiCol]);
        set(ah, 'xtick', []);
        set(ah, 'ytick', []);
        set(ah, 'ydir', 'reverse');
    end
    %     for c=1:nCols
    %         rr=((c-1)*50 + 1):min(c*50, length(all_h));
    %         if c==1
    %             ah(c) = gca;
    %             lh(c) = legend(all_h(rr), all_sym(rr), 'color', 'k', 'textcolor', 'w', 'fontsize', FontSize, 'location', 'east');
    %             set(lh(c), 'color', 'k');
    %             pos(c,:) = get(lh(c), 'position');
    %         else
    %             ah(c) = axes('position',get(gca,'position'), 'visible','off');
    %             lh(c) = legend(ah(c), all_h(rr), all_sym(rr), 'color', 'k', 'textcolor', 'w', 'fontsize', FontSize, 'location', 'east');
    %             set(lh(c), 'position', pos(c-1,:) + [1.1 0 0 0]*pos(c-1,3));
    %             uistack(lh(c), 'top');
    %         end
    %     end
    %     axes(ah(1));
    
    %    error('multicolumn not done yet!');
    % end
    %     for i=1:nCols
    %         first = 1+(i-1)*nCols;
    %         last = min(i*nCols,length(all_h));
    %         lh = legend(all_h(first:last), ...
    %             all_sym{first:last});%, ...
    %             %'color', 'k', 'textcolor', 'w', 'fontsize', FontSize);
    %         set(lh, 'color', 'k');
    %     end
    set(gcf, 'color', 'k');
    set(gcf, 'InvertHardcopy', 'off');
    
    axes(MainAxes)
%     uistack(ah, 'top');
end

end