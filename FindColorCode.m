% o.CodeFile = 'C:\Users\bugeon\Documents\\codebook_73g_ctx.txt';
o.CodeFile = 'G:\Documents\C_doc\codebook_7rounds.txt';
Colors = {'AF405', 'AF488', 'Dy485XL', 'Cy3', 'TexasRED', 'Cy5','Atto425'};


fp = fopen(o.CodeFile, 'r');
tmp = textscan(fp, '%s %s', inf);
GeneName=tmp{1};
CharCode=tmp{2};
fclose(fp);
%%
T = table();
T.GeneName = GeneName;
AllCodes =[];
for i=1:length(CharCode{1})
dd = str2num(cellfun(@(x) x(i),CharCode));
AllCodes = [AllCodes , dd];
GeneColor = Colors(dd+1);
roundname = ['GeneColor_r',num2str(i-1)];
T.(roundname) = GeneColor';
end
writetable(T,fullfile('C:\Users\bugeon\Documents\73gColorRound.xlsx'))
%% find bridge probes to pool to make IN mix bridges!
IN_list = {'Calb2','Gad1','Htr3a','Kit','Lhx6','Ndnf','Pvalb','Sst','Vip','Npy'}; 
B = readtable('C:\Users\bugeon\Documents\BridgeProbeList73g.xlsx');
Tb = table2array(T);
Bn={};ListB =table();
for i=1:length(Colors)
    gg = strcmp(Tb,Colors{i});
    for j=1:length(IN_list)
        ThisG = find(strcmp(IN_list{j},GeneName));
        dd = find(gg(ThisG,:));
        for k=1:length(dd)
            Bn{i}{j,k} = [GeneName{ThisG},'_r',num2str(dd(k)-2)];
            kk = find(strcmp(Bn{i}{j,k},B.OligoName));
            for b = 1:length(kk)
                NewT = B(kk(b),:);
                NewT.Color = Colors(i);
                ListB = [ListB;NewT];
            end
        end
    end

end
ListB = unique( ListB,'rows');
ListB = sortrows(ListB,'Color','ascend');
writetable(ListB,'C:\Users\bugeon\Documents\IN_bridgeMix.xlsx')