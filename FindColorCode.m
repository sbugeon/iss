o.CodeFile = 'C:\Users\bugeon\Documents\\codebook_73g_ctx.txt';
Colors = {'AF405', 'AF488', 'Dy485XL', 'Cy3', 'TexasRED', 'Cy5','Atto425'};


fp = fopen(o.CodeFile, 'r');
tmp = textscan(fp, '%s %s', inf);
GeneName=tmp{1};
CharCode=tmp{2};
fclose(fp);
%%
T = table();
T.GeneName = GeneName;

for i=1:length(CharCode{1})
dd = str2num(cellfun(@(x) x(i),CharCode));
GeneColor = Colors(dd+1);
roundname = ['GeneColor_r',num2str(i-1)];
T.(roundname) = GeneColor';
end
