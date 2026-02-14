%% David Ritz
% PhD candidate | Schultz lab
% Created: Fall 2023 - Summer 2025
% Description: Create collapsed tree from efflux pump non-synonymous 
% polymorphisms from Deb's excel sheet of P. aeruginosa isolates

% Read excel sheet from Deb

baseFolder = pwd;
inputFolder = [pwd,'\input_mex'];
if ~isfolder(inputFolder)
    mkdir(inputFolder);
end
cd(baseFolder)
addpath(genpath(baseFolder));

% input tree for ITOL export
Tree = phytreeread('hogan_tree_newick.txt');
% Get leaves of tree
leaves = get(Tree,'LeafNames');

% Assign clades
clade = ones(size(leaves));

cutoffs = {'T1_RUL_B2','T1_RUL_C4','T2_RUL_E10','T1_RML_H2','break'};
ccount = 1;
for c = 1:size(clade,1)
    if strcmp(leaves{c},cutoffs{ccount})
        ccount = ccount+1;
        clade(c:end) = ccount;
    end
end

% read polymorphisms of strains
dat = readtable('deb_onesheet.txt','ReadVariableNames',true,'VariableNamingRule','preserve');

syn = "synonymous_variant";

genes = ['PA2018';'PA2019';'PA5471';... % mexXY, armZ
    'PA0425';'PA0426';... % mexAB
    'PA2527';... % muxB
    'PA2493';'PA2491';... % mexE, mexS
    'PA3676';... % mexK
    'PA4598';'PA4597'... % mexD, oprJ 
    ];

% build matrix with all possible polymorphisms of interest by leaf.
% Size of nCols may not be size of pLeaf columns, because there may be multiple
% polymorphisms in same gene
nRows = size(leaves,1);
nCols = size(genes,1);

pLeaf = cell(nRows,nCols);
polyCount = zeros(nCols,1);
lobe = zeros(nRows,3);
trunk = zeros(19,1);

% red, blue, yellow, green, purple, orange: 2 at most
colors = ["#FF0000";"#0096FF";"#FDDA0D";"#008000";"#FF028D";'#CD7F32';...
    "#FF6347";"#5CB3FF";"#FFFF00";"#52D017";"#FF77FF";'#FFA500'];


% cd into where output metadata is stored
cd(inputFolder)

%% Step through Excel sheet

masterCount = 0;
countCat = cell(size(genes,1)+7,2);
% add this for counting category. Because, {} + 1 = {}
countCat(1:6,2) = {0};
colNames = dat(1,:).Properties.VariableNames;
poly=1;


for g = 1:size(genes,1)

    % get gene and user-defined color
    gene = genes(g,:);
    
    % mexXY, armZ
    if g < 4
        num=1;
        color = colors(num+6,:);
        id = 'mexXY';

    % mexAB
    elseif g < 6
        num=2;
        color = colors(num+6,:);
        id = 'mexAB';

    % muxB
    elseif g < 7
        num=3;
        color = colors(num+6,:);
        id = 'muxB';

    % mexE, mexS
    elseif g < 9
        num=4;
        color = colors(num+6,:);
        id = 'mexE';

    % mexK
    elseif g < 10
        num=5;
        color = colors(num+6,:);
        id = 'mexK';
    
    % mexD
    else
        num=6;
        color = colors(num+6,:);
        id = 'mexD';
    end
    
    for row = 1:size(dat,1)
        % if there's a polymorphism in a gene of interest and it's non-synonymous
        if strcmp(gene,cell2mat(table2array(dat(row,2)))) && ~strcmp(syn,cell2mat(table2array(dat(row,11))))
            
            polyCount(poly,1) = num;
            if strcmp('stop_gained',cell2mat(table2array(dat(row,11))))
                trunk(poly,1) = 1;
            end

            % Count number in each category
            countCat(num,1) = {id};
            countCat(num,2) = {cell2mat(countCat(num,2))+1};

            % Find what the nt changed to
            ntFin = checkNt(dat,row,colNames);
            
            % get polymorphism identification
            mut = cell2mat(table2array(dat(row,13)));
            mutSave = mut(3:end);
            mutSave = strrep(mutSave,'*','--');

            % Shorten 3 letter amino acid to 1 letter
            if strcmp(mut(end),'*')
                mut = [aminolookup(mut(3:5)),mut(6:end)];
            else
                mut = [aminolookup(mut(3:5)),mut(6:end-3),aminolookup(mut(end-2:end))];
            end
            

            % replace invalid txt file name character with hyphen (if
            % tracking nt, but I've swapped to AA)
            % mut = strrep(mut,'>','-');
            
            
            % If gene has name on Excel sheet, use that for file name. If not, use hypothetical gene ID
            if isempty(cell2mat(table2array(dat(row,3))))
                geneid = gene;
            else
                geneid = cell2mat(table2array(dat(row,3)));
            end


            if strcmp(geneid,'PA1797')
                geneid = 'mipB';
            end
                
            % define output text id for metadata
            outText = strcat(id,'_',geneid,'_',mutSave);

            % Initiate output metadata .txt file with 2 columns
            out = ["DATASET_COLORSTRIP" " ";
            "SEPARATOR" "TAB";
            "DATASET_LABEL"	outText;" " " ";
            "DATA" " "];
            
            % Initialize row step for out.txt
            count=0;
            
            % Record which strains have the polymorphism
            for col = 15:size(colNames,2)

                ntStrain = upper(cell2mat(table2array(dat(row,col))));
                if strcmp(ntFin,ntStrain)
                   
                    % strain name
                    strain = cell2mat(colNames(col));

                    % record entry
                    out(6+count,1) = strain;
                    out(6+count,2) = color;

                    % step to next row
                    count = count+1;

                    % add polymorphism to strain to pLeaf matrix
                    strain_num = find((string(leaves) == strain)==1);
                    pLeaf(strain_num,poly) = num2cell(1);

                end
            end

            % Count for next polymorphism for pLeaf
            poly=poly+1;
            % Output color-coded metadata for each mutation of interest
            writematrix(out,strcat(outText,'.txt'),'Delimiter','\t',"QuoteStrings","none")

            % output to other txt file too
            countCat(8+masterCount,1) = {strcat(geneid,'_',mut)};
            countCat(8+masterCount,2) = {count};

            % step to next row
            masterCount = masterCount+1;
        end
    end
end
writecell(countCat,'AA_count_categories.txt','Delimiter','\t',"QuoteStrings","none")
save = pLeaf;

%% Remove neighbor in pLeaf if they are the same.. keep track of lobe amount, number of strains during each collapse
% create temporary variable and replace empty entries with 0's
temp = pLeaf;
temp = ~cellfun('isempty', temp);
headerRow = countCat(8:end,1)';
startSize = size(pLeaf,1);

for s = 1:size(pLeaf,1)
    strain = string(leaves(s));
    % different lung locations 
    if extract(strain,5) == "U"
        locCat = 1;
    elseif extract(strain,5) == "M"
        locCat = 2;
    elseif extract(strain,5) == "L"
        locCat = 3;
    else
        locCat = 0;
    end
    
    % record lung location of poly, for each strain (U, M, L, respectively)
    if locCat > 0
        lobe(s,locCat) = lobe(s,locCat)+1;
    end
end

k=1;
for s = size(pLeaf,1)-1:-1:1
    if temp(s+1,:) == temp(s,:)
        lobe(s,:) = lobe(s,:) + lobe(s+1,:);
        delLeaves(k+1,:) = leaves(s+1,:);
        lobe(s+1,:) = [];
        pLeaf(s+1,:) = [];
        leaves(s+1,:) = [];
        temp(s+1,:) = [];
        clade(s+1) = [];
        
        k=k+1;
    end
end

% remove all blank columns (mutation not in any strain)
keep = any(~cellfun('isempty',pLeaf),1);
pLeaf = pLeaf(:,keep);
headerRow = headerRow(:,keep);
polyCount = polyCount(keep);

%% Prune and output phylogenetic tree

% Prune leaves I want to remove
ind = getbyname(Tree,delLeaves,'Exact', true);

smallTree = prune(Tree,any(ind,2));
cd(baseFolder)

% pruned tree orders leaves differently AFTER outputting from Matlab using 
% phytreeWrite due to Matlab phylogeny optimization.
% -> Can see differences in outputted .tree file
% -> Matlab reorder() function is not sufficient to solve this
% 
% Create phytreeOrdered function to do this

% input mexPruned.nwk into ITOL
phytreeOrdered(smallTree, leaves, 'mexPruned.nwk')

% Remove rooted origin strain
lobe(1,:) = [];
pLeaf(1,:) = [];
leaves(1,:) = [];
clade(1) = [];

%reTree = phytreeread('mexPruned.nwk');

%% Fill out color scheme for each collapsed strain. One box is one polymorphism
colorCols = cell(size(pLeaf,1),size(pLeaf,2));

for s = 1:size(pLeaf,1)
    for m = 1:size(pLeaf,2)
        cat = polyCount(m);
        if sum(cell2mat(pLeaf(s,m))) > 0
            colorCols(s,m) = cellstr(colors(cat));
        else
            colorCols(s,m) = cellstr('#FFFFFF');
        end
    end
end


%% Output subplot graph with each row being (1) triangle with num of strains
% next sections: (2) matrix of polymorphism colors (3) pie graph percent of
% lobe location (4) clade number & star

% Make triangle's size correspond to number of strains
% Should be on same x (2nd number), 2 units between each y (1st number)
len = 1;
start = 0;
stepY = 1.3;
figure;

% Rescale by dividing by max number of strains in clade
totStrains = sum(lobe,2);
facTot = rescale(totStrains,0.1,1);

subplot(1,2,1)
for k = 0:size(pLeaf,1)-1
    w = len;
    h = len*facTot(k+1)/2;
    yLoc = start-stepY*k-0.8;
    tri = makeTriangle(w,h,yLoc);

    hold on
    axis off
end
axLimY = [-stepY*(size(pLeaf,1)-1)-stepY start+len*2+stepY];
axLimX = [0 15];

ylim(axLimY)
xlim(axLimX)

aaa=subplot(1,2,1);
aaa.Position(2)=aaa.Position(2)-0.05; % y;

hold off


% Make matrix of polymorphism colors
% Each square represents one category.
% red, blue, yellow, green, purple, brown
cats = ["{\itmexXY}","{\itmexAB}","{\itmuxABC}","{\itmexEF}","{\itmexJK}","{\itmexCD}"];

starColor = ["#1171BE";"#3BAA32";"#8516D1"];
starCount = 1;

len = 1;
startX = 0;
startY = 0;
stepX = 1.5;
stepY = 1.3;
subplot(1,2,2)
for k = 0:size(pLeaf,1)-1
    for s = 1:size(colorCols,2)
        yLoc = startY-stepY*k-0.8;
        xLoc = startX+stepX*(s-1);

        % add color legend
        if k == 0 && s == 1
            for j = 1:size(cats,2)
                legX = xLoc+40;
                legY = yLoc+2-(j-1)*5;
                text(legX,legY-0.5,cats(:,j),'Color',[0 0 0],'FontSize', 13,'VerticalAlignment', 'bottom','HorizontalAlignment','right')
                pos = [legX+0.5 legY-1 len+2 len+2];
                squ = rectangle('Position',pos,'FaceColor',colors(j));
            end

            % legend entry for lobe pie graph
            legY = legY -5.5;
            text(legX-0.5,legY,'lobe','Color',[0 0 0],'FontSize', 13,'VerticalAlignment', 'bottom','HorizontalAlignment','right')
        end

        % bottom left x, bottom left y, width, height
        pos = [xLoc-1 yLoc len len];

        % Mark gene color
        squ = rectangle('Position',pos,'FaceColor',string(colorCols(k+1,s)));
        
        % Add allele info above first row
        if k == 0
            allele = string(headerRow(s));
            allele = allele{1};
            textOut = "{\it"+allele(1:4)+"}"+" "+allele(6:end);
            if trunk(s) == 1
                alleleText = text(xLoc-1+stepX-0.2,yLoc+1.2,textOut,'Color',[0 0 0],'FontSize', 12,'VerticalAlignment', 'bottom','HorizontalAlignment','left','FontWeight','bold');
            else
                alleleText = text(xLoc-1+stepX-0.2,yLoc+1.2,textOut,'Color',[0 0 0],'FontSize', 12,'VerticalAlignment', 'bottom','HorizontalAlignment','left');
            end            
            set(alleleText,'Rotation',90);
        end
    end
    axLimY = [-stepY*(size(pLeaf,1)-1)-stepY startY+len*2+stepY];
    axLimX = axLimY+35;
    axLimX(2) = axLimX(2)+8;
    %axLimY(2) = axLimY(2)+5.5;

    ylim(axLimY)
    xlim(axLimX)
    
    hold on
    axis off
end
hold off

% move subplot 2 closer to subplot 1
aa=subplot(1,2,2);
aa.Position(1)=0.152; % x;
%aa.Position(2)=0.1005; % y
aa.Position(2)=0.05; % y

% Make pie graph for each lobe and mark clade #
% Pie charts are hard to work with in Matlab.. can't use subplot nor
% nexttile to add to an existing figure. But can place each specific pie in
% the existing figure using coordinates

% Lobe order = Upper, middle, lower
len = 0.5;
start = 0;
stepY = 1.3;

hold on
colorsClade = {"#1171BE",...
        "#3BAA32"...
        "#8516D1",...
        "#FF0000",...
        "#d3af37"};
for q = 1:size(lobe,1)

    % 3 equal slices for pie graph
    shapePie = [1,1,1];
    pie_chart = pie(shapePie, repmat({''},size(shapePie)));

    % Shade of each slice of pie chart, proportional to number of strains found in each lobe
    X = lobe(q,:);
    shadeScale = rescale(X,0,1);
    
    % Tweek rescale so total sum of all slices is 1. 
    % E.g. 1 strain showing up in 2 different lobes should have 2 pie slices be gray, not black
    shadeScale = shadeScale/sum(shadeScale);
    
    shade = zeros(3,3);
    for s = 1:3
        shade(s,:) = [1-shadeScale(s),1-shadeScale(s),1-shadeScale(s)];
    end


    scale=0.55;
    j=1;
    for k = 1:length(pie_chart) 
        if pie_chart(k).Type=="patch" 
            XData = pie_chart(k).XData;
            YData = pie_chart(k).YData;
            
            % Make pie chart of certain size set by scale, and shift it
            % to the right of colored boxes using 6.5 and down using
            % stepY and loop counter
            set(pie_chart(k),'XData', XData*scale+28);
            set(pie_chart(k),'YData', YData*scale - (q-1)*stepY+0.5-0.8);
            set(pie_chart(k),'FaceColor', shade(j,:));
            set(pie_chart(k),'EdgeAlpha',1) 
            j=j+1;
            
            
            % Add pie chart to figure legend
            if q == size(lobe,1)
                move_x = 1.8;
                pie_chart2 = pie(shapePie, repmat({''},size(shapePie)));
                for g = 1:length(pie_chart) 
                    if pie_chart2(g).Type=="patch"
                        set(pie_chart2(g),'XData', XData*1.75+legX+move_x);
                        set(pie_chart2(g),'YData', YData*1.75 + legY+0.5);
                        set(pie_chart2(g),'FaceColor', [1 1 1]);
                        set(pie_chart2(g),'EdgeAlpha',1);
                        set(pie_chart2(g),'LineWidth',1);
                    end
                    text(legX+move_x-0.85,legY+0.5,'U','Color',[0 0 0],'FontSize', 10,'VerticalAlignment', 'bottom','HorizontalAlignment','center')
                    text(legX+move_x+0.85,legY+0.5,'L','Color',[0 0 0],'FontSize', 10,'VerticalAlignment', 'bottom','HorizontalAlignment','center')
                    text(legX+move_x,legY-1,'M','Color',[0 0 0],'FontSize', 10,'VerticalAlignment', 'bottom','HorizontalAlignment','center')
                end
            end
        end
    end

    % Mark star for those clades we tested resistance
    adj = 0.35;
    if q == 2 || q == 4 || q == 5
        % central x = 30.37
        plot(31,-(q-1)*stepY-adj+0.09,'Marker','p','MarkerSize',13,'MarkerFaceColor',starColor(starCount),'MarkerEdgeColor','k','LineWidth',1.25)
        starCount=starCount+1;
    end

    % Mark clade number adjacent to pie graph
    %rectangle('Position',[28.8,-(q-1)*stepY-adj-0.64,1,1.2],'FaceColor',hex2rgb(colorsClade(clade(q))),'EdgeColor','k');
    text(28.8,-(q-1)*stepY-adj,num2str(clade(q)),'FontSize',14);


end

% Ensure aspect ratio is good
ylim(axLimY)
xlim(axLimX)

axis off;

% high quality output of open figure
iptsetpref('ImshowBorder','tight');

% Slight tweek aspect ratio, because need to extend into the
% x-direction a bit farther
set(gcf, 'Position', [0,130,1540,580]);
%export_fig mexTree.tiff -m2.5 -q101;

