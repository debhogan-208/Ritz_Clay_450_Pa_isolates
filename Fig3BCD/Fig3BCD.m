%% David Ritz
% PhD candidate | Schultz lab
% Created: Dec 2025
% Description: Read and analyze mucoid data of Deb CF isolates


%% 

fn='Data on mucoidy.xlsx';
dat = table2array(readtable(fn))*100;

% Remove blank rows
dat = dat(~all(isnan(dat), 2), :);

colorsClade{1} = hex2rgb(['#994500';'#ff7100';'#ffa966']);
colorsClade{2} = hex2rgb(['#1f1782';'#877de8';'#cbc8ec']);
colorsClade{3} = hex2rgb(['#178262';'#7de8c8']);

%% Find mean and std of data

err = zeros(size(dat,1)/3,size(dat,2));
avg = zeros(size(dat,1)/3,size(dat,2));
kk=1;
for k=[1 4 7]
    for i = 1:size(dat,2)
        err(kk,i) = std(dat(k:k+2,i));
        avg(kk,i) = mean(dat(k:k+2,i));
    end
    kk=kk+1;
end

%% Plot bar graphs

for k = 1:size(avg,1)

    if k == size(avg,1)
        xLoc = [0.5 2];
        lab = {'Upper Lobe';'Lower Lobe'};
        upX = 3;
    else
        xLoc = [0.5 2 3.5];
        lab = {'Upper Lobe';'Middle Lobe';'Lower Lobe'};
        upX = 4.5;
    end

    figure(k);
    hold on
    for kk = 1:size(xLoc,2)
        b=bar(xLoc(kk)',avg(k,kk)','EdgeColor',[0 0 0],'LineWidth',3,'BarWidth',1,'FaceColor',colorsClade{k}(kk,:));
    end
    
    % built-in scatter jitter doesnt work for this, so create my own
    % jitter in x-direction
    startRow = [1 4 7];
    for i = 1:size(xLoc,2)
        for r = startRow(k):startRow(k)+2
            
            a=-0.25;
            aa=0.25;
            rr = a + (aa-a).*rand();

            sc = scatter(rr+xLoc(i),dat(r,i),100,'filled','MarkerEdgeColor',[0 0 0],'LineWidth',3,'MarkerFaceColor',colorsClade{k}(i,:));

        end
        % Put error bars on top of everything. Only need to do this during 1 loop
        if i == size(xLoc,2)
            idx = ~isnan(avg(k,:));
            errorbar(xLoc,avg(k,idx),-err(k,idx),err(k,idx),'Color',[0 0 0],'LineWidth',3,'LineStyle','none','CapSize',20);
        end
    end
    
    % easier to first (1) use numerical label, then (2) swap to non-numerical
    xlim([-0.5 upX])
    xticks(xLoc)
    
    set(gca,'XTickLabel',lab)
    ylim([0 60])
    yticks([0 20 40 60])

    ylabel('Percent mucoid (%)')
    hh=gca;
    hh.XAxis.TickLength = [0 0];


    set(findall(gcf,'-property','FontSize'),'FontSize',28)
    
    % Stats between clades
    if k==1
        newY = [45 55 20];
    elseif k == 2
        newY = [35 45 10];
    else
        newY = [10];
    end
        
    % check if normal assumption holds before ANOVA.. user will be alerted if not
    % h is binary yes/no if reject null of normal. p is p-value
    
    for kk = 1:size(xLoc,2)
        vec = dat(startRow(k):startRow(k)+2,kk);
        [h(kk,1),p(kk,1),~] = swtest(vec,0.05);
    end
    
    if any(h == 1)
        disp('Normal distribution assumption is rejected for at least one dataset within clades.')
    end
    
    if k < 3
        comparisons = [...
            1 2;...
            1 3;...
            2 3];
        cladeStat = dat(startRow(k):startRow(k)+2,:);
        cladeStat = cladeStat(:);
        cats = [1;1;1;2;2;2;3;3;3];
        [~,~,stats] = anova1(cladeStat,cats,'off');
        c = multcompare(stats,'Display','off','CriticalValueType','tukey-kramer');
        szComp = size(comparisons,1);
        pPairwise = zeros(szComp,1);
    
        % Only check user-inputted pairwise comparisons for significance
        for w = 1:szComp
            % find [1 1] vector -- where c row and comparisons row are
            % equal
            idx = find(sum((c(:,1:2) == comparisons(w,:)),2) == 2);
            pPairwise(w) = c(idx,6);
        end
        comp = [xLoc(1),xLoc(2);xLoc(1),xLoc(3);xLoc(2),xLoc(3)];
    else
        % Normal distribution rejected for timepoint 3
        [pPairwise,~] = ranksum(dat(startRow(k):startRow(k)+2,1),dat(startRow(k):startRow(k)+2,2));
        comp = [xLoc(1),xLoc(2)];
    end

    % STYLISHLY add stats to current plot using ritzStar
    if any(pPairwise<=0.05)
        ritzStar(comp,pPairwise,newY,0);
    end
            
    % high quality output of open figure
    iptsetpref('ImshowBorder','tight');
    set(gcf, 'Position', [1 49 1550/3 700]);
    saveas(gcf,['timepoint_',num2str(k),'.png'])
end
