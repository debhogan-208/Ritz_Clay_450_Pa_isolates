%% David Ritz
% PhD candidate | Schultz lab
% Created: Summer 2025
% Description: Read and analyze plate reader data of Deb CF isolates grown 
% in gradient of drug. Calculate IC50 and MIC.

load data.mat

%% Plot bar graphs of IC50, MIC
% Not enough info to interpolate to find WT resistance values 
% because concentrations chosen to fit CF strains

% resistance bar graphs of each strain within each clade
xLoc = [0.5 1 1.5 2.5 3 3.5 4.5 5 5.5];
for k = 1:size(avg,3)
    checkSig=0;
    for n = 1:size(avg,2)
    
        figure(n+k*10);
        hold on
        toPlot = squeeze(avg(1:9,n,k));
        z=1;
        for c = [1 4 7]
            b=bar(xLoc(c:c+2)',toPlot(c:c+2)','EdgeColor',[0 0 0],'LineWidth',3,'BarWidth',1,'FaceColor',hex2rgb(colorsClade(z)));
            z = z+1;
        end
        
        m=1;
        for i = 1:size(avg,1)-2
            if i>1 && mod(i-1,3) == 0
                m=m+1;
            end
            for r = 1:size(MIC{1},2)
                
                % built-in scatter jitter doesnt work for this, so create my own
                % jitter in x-direction
                a=-0.125;
                aa=0.125;
                rr = a + (aa-a).*rand();

                if k == 1
                    resistType = MIC{n}(i,r);
                else
                    resistType = IC50{n}(i,r);
                end
                
                
                sc = scatter(rr+xLoc(i),resistType,100,'filled','MarkerEdgeColor',[0 0 0],'LineWidth',3,'MarkerFaceColor',hex2rgb(colorsClade(m)));

            end
            % Put error bars on top of everything. Only need to do this during 1 loop
            if i == size(avg,1)-2
                errorbar(xLoc,avg(1:9,n,k),-err(1:9,n,k),err(1:9,n,k),'Color',[0 0 0],'LineWidth',2,'LineStyle','none','CapSize',10);
            end
        end
        
        % easier to first (1) use numerical label, then (2) swap to non-numerical
        xlim([0 6])
        xticks([1 2 3 4 5 6])
        
        set(gca,'XTickLabel',{'Clade 1';'';'Clade 2';'';'Clade 3';''})
        if k == 1
            if n == 1
                ylim([0 700])
                yticks([0 150 300 450 600 750])
            elseif n == 2
                ylim([0 90])
                yticks([0 25 50 75 100])
            elseif n == 4
                ylim([0 700])
                yticks([0 150 300 450 600 750])
            else
                ylim([0 3400])
                yticks([0 750 1500 2250 3000])
            end
        else
            if n == 1
                ylim([0 500])
                yticks([0 150 300 450])
            elseif n == 2
                ylim([0 80])
                yticks([0 25 50 75 100])
            elseif n == 4
                ylim([0 400])
                yticks([0 150 300 450])
            else
                ylim([0 2500])
                yticks([0 1250 2500])
            end
        end
        ylabel([drugs{n}, ' (\mug/mL)'])
        hh=gca;
        hh.XAxis.TickLength = [0 0];
    

        set(findall(gcf,'-property','FontSize'),'FontSize',28)
        
        % Stats within clades
                
        % check if normal assumption holds before ANOVA.. user will be alerted if not
        % h is binary yes/no if reject null of normal. p is p-value
        
        comparisons = [...
            1 2;...
            1 3;...
            2 3];
        
        if k == 1
            cladeStat = MIC{n}(1:9,:);
        else
            cladeStat = IC50{n}(1:9,:);
        end
        
        for kk = 1:size(cladeStat,1)
            vec = cladeStat(kk,:);
            vec = vec(:);
    
            [h(kk,1),p(kk,1),~] = swtest(vec,0.05);
        end
        
        if any(h == 1)
            disp('Normal distribution assumption is rejected for at least one dataset within clades.')
        end
        
        % can't assume normal, so use kruskalwallis: 1st input -> data, 2nd input -> grouping
        % assume alphas = 0.05
        cats = [1;2;3;1;2;3;1;2;3];
        for g = [1 4 7]
            cclade = cladeStat(g:g+2,:);
            cclade = cclade(:);
            [~,~,stats] = kruskalwallis(cclade,cats,'off');
            c = multcompare(stats,'Display','off','CriticalValueType','bonferroni');
            szComp = size(comparisons,1);
            pPairwise = zeros(szComp,1);
    
    
            % Only check user-inputted pairwise comparisons for significance
            for w = 1:szComp
                % find [1 1] vector -- where c row and comparisons row are
                % equal
                idx = find(sum((c(:,1:2) == comparisons(w,:)),2) == 2);
                pPairwise(w) = c(idx,6);
            end
            if any(pPairwise<=0.05)
                % STYLISHLY add stats to current plot using ritzStar
                newY = ritzStar(comparisons,pPairwise,newY,offset(k));
                checkSig=1;
            end
        end

        clear h p

        % Stats between clades
        if k==1
            newY = [350 45 1800 325];
        else
            newY = [250 45 1300 225];
        end
            
        % check if normal assumption holds before ANOVA.. user will be alerted if not
        % h is binary yes/no if reject null of normal. p is p-value
        
        comparisons = [...
            1 2;...
            1 3;...
            2 3];
        
        cladeStat = avg(1:9,n,k);
        
        for kk = 1:3:size(cladeStat,1)
            vec = cladeStat(kk:kk+2);
            [h(kk,1),p(kk,1),~] = swtest(vec,0.05);
        end
        
        if any(h == 1)
            disp('Normal distribution assumption is rejected for at least one dataset between clades.')
        end
        
        % ANOVA: 1st input -> data, 2nd input -> grouping
        % assume alphas = 0.05
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
    
        % STYLISHLY add stats to current plot using ritzStar
        if any(pPairwise<=0.05)
            ritzStar(([xLoc(2),xLoc(5);xLoc(2),xLoc(8);xLoc(5),xLoc(8)]),pPairwise,newY(n),0);
        end
                
        % high quality output of open figure
        iptsetpref('ImshowBorder','tight');
        set(gcf, 'Position', [1 49 1550/4 700]);
        if k == 1
            saveas(gcf,[drugs{n},'_mic.png'])
        else
            saveas(gcf,[drugs{n},'_ic50.png'])
        end
    end
    if checkSig == 0 && k == 1
        disp('No significance found within clades for MIC.')
    elseif checkSig == 0 && k == 2
        disp('No significance found within clades for IC50. ')
    end
end
