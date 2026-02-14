function [MIC,IC50] = plotTrend(resist,conc,drugType,avg,err,numbers,colors)

    MIC = zeros(11,3);
    IC50 = zeros(11,3);
    pOrder = [1,4,7,10,12];

    for k=1:4
        figure
        
        hold on
    
        % Plot trendline (MIC)
        for i = pOrder(k):pOrder(k+1)-1
            for j = 1:size(resist,4)
                f1=fittype( 'linearinterp' );
                drug_range=linspace(conc(1),conc(end),1000);
                pd = fit( conc' , resist(i,:,1,j)' , f1 ,'Normalize','off');
        
                vx = feval( pd, drug_range );

                % MIC threshold = ~90% reduction in final OD
                t = 0.1;
        
                ix = find( vx < t , 1);
                if ~(isempty(ix))
                    MIC(i,j) = drug_range(ix);
                else
                    MIC(i,j) = conc(end);
                end
    
                % plot line, error bar, scatterpoint

                % Uncomment to plot each individual replicate
                % p(i) = plot(drug_range,vx,'Color',colors{i},'LineWidth',6);
                % s(i) = scatter(conc,resist(i,:,1,j),200,'filled','MarkerFaceColor',colors{i},'MarkerEdgeColor','black','LineWidth',3);


                % Plot average of each
                plot(conc,avg(i,:,1),'Color','k','LineWidth',8);
                p(i) = plot(conc,avg(i,:,1),'Color',colors{i},'LineWidth',6);
                e(i) = errorbar(conc,avg(i,:,1),err(i,:,1),'LineStyle','none','Color','black','LineWidth', 2,'CapSize',14);
                s(i) = scatter(conc,avg(i,:,1),400,'filled','MarkerFaceColor',colors{i},'MarkerEdgeColor','black','LineWidth',6);

                % Plot MIC threshold, which is around 0.1 for all
                thresh = plot([-100 3000],[0.1 0.1],'--k','LineWidth',3);
            end
        end

        axis([conc(1)-conc(2)*0.5 conc(end) -0.05 1.5])
        xlabel([drugType,' (\mug/mL)'],'FontSize',20)
        ylabel('OD','FontSize',20)
        yticks([0 0.5 1 1.5])
        l=legend(p(pOrder(k):pOrder(k+1)-1),numbers{pOrder(k):pOrder(k+1)-1},'Location','northeast');
        l.FontSize=20;
        set(findall(gcf,'-property','FontSize'),'FontSize',28)
        hold off

        iptsetpref('ImshowBorder','tight');
        set(gcf, 'Position', [1 161 1550/2 600]);
        saveas(gcf,[drugType,'_mic_clade',num2str(k),'.png'])
        %close all
        

        clear p e s

        % IC50
        figure
        
        hold on
    
        % Plot trendline (IC50)
        for i = pOrder(k):pOrder(k+1)-1
            for j = 1:size(resist,4)
                f1=fittype( 'linearinterp' );
                drug_range=linspace(conc(1),conc(end),1000);
                pd = fit( conc' , resist(i,:,2,j)' , f1 ,'Normalize','off');
        
                vx = feval( pd, drug_range );

                % IC50 threshold = half growth
                t = resist(i,1,2,j)*0.5;
        
                ix = find( vx < t , 1);
                if ~(isempty(ix))
                    IC50(i,j) = drug_range(ix);
                else
                    IC50(i,j) = conc(end);
                end
    
                % plot line, error bar, scatterpoint

                % Uncomment to plot each individual replicate
                % p(i) = plot(drug_range,vx,'Color',colors{i},'LineWidth',6);
                % s(i) = scatter(conc,resist(i,:,2,j),200,'filled','MarkerFaceColor',colors{i},'MarkerEdgeColor','black','LineWidth',3);


                % Plot average of each
                plot(conc,avg(i,:,2),'Color','k','LineWidth',8);
                p(i) = plot(conc,avg(i,:,2),'Color',colors{i},'LineWidth',6);
                e(i) = errorbar(conc,avg(i,:,2),err(i,:,2),'LineStyle','none','Color','black','LineWidth', 2,'CapSize',14);
                s(i) = scatter(conc,avg(i,:,2),400,'filled','MarkerFaceColor',colors{i},'MarkerEdgeColor','black','LineWidth',6);

                % Plot IC50 threshold, which is around 0.12 for CF, 0.4 for
                % WT
                if k < 4
                    plot([-100 3000],[0.2 0.2],'--k','LineWidth',3);
                else
                    plot([-100 3000],[0.4 0.4],'--k','LineWidth',3);
                end
            end
        end

        axis([conc(1)-conc(2)*0.5 conc(end) -0.05 1.2])
        xlabel([drugType,' (\mug/mL)'],'FontSize',20)
        ylabel({'Growth';'(doub./hour)'},'FontSize',20)
        yticks([0 0.5 1])
        l=legend(p(pOrder(k):pOrder(k+1)-1),numbers{pOrder(k):pOrder(k+1)-1},'Location','northeast');
        l.FontSize=20;
        set(findall(gcf,'-property','FontSize'),'FontSize',28)
        hold off

        iptsetpref('ImshowBorder','tight');
        set(gcf, 'Position', [1 161 1550/2 600])
        saveas(gcf,[drugType,'_ic50_clade',num2str(k),'.png'])
        %close all

    end
    
end