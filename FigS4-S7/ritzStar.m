function newY = ritzStar(comparisons,pPairwise,upperY,offset)
    yInc1 = upperY/10;
    yInc2 = yInc1/2.5;
    loop=1;
    for q = 1:size(comparisons,1)
        if pPairwise(q) <= 0.05
            % Add line for each pairwise comparison
            newY = upperY+yInc1*loop*2.8-upperY/40;
    
            % x-values slightly off
            line([comparisons(q,1) comparisons(q,2)]+offset,[newY newY],'Color','black','LineWidth',2)
            halfPt = sum(comparisons(q,:))/2+offset;
    
            % Draw significance marking for pairwise comparison
            if pPairwise(q) <= 0.00001
                text(halfPt,newY+yInc2-upperY/40,'*****','FontSize',50,'FontWeight','bold','HorizontalAlignment','center');
            elseif pPairwise(q) <= 0.0001
                text(halfPt,newY+yInc2-upperY/40,'****','FontSize',50,'FontWeight','bold','HorizontalAlignment','center');
            elseif pPairwise(q) <= 0.001
                text(halfPt,newY+yInc2-upperY/40,'***','FontSize',50,'FontWeight','bold','HorizontalAlignment','center');
            elseif pPairwise(q) <= 0.01
                text(halfPt,newY+yInc2-upperY/40,'**','FontSize',50,'FontWeight','bold','HorizontalAlignment','center');
            elseif pPairwise(q) <= 0.05
                text(halfPt,newY+yInc2-upperY/40,'*','FontSize',50,'FontWeight','bold','HorizontalAlignment','center');
            else
                %text(halfPt,newY+yInc2*1.5,'n.s.','FontSize',30,'FontWeight','bold','HorizontalAlignment','center');
            end
            loop=loop+1;
        end
    end
    
end
