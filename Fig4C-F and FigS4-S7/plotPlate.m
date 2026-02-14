function resist=plotPlate(dataAll,numbers,drugType,conc,colors)

    for plate = 1:size(dataAll,2)
        data=dataAll{plate};
        % convert time to hours
        t0=[0:1/12:(size(data,3)-1)/12];
    
        d1=size(data,1);
        d2=size(data,2);
        
        if plate == 1
            resist=NaN(d1,d2,1,size(dataAll,2)); % max OD, growth rate
            resist(:,:,2,:)=0;
        end
    
        figure;
        for i=1:d1
            % time drug is added
            i1=1;
            for j=1:d2
        
                v=data(i,j,:);
                v=squeeze(v(~isnan(v)));
                t=squeeze(t0(~isnan(v)));
                
                % Account for cell stuff accumulating
                if strcmp(drugType,'Carbenicillin') || strcmp(drugType,'Ceftazidime')
                    v=filter_spikes(v);
                    v = sgolayfilt(v,3,11);
                end
                nonlogV = v;
                v=log2(v);

                % OD at end of 20 hrs minus background-> for MIC
                resist(i,j,1,plate) = mean(nonlogV(end-2:end))-mean(nonlogV(4:6));
                
                % log OD at time drug added and 1/2 doubling after
                o1=v(i1);
                o2=o1+0.5;

                % index at which log OD is 1/2 doubling from drug added OD
                i2=find(v(i1:end)>o2,1);

                if(~isempty(i2)) 
            
                    i2=i2+i1-1;
        
                    if(j==1) 
                        t1=t;
                        v1=v;
                    end
        
                    o3=o2;
                    mm=max(v);
                    %o4=max(o3,mm-0.3);
                    o4 = o1+1.5;
                    i3=find(v(i2:end)>o3,1);
                    
                    if(~isempty(i3)) 
                        
                        i3=i3+i2-1;
                        i4=find(v(i3:end)>=o4,1);

                        if(~isempty(i4))
                            i4=i4+i3-1;
                        
                        % Reduce search area to 1 doubling if 1.5 doublings
                        % is too much
                        else
                            % o4 = o1+1.25;
                            % i4=find(v(i3:end)>=o4,1);
                            % if(~isempty(i4))
                            %     i4=i4+i3-1;
                            % end
                            i4 = i3+40;
                            if i4 > length(v)
                                i4 = length(v);
                            end
                        end

                    end
        
                    % growth->for IC50
                    c=polyfit(t(i3:i4),v(i3:i4)',1);
                    resist(i,j,2,plate)=c(1);
                       
                end 
       
                if(j==1) 
                    t1=t;
                    v1=v;
                end
                
                subplot(d2,d1,(j-1)*d1+i)
                
                plot(t1,v1,'--k','LineWidth',2)
                hold on
                plot(t,v,'Color',colors{i},'LineWidth',2)
                if(~isempty(i2)) 
                    plot(t(i3:i4),polyval(c,t(i3:i4)),'Color','r','LineWidth',5)
                end
                axis([0 20 -3.85 0.6])
                if j==1
                    n=numbers{i};
                    if strcmp(n,'PAO1')
                        text(3.5,3,'PAO1','FontSize',20)
                    elseif strcmp(n,'PA14')
                        text(3.5,3,'PA14','FontSize',20)
                    else
                        text(-1,3,[num2str(n)],'FontSize',20)
                    end
                end
        
                if(i==1 && j==8) 
                    set(gca,'box','off')
                    ylabel('doub.','FontSize',15,'FontWeight', 'bold')
                    xlabel('Time (hours)','FontSize',15,'FontWeight', 'bold')
                    xticks([0,20]);
                    yticks([-3,0]);
                    ax = gca;
                    ax.FontSize = 15; 
                    ax.FontWeight = 'bold';
                else
                    axis off
                end

                if i == 1
                    text(-50,-2,conc{j},'FontSize',20)
                end
        
            end
         
        end
        iptsetpref('ImshowBorder','tight');
        set(gcf, 'Position', [1 161 1550 550]);
        saveas(gcf,[drugType,num2str(plate),'.png'])
        %close all
        
    end
end