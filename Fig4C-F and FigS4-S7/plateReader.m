function [MIC,IC50] = plateReader(fn,drugType,conc,colors2,colorsHold)
    
    %% Read data
    shts=sheetnames(fn);
    dat = {};
        
    for s = 1:size(shts,1)
        tab = table2array(readtable(fn,'Sheet',shts(s),'ReadVariableNames',true,'NumHeaderLines',0,'ReadRowNames',true));
    
        % remove all rows with any NaN
        tab( any( isnan( tab ), 2 ), : ) = [];
    
        tab(:,1) = [];
        tab(:,1) = [];

        dat{s} = tab;
    end
        
    %% change dimensions for existing software pipeline
    
    data = {};
    
    for plate = 1:size(dat,2)
        data{plate} = zeros(12,8,size(dat{plate},1));
        for strain=1:12
            k=strain;
            for drug=1:8
                data{plate}(strain,drug,:) = dat{plate}(:,k);
                k=k+12;
            end
        end
        % take first 20 hrs
        data{plate} = data{plate}(:,:,1:240);

        % Delete blank column. 1st column is blank for last rep of carb
        if strcmp(drugType,'Carbenicillin') && plate == 3
            data{plate}(1,:,:) = [];
        else
            data{plate}(12,:,:) = [];
        end
        
        % Match clades in mexTree
        dataHold = data{plate}(1:3,:,:);
        data{plate}(1:3,:,:) = data{plate}(4:6,:,:);
        data{plate}(4:6,:,:) = dataHold;
    end
    
    
    clades={'','Clade 1','','','Clade 2','','','Clade 3','','PA14','PAO1'};
    numbers={'C1S1','C1S2','C1S3','C2S1','C2S2','C2S3','C3S1','C3S2','C3S3','PA14','PAO1'};
    for j = 1:size(conc,2)
        concLab{j} = {[num2str(conc(j)),' ',char(181),'g/mL']};
    end

    %% plot plate reader data growth curves

    % multiple assignments in Matlab is annoying..
    j=1;
    for k=1:11
        if k == 4
            j=2;
        elseif k == 7
            j=3;
        elseif k == 10
            j=4;
        elseif k == 11
            j=5;
        end
        colors{k} = colorsHold{j};
    end

    resist = plotPlate(data,clades,drugType,concLab,colors);

    %% plot plate reader growth and MIC trendlines
    [err,avg] = errAvg(resist);
    [MIC,IC50] = plotTrend(resist,conc,drugType,avg,err,numbers,colors2);
    
    
end