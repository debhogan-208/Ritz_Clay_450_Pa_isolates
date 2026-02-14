function [err,avg] = errAvg(resist)
    % MIC, IC50
    err = zeros(size(resist,1),size(resist,2),2);
    avg = zeros(size(resist,1),size(resist,2),2);
    
    for i = 1:size(resist,1)
        for j = 1:size(resist,2)
            err(i,j,1) = std(resist(i,j,1,:));
            avg(i,j,1) = mean(resist(i,j,1,:));
    
            err(i,j,2) = std(resist(i,j,2,:));
            avg(i,j,2) = mean(resist(i,j,2,:));
        end  
    end
end