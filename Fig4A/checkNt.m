function ntSend = checkNt(dat,row,colNames)

    % extract nt polymorphism info
    mutString = reverse(cell2mat(table2array(dat(row,12))));
    ntStart = mutString(3);
    ntEnd = mutString(1);

    
    % If polymorphism is ambiguous, don't know which strand polymorphism call is on. In this
    % case, go with simple majority to determine ancestor nt. Safe
    % assumption for everything in mex/alg trees-> verified transversions
    % by checking strand direction on pseudomonas.com
    if (strcmp(ntStart,'A') || strcmp(ntStart,'T')) && (strcmp(ntEnd,'A') || strcmp(ntEnd,'T'))
        majority = 1;
    elseif (strcmp(ntStart,'G') || strcmp(ntStart,'C')) && (strcmp(ntEnd,'G') || strcmp(ntEnd,'C'))
        majority = 1;
    else
        majority = 0;
    end
    
    if majority == 0
        if strcmp(ntEnd,'A')
            col2 = 7;
            ntEnd2 = 'T';
        elseif strcmp(ntEnd,'T')
            col2 = 8;
            ntEnd2 = 'A';
        elseif strcmp(ntEnd,'C')
            col2 = 9;
            ntEnd2 = 'G';
        elseif strcmp(ntEnd,'G')
            col2 = 10;
            ntEnd2 = 'C';
        end
        
        
        % Check which strand polymorphism is on, so then we can check which nt to
        % search for in strains
        if table2array(dat(row,col2)) > 0
            ntSend = ntEnd;
        else
            ntSend = ntEnd2;
        end
    else
        check = 0;
        for c = 7:10
            if table2array(dat(row,c)) > 0 && check == 0
                col1 = c;
                check = check+1;
            elseif table2array(dat(row,c)) > 0 && check > 0
                col2 = c;
            end
        end
    
        % send back one with less entries
        if table2array(dat(row,col1)) > table2array(dat(row,col2))
            ntSend = cell2mat(colNames(col2));
        else
            ntSend = cell2mat(colNames(col1));
        end
    end

end