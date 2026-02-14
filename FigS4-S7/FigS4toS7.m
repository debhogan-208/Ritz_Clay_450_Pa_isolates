%% David Ritz
% PhD candidate | Schultz lab
% Created: Summer 2025
% Description: Read and analyze plate reader data of Deb CF isolates grown 
% in gradient of drug. Calculate IC50 and MIC.

%% MIC and IC50 of carb, ceft, gent, tob
colorsStrain = {'#8AA3FE','#012ED0','#011E87',...
        '#00FF00','#00A300','#004A00',...
        '#D52EFF','#A700D1','#5E0075',...
        '#000000','#808080'};

colorsClade = {"#1171BE",...
        "#3BAA32"...
        "#8516D1",...
        '#000000','#808080'};

drugs = {'Carbenicillin','Ceftazidime','Gentamicin','Tobramycin'};


fn='carb_08.04.2025.xlsx';
conc = [0,10,25,50,75,100,200,400];
[MIC{1},IC50{1}] = plateReader(fn,drugs{1},conc,colorsStrain,colorsClade);

fn='ceft_08.19.2025.xlsx';
conc = [0,5,10,15,20,30,40,50];
[MIC{2},IC50{2}] = plateReader(fn,drugs{2},conc,colorsStrain,colorsClade);

fn='gent_08.04.2025.xlsx';
conc = [0,8,250,500,750,1000,1500,2000];
[MIC{3},IC50{3}] = plateReader(fn,drugs{3},conc,colorsStrain,colorsClade);

fn='tob_08.18.2025.xlsx';
conc = [0,8,50,100,200,250,300,350];
[MIC{4},IC50{4}] = plateReader(fn,drugs{4},conc,colorsStrain,colorsClade);

close all

%% Find mean and std of biological replicates

% [strain,drug,resistance] -> (1) resistance = MIC; (2) resistance = IC50
err = zeros(size(MIC{1},1),size(MIC,2),2);
avg = zeros(size(MIC{1},1),size(MIC,2),2);
for k=1:size(MIC,2) 
    for i = 1:size(MIC{k},1)
        err(i,k,1) = std(MIC{k}(i,:));
        avg(i,k,1) = mean(MIC{k}(i,:));
    
        err(i,k,2) = std(IC50{k}(i,:));
        avg(i,k,2) = mean(IC50{k}(i,:));
    end
end

% save data MIC IC50 err avg colorsClade