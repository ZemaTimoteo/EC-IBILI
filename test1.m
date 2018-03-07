% clear all
% clc
% load('data_10a80_SNR8.mat')

% ------- Allowed to add data  ------


load('data_100a200_SNR8.mat')
dat = data;


pointsVect  = [10 20 30 40 60 80 100 160 200]; % [3000]; % [10 20 30 40 60 80 100 160 200]; % [40 80]; % [40:5:80]; % 
blocksVect  = [ 2  4  8 13 15 17  25  30  40]; % [ 20]; % [ 2  4  8 13 15 17  25  30  40]; % [4 25]; % 40; %

%%
% point=1;
points = 1;
n_cond= 3; cond_names = {'simulate05','simulate250','simulateOriginal'};

while point<=size(pointsVect,2) % itterates on the number of points
    block=1;
    while block<=size(blocksVect,2)  % itterates on the number of blocks
        for cond = 1:n_cond
            datas.results.sensi{point,block}.(cond_names{cond}) = dat.results.sensi{points,block}.(cond_names{cond});
            datas.results.speci{point,block}.(cond_names{cond}) = dat.results.speci{points,block}.(cond_names{cond});
            datas.results.Fone{point,block}.(cond_names{cond})  = dat.results.Fone{points,block}.(cond_names{cond});
            datas.results.bACC{point,block}.(cond_names{cond})  = dat.results.ACC{points,block}.(cond_names{cond});
            datas.results.ACC{point,block}.(cond_names{cond})   = dat.results.bACC{points,block}.(cond_names{cond});
            for blo=1:blocksVect(block)
                datas.cons.(cond_names{cond})(block,blo,point) = dat.cons.(cond_names{cond})(block,blo,points);
            end    
        end
        block = block + 1; 
    end
    point = point + 1;
    points = points + 1;

end

data = datas;