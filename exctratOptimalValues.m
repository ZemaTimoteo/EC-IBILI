load('C:\Users\admin\Documents\Tiago\IBILI\Results\Testes\MatrixResilts_Test\5var_LowCorr_SSGC_timeRev_SNR200_ChannelGap10-200.mat')

%%
point = 1;
n_cond= 3; cond_names = {'simulate05','simulate250','simulateOriginal'};
maxFone = 0;
i=0;
indSensi_point = NaN(64,1);
indSpeci_block = NaN(64,1);
n_cond= 3; cond_names = {'simulate05','simulate250','simulateOriginal'};

while point<=size(pointsVect,2) % itterates on the number of points
    block=1;
    
    while block<=size(blocksVect,2)  % itterates on the number of blocks
        for cond = 1:n_cond            
            if data.results.Fone{point,block}.(cond_names{cond})>maxFone
                maxFone = data.results.sensi{point,block}.(cond_names{cond});
                indMaxFone_point = point;
                indMaxFone_block = block;
            end
            if data.results.speci{point,block}.(cond_names{cond})>0.55 && data.results.sensi{point,block}.(cond_names{cond})>0.55
                i=i+1;
                indSensi_point(i) = point;
                indSpeci_block(i) = block;
            end
            
        end
        
        block = block + 1;
    end
    
    point = point + 1;
end
