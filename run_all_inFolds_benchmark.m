% Image alignment + linear optimization
%
% Wu Shiqian. 11 Sep 2010

%% specify the directory that contains your range of differently exposed
clear all; close all; clc
%pathname0 = 'D:\LDR Images\Our data_large movement';
pathname0 = 'D:\LDR Images\Still Cam Still Scene\Benchmark_synthesized';
thumbs_file = fullfile(pathname0,'thumbs.db');
delete(thumbs_file)
dir_struct = dir(pathname0);
%%% (1)name, (2)date, (3)bytes and (4) isdir
[sorted_names,sorted_index]=sortrows({dir_struct.name}');
a=[dir_struct.isdir];
a(1:2)=[];             %%% Delete the dot dir
sorted_names(1:2)=[]; 
[nFold,m]=size(sorted_names);
TH1 = 2; TH2 = 1;  
nLevel = 1;
iteration = 1;
 
MI_imfLTP = zeros(nLevel,16,nFold,iteration);  alfa_imfLTP = zeros(nLevel,16,nFold,iteration);NMI_imfLTP = zeros(nLevel,16,nFold,iteration);  
shifts_imfLTP = zeros(nLevel,16,2,nFold,iteration);

MI_unified = zeros(nLevel,16,nFold,iteration); alfa_unified = zeros(nLevel,16,nFold,iteration);NMI_unified = zeros(nLevel,16,nFold,iteration); 
shifts_unified = zeros(nLevel,16,2,nFold,iteration);
 
MI_unified2F = zeros(nLevel,16,nFold,iteration); alfa_unified2F = zeros(nLevel,16,nFold,iteration);NMI_unified2F = zeros(nLevel,16,nFold,iteration);  
shifts_unified2F = zeros(nLevel,16,2,nFold,iteration);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MI_imfLTP9 = zeros(nLevel,16);  alfa_imfLTP9 = zeros(nLevel,16);  shifts_imfLTP9 = zeros(nLevel,16);
NMI_imfLTP9 = zeros(nLevel,16); SSIM_imfLTP9 = zeros(nLevel,16); 

MI_unified9 = zeros(nLevel,16); alfa_unified9 = zeros(nLevel,16); shifts_unified9 = zeros(nLevel,16);
NMI_unified9 = zeros(nLevel,16); SSIM_unified9 = zeros(nLevel,16); 

MI_unified2F9 = zeros(nLevel,16); alfa_unified2F9 = zeros(nLevel,16); shifts_unified2F9 = zeros(nLevel,16);
NMI_unified2F9 = zeros(nLevel,16); SSIM_unified2F9 = zeros(nLevel,16); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
for iter = 1:iteration
    for jfold = 1:nFold
        foldName = sorted_names{jfold};
        select_fold = fullfile(pathname0,sorted_names{jfold});
        [path_names, exposures, numExposures] = readDir_Revised(select_fold);
        a = imread(path_names{1});
        [row,col,h] = size(a);
        I = cell(1,numExposures);
        MI = zeros(nLevel,16);
        NMI = zeros(nLevel,16);
        SSIM = zeros(nLevel,16);
        Dshifts = zeros(nLevel,16); Dalfa = zeros(nLevel,16);
        for jnoise = 1:nLevel
            noiseLevel = 10*(jnoise-1);
            for i=1:numExposures
                a = imread(path_names{i});
                if size(a,3)==3
                    a = rgb2gray(a);
                end
                a = uint8(double(a)+ noiseLevel*randn(size(a)));
                I{i} = a;    
            end            
            PM = cell(1,2);
            PM{1} = I{1};  
            reference = I{1};
            reference_50 = reference(51:end-50,51:end-50);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TH_imfLTP
            for i = 2:numExposures
                PM{2} = I{i};
                [sft, beta] = TH_imfLTPalignment(PM,2,1);
                %%% compute mutual information    
                tmp = PM{2};
                if any(abs(sft)>0.25)
                    xx = floor(sft(2)*4+0.5)/4;
                    yy = floor(sft(1)*4+0.5)/4;
                    tmp = shift(tmp,xx,yy);
                end
                if abs(beta)>0.01
                    tmp = imrotate(tmp,-beta,'bicubic','crop');
                end
                jmi = mutualInformation(reference,tmp,50);
                jnmi = nmi(reference,tmp,50);
                MI(jnoise,i) = jmi;
                NMI(jnoise,i) = jnmi;
                tmp1 = tmp(51:end-50,51:end-50);
                jssim = ssim_index(reference_50,tmp1);
                SSIM(jnoise,i) = jssim;
                a1 = abs(sft(1)-30);
                a2 = abs(sft(2)-10);
                a12 = [a1,a2];
                Dshifts(jnoise,i)= max(a12);
                Dalfa(jnoise,i) = abs(beta+5);
                MI_imfLTP(jnoise,i,jfold,iter) = jmi; 
                NMI_imfLTP(jnoise,i,jfold,iter) = jnmi; 
                alfa_imfLTP(jnoise,i,jfold,iter)= beta;
                shifts_imfLTP(jnoise,i,:,jfold,iter) = sft;
            end
        end
        MI_imfLTP9 = MI_imfLTP9 + MI;
        NMI_imfLTP9 = NMI_imfLTP9 + NMI;
        SSIM_imfLTP9 = SSIM_imfLTP9 + SSIM;
        alfa_imfLTP9 = alfa_imfLTP9 + Dalfa;
        shifts_imfLTP9 = shifts_imfLTP9 + Dshifts;
        savefile = ['imfLTP_',foldName,'.mat'];
        save(savefile, 'MI', 'NMI','SSIM','Dalfa','Dshifts')
    end
end
MI_imfLTP9N = sum(MI_imfLTP9,2)/56/iteration;
NMI_imfLTP9N = sum(NMI_imfLTP9,2)/56/iteration;
SSIM_imfLTP9N = sum(SSIM_imfLTP9,2)/56/iteration;
alfa_imfLTP9N = sum(alfa_imfLTP9,2)/56/iteration
shifts_imfLTP9N = sum(shifts_imfLTP9,2)/56/iteration
savefile = ['imfLTP_All_results','.mat'];
save(savefile, 'MI_imfLTP', 'NMI_imfLTP','alfa_imfLTP','shifts_imfLTP')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
for iter = 1:iteration
    for jfold = 1:nFold
        foldName = sorted_names{jfold};
        select_fold = fullfile(pathname0,sorted_names{jfold});
        [path_names, exposures, numExposures] = readDir_Revised(select_fold);
        a = imread(path_names{1});
        [row,col,h] = size(a);
        I = cell(1,numExposures);
        MI = zeros(nLevel,16);
        NMI = zeros(nLevel,16);
        SSIM = zeros(nLevel,16);
        Dshifts = zeros(nLevel,16); Dalfa = zeros(nLevel,16);
        for jnoise = 1:nLevel
            noiseLevel = 10*(jnoise-1);
            for i=1:numExposures
                a = imread(path_names{i});
                if size(a,3)==3
                    a = rgb2gray(a);
                end
                a = uint8(double(a)+ noiseLevel*randn(size(a)));
                I{i} = a;    
            end            
            PM = cell(1,2);
            PM{1} = I{1};  
            reference = I{1};
            reference_50 = reference(51:end-50,51:end-50);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% unified
            for i = 2:numExposures
                PM{2} = I{i};
                [sft, beta] = unifiedAlignment(PM,0);
                %%% compute mutual information    
                tmp = PM{2};
                if any(abs(sft)>0.25)
                    xx = floor(sft(2)*4+0.5)/4;
                    yy = floor(sft(1)*4+0.5)/4;
                    tmp = shift(tmp,xx,yy);
                end
                if abs(beta)>0.01
                    tmp = imrotate(tmp,-beta,'bicubic','crop');
                end
                jmi = mutualInformation(reference,tmp,50);
                jnmi = nmi(reference,tmp,50);
                MI(jnoise,i) = jmi;
                NMI(jnoise,i) = jnmi;
                tmp1 = tmp(51:end-50,51:end-50);
                jssim = ssim_index(reference_50,tmp1);
                SSIM(jnoise,i) = jssim;
                a1 = abs(sft(1)-30);
                a2 = abs(sft(2)-10);
                a12 = [a1,a2];
                Dshifts(jnoise,i)= max(a12);
                Dalfa(jnoise,i) = abs(beta+5);
                MI_unified(jnoise,i,jfold,iter) = jmi; 
                NMI_unified(jnoise,i,jfold,iter) = jnmi; 
                alfa_unified(jnoise,i,jfold,iter)= beta;
                shifts_unified(jnoise,i,:,jfold,iter) = sft;
            end
        end
        MI_unified9 = MI_unified9 + MI;
        NMI_unified9 = NMI_unified9 + NMI;
        SSIM_unified9 = SSIM_unified9 + SSIM;
        alfa_unified9 = alfa_unified9 + Dalfa;
        shifts_unified9 = shifts_unified9 + Dshifts;
        savefile = ['unified_',foldName,'.mat'];
        save(savefile, 'MI', 'NMI','SSIM','Dalfa','Dshifts')
    end
end
MI_unified9N = sum(MI_unified9,2)/56/iteration;
NMI_unified9N = sum(NMI_unified9,2)/56/iteration;
SSIM_unified9N = sum(SSIM_unified9,2)/56/iteration;
alfa_unified9N = sum(alfa_unified9,2)/56/iteration
shifts_unified9N = sum(shifts_unified9,2)/56/iteration
savefile = ['unified_All_results','.mat'];
save(savefile, 'MI_unified', 'NMI_unified','alfa_unified','shifts_unified')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
for iter = 1:iteration
    for jfold = 1:nFold
        foldName = sorted_names{jfold};
        select_fold = fullfile(pathname0,sorted_names{jfold});
        [path_names, exposures, numExposures] = readDir_Revised(select_fold);
        a = imread(path_names{1});
        [row,col,h] = size(a);
        I = cell(1,numExposures);
        MI = zeros(nLevel,16);
        NMI = zeros(nLevel,16);
        SSIM = zeros(nLevel,16);
        Dshifts = zeros(nLevel,16); Dalfa = zeros(nLevel,16);
        for jnoise = 1:nLevel
            noiseLevel = 10*(jnoise-1);
            for i=1:numExposures
                a = imread(path_names{i});
                if size(a,3)==3
                    a = rgb2gray(a);
                end
                a = uint8(double(a)+ noiseLevel*randn(size(a)));
                I{i} = a;    
            end            
            PM = cell(1,2);
            PM{1} = I{1};  
            reference = I{1};
            reference_50 = reference(51:end-50,51:end-50);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% unified2F
            for i = 2:numExposures
                PM{2} = I{i};
                [sft, beta] = unified2FeatureAlignment(PM,2,1);
                %%% compute mutual information    
                tmp = PM{2};
                if any(abs(sft)>0.25)
                    xx = floor(sft(2)*4+0.5)/4;
                    yy = floor(sft(1)*4+0.5)/4;
                    tmp = shift(tmp,xx,yy);
                end
                if abs(beta)>0.01
                    tmp = imrotate(tmp,-beta,'bicubic','crop');
                end
                jmi = mutualInformation(reference,tmp,50);
                jnmi = nmi(reference,tmp,50);
                MI(jnoise,i) = jmi;
                NMI(jnoise,i) = jnmi;
                tmp1 = tmp(51:end-50,51:end-50);
                jssim = ssim_index(reference_50,tmp1);
                SSIM(jnoise,i) = jssim;
                a1 = abs(sft(1)-30);
                a2 = abs(sft(2)-10);
                a12 = [a1,a2];
                Dshifts(jnoise,i)= max(a12);
                Dalfa(jnoise,i) = abs(beta+5);
                MI_unified2F(jnoise,i,jfold,iter) = jmi; 
                NMI_unified2F(jnoise,i,jfold,iter) = jnmi; 
                alfa_unified2F(jnoise,i,jfold,iter)= beta;
                shifts_unified2F(jnoise,i,:,jfold,iter) = sft;
            end
        end
        MI_unified2F9 = MI_unified2F9 + MI;
        NMI_unified2F9 = NMI_unified2F9 + NMI;
        SSIM_unified2F9 = SSIM_unified2F9 + SSIM;
        alfa_unified2F9 = alfa_unified2F9 + Dalfa;
        shifts_unified2F9 = shifts_unified2F9 + Dshifts;
        savefile = ['unified2F_',foldName,'.mat'];
        save(savefile, 'MI', 'NMI','SSIM','Dalfa','Dshifts')
    end
end
MI_unified2F9N = sum(MI_unified2F9,2)/56/iteration;
NMI_unified2F9N = sum(NMI_unified2F9,2)/56/iteration;
SSIM_unified2F9N = sum(SSIM_unified2F9,2)/56/iteration;
alfa_unified2F9N = sum(alfa_unified2F9,2)/56/iteration
shifts_unified2F9N = sum(shifts_unified2F9,2)/56/iteration
savefile = ['unified2F_All_results','.mat'];
save(savefile, 'MI_unified2F', 'NMI_unified2F','alfa_unified2F','shifts_unified2F')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Finish image alignment\n')

