
clc
clear

%% add toolboxes to path
addpath(genpath('/data/anw/anw-gold/NP/doorgeefluik/toolboxes/BCT'))
addpath(genpath('/data/anw/anw-gold/NP/projects/data_chris/CORE/topological_filtering_networks'))


%% input variables
%% input variables
headdir='/data/anw/anw-gold/NP/projects/data_chris/CORE/func';
session='ses-T1';
%headdir='/s4ever/anwG/NP/yvdwerf/data_MOTAR/analysis'
timedir=strcat(headdir,filesep,session);
atlas='300P7N';

outputdir=strcat(timedir,filesep,'graph');
samplename='CORE';
Nnodes=314;
Nvols=269;
mkdir(outputdir)



%% load and rearrange data
cd(timedir)

files=dir('sub-*_timeseries.csv');


conn_matrices=nan(Nnodes,Nnodes,length(files));
subjects=cell(length(files),1);

for a=1:length(files)

    file=files(a).name;
    subjects(a,1)=strcat('sub-',extractBetween(file,'sub-','_atlas'));
    subj=subjects{a,1};
    if exist(strcat(outputdir,filesep,subj,'_acq-func_atlas-',atlas,'_OMST.mat'),'file')~=2
    disp(subj)

        % timeseries file
        tseries = readtable(file);
        epoch=table2array(tseries);

        % check size
        if size(epoch,1)~=Nvols || size(epoch,2)~=Nnodes
            error(strcat(file, ' does not have the right dimensions'))

        end



        %% perform pearson correlation
        temp=corr(epoch);
        temp(logical(eye(size(temp)))) = 0;


    %    temp =.5*log((1+temp)./(1-temp)); % convert to Fisher z-scores
        %temp=abs(temp); % only take absolute into account (might also want to check for other ways, such as (prop) thresh)

        % normalize
        pmatrix_temp=weight_conversion(temp,'autofix'); % autofix to ensure that matrix is truly symmetrical,
     %   pmatrix_temp=weight_conversion(pmatrix_temp,'normalize');

        %         % find zero entries in matrix to add back after Z-transformation and
        %         % rescaling
        %         logic0=pmatrix_temp==0;
        %         disp(strcat('there are ',' ', num2str(sum(sum(logic0))-Nnodes), ' off-diagonal zero entries in the matrix'))
        %
        %         Mu=mean(pmatrix_temp(~logic0));
        %         Sig=std(pmatrix_temp(~logic0));
        %
        %         pmatrix_temp2=(pmatrix_temp - Mu) / Sig;
        %         pmatrix_temp2(logic0)=0;
        %
        %         pmatrix_temp3=pmatrix_temp2 - min(pmatrix_temp2(:));
        %
        %         pears_matrix=pmatrix_temp3 ./ max(pmatrix_temp3(:));
        %         pears_matrix(logic0)=0;

        W=abs(pmatrix_temp);

        [~, W_OMST,~]=threshold_omst_gce_wu(squeeze(W),1);
        W_OMST_MST=kruskal_algorithm(squeeze(W_OMST));
        W_MST=kruskal_algorithm(squeeze(W));

        disp(['similarity between OMST - MST =' num2str(dice(W_MST,W_OMST_MST))])


        disp('density =')
        disp(num2str(density_und(W_OMST)));

        %save 2D coherence matrix (MXN) where M and N are number of nodes
        save(strcat(outputdir,filesep,subj,'_acq-func','_atlas-',atlas,'_OMST.mat'),'W_OMST', '-v7.3')
        save(strcat(outputdir,filesep,subj,'_acq-func','_atlas-',atlas,'_OMST_MST.mat'),'W_OMST_MST', '-v7.3')

        %% Partial Correlations

%         C=cov(epoch);
%         invC=inv(C);
%         pCorr=zeros(size(C));
% 
% 
%         % Step 4: Calculate partial correlations
%         for i = 1:size(C,1)
%             for j = 1:size(C,2)
%                 if i ~= j
%                     pCorr(i,j) = -invC(i,j) / sqrt(invC(i,i) * invC(j,j));
%                 else
%                     pCorr(i,j) = 1;  % Diagonal elements are 1 (correlation with itself)
%                 end
%             end
%         end
%         pCorr=real(pCorr);
%         pCorr(logical(eye(size(pCorr)))) = 0;
% 
%         pCorr =.5*log((1+pCorr)./(1-pCorr)); % convert to Fisher z-scores
%         %temp=abs(temp); % only take absolute into account (might also want to check for other ways, such as (prop) thresh)
% 
%         % normalize
%         pmatrix_temp=weight_conversion(pCorr,'autofix'); % autofix to ensure that matrix is truly symmetrical,
%         pmatrix_temp=weight_conversion(pmatrix_temp,'normalize');
% 


        clear temp pmatrix* epoch*

    end

end


