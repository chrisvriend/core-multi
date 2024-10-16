
clc
clear

%% source toolboxes
addpath(genpath('/data/anw/anw-work/NP/doorgeefluik/toolboxes/BCT')) % Brain Connectivity toolbox - https://sites.google.com/site/bctnet/
addpath(genpath('/data/anw/anw-work/NP/doorgeefluik/toolboxes/network_community')) % Network Community toolbox - https://commdetect.weebly.com/
addpath(genpath('/data/anw/anw-work/NP/doorgeefluik/toolboxes/GenLouvain2.2')) % https://github.com/GenLouvain/GenLouvain
addpath(genpath('/data/anw/anw-work/NP/projects/data_chris/CORE/topological_filtering_networks')) % https://github.com/stdimitr/topological_filtering_networks



%% input variables
%% input variables
headdir='/data/anw/anw-work/NP/projects/data_chris/CORE/func';
session='ses-T1';
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
        pmatrix_temp=weight_conversion(pmatrix_temp,'normalize');
     
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


        clear pmatrix* W* temp

    end

end


