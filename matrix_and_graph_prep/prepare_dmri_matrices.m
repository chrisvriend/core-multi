%% C. Vriend - Amsterdam UMC - Jul '24
%% construct dMRI matrices from MRTRIX3 derived connectivity matrices saved as csv file
%% apply OMST thresholding

clc 
clear

%% source toolboxes
addpath(genpath('/data/anw/anw-work/NP/doorgeefluik/toolboxes/BCT')) % Brain Connectivity toolbox - https://sites.google.com/site/bctnet/
addpath(genpath('/data/anw/anw-work/NP/doorgeefluik/toolboxes/network_community')) % Network Community toolbox - https://commdetect.weebly.com/
addpath(genpath('/data/anw/anw-work/NP/doorgeefluik/toolboxes/GenLouvain2.2')) % https://github.com/GenLouvain/GenLouvain
addpath(genpath('/data/anw/anw-work/NP/projects/data_chris/CORE/topological_filtering_networks')) % https://github.com/stdimitr/topological_filtering_networks




% kruskal_algorithm for multilayer MST in topological_filtering_networks folder

%% input variables
headdir='/data/anw/anw-work/NP/projects/data_chris/CORE/dwi';
session='ses-T0';
condir=strcat(headdir,filesep,session);
outputdir=strcat(condir,filesep,'graph');
mkdir(outputdir)
atlas='300P7N';
samplename='CORE';
Nnodes=314;

niter_louv = 100;           % number of iterations community detection (default 100)
G = 1.00;                   % Gamma variable for community detection
n_try = 5;
n_iter = 100;

atlasparcels='/data/anw/anw-work/NP/projects/data_chris/CORE/func/300P7N-to-network.legend';
parcels=readtable(atlasparcels,'FileType','text');

subnetworks=unique(parcels{:,3});
subnetID=cell(length(subnetworks),1);

for j = 1:length(subnetworks)

    % Create a logical index for rows where the third column equals 'DAN'
    idx = strcmp(parcels{:, 3}, subnetworks{j});

    % Extract the rows that match the condition
    subnetID{j,1}= table2array(parcels(idx, 1))';


end
% convert atlas txt labels to numbers
[~, ~, atlasIDS4PC] = unique(parcels{:,3});

%-----------------------

% find matrices
cd(condir)
matrices=dir('sub-*matrix.csv');


% preallocation

subjects=cell(length(matrices),1);


%% calc graph measures
for i = 1:length(matrices)

     subjects{i,1}=extractBefore(matrices(i).name,strcat('_atlas-',atlas,'_desc-streams_connmatrix.csv'));


    if exist(strcat(outputdir,filesep,subjects{i,1},'_acq-dwi','_atlas-',atlas,'_OMST.mat'),'file')~=2
     disp(['working on ' subjects{i,1}])

    file=matrices(i).name;
    W=readmatrix(file);
    W=weight_conversion(weight_conversion(W,'normalize'),'autofix');
    
    if size(W,1)~=Nnodes
        error('matrix dimension not consistent with number of nodes')

    end
    disp(['density before ' num2str(density_und(W)) ])
    tic
    [temp,W_OMST,~]=threshold_omst_gce_wu(squeeze(W),1);
    toc
    disp(['density after ' num2str(density_und(W_OMST)) ])
    % MST for multilayer
    W_OMST_MST=kruskal_algorithm(squeeze(W_OMST));
    W_MST=kruskal_algorithm(squeeze(W));
    disp(['similarity between OMST - MST =' num2str(dice(W_MST,W_OMST_MST))])

% 
%     % construct random network
%     tic
%     disp('randomise network')
%     W_OMST_rand = null_model_und_sign(squeeze(W_OMST),100);
%     
%     W_OMST_rand_MST=kruskal_algorithm(squeeze(W_OMST_rand));
%     toc




    save(strcat(outputdir,filesep,subjects{i,1},'_acq-dwi','_atlas-',atlas,'_OMST.mat'),'W_OMST')
    save(strcat(outputdir,filesep,subjects{i,1},'_acq-dwi','_atlas-',atlas,'_OMST_MST.mat'),'W_OMST_MST')
%     save(strcat(outputdir,filesep,subjects{i,1},'_acq-dwi','_mod-rand','_atlas-',atlas,'_OMST_MST.mat'),'W_OMST_rand_MST')

     else
     disp('already converted')
     end


end
