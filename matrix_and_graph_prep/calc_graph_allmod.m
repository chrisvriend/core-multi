
clc
clear
%% source toolboxes
addpath(genpath('/data/anw/anw-gold/NP/doorgeefluik/toolboxes/BCT'))
addpath(genpath('/data/anw/anw-gold/NP/doorgeefluik/toolboxes/network_community'))
addpath(genpath('/data/anw/anw-gold/NP/doorgeefluik/toolboxes/GenLouvain2.2'))
addpath(genpath('/data/anw/anw-gold/NP/projects/data_chris/CORE/topological_filtering_networks'))

% kruskal_algorithm for multilayer MST in topological_filtering_networks folder

modality={'dwi'};
atlas='300P7N';
samplename='CORE';
Nnodes=314;

niter_louv = 100;           % number of iterations community detection (default 100)
n_iter=100;
n_try=5;

% reduce to 300x300 300P7N matrix (minus subcortical areas) for MIND cf
% analyses
MIND=0;
% ------

for jj = 1:length(modality)
    mod=modality{jj};
    disp(mod)


    %% input variables
    headdir=strcat('/data/anw/anw-gold/NP/projects/data_chris/CORE',filesep,mod);
    session='ses-T1';
    condir=strcat(headdir,filesep,session,filesep,'graph');
    outputdir=condir;
    mkdir(outputdir)


    if strcmp(mod,'func')
        desc='func';
        G=1.00;
    elseif strcmp(mod,'dwi')
        desc='dwi';
        G=0.48;

    end
    % Gamma variable for community detection | CV 210624 with Gamma = 1 |
    % N comm = 8.1 +/- 0.08 for func
    % Gamma = 0.48 | N comm = 8.00 +/- 0.07


    atlasparcels='/data/anw/anw-gold/NP/projects/data_chris/CORE/func/300P7N-to-network.legend';
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
    %matrices=dir('*matrix.csv');

    matrices=dir('*300P7N_OMST.mat');

    % preallocation

    subjects=cell(length(matrices),1);
    Eglob=nan(length(matrices),1);
    Eglob_norm=Eglob;
    lambda=nan(length(matrices),1);
    lambda_norm=lambda;
    avclust=nan(length(matrices),1);
    avclust_norm=avclust;
    consensus_simm=nan(length(matrices),1);
    modules=nan(Nnodes,length(matrices));
    modidx=nan(length(matrices),1);
    PC_glob_ind=nan(length(matrices),1);
    PC_glob_Schaef=nan(length(matrices),1);
    PC_Schaef=nan(length(matrices),length(subnetworks));
    S=nan(length(matrices),1);
    clustering=nan(length(matrices),Nnodes);
    btwneeness=nan(length(matrices),Nnodes);
    Qyeo=nan(length(matrices),1);
    % Preallocate cell arrays for tables
    resultTables = cell(length(matrices), 1);

    %% calc graph measures
    for i = 1:length(matrices)


        %    subjects{i,1}=extractBefore(matrices(i).name,strcat('_atlas-',atlas,'_desc-streams_connmatrix.csv'));

        subjects{i,1}=extractBefore(matrices(i).name,strcat('_acq-',desc,'_atlas-',atlas,'_OMST.mat'));

        disp(['working on ' subjects{i,1}])

        file=matfile(matrices(i).name);
        W_OMST=file.W_OMST;

        if exist(strcat(subjects{i,1},'_acq-',desc,'_atlas-',atlas,'_rand_OMST.mat'),'file')~=2
            error('no random matrix found')
        end

        if size(W_OMST,1)~=Nnodes
            error('matrix dimension not consistent with number of nodes')

        end

        file_rnd=matfile(strcat(subjects{i,1},'_acq-',desc,'_atlas-',atlas,'_rand_OMST.mat'));
        W_rnd_OMST=file_rnd.W_rnd_OMST;


        %% Exclude subcortical regions (optional for comparison with MIND analyses)
        if MIND==1
            disp('modify size matrix too 300x300')

            subcnodes=cell2mat(subnetID(strcmp(subnetworks,'SUBC')));
            W_OMST(subcnodes,:)=[];
            W_OMST(:,subcnodes)=[];
            W_rnd_OMST(subcnodes,:)=[];
            W_rnd_OMST(:,subcnodes)=[];

        else
            disp(['using orig size matrix = ' num2str(Nnodes)])

        end


        %% -----

        % Compute basic graph indices
        N=length(W_OMST);                          % number of vertices
        k = full(sum(W_OMST));                     % degree
        twom = sum(k);                        % number of edges (each undirected edge is counted twice)
        B =full(W_OMST - G * (k.' * k)/twom);      % modularity matrix

        communities = zeros(size(W_OMST,1), niter_louv);

        disp('running community detection')
        tic
        for x = 1:niter_louv
            [communities(:,x),Qual,n_it]=iterated_genlouvain(B,10000,0);
            %disp(strcat('Quality of parcellation =  ',num2str(Qual)))
            %disp(strcat('number of iterations = ', num2str(n_it)))
        end
        toc

        % perform consensus similarity across N=niter_louv
        [consensus, consensus_simm(i,1), ~] = consensus_similarity(transpose(communities));
        disp(['there are ' num2str(length(unique(consensus))) ' modules']);

        disp(strcat('similarity of consensus to all others = ', num2str(consensus_simm(i,1))))

        % adapted from modularity_und.m
        s=consensus(ones(1,N),:);                      %compute modularity (transpose of original script)
        Q=~(s-s.').*B/twom;
        Q=sum(Q(:));

        % assign to variable
        modules(:,i)=consensus;



        %% nodal measures

        % clustering coefficient
        clustering(i, :) = clustering_coef_wu(W_OMST)';
        %avclust=mean(clustering(i,:));

        % betweenness centrality | normalised to the range [0,1] using
        % BC/[(N-1)(N-2)] with N = #nodes
        btwneeness(i,:) = betweenness_wei(weight_conversion(W_OMST,'lengths'))/((Nnodes-1)*(Nnodes-2));

        % Participation Coefficient

        PC_individ=participation_coef(W_OMST,modules(:,i));
         PC_atlas=participation_coef(W_OMST,atlasIDS4PC);

        % normalized

        Eglob_rnd=nan(length(n_iter),1);
        avclust_rnd=nan(length(n_iter),1);
        lambda_rnd=nan(length(n_iter),1);
        disp('calculating randomized measures')
        tic
        for kk = 1:n_iter

            % Random network GE and CC
            Eglob_rnd(kk,1) = efficiency_wei(W_rnd_OMST(:,:,kk)); % GE of randomized networks
            avclust_rnd(kk,1) = mean(clustering_coef_wu(W_rnd_OMST(:,:,kk))); % CC of randomized networks
            lambda_rnd(kk,1)=charpath(distance_wei(weight_conversion(W_rnd_OMST(:,:,kk),'lengths')));
        end
        toc


        %% global measures
        % global efficiency
        Eglob(i,1) = efficiency_wei(W_OMST);
        Eglob_norm(i,1) = Eglob(i,1)/mean(Eglob_rnd);

        % Weighted clustering coefficient
        avclust(i,1) = mean(clustering(i, :));
        avclust_norm(i,1) = avclust(i,1)/mean(avclust_rnd);

        % Calculate the small-worldness coefficient

        lambda(i,1)=charpath(distance_wei(weight_conversion(W_OMST,'lengths')));
        lambda_norm(i,1)=lambda(i,1)/mean(lambda_rnd);
        % small worldness - sigma
        S(i,1) = avclust_norm(i,1) / lambda_norm(i,1);

        % Modularity
        modidx(i,1)=Q;
        clear Q
        % modularity based on Yeo + subcort parcellation
        Qyeo(i,1)=calc_Q_Schaefer(W_OMST,atlasIDS4PC);

        PC_glob_ind(i,1)=mean(PC_individ);
        PC_glob_Schaef(i,1)=mean(PC_atlas);


        %% network-network matrix


        % pre-allocation
        bnconn=zeros(length(subnetworks),length(subnetworks));

        % between system connections / Pc
        for k = 1:length(subnetworks)-1

            for l = k+1:length(subnetworks)
                bnconn(k,l) = mean2(W_OMST(subnetID{k}, subnetID{l}));
            end


        end

%         
%         % Calculate the number of non-zero elements
%         num_nonzero_elements = nnz(W_OMST);
%         
%         % Calculate the mean of non-zero elements
%         if num_nonzero_elements > 0
%             mean_connectivity_excluding_zeros = total_connectivity_strength_nonzero / num_nonzero_elements;
%         else
%             mean_connectivity_excluding_zeros = 0; % handle the case where there are no connections
%         end



        % symmetrize matrix
        bnconn=bnconn+bnconn';


        % add within system connections and calc PC
        for k = 1:length(subnetworks)
          %  bnconn(k,k) = mean2(W_OMST(subnetID{k}, subnetID{k}));

            PC_Schaef(i,k)=mean(PC_atlas(subnetID{k}));

        end

        bnconn=weight_conversion(bnconn,'normalize');

        %sysconn=array2table(bnconn,'VariableNames', subnetworks, 'RowNames', subnetworks);
        valueslong=reshape(bnconn,[],1);

        % Initialize arrays for the long format
        % Generate all combinations of these entries including themselves
        numEntries = numel(subnetworks);
        [combIdx1, combIdx2] = ndgrid(1:numEntries, 1:numEntries);

        % Flatten the indices
        combIdx1 = combIdx1(:);
        combIdx2 = combIdx2(:);

        % Preallocate cell arrays for ROI1 and ROI2
        ROI1 = cell(numel(combIdx1), 1);
        ROI2 = cell(numel(combIdx2), 1);
        subjid=cell(numel(combIdx1),1); 

        % Fill ROI1 and ROI2 with the corresponding subnetworks entries
        for j = 1:numel(combIdx1)
            ROI1{j} = subnetworks{combIdx1(j)};
            ROI2{j} = subnetworks{combIdx2(j)};
            subjid{j}=subjects{i,1};
        end


        % Create a table
        resultTables{i} = table(subjid,ROI1, ROI2,valueslong,'VariableNames', {'Subj','ROI1', 'ROI2', 'Y'});


        clear bnconn ROI1 ROI2 valueslong W_rnd_OMST W_OMST Eglob_rnd avclust_rnd lambda_rnd


    end



    %% save output

    % global measures
    globalheads={'GE','GE_norm','GCC','GCC_norm','lambda', 'lambda_norm', 'Sigma','Q','Qyeo','PC_glob_ind','PC_glob_Schaef'};
    globalX=array2table([Eglob Eglob_norm avclust avclust_norm lambda lambda_norm S modidx Qyeo PC_glob_ind PC_glob_Schaef],'VariableNames',globalheads);

    globalX.('subjID')=subjects;
    globalX=[globalX(:,end) globalX(:,1:end-1)];
    writetable(globalX,strcat(outputdir, filesep,'graph_global_acq-',desc,'_',samplename, '.csv'),'Delimiter',',','WriteVariableNames', 1);

    PCTable=array2table(PC_Schaef,'VariableNames',subnetworks);
    PCTable.('subjID')=subjects;
    writetable(PCTable,strcat(outputdir,filesep,'graph_PC_subnetworks_acq-',desc,'_',samplename,'.csv'),'Delimiter',',','WriteVariableNames', 1);

    MBAinput=vertcat(resultTables{:});
    writetable(MBAinput,strcat(outputdir,filesep,'MBA_input_acq-',desc,'_',samplename,'.txt'),'Delimiter',' ','WriteVariableNames',1)


end
