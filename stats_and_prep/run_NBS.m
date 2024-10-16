

clc
clear

%% source software
addpath(genpath('/data/anw/anw-gold/NP/doorgeefluik/toolboxes/NBS1.2'))



modality={'func','dwi'};
atlas='300P7N';
samplename='CORE';
Nnodes=314;


atlasparcels='/data/anw/anw-gold/NP/projects/data_chris/CORE/func/300P7N-to-network.legend';
parcels=readtable(atlasparcels,'FileType','text');
labels=parcels{:,2};
% select only OCD trials
OCD=0

% T-threshold for NBS
Tthresh=3.1
compsize='extent'
% change this
analysisname='responder';
NBS

for jj = 1:length(modality)
    mod=modality{jj};
    disp(mod)


    %% input variables
    clindir=strcat('/data/anw/anw-gold/NP/projects/data_chris/CORE/clin');
    headdir=strcat('/data/anw/anw-gold/NP/projects/data_chris/CORE',filesep,mod);
    statsdir=strcat('/data/anw/anw-gold/NP/projects/data_chris/CORE/stats');
    session='ses-T0';
    condir=strcat(headdir,filesep,session,filesep,'graph');

   

    save(strcat(condir,filesep,'NBS_300P7N_labels_temp.mat'),'labels');


    % find matrices
    cd(condir)
    %matrices=dir('*matrix.csv');
    NBStempdir=strcat(condir,filesep,'NBS_temp');
    [~,~,~]=mkdir(strcat(NBStempdir,filesep,'conn'));
    [~,~,~]=mkdir(strcat(NBStempdir,filesep,'design'));
    [~,~,~]=mkdir(strcat(statsdir,filesep,'NBS'));



    if OCD == 1
        matricesA=dir('*ARRIBA*300P7N_OMST.mat');
        matricesB=dir('*TIPICCO*300P7N_OMST.mat');
        matrices=vertcat(matricesA,matricesB);
        Tclin=readtable(strcat(clindir,filesep,'OCD_trials4analyses.xlsx'));

    else
        matrices=dir('*300P7N_OMST.mat');
        Tclin=readtable(strcat(clindir,filesep,'OCDPTSD_trials4analyses.xlsx'));

    end

    connmatrices=nan(Nnodes,Nnodes,length(matrices));
    subjects=cell(length(matrices),1);

    for i = 1:length(matrices)

        subjects{i,1}=extractBefore(matrices(i).name,strcat('_acq-',mod,'_atlas-',atlas,'_OMST.mat'));

        disp(['working on ' subjects{i,1}])

        file=matfile(matrices(i).name);
        W_OMST=file.W_OMST;
        connmatrices(:,:,i)=W_OMST;
        clear W_OMST

    end

    save(strcat(NBStempdir,filesep,'conn',filesep,'NBS_acq-',mod,'_connmatrices_temp.mat'),'connmatrices')

    subjs=strrep(subjects,strcat('_',session),'');


    % subset table
    Tclin2=Tclin(ismember(Tclin.Subj, subjs),:);

    Tclin2=sortrows(Tclin2,'Subj');


    if ~isequal(Tclin2.Subj, subjs)
        error('Tclin_subset.Subj is NOT the same as subjs, or the order is different.');
    end
    % convert sex column to numbers
    Tclin2.seks = double(strcmp(Tclin2.seks, 'm'));

    % define design matrix and contrast
    % w/ dYBOCS
    %   design=horzcat(ones(length(subjs),1),Tclin2.dYBOCS,Tclin2.baseline_sev,Tclin2.age,Tclin2.seks);

    if strcmp(analysisname,'perc_improv')

        design=horzcat(ones(length(subjs),1),Tclin2.perc_improv,Tclin2.age,Tclin2.seks,Tclin2.responder);
        con1=[0 1 0 0 0];
        con2=[0 -1 0 0 0];
    elseif strcmp(analysisname,'responder')

        design=horzcat(Tclin2.responder, ...
            ~Tclin2.responder,Tclin2.age,Tclin2.seks);
        con1=[1 -1 0 0 ];
        con2=[-1 1 0 0 ];
    end

    save(strcat(NBStempdir,filesep,'design',filesep,'NBS_acq-',mod,'_design_temp.mat'),'design')
    save(strcat(NBStempdir,filesep,'design',filesep,'NBS_acq-',mod,'_con1_temp.mat'),'con1')
    save(strcat(NBStempdir,filesep,'design',filesep,'NBS_acq-',mod,'_con2_temp.mat'),'con2')

    % for both contrasts
    for i=[1 2]
        UI.method.ui='Run NBS';
        UI.test.ui='F-test';
        UI.size.ui=compsize;
        UI.thresh.ui=num2str(Tthresh);
        UI.perms.ui='5000';
        UI.alpha.ui='0.05';
        UI.design.ui=strcat(NBStempdir,filesep,'design',filesep,'NBS_acq-',mod,'_design_temp.mat');
        UI.exchange.ui='';
        UI.matrices.ui=strcat(NBStempdir,filesep,'conn',filesep,'NBS_acq-',mod,'_connmatrices_temp.mat');
        UI.contrast.ui=strcat(NBStempdir,filesep,'design',filesep,'NBS_acq-',mod,'_con',num2str(i),'_temp.mat');
        UI.node_coor.ui='';
        UI.node_label.ui=strcat(condir,filesep,'NBS_300P7N_labels_temp.mat');


        NBSrun(UI,[])

        global nbs

        output_file = strcat(statsdir,filesep,'NBS',filesep,'NBS_acq-',mod,'_',analysisname,'_con', num2str(i), '_',num2str(Tthresh),'_',compsize,'_perm5000_extent.mat');

        save(output_file, 'nbs')
    end
end
