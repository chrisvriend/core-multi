

clc
clear
% 
% conda activate /scratch/anw/share-np/multilayer
% conda activate python env with Spyder
% 
% Open Spyder > go to Tools > Preferences > Python Interpreter > Change the 'use the following Python interpreter' to your environment's one. 
% /scratch/anw/share-np/multilayer/bin/python > Apply and close tools >Restart Kernel console in spyder



session='ses-T1';
basedir='/data/anw/anw-gold/NP/projects/data_chris/CORE';

fmridir=strcat(basedir,filesep,'func',filesep,session,filesep,'graph');
dmridir=strcat(basedir,filesep,'dwi',filesep,session,filesep,'graph');
multidir=strcat(basedir,filesep,'multi',filesep,session,filesep,'graph');
subjectfile=strcat(basedir,filesep,'multi',filesep,session,filesep,'subjects_',session,'_multilayer.txt');
%subjectfile=strcat(basedir,filesep,'multi',filesep,session,filesep,'HC_subjects_','ses-T0','_multilayer.txt');


subjects=readtable(subjectfile,'ReadVariableNames', false); % col 1 = subj | col 2 = ses

Nnodes=314;
% N modalities
nlrs=2;
id = Nnodes:Nnodes:Nnodes*nlrs;


% construct weighted supra-adjacency matrix
for c = 1:height(subjects)

   subjbase=cell2mat(strcat(subjects{c,1},'_',subjects{c,2}));
% for HC data
 %   subjbase=cell2mat(strcat(subjects{c,1}));

    fprintf(1, 'Now constructing supra_weighted for sub %s!\n', num2str(c))

    dwimat=strcat(dmridir,filesep,subjbase,'_acq-dwi_atlas-300P7N_OMST_MST.mat');
    funcmat=strcat(fmridir,filesep,subjbase,'_acq-func_atlas-300P7N_OMST_MST.mat');

    dwi=matfile(dwimat);
    func=matfile(funcmat);
    
    supramat = blkdiag(squeeze(func.W_OMST_MST),squeeze(dwi.W_OMST_MST));
    for i = 1:length(id)
        supramat(id(i)*size(supramat,1)+1:size(supramat,1)+1:end) = 1;
        supramat(id(i)+1:size(supramat, 1)+1:1+size(supramat, 1)*min(size(supramat, 1)-id(i),size(supramat, 2))) = 1;
    end
  


    save(strcat(multidir,filesep,subjbase,'_atlas-300P7N_multiplex.mat'),'supramat')
    clear supramat dwi* func* subjbase
end


cd(multidir)


files=dir('*_multiplex.mat');

multi_matrices=nan(Nnodes*nlrs,Nnodes*nlrs,length(files));

for i=1:length(files)

    multimat=matfile(files(i).name);
multi_matrices(:,:,i)=multimat.supramat;

end
save(strcat(multidir,filesep,'CORE_multiplexes.mat'),'multi_matrices')


