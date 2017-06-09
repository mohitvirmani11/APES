function [SourceWaveform]=ComputeSourceWaveform(ProtocolName, nSubjects, nConditions)
% Please read the DOCUMETATION carefully to utilize this code. 

% begin DOCUMENTATION

% Author: Mohit Virmani 
% March 2017

% DESCRIPTION : This code computes the source waveform at all frames using
% the Linearly Constrained Minimum Variance (LCMV) Method Beamforming
% Algorithm. The script automates the source reconstruction for multiple
% subjects and multiple conditions. It uses the brainstorm[1] toolbox. 



% INPUT:
%   - ProtocolName: ProtocolName is a valid Brainstorm Protocol Name. Ex: 'Protocol1' 
%   - nSubjects   : No. of subjects. Numeric input for subjects. Ex: 20
%   - nConditions : No. of conditions. Numeric input for conditions. Ex: 3     

% In addition, the following inputs are asked from the user when the
% program is run. All the below mentioned input folder/ files should in the same
% placed in the same folder. 

%     - 'Choose the brainstorm toolbox directory' : Select the brainstorm folder
%     containing the brainstorm.m file. 

%     - 'Choose Data Directory which contains all subjects data, channel
%     file and headmodel file' : This directory should contain all subject
%     folders (explained below) which contains the .set preprocessed files,
%     channel file, and the headmodel file. 

%     - 'Choose warped anatomy template': This folder contains all the
%     warped anatomy files. The warped anatomy is obtained by warping the
%     channel file over the headmodel, thus representing the fitting of the
%     net over the headmodel. 

%     - Çhoose the channel file : (ex. GSN129.sfp) which contains the coordinates for the
%     all the electrode locations

%     - Choose the headmodel file: which contains the leadfield matrix (also called gain matrix or the forward model). 

%   Subjects folder - The Data Directory should contain one folder each
%   Subject. And each Subject's folder should contain one folder for each
%   condition. Each condition folder contains a .set and a .fdt file (these files are generated with EEGLAB[2], which
%   are the preprocessed files containing the corresponding data. 
%   Ex: for 2 subjects, the folder structures would be: 
%   ...Data\Subject1\EyesOpen\Subject1_EyesOpen.set
%   ...Data\Subject1\EyesOpen\Subject1_EyesOpen.fdt
%   ...Data\Subject1\Sternberg\Subject1_Sternberg.set
%   ...Data\Subject1\Sternberg\Subject1_Sternberg.fdt
        
%   ...Data\Subject2\EyesOpen\Subject2_EyesOpen.set
%   ...Data\Subject2\EyesOpen\Subject2_EyesOpen.fdt
%   ...Data\Subject2\Sternberg\Subject2_Sternberg.set
%   ...Data\Subject2\Sternberg\Subject2_Sternberg.fdt

% OUTPUT: 
%      - SourceWaveform 
%      (to output Averaged SourceWaveform uncomment the code segment
%      "Averaging all data points")


% Credits: Thanks to Prof. Ratna Sharma for the support in this work, 
% Dr. Navdeep Ahuja for providing the useful insights while developing this pipeline and Nishi Pegwal and other  
% members of Stress and Cognitive Electro Imaging Lab (SCEL),
% Department of Physiology, All India Institute of Medical Science (AIIMS), New Delhi for the their inputs. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References:
% [1] Brainstorm toolbox is documented and freely available 
% for download online under the GNU general public license (http://neuroimage.usc.edu/brainstorm).
% Tadel F, Baillet S, Mosher JC, Pantazis D, Leahy RM, “Brainstorm: A User-Friendly Application for MEG/EEG Analysis,” 
% Computational Intelligence and Neuroscience, vol. 2011, Article ID 879716, 13 pages, 2011. doi:10.1155/2011/879716

% [2] EEGLAB Toolbox 
% A Delorme & S Makeig (2004) EEGLAB: an open source toolbox for analysis of single-trial EEG dynamics (pdf, 0.7 MB) Journal of Neuroscience Methods 134:9-21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% end DOCUMENTATION

%% ============= DIRECTORIES ========== 

% Select brainstorm toolbox directory
disp('Choose the brainstorm toolbox directory');
Brainstorm_Dir= uigetdir;

% Select data directory
disp('Choose Data Directory which contains all subjects data, channel file and headmodel file');
DataDir=uigetdir;

% Warped anatomy folder      
disp('Choose warped anatomy template');
WarpedAnatomyTemplate = uigetdir; 

% Select Channel File
disp('mChoose channel file');
ChannelFileName = uigetfile('*.sfp');
ChannelFile=fullfile(DataDir,ChannelFileName);

% Select headmodel file
disp('Choose headmodel file');
HeadmodelName= uigetfile;
HeadmodelFile=fullfile(DataDir,HeadmodelName);

% Start timer
tic;
%% Start brainstorms
cd(Brainstorm_Dir);

if ~brainstorm('status')
    brainstorm nogui
end


%% ================== DECLARE OUTPUT FILE SIZES ============ 

SourceWaveform= struct([]); % Solution structure containing all the structures from Subject1 to Subject21.
EEGdata=zeros(129,1000,1,27);    % EEGdata output, for all epochs(3) and all subjects(4)
%% Create a new report 
% Start a new report
% bst_report('Start', sFiles);


%% 1. Create the protocol.

%Delete the protocol, if exists
gui_brainstorm('DeleteProtocol',ProtocolName);

% Create a new protocol 
gui_brainstorm('CreateProtocol',ProtocolName,0,0); % 0,0 for not using default anatomy and for not using default channel;

%Get Protocol Info
ProtocolInfo=bst_get('ProtocolInfo');

%nSubjects=20;
%nConditions=4;
SubjectNames=cell(nSubjects,1);
RawDataFiles=cell(nSubjects,nConditions);

cnt1=1; %counter for EO 
cnt2=1; %counter for EC
%% ================= CREATE SUBJECTS AND CONDITIONS ==========

for i=1:nSubjects
    SubjectNames{i}=strcat('Subject',num2str(i));
    [sSubject{i},iSubject{i}]=db_add_subject(SubjectNames{i});
    
    for k=1:nConditions
        switch k
               case 1
                    Condition{k}='EO';
%               case 2
%                    Condition{k}='EC';
%            case 1
%                 Condition{k}='Sternberg1';
%            case 2
%                Condition{k}='Sternberg2';
        end
        
        %% ================ LINK RAW FILE ==============
        FileName=strcat(SubjectNames{i},{'_'},Condition{k},{'.set'});
        RawFiles = fullfile(DataDir, SubjectNames{i},Condition{k},FileName{1});
        
        % Process: Create link to raw file
        RawDataFiles{i}{k} = bst_process('CallProcess', 'process_import_data_raw', [], [], ...
            'subjectname',    SubjectNames{i}, ...
            'datafile',       {RawFiles, 'EEG-EEGLAB'}, ...
            'channelreplace', 1, ...
            'channelalign',   1, ...
            'evtmode',        'value');
       
        %% ============== EXTRACT EPOCHS TO BRAINSTORM DATABASE =============== 
        %%%% CAUTION: Dont perform the below mentioned step before Setting channel
        %%%% file. 
        EpochDataFiles = bst_process('CallProcess', 'process_import_data_epoch', RawDataFiles{i}{k}.FileName, [], ...
            'subjectname', SubjectNames{i}, ...
            'condition',   '', ...
            'iepochs',     [], ...
            'eventtypes',  '', ...  
            'createcond',  0, ...
            'usectfcomp',  1, ...
            'usessp',      1, ...
            'freq',        [], ...
            'baseline',    []);

        % Number of Epochs
        nEpochs=size(EpochDataFiles,2); 
        EpochData=cell(nEpochs,1);
        InverseSolution=cell(nEpochs,1);
        
        %% ================ SET CHANNEL FILE ============

        % Process: Set channel file
        ChannelDataFile = bst_process('CallProcess', 'process_import_channel', RawDataFiles{i}{k}.FileName, [], ...
        'channelfile',  {ChannelFile, 'EGI'}, ...
        'usedefault',   1, ... 
        'channelalign', 1);

        %% ================ COPY WARP DEFAULT ANATOMY ===========
       
        %get the subject's folder
        AnatDir = fullfile(ProtocolInfo.SUBJECTS, SubjectNames{i});
        
        % get the names of all files. 
        dirListing = dir(WarpedAnatomyTemplate);

        % loop through the files and open. Note that dir also lists the directories, so you have to check for them.
        for d = 1:length(dirListing)
            if ~dirListing(d).isdir
             fileName = fullfile(WarpedAnatomyTemplate,dirListing(d).name); % use full path because the folder may not be the active path

                %copy file to the subject's directory
                copyfile(fileName,AnatDir);

            end % if-clause
        end % for-loop
        
        % Reload the anatomy of the selected subject
        db_reload_subjects(iSubject{i});

        %% ============== COPY HEAD MODEL =====================
       
            sStudy=bst_get('Study');
            iStudy=EpochDataFiles(1).iStudy;
            db_add(iStudy,HeadmodelFile);

       
        %% ================== DATA FOR EACH EPOCH ===============
        % Import data of the first epoch: Epochdata_1.Data has the 129X500 matrix data.
        for j=1:nEpochs
            EpochData{j}= bst_process('LoadInputFile',EpochDataFiles(j).FileName);

          %% =================== COMPUTE IDENTITY NOISE COVARIANCE MATRIX ============= 

         % Process: Compute identity noise covariance matrix
         bst_process('CallProcess', 'process_noisecov', EpochDataFiles(j).FileName, [], ...
            'baseline',       [0, 0.998], ...
            'datatimewindow', [0, 0.998], ...
            'sensortypes',    'EEG', ...
            'target',         1, ...  % Data covariance      (covariance over data time window)
            'dcoffset',       1, ...  % Block by block, to avoid effects of slow shifts in data
            'identity',       1, ...
            'copycond',       0, ...
            'copysubj',       0, ...
            'replacefile',    3);  % Replace

            
            
            
            %% ================= COMPUTE DATA COVARIANCE ============
            % Process: Compute data covariance
            bst_process('CallProcess', 'process_noisecov', EpochDataFiles(j).FileName , [], ...
                'baseline',       [0, 0.998], ...
                'datatimewindow', [0, 0.998], ...
                'sensortypes',    'EEG', ...
                'target',         2, ...  % Data covariance (covariance over data time window)
                'dcoffset',       1, ...  % Block by block, to avoid effects of slow shifts in data
                'identity',       0, ...
                'copycond',       0, ...
                'copysubj',       0, ...
                'replacefile',    1); %Replace

            % For each iteration, import the data covariance structure in a cell of
            % DataCovariance. This is done as data covariance for each epoch overwrites
            % the data covariance of the previious epoch. 

            %% ================ COMPUTE INVERSE SOLUTION ===========
            %Compute inverse solution 
            bst_process('CallProcess', 'process_inverse_2016', EpochDataFiles(j).FileName, [], ...
            'output',  3, ...  % Full results: one per file
            'inverse', struct(...
                 'Comment',        'PNAI: EEG', ...
                 'InverseMethod',  'lcmv', ...
                 'InverseMeasure', 'nai', ...
                 'SourceOrient',   {{'fixed'}}, ...
                 'Loose',          0.2, ...
                 'UseDepth',       1, ...
                 'WeightExp',      0.5, ...
                 'WeightLimit',    10, ...
                 'NoiseMethod',    'median', ...
                 'NoiseReg',       0.1, ...
                 'SnrMethod',      'rms', ...
                 'SnrRms',         1e-06, ...
                 'SnrFixed',       3, ...
                 'ComputeKernel',  0, ...
                 'DataTypes',      {{'EEG'}}));
               
           

            % Save the output
            sStudy=bst_get('Study'); %Get the details of the current study. 
             
            %Source Waveform
            Source_FullPath=file_fullpath(sStudy.Result(j).FileName);
            SourceSolution{j}= load(Source_FullPath);
            
            %Computing Source WaverForm 
            if j==1
                SourceWaveform(1).(SubjectNames{i}).(Condition{k})=SourceSolution{j}.ImageGridAmp;  
            else 
                SourceWaveform(1).(SubjectNames{i}).(Condition{k})(:,:,j)=SourceSolution{j}.ImageGridAmp;
            end 
            
          
        end % Epochs loop
        
        %% Averaging all data points
       % SourceWaveform.(SubjectNames{i}).(Condition{k})=mean(SourceWaveform(1).(SubjectNames{i}).(Condition{k})(:,1:nframes,:),2);

%         
    end % Conditions loop
    
end % Subjects loop
toc;
end
