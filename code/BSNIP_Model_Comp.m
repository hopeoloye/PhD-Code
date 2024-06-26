% -------------------------------------------------------------------------
% Script: SPM 
% Dataset: BSNIP
% Neuroimaging: EEG/ERP
% Paradigmn: P50 & P300
% Author: Lioba Bendt
% Date: 04/06/2024
% Topic: DCM Model Comparison
% Notes: Edited by Hope, updated R2 & DistCorr only 
% -------------------------------------------------------------------------
  clearvars -except subj_sheet R2_comp R2_mean
  dbstop if error
% -------------------------------------------------------------------------
% Set paths
% -------------------------------------------------------------------------
%   data_path       = '/Volumes/T7/NAPLS2/MMN/DCM_v08_PEB3/'; % folder containing DCM files
%   NAPLS_path      = '/Volumes/T7/NAPLS2/NAPLS2_ERP_Data/';
%   code_path       = '/Users/liobaberndt/Desktop/spm_materials/';
%   save_path       = '/Volumes/T7/NAPLS2/MMN/DCM_v08_PEB3/';
%   spm_path        = '/Users/liobaberndt/Desktop/spm_materials/';
%   fun_path        = '/Users/liobaberndt/Desktop/spm_materials/';
%   orig_path       = '/Volumes/T7/NAPLS2/MMN/SPMdata_epoch/';  

% data_path       = '/Volumes/Large Harddrive Backup/dcm_ei_sim/Hope_results/P50_preprocessed_final_dcm_inv_v8/dcms';  
%    data_path       = '/Volumes/Large Harddrive Backup/dcm_ei_sim/Hope_results/P300_preprocessed_final_dcm_inv_v16_bpa/dcms'; 
   data_path       = '/Volumes/Large Harddrive Backup/dcm_ei_sim/Hope_results/P300_preprocessed_final_dcm_inv_v16/dcms';    
  BSNIP_path      = '/Users/hopeoloye/Dropbox/Hope_PhD_New_22/BSNIP/';
  code_path       = '/Volumes/Large Harddrive Backup/dcm_ei_sim/code/';
%     save_path       = '/Volumes/Large Harddrive Backup/dcm_ei_sim/Hope_results/P50_preprocessed_final_dcm_inv_v8/diagnostics';
%    save_path       = '/Volumes/Large Harddrive Backup/dcm_ei_sim/Hope_results/P300_preprocessed_final_dcm_inv_v16_bpa/diagnostics';
   save_path       = '/Volumes/Large Harddrive Backup/dcm_ei_sim/Hope_results/P300_preprocessed_final_dcm_inv_v16/diagnostics';
  spm_path        = '/Users/hopeoloye/Applications/spm12'; %'/Users/hopeoloye/spm12/';
  fun_path        = '/Users/hopeoloye/Applications/spm12'; %'/Users/hopeoloye/spm12/';
%   orig_path       = '/Volumes/Large Harddrive Backup/BSNIP/non_DSS/P50/P50_preprocessed_final/';
  orig_path       = '/Volumes/Large Harddrive Backup/BSNIP/non_DSS/P300/P300_preprocessed_final/';

  addpath(spm_path, data_path, BSNIP_path, orig_path)
  
  cd(BSNIP_path)
  subj_sheet = xlsread('Subjects_BSNIP_Master.xlsx');

  % Get list of subjects from folders
cd(orig_path)
subj_struct = dir(orig_path);
for sub = length(subj_struct):-1:1
    % Remove Folders starting with '.' or other non-subject files/folders
    fname = subj_struct(sub).name;
    if fname(1) == '.'
        subj_struct(sub) = [ ];
    elseif fname(1) == 'E' || fname(1) == 'c' || fname(1) == 'l'
        subj_struct(sub) = [ ];
    end
end
subj = {subj_struct.name}'; Subjects = str2double(subj);

  
%  rmpath([fun_path filesep 'gifti-1.6']) % one copy in SPM12 already
% -------------------------------------------------------------------------
% Set model comparison
% -------------------------------------------------------------------------
  doBMS           = 0;      % do Bayesian model comparison (multiple models)
  doR2            = 1;      % look at R2 for one model across conditions
  doEp            = 0;      % get priors from best-fitting models (run doR2 first)
  compareR2       = 0;      % compare R2s from different inversions (rpt doR2 each time)
  doGrdMeans      = 0;      % compare grandmeans
% -------------------------------------------------------------------------
% R2 parameters
% -------------------------------------------------------------------------
  no_modes        = 4;      % number of data modes of interest (total = size(DCM.H{1},2); 4 probably enough for ERP, ?have to do all for CSD
  no_cond         = 2;      % ERP with 2 conditions or CSD with 1 condition?
% -------------------------------------------------------------------------
% Further parameters
% -------------------------------------------------------------------------  
%   EEG_strg = 'DCM_v0'; %'P50_DCM_v'; % title of DCM folder (without version number)
%   EEG_strg = 'P50_preprocessed_final_dcm_inv_v'; %'P50_DCM_v'; % title of DCM folder (without version number)
  EEG_strg = 'P300_preprocessed_final_dcm_inv_v'; %'P50_DCM_v'; % title of DCM folder (without version number)  
  %  DCM_name_length = 4;       % 4 if BSNIP, 3 if Uhlhaas ASSR data
  versions = 16; %[1 2] ;      % DCM version(s) to use - also ensure correct no. columns in line 127 if doBMS
% -------------------------------------------------------------------------
% Start SPM
% -------------------------------------------------------------------------    
  setup_paths  
  spm('defaults', 'eeg');

% -------------------------------------------------------------------------
% Load group info
% -------------------------------------------------------------------------
  cd(BSNIP_path)
%   subj_sheet = readtable('MMN Database March 2022 share.csv');
%   subj_sheet = xlsread('Subjects_BSNIP_Master.xlsx');
%   subj_sheet = subj_sheet(:,[3,4,6]);
  subj_sheet = subj_sheet(:,[5 6]);
%   filePattern = fullfile(orig_path,'bmdspmeeg_BSNIP*1vodfaster.mat');
%   filePattern = fullfile(orig_path, '*', '60-2avMrejdbMspmeeg_*.mat');
  filePattern = fullfile(orig_path, '*', 'f60-2avMrejdbMspmeeg_*.mat');
  
%   for sub = 1:length(subj)
%   
%       subj_path = [orig_path, subj{sub}, filesep];
%       cd(subj_path)
%       filename = dir(filePattern);
%       filename = filename.name;
%   end
%   
  
  matFiles = dir(filePattern);
%   DCM_files = fullfile(data_path,'DCM*.mat');
%   DCM_files = fullfile(data_path,'*_dcm_p50_cmc_peb_iter_1_v8.mat');
  DCM_files = fullfile(data_path,'*_dcm_p300_cmc_all_groups_inv1_v16.mat');  
  DCM_struct = dir(DCM_files);
  DCM_names  = {DCM_struct.name}';
  c = 1; n = 1;
% =========================================================================
% BMS 
% =========================================================================
  if doBMS == 1
     DCM_subj = [];
     DCM_subj_01 = [];
     DCM_subj_02 = [];
     model1 = '/Volumes/T7/BSNIP2/MMN/DCM_v04_noPEB/';
     model2 = '/Volumes/T7/BSNIP2/MMN/DCM_v08_noPEB/';
     filePattern = fullfile(model1,'DCM_*.mat');
     matFiles = dir(filePattern);
     matlabbatch{1}.spm.dcm.bms.inference.dir = {save_path};
     for r = 1 :length(matFiles)
         matFilename = fullfile(model1, matFiles(r).name);
         subject = extractBetween(matFilename, "/Volumes/T7/NAPLS2/MMN/DCM_v04_noPEB/DCM_", ".mat");
         subj_model1 = [model1, 'DCM_', subject,'.mat'];
         subj_model1 = strjoin(subj_model1);
         subj_model1 = num2str(subj_model1);
         subj_model1 = erase(subj_model1, " ");
         subj_model2 = [model2, 'DCM_', subject,'.mat'];
         subj_model2 = strjoin(subj_model2);
         subj_model2 = num2str(subj_model2);
         subj_model2 = erase(subj_model2, " ");
         matlabbatch{1}.spm.dcm.bms.inference.sess_dcm{r}.dcmmat = {subj_model1; subj_model2};
         clear subj_model1
         clear subj_model2
         clear subject
         continue
     end
     matlabbatch{1}.spm.dcm.bms.inference.model_sp = {''};
     matlabbatch{1}.spm.dcm.bms.inference.load_f = {''};
     matlabbatch{1}.spm.dcm.bms.inference.method = 'FFX';
     matlabbatch{1}.spm.dcm.bms.inference.family_level.family_file = {''};
     matlabbatch{1}.spm.dcm.bms.inference.bma.bma_no = 0;
     matlabbatch{1}.spm.dcm.bms.inference.verify_id = 1;
     spm_jobman('run',matlabbatch); 
  end 
% =========================================================================
% R2
% =========================================================================
  if doR2
%% -------------------------------------------------------------------------
% Load DCM results
% -------------------------------------------------------------------------   
%      filePattern = fullfile(data_path,'DCM*.mat');
    cd(data_path)     
%     filePattern = fullfile(data_path,'*_dcm_p50_cmc_peb_iter_1_v8.mat');
    filePattern = fullfile(data_path,'*_dcm_p300_cmc_all_groups_inv1_v16.mat');     
     matFiles = dir(filePattern);
     R2 = NaN(length(matFiles),no_cond+1);
     dcor1 = NaN(length(matFiles),no_modes);
     dcor2 = NaN(length(matFiles),no_modes);
     dcor3 = NaN(length(matFiles),no_modes);
    
     for f = 1:length(matFiles) 
         load(matFiles(f).name) 
% -------------------------------------------------------------------------
% Subjects with enough modes in the DCM 
% -------------------------------------------------------------------------
         try 
            if no_cond == 2
               % for m = 1:no_modes
% -------------------------------------------------------------------------
% Compute and store R2 for ERP
% -------------------------------------------------------------------------
                  % for m = 1:no_modes
                  %   % compute and store R2 for ERP
                  %   MSE_1 = sum(DCM.R{1}(:,m).^2); % for condition 1
                  %   MSE_2 = sum(DCM.R{2}(:,m).^2); % for condition 2
                  %   SS_1  = sum((DCM.H{1}(:,m)+DCM.R{1}(:,m)).^2);
                  %   SS_2  = sum((DCM.H{2}(:,m)+DCM.R{2}(:,m)).^2);
                  %   R2(f,m,1) = 1-(MSE_1/SS_1);                 % standard
                  %   R2(f,m,2) = 1-(MSE_2/SS_2);                 % deviant
                  %   R2(f,m,3) = 1-((MSE_1+MSE_2)/(SS_1+SS_2));  % both
                  %   clear MSE_1 MSE_2 SS_1 SS_2
                  % end
                  for m = 1:no_modes
                    % compute and store R2 for ERP
                    MSE_1 = sum(DCM.R{1}(:,m).^2); % for condition 1
                    MSE_2 = sum(DCM.R{2}(:,m).^2); % for condition 2
                    SS_1  = sum((DCM.H{1}(:,m)+DCM.R{1}(:,m)).^2);
                    SS_2  = sum((DCM.H{2}(:,m)+DCM.R{2}(:,m)).^2);
                    R2(f,m,1) = 1-(MSE_1/SS_1);                 % standard
                    R2(f,m,2) = 1-(MSE_2/SS_2);                 % deviant
                    R2(f,m,3) = 1-((MSE_1+MSE_2)/(SS_1+SS_2));  % both
                    clear MSE_1 MSE_2 SS_1 SS_2
                  end
                  
                   for m = 1:no_modes
                    %distance correlation
                    y1 = DCM.H{1}(:,m)+DCM.R{1}(:,m);
                    y2 = DCM.H{2}(:,m)+DCM.R{2}(:,m);
                    y3 = y2 - y1;
                 
%                   yhat = DCM.H{2};​
                    yhat1 = DCM.H{1}(:,m);  
                    yhat2 = DCM.H{2}(:,m);                  
                    yhat3 = yhat2 - yhat1;
                  
                    dcor1(f,m,1) = distcorr(y1(:,1),yhat1(:,1));
                    dcor2(f,m,1) = distcorr(y2(:,1),yhat2(:,1));
                    dcor3(f,m,1) = distcorr(y3(:,1),yhat3(:,1));
                    dcor(f,m,1) = distcorr(y1(:,1),yhat1(:,1));
                    dcor(f,m,2) = distcorr(y2(:,1),yhat2(:,1));
                    dcor(f,m,3) = distcorr(y3(:,1),yhat3(:,1));
                    
                    dcor(f,m,1) = distcorr(y1(:,1),yhat1(:,1));
                    dcor(f,m,2) = distcorr(y2(:,1),yhat2(:,1));
                    dcor(f,m,3) = distcorr(y3(:,1),yhat3(:,1));
                   end
                  
                % end  
            elseif no_cond == 1
                   peak_ind = DCM.Hz>40 & DCM.Hz<80; % indices of 38-42 peak
% -------------------------------------------------------------------------
% Compute and store R2 for CSD (without imaginary part)
% -------------------------------------------------------------------------              
                   MSE_1 = real(sum(sum(sum(DCM.Rc{1}(~peak_ind,:,:).^2)))); % Rc = channel space
                   MSE_2 = real(sum(sum(sum(DCM.Rc{1}(peak_ind,:,:).^2)))); 
                   SS_1  = real(sum(sum(sum((DCM.Hc{1}(~peak_ind,:,:)+DCM.Rc{1}(~peak_ind,:,:)).^2))));
                   SS_2  = real(sum(sum(sum((DCM.Hc{1}(peak_ind,:,:)+DCM.Rc{1}(peak_ind,:,:)).^2))));
                   R2(f,1) = 1-(MSE_1/SS_1); % not including peak
                   R2(f,2) = 1-(MSE_2/SS_2); % only peak
                   R2(f,3) = 1-((MSE_1+MSE_2)/(SS_1+SS_2)); % all
                   clear MSE_1 MSE_2 SS_1 SS_2           
            end
        end
    end
%% -------------------------------------------------------------------------
% Remove zeros from R2
% ------------------------------------------------------------------------- 
    R2(R2==0) = NaN;
    figure; pl = 1;
% -------------------------------------------------------------------------
% For ERPs
% -------------------------------------------------------------------------
    if no_cond == 2 
       for c = 1:no_cond+1
           for m = 1:no_modes

               subplot(no_cond+1,no_modes,pl)
%                hist(R2(:,m,c),-0.45:0.1:0.95); xlim([0 1]); % ylim([0 ceil(length(R2)/40)*10])
               histogram(R2(:,m,c),-0.45:0.1:0.95); ylim([0 200]); %xlim([0 1]); % ylim([0 ceil(length(R2)/40)*10])
               xlabel('R2'); %ylabel(['mode ' num2str(m)]);
               pl = pl+1;
               
           end
           
       end
       R2_1 = num2str(nanmean(R2(:,1,1))); R2_1 = R2_1(1:4);
       R2_2 = num2str(nanmean(R2(:,1,2))); R2_2 = R2_2(1:4);
       R2_3 = num2str(nanmean(R2(:,1,3))); R2_3 = R2_3(1:4);
       R2_1NaN = sum(isnan(R2(:,1,1)));
       
       subplot(no_cond+1,no_modes,2)  
       sgtitle(['n=' num2str(length(R2)) ', 1st mode: Cond 1 R2 = ' R2_1 ', Cond 2 R2 = ' R2_2 ', Total R2 = ' R2_3 ', NaN = ' num2str(R2_1NaN)])
%        subtitle(['n=' num2str(length(R2)) ', 1st mode: Cond 1 R2 = ' R2_1 ', Cond 2 R2 = ' R2_2 ', Total R2 = ' R2_3 ', NaN = ' num2str(R2_1NaN)])

%% -------------------------------------------------------------------------
% Plot Distance Correlation
% -------------------------------------------------------------------------
    figure; pl = 1;
           for c = 1:no_cond+1
           for m = 1:no_modes

               subplot(no_cond+1,no_modes,pl)
%                hist(R2(:,m,c),-0.45:0.1:0.95); xlim([0 1]); % ylim([0 ceil(length(R2)/40)*10])
               histogram(dcor(:,m,c),-0.45:0.1:0.95); ylim([0 250]); %xlim([0 1]); % ylim([0 ceil(length(R2)/40)*10])
               xlabel('dcor'); %ylabel(['mode ' num2str(m)]);
               pl = pl+1;
               
           end
           
           end
           
       
% %        bar(1:8,dcor1)
%        title('distance correlation')

       
       dcor_1 = num2str(nanmean(dcor(:,1,1))); dcor_1 = dcor_1(1:4);
       dcor_2 = num2str(nanmean(dcor(:,1,2))); dcor_2 = dcor_2(1:4);
       dcor_3 = num2str(nanmean(dcor(:,1,3))); dcor_3 = dcor_3(1:4);
       dcor_1NaN = sum(isnan(dcor(:,1,1)));
       
%        subplot(no_cond+1,no_modes,2)  
%        sgtitle(['n=' num2str(length(R2)) ', 1st mode: Cond 1 R2 = ' R2_1 ', Cond 2 R2 = ' R2_2 ', Total R2 = ' R2_3 ', NaN = ' num2str(R2_1NaN)])
       sgtitle(['n=' num2str(length(dcor)) ', 1st mode: Std dcor = ' dcor_1 ', Dev dcor = ' dcor_2 ', Dev-Std dcor = ' dcor_3 ', NaN = ' num2str(dcor_1NaN)])

% -------------------------------------------------------------------------
% Find subjects with best and worst R2 (across both conditions)
% -------------------------------------------------------------------------        
       best_fits = find(R2(:,1,3)>0.8);
       worst_fits = find(R2(:,1,3)<0.2);
% -------------------------------------------------------------------------
% For CSDs
% -------------------------------------------------------------------------   
    elseif no_cond == 1 
           for c = 1:3
               subplot(1,3,pl)
               hist(R2(:,m,c),-0.45:0.1:0.95); xlim([-1 1]); % ylim([0 60])%ylim([0 ceil(length(R2)/40)*10])
               xlabel('R2'); ylabel(['Version ' num2str(c)]);
               pl = pl+1;
           end
           R2_1 = num2str(nanmean(R2(:,1))); 
           R2_2 = num2str(nanmean(R2(:,2))); 
           R2_3 = num2str(nanmean(R2(:,3))); 
           R2_3NaN = sum(isnan(R2(:,3)));
           suptitle(['DCM v' num2str(versions) ', n=' num2str(length(R2)) ', Excl peak: R2 = ' R2_1 ', Peak only: R2 = ' R2_2 ', Total R2 = ' R2_3 ', NaN = ' num2str(R2_3NaN)])
% -------------------------------------------------------------------------
% Find subjects with best R2 (index 2 = peak only, 3 = overall)
% -------------------------------------------------------------------------        
           index = 2;
           criterion = 0.9999;
           best_fits = find(R2(:,index)>criterion);
           worst_fits = find(R2(:,index)<criterion);
    end
% -------------------------------------------------------------------------
% Which meanR2 measure to use? Of mode 1? modes 1-4? etc
% -------------------------------------------------------------------------    
%    meanR2 = nanmean(R2(:,[1:4]),2); %R2(:,1); %nanmean(R2,2);
%    figure('DefaultAxesFontSize',13)
%    hist(meanR2,10); box off
%    h(1) = findobj(gca,'Type','patch');
%    h(1).FaceColor = [1 1 1]*.8;
%    h(1).EdgeColor = [1 1 1]/2;
%    xlim([0 1]); %ylim([0 45]) 
%    title(['DCM v32 Mode 1-4 R2s; mean = ' num2str(nanmean(meanR2))])
%    xlabel('R2'); ylabel('Frequency')    
  end
% =========================================================================
% Ep
% =========================================================================
  if doEp
     if isequal(data_path(11:end),'/MMN/')
         cd([data_path EEG_strg '0' num2str(versions)])
     else
         cd([data_path EEG_strg num2str(versions)])
     end
     for b = 1:length(best_fits)
         fname = files(best_fits(b)).name;
         load(fname)
% -------------------------------------------------------------------------
% Vectorise posteriors
% ------------------------------------------------------------------------- 
          Ep_all(:,b) = spm_vec(DCM.Ep); 
      end    
% -------------------------------------------------------------------------
% Average all posteriors and put back into Ep structure
% -------------------------------------------------------------------------    
      Ep = mean(Ep_all,2);    
      Ep = spm_unvec(Ep_all,DCM.Ep);
% -------------------------------------------------------------------------
% Get covariances from initial DCM and save - CHECK CORRECT
% -------------------------------------------------------------------------    
      load('/Volumes/T7/NAPLS2/MMN/DCM/DCM_01-0003-1.mat')
      pC = DCM.M.pC;
      pE = Ep;
      cd(data_path)
      save(num2str(versions,'Priors_BSNIP%.f.mat'), 'pE', 'pC')
      cd(code_path)    
  end
% =========================================================================
% Compare R2
% =========================================================================
  if compareR2
     if ~exist('R2_comp','var') 
        R2_comp = []; % initialise
     end
     if size(R2_comp,2) == 0        
        R2_comp(:,1)     = R2(:,1,3); 
        R2_mean          = nanmean(R2(:,3));
     else
         try
         R2_comp(:,end+1) = R2(:,1,3);
         end
         R2_mean(1,end+1) = nanmean(R2(:,3));
     end    
     figure; 
     bar(R2_mean)
     xticklabels({'2','5','6','7','8','9','10','11','12','13'}); xlabel('DCM version'); ylabel('R2')
     suptitle(['MMN model R2s, ' num2str(no_modes) ' modes used'])    
  end
%