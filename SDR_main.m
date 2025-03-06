%Purpose:
%   Main function of the GPS software-defined receiver (SDR) based on 
%   vector tracking
%
%--------------------------------------------------------------------------
%                           GPSSDR_vt v1.0
% 
% Written by B. XU and L. T. HSU


% 
clear; 
format long g;
addpath geo             %  
addpath acqtckpos       % Acquisition, tracking, and postiong calculation functions


%% Parameter initialization 
%[file, signal, acq, track, solu, dyna, cmn] = initParametersOpenSky();
[file, signal, acq, track, solu, dyna, cmn] = initParametersUrban();

 
%% Acquisition 
if ~exist(['Acquired_',file.fileName,'_',num2str(file.skip),'.mat'])
    Acquired = acquisition(file,signal,acq); 
    save(['Acquired_',file.fileName,'_',num2str(file.skip)],'Acquired');    
else
    load(['Acquired_',file.fileName,'_',num2str(file.skip),'.mat']);
end 
fprintf(['Acquisi' ...
    'tion Completed. \n\n']);
 
% return


%% Do conventional signal tracking and obtain satellites ephemeris
fprintf('Tracking ... \n\n');
if ~exist(['eph_',file.fileName,'_',num2str(track.msToProcessCT/1000),'.mat'])
    % tracking using conventional DLL and PLL
    if ~exist(['TckResult_Eph',file.fileName,'_',num2str(track.msToProcessCT/1000),'.mat']) %
        [TckResultCT, CN0_Eph] =  trackingCT(file,signal,track,Acquired); 
        TckResult_Eph = TckResultCT;
        save(['TckResult_Eph',file.fileName,'_',num2str(track.msToProcessCT/1000)], 'TckResult_Eph','CN0_Eph');        
    else   
        load(['TckResult_Eph',file.fileName,'_',num2str(track.msToProcessCT/1000),'.mat']);
    end 
    
    % navigaion data decode
    fprintf('Navigation data decoding ... \n\n');
    [eph, ~, sbf] = naviDecode(Acquired, TckResult_Eph);
    save(['eph_',file.fileName,'_',num2str(track.msToProcessCT/1000)], 'eph');
    save(['sbf_',file.fileName,'_',num2str(track.msToProcessCT/1000)], 'sbf');
%     save(['TckRstct_',file.fileName,'_',num2str(track.msToProcessCT/1000)], 'TckResultCT'); % Track results are revised in function naviDecode for 20 ms T_coh
else
    load(['eph_',file.fileName,'_',num2str(track.msToProcessCT/1000),'.mat']);
    load(['sbf_',file.fileName,'_',num2str(track.msToProcessCT/1000),'.mat']);
    load(['TckResult_Eph',file.fileName,'_',num2str(track.msToProcessCT/1000),'.mat']);
end 
 
  
%% Find satellites that can be used to calculate user position
posSV  = findPosSV(file,Acquired,eph);
 
%% Do positiong in conventional or vector tracking mode
cnslxyz = llh2xyz(solu.iniPos);

if cmn.vtEnable == 1    
    fprintf('Positioning (VTL) ... \n\n');
  
    % load data to initilize VT
    load(['nAcquired_',file.fileName,'_',num2str(file.skip),'.mat']); % load acquired satellites that can be used to calculate position  
    Acquired = nAcquired;
    
    load(['eph_',file.fileName,'_',num2str(track.msToProcessCT/1000),'.mat']); % load eph
    load(['sbf_',file.fileName,'_',num2str(track.msToProcessCT/1000),'.mat']); % 
    
    load(['tckRstCT_1ms_',file.fileName,'.mat']);%,'_Grid'
    load(['navSolCT_1ms_',file.fileName,'.mat']); 
     
    [TckResultVT, navSolutionsVT] = ...
                  trackingVT_POS(file,signal,track,cmn,solu,Acquired,cnslxyz,eph,sbf,TckResult_Eph, TckResultCT_pos,navSolutionsCT);
else 
    load(['nAcquired_',file.fileName,'_',num2str(file.skip),'.mat']); % load acquired satellites that can be used to calculate position  
    Acquired = nAcquired;
    
    [TckResultCT_pos, navSolutionsCT] = ...
           trackingCT_POS(file,signal,track,cmn,Acquired,TckResult_Eph, cnslxyz,eph,sbf,solu); %trackingCT_POS_multiCorr_1ms
                 
end 

fprintf('Tracking and Positioing Completed.\n\n');
