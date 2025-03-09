function [file, signal, acq, track, solu, dyna, cmn] = initParameters()
%Purpose:
%   Parameter initialization
%Inputs: 
%	None
%Outputs:
%	file        - parameters related to the data file to be processed,a structure
%	signal      - parameters related to signals,a structure
%	acq         - parameters related to signal acquisition,a structure
%	track       - parameters related to signal tracking,a structure
%	solu        - parameters related to navigation solution,a structure
%	cmn         - parameters commmonly used,a structure
%--------------------------------------------------------------------------
%                           GPSSDR_vt v1.0
% 
% Written by B. XU and L. T. HSU

%% File parameters
file.fileName       = 'Urban';
file.fileRoute      = ['Data\',file.fileName,'.dat'];  
file.skip        	= 0; % in unit of ms

global ALPHA BETA % Parameters for iono. and trop. correction; From RINEX file
ALPHA = [0.1118e-07  0.7451e-08 -0.5960e-07 -0.5960e-07];
BETA  = [0.9011e+05  0.1638e+05 -0.1966e+06 -0.6554e+05];  
cmn.doy = 90; % Day of year 

%% File parameters
file.fid           	= fopen(file.fileRoute,'r','ieee-le');
file.skiptimeVT     = 0; % skip time from the first measurement epoch of CT, in uint of msec
file.dataType       = 2;    %1:I; 2:IQ
file.dataPrecision  = 1;    %1:int8 or byte; 2; int16 

%% Signal parameters
signal.IF               = 0; % unit: Hz 
signal.Fs               = 26e6;	% unit: Hz
signal.Fc               = 1575.42e6; % unit: Hz	
signal.codeFreqBasis	= 1.023e6; % unit: Hz 	
signal.ms               = 1e-3; % unit: s
signal.Sample           = ceil(signal.Fs*signal.ms);	
signal.codelength       = signal.codeFreqBasis * signal.ms;

%% Acquisition parameters
acq.prnList     = 1:32;	% PRN list
acq.freqStep    = 500;	% unit: Hz
acq.freqMin     = -10000;   % Minimum Doppler frequency
acq.freqNum     = 2*abs(acq.freqMin)/acq.freqStep+1;    % number of frequency bins
acq.L           = 10;   % number of ms to perform FFT

%% Tracking parameters
track.mode                  = 0;    % 0:conventional tracking; 1:vector tracking
track.CorrelatorSpacing  	= 0.5;  % unit: chip
track.DLLBW               	= 5;	% unit: Hz
track.DLLDamp           	= 0.7; 
track.DLLGain            	= 0.1;	% OriginalL: 0.1
track.PLLBW              	= 15;
track.PLLDamp             	= 0.7;
track.PLLGain              	= 0.25; 	
track.msToProcessCT       	= 90000; % unit: ms 40000
track.msPosCT               = 90000; % unit: ms
track.msToProcessVT         = 5000; %track.msPosCT - file.skiptimeVT; %
track.pdi                   = 1; %


%% Navigation solution parameters
solu.navSolPeriod = 10; % unit: ms 
solu.mode  	= 0;    % 0:conventional LS/WLS; 1:conventionalKF; 2:VT
solu.flag_spaceUsr = 0; % 0:ground user; 1:spaceborne user, for tropo correction switching
solu.accel_navStateVT = 0; % 0:no accelerations as navigation states; 1:accelerations as navigation states in VTL EKF
solu.sat_dynamic = 0; % 0: no dynamic model involved; 1: satellite dynamic model aiding in VTL EKF
solu.iniCoeffAlpha_pos = 0.5; % Position weighting of VTL EKF in VTL+DM integration
solu.iniCoeffBeta_pos = 0.5; % Position weighting of DM in VTL+DM integration
solu.iniCoeffAlpha_vel = 0.5; % Velocity weighting of VTL EKF in VTL+DM integration
solu.iniCoeffBeta_vel = 0.5; % Velocity weighting of DM in VTL+DM integration

%% Navigation solution initial conditions
% solu.iniPos	= [3091670.971904	376405.140683	6150286.705625]; % Ground truth location at HD_GPSL1Exc1 localTime 264894.0
% solu.iniPos	= [2599519.854832	1497566.305129	6207014.393900]; % Ground truth location at HD_GPSL1 (FULL) localTime 259200.0
solu.GPSWN = 1046; % Ground truth GPS week number
solu.GPSdSOW = 264894.0; % Ground truth localTime (seconds of week GPST)
solu.iniPos = [22.3198722 * pi/180, 114.209101777778 * pi/180, 0]; % Ground truth location at HD_GPSL1 (FULL) localTime 264894.0
solu.iniVel = [0 0 0];% Ground truth velocity at HD_GPSL1 (FULL) localTime 264894.0
%% Dynamic model parameters
dyna.updInterval = 0.2; % unit: s

%% commonly used parameters
cmn.vtEnable  	= 0;%   % 0: disable vector tracking; 1:enable vector tracking
cmn.cSpeed      = 299792458;    % speed of light, [m/s]

%% end

