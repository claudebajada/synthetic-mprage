function [T1w_image , T1w_header] = synthetic_mprage(R1, PD, R2_star)

%% Function to create a synthetic MPRAGE based on Lorio et al. (2015)
% Neurobiological origin of spurious brain morphological changes: A quantitative MRI study
% The synthetic T1w parameters were obtained from this paper. These were in
% turn obtained from Tardif et al. (2009).
% 
% [T1w_image , T1w_header] = synthetic_mprage(R1, PD, R2_star)
% 

% probably should insert some check to see if R1, PD and R2_star seem to be
% in seconds (basically check that they have the same order of magnitude. 
% 
% would be nice if i could give the option to insert own parameters.

%% read in the header
header = spm_vol(R1);

%% read in the R1 PD and R2_star maps (assumes that they are in seconds)
R1 = spm_read_vols(header);
R1(R1<0)=0;% remove negatives

PD = spm_read_vols(spm_vol(PD));
PD(PD<0) = 0; % remove negatives

R2_star = spm_read_vols(spm_vol(R2_star));
R2_star(R2_star<0) = 0; % remove negatives

%% Set parameters
flip_angle = 9; % alpha - in degrees
TE = 0.00419;
ES = 0.0099; % echo spacing in seconds (this is because R1, PD and R2_star are in seconds in our data)
TI = .960; % in seconds
readout_duration = 176 * ES; % in seconds
TR = 2.42;
TD = TR - (TI + (readout_duration/2)); % delay time

%% MPRAGE signal equation
T1_star = (R1 - ((1/ES) .* (log(cosd(flip_angle))))).^-1;
E1 = exp(-TI .* R1);
E2 = exp(-TD*R1);
E3 = exp(-(readout_duration ./  T1_star));
E4 = exp(-(readout_duration ./ (2*T1_star)));

x1 = (1 - 2.*E1 + E1.*E2) ./ (1 + E1.*E2.*E3);
x2 = (1 + E1.*E2.*E3 - E4 - E4.*E1.*E2) ./ (1 + E1.*E2.*E3);

fPD = PD;
fR2_star = exp(-TE * R2_star);
fR1 = E4 .* x1 + T1_star .* R1 .* x2; 

% T1w_image = fR1 .* sind(flip_angle); % this is here just in case someone wants to get an MPRAGE from only R1 data
% T1w_image = fR1 .* fPD .* sind(flip_angle); 
T1w_image = fR1 .* fPD .* fR2_star .* sind(flip_angle); 

%% Convert T1w image to milliseconds
T1w_image = T1w_image * 1000;

%% Set up descriptors
T1w_header = header;
T1w_header.fname = 'synthetic_mprage.nii';
T1w_header.descrip = ['synthetic_mprage ' , ...
    'TE=' , num2str(TE*1000) , 'ms ' , ...
    'TR=' , num2str(TR*1000) , 'ms ' , ...
    'ES=' , num2str(ES*1000) , 'ms ' , ...
    'TD=' , num2str(TD*1000) , 'ms ' , ...
    'TI=' , num2str(TI*1000) , 'ms ' , ...
    'Readout Duration (theta)=' , num2str(readout_duration*1000) , 'ms'];
