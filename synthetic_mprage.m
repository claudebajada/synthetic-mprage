function [T1w_image , T1w_header] = synthetic_mprage(R1, PD, R2_star)

%% Function to create a synthetic MPRAGE based on Lor et al. (2015)
% Neurobiological origin of spurious brain morphological changes: A quantitative MRI study
% The synthetic T1w parameters were obtained from this paper. These were in
% turn obtained from Tardif et al. (2009).
% 
% T1w = synthetic_mprage(R1)
% 

% probably should insert some check to see if R1, PD and R2_star seem to be
% in seconds (basically check that they have the same order of magnitude. 
% 
% would be nice if i could give the option to insert own parameters.

header = spm_vol(R1);
R1_image = spm_read_vols(header);

PD = spm_read_vols(spm_vol(PD));
R2_star = spm_read_vols(spm_vol(R2_star));

flip_angle = 9; % in degrees
echo_spacing = 0.0099; % in seconds (this is because R1, PD and R2_star are in milliseconds in our data)
TI = .960; % in seconds
readout_duration = 176 * echo_spacing; % in seconds
TR = 2.42;
delay_time = TR - (TI + (readout_duration/2));

T1_star = (R1_image - ((1/echo_spacing) .* (log(cosd(flip_angle))))).^-1;
E1 = exp(-TI .* R1_image);
E2 = exp(-delay_time*R1_image);
E3 = exp(-(readout_duration ./  T1_star));
E4 = exp(-(readout_duration ./ (2*T1_star)));

t1 = (1 - 2.*E1 + E1.*E2) ./ (1 + E1.*E2.*E3);
t2 = (1 + E1.*E2.*E3 - E4 - E4.*E1.*E2) ./ (1 + E1.*E2.*E3);
fR1 = E4 .* t1 + T1_star .* R1_image .* t2; % this bit needs more thourough checking
fPD = PD;
fR2_star = exp(-echo_spacing * R2_star);

% T1w_image = real(fR1 .* sind(flip_angle)); % this is here just in case someone wants to get an MPRAGE from only R1 data
T1w_image = fR1 .* fPD .* fR2_star .* sind(flip_angle); % I am unsure about this bit where I only take the real component of the complex number ... I dont know if i should be taking the magnitude

T1w_image = T1w_image * 1000; % to convert to milliseconds

T1w_header = header;
T1w_header.fname = 'synthetic_mprage.nii';
T1w_header.descrip = 'synthetic_mprage'; % perhaps should put more data in the descriptor eg: the parameters used for the generation of the mprage
