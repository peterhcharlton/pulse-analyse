function [pw_inds, fid_pts, pulses, sigs] = PulseAnalyse(S, options)
% PULSEANALYSE  Extracts pulse wave (PW) indices from pulse waves.
%   
%  Inputs:
%
%    S        -  a pulsatile signal, consisting of either a single pulse, or a
%                 signal containing several pulses. S should be a structure, containing
%                 a vector of amplitudes, S.v, and the sampling frequency (in Hz), S.fs,
%                 and optionally the subject's height in metres, S.ht.
%    options  -  (optional) a structure of options which determine the settings used for the analysis:
%                    options.do_plot                    - (1 or 0, default value of 1) A logical indicating whether or not to make plots
%                    options.exclude_low_quality_data   - (default of 1) A logical indicating whether or not to exclude low quality pulse waves from the calculation of median values of pulse wave indices
%                    options.plot_third_deriv           - (default of 0) A logical indicating whether or not to plot the third derivative
%                    options.close_figures              - (default of 1) A logical indicating whether or not to close all Matlab figures before running the analysis
%                    options.do_filter                  - (default of 1) A logical indicating whether or not to filter pulse waves prior to analysis
%                    options.save_folder                - (default of '') A string containing the path of a directory in which to save plots
%                    options.save_file                  - (default of '') A string containing the filename under which to save plots (without an extension)
%                    options.beat_detector              - (default of 'IMS') The type of beat detector to use when analysing a prolonged recording of multiple pulse waves. Definitions: IMS - [incremental merge segmentation algorithm](https://doi.org/10.1109/EMBC.2012.6346628) (implemented by M.A.F Pimentel as part of the [RRest](https://github.com/peterhcharlton/RRest) toolbox).
%                    options.verbose                    - (default of 0) A logical indicating whether or not to provide text outputs for any warnings
%                    options.plot_areas                 - (default of 0) A logical indicating whether or not to plot the systolic and diastolic areas on the pulse wave
%                    options.plot_pw_only               - (default of 0) A logical indicating whether or not to plot only the pulse wave and not its derivatives
%                    options.normalise_pw               - (default of 1) A logical indicating whether or not to normalise the pulse wave to occupy a range of 0 to 1
%
%  Outputs:
%
%    pw_inds  -  a structure containing the calculate pulse wave indices. 
%                   For instance, pw_inds.AI.v is the value of the
%                   augmentation index (AI). If the input S contains
%                   several pulses, then  pw_inds.AI.v is the median AI
%                   value, and pw_inds.AI.raw is a vector of raw values for
%                   each pulse.
%    fid_pts  -  a structure containing:
%                   fid_pts.ind: the indices of the fiducial points
%                   fid_pts.amp: the ampltiudes of the fiducial points
%                   fid_pts.amp_norm: the amplitudes of the fiducial points (when the pulse wave is normalised to occupy a range of 0 to 1) 
%                   fid_pts.t: the times of the fiducial points
%                   For instance, fid_pts.ind.dic provides the indices of the
%                   dicrotic notch(es).
%    pulses   -  a structure containing information on the pulses:
%                   pulses.peaks    - the indices of systolic peaks
%                   pulses.onsets   - the indices of pulse onsets (i.e. beginning of systole)
%                   pulses.quality  - a logical indicating the signal quality of each pulse (true indicates high quality)
%    sigs     -  the following signals are provided:
%                   sigs.orig - the original signal
%                   sigs.filt - a filtered signal (if asked for)
%                   sigs.first_d - first derivative
%                   sigs.second_d - second derivative
%                   sigs.third_d - third derivative
%
%  Exemplary usage:
%
%    [pw_inds, fid_pts, pulses, sigs] = PulseAnalyse            performs a demo pulse wave analysis (of a photoplethysmogram, PPG, pulse wave)
%    [pw_inds, fid_pts, pulses, sigs] = PulseAnalyse([],'pressure_example')     performs a demo pressure pulse wave analysis
%    [pw_inds, fid_pts, pulses, sigs] = PulseAnalyse(S)         analyses the pulsatile signal specified in S
%    [pw_inds, fid_pts, pulses, sigs] = PulseAnalyse(S, options)                uses options to specify the settings for the analysis
%
%  For further information please see the accompanying manual: https://github.com/peterhcharlton/pulse-analyse/wiki/Examples
%
%   Licence:
%       Available under the GNU public license - please see the accompanying
%       file named "LICENSE"
%
%   Sources:
%       This script contains items either copied or modified from the RRest
%       toolbox which is covered by the GNU public licence (<a href="http://github.com/peterhcharlton/RRest/">link</a>).
%
% Author: Peter H. Charlton, King's College London, June 2019.
% Version 1.0

%% Setup

% - check to see whether a signal has been provided
if nargin == 2 && strcmp(options, 'pressure_example'), fprintf('\n - Running basic pressure example'), S = basic_example_S('pressure'); clear options, options.normalise_pw = 0; options.do_filter = 0; end
if nargin == 0 || isempty(S) || ~sum(strcmp(fieldnames(S),'v')), fprintf('\n - No input signal provided, so running basic PPG example'), S = basic_example_S('ppg'); end
% - create options structure if it hasn't been provided
if nargin < 2, options = struct; end
% - insert nan for height if it's not provided
if ~sum(strcmp(fieldnames(S), 'ht')), S.ht = nan; end
% - convert height into metres if it's in cm
if S.ht > 20, S.ht = S.ht/100; end
% - make sure S.v is a column vector
S.v = S.v(:);
% - setup universal parameters
up = setup_up(S, options);

%% Beat Detection
pulses = beat_detection(S, up);

%% Quality Assessment
pulses.quality = assess_signal_quality(S, pulses.peaks);

%% Filter Signal
sigs = filter_signal(S, up);

%% Calculate Derivatives
sigs = calculate_derivs(sigs, up);

%% Identify Fiducial Points, and Calculate Fiducial Point Timings and Amplitudes
fid_pts = identify_fiducial_point_indices(sigs, pulses, up);

%% Calculate Pulse Wave Indices
pw_inds = calculate_pw_inds(fid_pts, sigs, pulses, S.ht, up);

%% Plot Results
plot_results(sigs, pulses, fid_pts, pw_inds, up);

end

function options = setup_options(options)

if isempty(fieldnames(options)) || ~sum(strcmp(fieldnames(options), 'do_plot'))
    options.do_plot = 1;
end

if ~strcmp(fieldnames(options), 'exclude_low_quality_data')
    options.exclude_low_quality_data = 1;
end
if ~strcmp(fieldnames(options), 'plot_third_deriv')
    options.plot_third_deriv = 0;
end
if ~strcmp(fieldnames(options), 'close_figures')
    options.close_figures = 1;
end
if ~strcmp(fieldnames(options), 'do_filter')
    options.do_filter = 1;
end
if ~strcmp(fieldnames(options), 'save_folder')
    options.save_folder = '';
end
if ~strcmp(fieldnames(options), 'save_file')
    options.save_file = '';
end
if ~strcmp(fieldnames(options), 'beat_detector')
    options.beat_detector = 'IMS';
end
if ~strcmp(fieldnames(options), 'verbose')
    options.verbose = false;
end
if ~strcmp(fieldnames(options), 'plot_areas')
    options.plot_areas = false;
end
if ~strcmp(fieldnames(options), 'plot_pw_only')
    options.plot_pw_only = false;
end
if ~strcmp(fieldnames(options), 'normalise_pw')
    options.normalise_pw = true;
end

end

function S = basic_example_S(signal_type)

if strcmp(signal_type, 'ppg')
    % Data from the PPGDiary Pilot Study (a single pulse wave acquired using an infrared PPG sensor on the thumb)
    S.v = [-14880;-14508;-13390;-11384;-8573;-5236;-1745;1569;4502;6979;8992;10536;11600;12196;12396;12328;12123;11860;11547;11160;10700;10218;9789;9458;9211;8990;8742;8465;8206;8016;7894;7778;7573;7223;6747;6239;5814;5554;5478;5556;5744;6010;6339;6720;7137;7561;7961;8309;8580;8760;8840;8825;8732;8580;8384;8151;7884;7595;7306;7042;6817;6629;6460;6286;6093;5874;5630;5364;5077;4771;4453;4132;3822;3534;3273;3041;2828;2622;2410;2187;1967;1771;1617;1501;1391;1243;1029;760;486;269;140;79;28;-69;-231;-425;-599;-716;-784;-847;-949;-1104;-1291;-1475;-1639;-1800;-1993;-2239;-2524;-2798;-3011;-3150;-3256;-3395;-3613;-3898;-4189;-4417;-4562;-4658;-4771;-4942;-5167;-5410;-5645;-5865;-6077;-6270;-6414;-6491;-6521;-6568;-6694;-6911;-7167;-7385;-7519;-7588;-7653;-7778;-7982;-8245;-8523;-8777;-8992;-9169;-9323;-9475;-9648;-9855;-10093;-10344;-10577;-10768;-10908;-11012;-11116;-11267;-11502;-11815;-12151;-12427;-12585;-12644;-12698;-12845;-13082;-13241];
    S.fs = 100; % in Hz
elseif strcmp(signal_type, 'pressure')
    % A simulated carotid pressure pulse wave, from the pulse wave database described at: https://peterhcharlton.github.io/pwdb/'
    S.v = [74.83;74.86;75.02;75.35;75.89;76.6;77.44;78.36;79.31;80.27;81.2;82.12;82.99;83.83;84.62;85.36;86.05;86.71;87.32;87.92;88.49;89.06;89.61;90.16;90.71;91.26;91.8;92.34;92.88;93.41;93.95;94.48;95.02;95.56;96.11;96.67;97.22;97.77;98.32;98.86;99.39;99.9;100.4;100.9;101.3;101.8;102.2;102.5;102.9;103.2;103.5;103.7;103.9;104.1;104.2;104.3;104.4;104.4;104.4;104.4;104.4;104.3;104.3;104.2;104.1;103.9;103.8;103.7;103.5;103.4;103.2;103;102.9;102.7;102.5;102.3;102.1;102;101.8;101.6;101.4;101.2;101.1;100.9;100.7;100.5;100.4;100.2;100.1;99.94;99.81;99.69;99.59;99.49;99.4;99.32;99.25;99.19;99.15;99.11;99.08;99.06;99.05;99.05;99.06;99.08;99.11;99.14;99.18;99.23;99.28;99.34;99.4;99.46;99.53;99.6;99.67;99.74;99.81;99.88;99.95;100;100.1;100.1;100.2;100.2;100.2;100.2;100.2;100.2;100.2;100.1;100.1;100;99.9;99.78;99.63;99.46;99.27;99.04;98.79;98.52;98.23;97.93;97.61;97.29;96.97;96.65;96.35;96.11;95.96;95.92;95.99;96.14;96.29;96.42;96.49;96.5;96.46;96.37;96.25;96.11;95.95;95.77;95.59;95.41;95.23;95.06;94.89;94.72;94.55;94.39;94.24;94.09;93.94;93.81;93.68;93.57;93.47;93.37;93.3;93.23;93.17;93.12;93.07;93.03;92.99;92.96;92.92;92.88;92.84;92.79;92.75;92.7;92.64;92.58;92.52;92.46;92.39;92.32;92.25;92.17;92.1;92.02;91.94;91.85;91.77;91.69;91.61;91.53;91.45;91.36;91.28;91.2;91.12;91.04;90.96;90.88;90.8;90.72;90.63;90.54;90.45;90.36;90.27;90.17;90.08;89.98;89.88;89.78;89.69;89.59;89.49;89.39;89.3;89.2;89.11;89.02;88.93;88.84;88.75;88.67;88.58;88.5;88.42;88.34;88.26;88.19;88.12;88.04;87.97;87.9;87.83;87.77;87.7;87.64;87.57;87.51;87.45;87.38;87.32;87.26;87.2;87.14;87.08;87.02;86.96;86.9;86.84;86.78;86.72;86.66;86.6;86.53;86.47;86.41;86.34;86.28;86.21;86.15;86.08;86.01;85.94;85.87;85.8;85.73;85.66;85.58;85.51;85.43;85.36;85.28;85.2;85.12;85.03;84.95;84.87;84.78;84.69;84.61;84.52;84.43;84.34;84.25;84.15;84.06;83.96;83.87;83.77;83.68;83.58;83.48;83.39;83.29;83.19;83.09;82.99;82.9;82.8;82.7;82.6;82.51;82.41;82.31;82.21;82.12;82.02;81.92;81.83;81.73;81.64;81.54;81.45;81.35;81.26;81.16;81.07;80.98;80.88;80.79;80.7;80.61;80.52;80.42;80.33;80.24;80.15;80.06;79.97;79.88;79.79;79.71;79.62;79.53;79.44;79.35;79.27;79.18;79.09;79.01;78.92;78.84;78.75;78.66;78.58;78.49;78.41;78.33;78.24;78.16;78.07;77.99;77.91;77.82;77.74;77.66;77.57;77.49;77.41;77.33;77.24;77.16;77.08;77;76.92;76.84;76.76;76.68;76.6;76.52;76.44;76.36;76.28;76.2;76.12;76.04;75.96;75.88;75.81;75.73;75.65;75.57;75.5;75.42;75.35;75.27;75.19;75.12;75.04;74.97;74.9];
    S.fs = 500; % in Hz
end

end

function pulses = beat_detection(S, up)

if strcmp(up.analysis.no_of_pulses, 'multiple')
    
    %% Filter to remove high frequencies
    filt_characteristics = up.filtering.beat_detection.elim_high_freqs;
    s_filt = elim_vhfs3(S, filt_characteristics);
    
    %% Filter to remove low frequencies
    filt_characteristics = up.filtering.beat_detection.elim_low_freqs;
    s_filt = elim_vlfs(s_filt, filt_characteristics, up);
    
    %% Detect beats
    pulses = detect_beats(s_filt, up);
    
else

    %% Locate onsets and peak in single pulse wave
    pulses.onsets = [1,length(S.v)];
    [~, pulses.peaks] = max(S.v);
    
end

end

function up = setup_up(S, options)

% setup options (using defaults if any of the options aren't specified)
up.options = setup_options(options);

if up.options.close_figures
    close all
end

%% Analysis settings

% Filter characteristics: Eliminate VLFs (below resp freqs): For 4bpm cutoff
up.paramSet.elim_vlf.Fpass = 0.157;  % in Hz
up.paramSet.elim_vlf.Fstop = 0.02;   % in Hz     (0.157 and 0.02 provide a - 3dB cutoff of 0.0665 Hz)
up.paramSet.elim_vlf.Dpass = 0.05;
up.paramSet.elim_vlf.Dstop = 0.01;

% Filter characteristics: Eliminate VHFs (above frequency content of signals)
up.paramSet.elim_vhf.Fpass = 38.5;  % in HZ
up.paramSet.elim_vhf.Fstop = 33.12;  % in HZ   (33.12 and 38.5 provide a -3 dB cutoff of 35 Hz)
% up.paramSet.elim_vhf.Fpass = 20;  % in HZ
% up.paramSet.elim_vhf.Fstop = 15;  % in HZ
up.paramSet.elim_vhf.Dpass = 0.05;
up.paramSet.elim_vhf.Dstop = 0.01;

% No of times to repeat a single pulse to perform VHF filtering
up.paramSet.no_pulse_repeats = 5;

%% Filter Characteristics

% - Beat detection: eliminating high freqs
up.filtering.beat_detection.elim_high_freqs.Fpass = 38.5;  % in HZ
up.filtering.beat_detection.elim_high_freqs.Fstop = 33.12;  % in HZ   (33.12 and 38.5 provide a -3 dB cutoff of 35 Hz)
up.filtering.beat_detection.elim_high_freqs.Dpass = 0.05;
up.filtering.beat_detection.elim_high_freqs.Dstop = 0.01;

% - Beat detection: eliminating low freqs
up.filtering.beat_detection.elim_low_freqs.Fpass = 0.3;  % in Hz
up.filtering.beat_detection.elim_low_freqs.Fstop = 0.1;  % in Hz
up.filtering.beat_detection.elim_low_freqs.Dpass = 0.05;
up.filtering.beat_detection.elim_low_freqs.Dstop = 0.01; 

% - Fiducial Points: eliminating high freqs
up.filtering.fiducial_points.elim_high_freqs.Fpass = 20;  % in HZ
up.filtering.fiducial_points.elim_high_freqs.Fstop = 15;  % in HZ
up.filtering.fiducial_points.elim_high_freqs.Dpass = 0.05;
up.filtering.fiducial_points.elim_high_freqs.Dstop = 0.01;

% - Fiducial Points: eliminating low freqs
up.filtering.fiducial_points.elim_low_freqs.Fpass = 0.3;  % in Hz
up.filtering.fiducial_points.elim_low_freqs.Fstop = 0.1;  % in Hz
up.filtering.fiducial_points.elim_low_freqs.Dpass = 0.05;
up.filtering.fiducial_points.elim_low_freqs.Dstop = 0.01; 

% - Derivative: number of points to use in Savitzky-Golay filter
up.filtering.derivatives.s_g_filt_len_no_filtering = 5;
up.filtering.derivatives.s_g_filt_len_no_filtered = 9;

%% Number of pulses

% Threshold signal duration to distinguish between a single pulse or multiple pulses: 
up.analysis.max_duration_of_single_pulse = 2.5;   % in secs

% Determine how many pulses there are:
up.analysis.no_of_pulses = determine_no_of_pulses(S, up);

end

function no_of_pulses = determine_no_of_pulses(S, up)
% DETERMINE_NO_OF_PULSES  Determines whether this is a single pulse wave,
% of a pulsatile signal containing multiple pulse waves.

signal_duration = (length(S.v)-1)/S.fs;

% If duration of signal is greater than a threshold then assume this is
% multiple pulse waves:
if signal_duration > up.analysis.max_duration_of_single_pulse
    no_of_pulses = 'multiple';
else
    no_of_pulses = 'single';
end

end

function s_filt = elim_vlfs(s, filt_characteristics, up)
%% Filter pre-processed signal to remove frequencies below resp
% Adapted from RRest

%% Eliminate nans
s.v(isnan(s.v)) = mean(s.v(~isnan(s.v)));

%% Make filter
flag  = 'scale';
[N,Wn,BETA,TYPE] = kaiserord([filt_characteristics.Fstop filt_characteristics.Fpass]/(s.fs/2), [1 0], [filt_characteristics.Dstop filt_characteristics.Dpass]);
b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag);
AMfilter = dfilt.dffir(b);

%% Check frequency response
% % Gives a -3 dB cutoff at ? Hz, using:
% freqz(AMfilter.Numerator)
% norm_cutoff_freq = 0.0266;    % insert freq here from plot
% cutoff_freq = norm_cutoff_freq*(fs/2);

if length(s.v) > (length(AMfilter.numerator)-1)*3
    % - Tukey Window to avoid edge effects
    win_edge_durn = 0.5; % in secs
    prop_win = win_edge_durn*2/((length(s.v)-1)/s.fs);
    tw = tukeywin(length(s.v),prop_win); 
    s_filt.v = filtfilt(AMfilter.numerator, 1, s.v.*tw);
    s_filt.v = s.v-s_filt.v;
else
    s_filt.v = s.v;
end
s_filt.fs = s.fs;
end

function pulses = detect_beats(s, up)

%% Detect peaks in pulse signal
switch up.options.beat_detector
    case 'IMS'
        [pulses.peaks, pulses.onsets, ~] = adaptPulseSegment(s.v,s.fs);
        for wave_no = 1 : length(pulses.peaks)-1
            [~, temp] = max(s.v(pulses.onsets(wave_no):pulses.onsets(wave_no+1)));
            pulses.peaks(wave_no) = temp+pulses.onsets(wave_no)-1;
        end
    case 'Aboy'
        pulses.peaks = abd_algorithm(s);
        pulses.onsets = nan(length(pulses.peaks),1);
        for wave_no = 1 : length(pulses.peaks)-1
            [~, temp] = min(s.v(pulses.peaks(wave_no):pulses.peaks(wave_no+1)));
            pulses.onsets(wave_no) = temp + pulses.peaks(wave_no) - 1;
        end
end

end

function s_filt = elim_vhfs3(s, filt_characteristics)
%% Filter signal to remove VHFs
% Adapted from RRest

s_filt.fs = s.fs;

%% Eliminate nans
s.v(isnan(s.v)) = mean(s.v(~isnan(s.v)));

%% Check to see if sampling freq is at least twice the freq of interest
if (filt_characteristics.Fpass/(s.fs/2)) >= 1
    % then the fs is too low to perform this filtering
    s_filt.v = s.v;
    return
end

%% Create filter
% parameters for the low-pass filter to be used
AMfilter = create_lpf(filt_characteristics, s);

%% Re-make filter if it requires too many samples for this signal
% check to see if it requires too many samples
req_sig_length = 3*(length(AMfilter.numerator)-1);
no_attempts = 0;
while no_attempts<4 && length(s.v)<=req_sig_length
    % change Fpass (i.e. the high frequency end of the filter band)
    filt_characteristics.Fpass = filt_characteristics.Fpass+0.75*(filt_characteristics.Fpass-filt_characteristics.Fstop);
    % re-make filter
    AMfilter = create_lpf(filt_characteristics, s);
    % update criterion
    req_sig_length = 3*(length(AMfilter.numerator)-1);
    % increment number of attempts
    no_attempts = no_attempts+1;
end
if length(s.v)<=req_sig_length
    fprintf('\n - Couldn''t perform high frequency filtering')
end

%% Check frequency response
% Gives a -3 dB cutoff at cutoff_freq Hz, using:
% freqz(AMfilter.Numerator)
% norm_cutoff_freq = INSERT_HERE;    % insert freq here from plot
% cutoff_freq = norm_cutoff_freq*(s.fs/2);

%% Remove VHFs
if length(s.v)>req_sig_length
    s_filt.v = filtfilt(AMfilter.numerator, 1, s.v);
else
    s_filt.v = s.v;
end

end

function AMfilter = create_lpf(filt_characteristics, s)

[N,Wn,BETA,TYPE] = kaiserord([filt_characteristics.Fstop filt_characteristics.Fpass]/(s.fs/2), [1 0], [filt_characteristics.Dstop filt_characteristics.Dpass]);
b  = fir1(N, Wn, TYPE, kaiser(N+1, BETA), 'scale');
AMfilter = dfilt.dffir(b);

end

function [peaks,onsets,clipp] = adaptPulseSegment(y,Fs,annot)
%ADAPTPULSESEGMENT perform adaptive pulse segmentation and artifact detection 
%in ppg signals
%   [peaks,onsets,artif] = adaptPulseSegment(y,annot)
%
% Inputs:
%       y      vector, ppg signal [Lx1] or [1xL], in which L = length(signal)
%       Fs      scalar, sampling frequency in Hz
%       annot   vector, timestamps (in samples) with location of the peaks
%
% Outputs:
%       peaks   vector, locations of peaks (in samples)
%       onsets  vector, locations of onsets of the beats (in samples)
%       artif   vector, locations of peaks classified as artefacts
% 
% References:
%       Karlen et al., Adaptive Pulse Segmentation and Artifact Detection in 
%       Photoplethysmography for Mobile Applications, 34th Internat. Conf. IEEE-EMBS 2012
%       
% Written by Marco A. Pimentel

doOptimise = 1;
doPlot = 0;
if nargin < 3
    % no annotations are provided, therefore, no optimisation will take
    % place
    doOptimise = 0; 
end


% The algorithm in the paper is applied to signals sampled at 125 Hz...
% We do not resample our signal
%ys = resample(y,125,Fs);
%Fs = 125;
% if Fs ~= 300
%     ys = resample(y,300,Fs);
%     Fs = 300;
% else
    ys = y;
% end

% The paper is not clear about the selection of the range of search for m
% We define the range of m to be between [0.005 - 0.100] secs (5ms to 100ms)
% We define "m" in terms of samples
opt.bounds = 0.005:0.005:0.100;
opt.m = unique(ceil(opt.bounds*Fs));

opt.perf = zeros(length(opt.m),4); % store results of performance for each m

if doOptimise
    % Perform optimisation
    for i = 1 : length(opt.m)
        % Determine peaks and beat onsets
        [linez,linezSig] = pulseSegment(ys,Fs,opt.m(i));
        % Calculate performance of the peak detection
        opt.perf(i,:) = evalPerf(annot,linez(:,2));
    end
    
else
    % Do not perform optimization; fix m
    opt.m = 10;
    [peaks,artifs,clipp,linez,linezSig] = pulseSegment(ys,Fs,opt.m);
end

if doPlot
    colData = {'g','y','r'};
    figure; 
    h(1) = subplot(211);
    plot(ys); hold on;
    for i = 1 : size(linez,1)
        %if linezSig(i) > -1
        plot(linez(i,:),ys(linez(i,:)),'-x','Color',colData{linezSig(i)+2});
        %end
    end
    
    h(2) = subplot(212);
    plot(ys,'g'); hold on;
    for i = 1 : size(peaks,1)
        plot(peaks(i,:),ys(peaks(i,:)),'-xr');
    end
    if ~isempty(artifs)
    for i = 1 : size(artifs,1)
        plot(artifs(i,:),ys(artifs(i,:)),'--^b');
    end
    end
    if ~isempty(clipp)
    for i = 1 : size(clipp,1)
        plot(clipp(i,:),ys(clipp(i,:)),'-om');
    end
    end
    linkaxes(h,'x');
    
end

% Correct for the downsmapling performed during the peak detection
onsets = peaks(:,1);
peaks  = peaks(:,2);
for i = 1 : size(peaks,1)
    [~,ind]  = min(ys(max([1 onsets(i)-opt.m]):min([length(ys) onsets(i)+opt.m])));
    onsets(i) = max([1 onsets(i)-opt.m]) + ind(1) - 1;
    [~,ind]  = max(ys(max([1 peaks(i)-opt.m]):min([length(ys) peaks(i)+opt.m])));
    peaks(i) = max([1 peaks(i)-opt.m]) + median(ind) - 1;
end

% Correct minimum value of onset of the beat
for i = 2 : length(onsets)
    [~,ind]   = min(ys(peaks(i-1):peaks(i)));
    onsets(i) = peaks(i-1) + ind - 1;
end

end

function [peaks,artifs,clipp,linez,linezSig] = pulseSegment(ys,Fs,m)
% Perform pulse segmentation in ys given m
% Inputs:
%       ys      vector, with ppg signal
%       m       scalar, with the length of each line (read paper for details)
% 
% Outputs:
%       linez      2-column vector, with location of beat onsets and peaks
%       linezSig   vector, with label of the slope of each line segment
%                  1 - positive slope; -1 - negative slope; 0 - constant
% 

% split signal in different segments
nseg = floor(length(ys)/m);    % number of segments
% nseg = round(length(ys)/m);    % number of segments
% intialize loop variables
seg = 1;    % segment counter
z = 1;      % line counter 
segInLine = 1;  % line controler
linez = zeros(nseg,2); linez(1,:) = [1,m];
% slope of segment/line
a = zeros(nseg,1); a(1) = slope(ys,linez(1,:));
% "classify" line segment according to the slope
linezSig = zeros(nseg,1); linezSig(1) = sign(a(1));
% Start loop over segments
z = z + 1; seg = seg + 1;
for i = 2 : nseg    % loop over segments
    linez(z,:) = [(seg-1)*m+1 seg*m];
    try
        a(z) = slope(ys,linez(z,:));
    catch
        a = 1;
    end
    linezSig(z) = sign(a(z));
    if sign(a(z)) == sign(a(z-1))
        linez(z-1,:) = [(seg-1-segInLine)*m+1 seg*m];
        seg = seg + 1;
        segInLine = segInLine + 1;
    else
        z = z + 1;
        seg = seg + 1;
        segInLine = 1;
    end
end

% remove extra spaces created in output variables
linezSig(sum(linez,2)==0,:) = [];
linez(sum(linez,2)==0,:) = [];

% Apply adaptive threshold algorithm
% For this algorithm to work, we need to first find a valide line segment 
% in order to intialize the thresholds! In order to this, we define a flag
% to control the intialization in the main loop
FOUND_L1 = 0;

% The algorithm includes the definition of 4 adaptation parameters
% We define the following adaptation parameters
% a =     | a_fast_low    a_fast_high |
%         | a_slow_low    a_slow_high |
% 
a = [0.5 1.6; ...
     0.6 2.0];
 
% Define fixed thresholds described in the paper
ThT  = 0.03 * Fs;    % Duration of the line
ThIB = 0.24 * Fs;    % Interbeat invertal (240 ms) 

% Define parameters used in the main loop
alpha = zeros(size(linez,1),1);
for i = 1 : size(linez,1)
    alpha(i) = slope(ys,linez(i,:));   % slopes of line segments
end
theta = diff(ys(linez),[],2);
durat = diff(linez,[],2);       % duration of line segments (in samples)

% remove lines that do not have the necessary duration
linez(durat<ThT,:) = [];
theta(durat<ThT,:) = [];
alpha(durat<ThT,:) = [];
horiz = horizontalLine(ys,linez,Fs);

FLAG = 0;
artifs = []; clipp = [];
% Select window for detect firs peaks!
wind = theta(theta>0);
try 
    wind = wind(1:10);
catch
    wind = wind;
end
ThAlow  = prctile(wind,95)*0.6;
ThAhigh = prctile(wind,95)*1.8;
peaks = [];
for z = 1 : size(linez,1)-1   % loop over line segments
    if FOUND_L1
        if alpha(z) > 0 && ... % slope must be positive
                alpha(z-1) ~= 0 && ...  % peaks before or after clipping are artefactual
                alpha(z+1) ~= 0
            if theta(z) >= ThAlow && theta(z) <= ThAhigh && ...
                    linez(z,2) >= peaks(end,2) + ThIB
                ThAlow  = (ThAlow + theta(z)*a(2,1))/2;
                ThAhigh = theta(z) * a(2,2);
                FLAG = 0;
                currTheta = [currTheta; theta(z)];
                peaks = [peaks; linez(z,:)];
            else
                if FLAG > 0
                    ThAlow  = (ThAlow + min(currTheta(max([1 end-4]):end))*a(1,1))/2;
                    ThAhigh = max(currTheta(max([1 end-4]):end)) * a(1,2);
                    %ThAlow  = (ThAlow + theta(z)*a(1,1))/2;
                    %ThAhigh = theta(z) * a(1,2);
                end
                FLAG = FLAG + 1;
                artifs = [artifs; linez(z,:)];
            end
        elseif theta(z) > 0 && ... 
                ((theta(z-1) ~= 0 || horiz(z-1) ~= 0) && ...
                (theta(z+1) ~= 0 || horiz(z+1) ~= 0))
            if theta(z) >= ThAlow && theta(z) <= ThAhigh && ...
                    linez(z,2) >= peaks(end,2) + ThIB
                ThAlow  = (ThAlow + theta(z)*a(2,1))/2;
                ThAhigh = theta(z) * a(2,2);
                FLAG = 0;
                currTheta = [currTheta; theta(z)];
                peaks = [peaks; linez(z,:)];
            else
                if FLAG > 0
                    %ThAlow  = (ThAlow + currTheta*a(1,1))/2;
                    %ThAhigh = currTheta * a(1,2);
                    ThAlow  = (ThAlow + min(currTheta(max([1 end-4]):end))*a(1,1))/2;
                    ThAhigh = max(currTheta(max([1 end-4]):end)) * a(1,2);
                    %ThAlow  = (ThAlow + theta(z)*a(1,1))/2;
                    %ThAhigh = theta(z) * a(1,2);
                end
                FLAG = FLAG + 1;
                artifs = [artifs; linez(z,:)];
            end
        elseif theta(z) == 0 && horiz(z) == 0
            artifs  = [artifs; linez(z,:)];
            clipp   = [clipp; linez(z,:)];
        end 
    else
        if alpha(z) > 0 && durat(z) >= ThT && ...
                theta(z) >= ThAlow && theta(z) <= ThAhigh 
            FOUND_L1 = 1;
            ThAlow  = theta(z)*0.5;
            ThAhigh = theta(z)*2.0;
            peaks = linez(z,:);    % loaction of onsets and peaks
            currTheta = theta(z);
        end
    end
end

end

function out = horizontalLine(ys,linez,Fs)
% Get horizontal lines from signal given linez
out = zeros(size(linez,1),1);
for i = 1 : size(linez,1)
    out(i) = median(abs(diff(ys(linez(i,1):linez(i,2)))));
    % check duration of the peaks
    if out(i) == 0 && diff(linez(i,:)) <= 0.200*Fs
        out(i) = 0.1;
    end
end

end

function out = slope(ys,interv)
start = interv(1); stop = interv(2);
out = sum(diff(ys([start:stop])))/(stop-start);
%out = median(gradient(ys(start:stop)));
end

function quality = assess_signal_quality(s, pulse_inds)
% ASSESS_SIGNAL_QUALITY  Assesses the signal quality of each beat of the
% pulsatile signal.
% Inputs:
%       s           -  pulsatile signal, a structure containing s.v (a
%                       vector of values), and s.fs (sampling frequency in Hz).
%       pulse_inds  -  indices of the pulse peaks
%
% Outputs:
%       quality     -  the signal quality of each beat (1 indicates high
%                       quality, 0 low quality).
%
% Adapted from RRest
%
% Reference: This function uses an adaptation of the signal quality index
% for the photoplethysmogram described in:
%     Orphanidou, C. et al., 2015. Signal-quality indices for the electrocardiogram and photoplethysmogram: derivation and applications to wireless monitoring. IEEE Journal of Biomedical and Health Informatics, 19(3), pp.8328. Available at: http://www.ncbi.nlm.nih.gov/pubmed/25069129.

%% Setup
s.t = [0:length(s.v)-1]/s.fs;

%% Segment into windows of 10s duration
win_durn = 10;   % in secs
win.deb = s.t(1):(win_durn-2):s.t(end);
win.fin = win.deb + win_durn;

high_quality_pulse_inds = [];
for win_no = 1 : length(win.deb)
    
    % identify data for this window
    
    rel_els = s.t >= win.deb(win_no) & s.t <= win.fin(win_no);
    first_el = find(rel_els,1);
    curr_sig.v = s.v(rel_els); 
    curr_sig.t = s.t(rel_els); clear rel_els
    curr_sig.t = curr_sig.t - curr_sig.t(1);
    curr_sig.pulse_ind_inds = find(s.t(pulse_inds) >= win.deb(win_no) & s.t(pulse_inds) <= win.fin(win_no));
    curr_pulse_inds = pulse_inds(curr_sig.pulse_ind_inds) - first_el + 1;
    
    % find beat-to-beat interval
    
    beat_to_beat_interval = median(diff(curr_sig.t(curr_pulse_inds)));
    beat_to_beat_samples = round(beat_to_beat_interval*s.fs); clear beat_to_beat_interval
    
    % find a template beat
    ts = [];
    rel_els = curr_pulse_inds>beat_to_beat_samples/2 & ...
        curr_pulse_inds+floor(beat_to_beat_samples/2)<length(curr_sig.v);
    rel_pulse_inds = curr_pulse_inds(rel_els);
    curr_sig.used_pulse_ind_inds = curr_sig.pulse_ind_inds(rel_els);
    % find beat morphologies
    for rel_pulse_no = 1 : length(rel_pulse_inds)
        this_pulse = rel_pulse_inds(rel_pulse_no);
        t = curr_sig.v(this_pulse-floor(beat_to_beat_samples/2):this_pulse+floor(beat_to_beat_samples/2));
        tt = t/norm(t); tt = tt(:)';
        ts = [ts; tt]; clear tt t
    end
    clear k l j
    
    % find all templates in current window
    avtempl = mean(ts,1);
    
    % now calculate correlation for every beat in this window
    r2 = nan(size(ts,1),1);
    for k = 1:size(ts,1)
        r2(k) = corr2(avtempl,ts(k,:));
    end
    clear k
    
    high_quality_beats = r2 > 0.86;
    
    high_quality_pulse_inds = [high_quality_pulse_inds, curr_sig.used_pulse_ind_inds(high_quality_beats)];
    
    % % calculate template using only high-quality beats
    % templ = mean(ts(r2>0.86,:),1);  %threshold = 0.66 for ECG, 0.86 for PPG;    % ecg cross-correlation threshold value (for sqi)
    
end

high_quality_pulse_inds = unique(high_quality_pulse_inds);
quality = false(length(pulse_inds),1);
for ind = 1 : length(quality)
    if intersect(high_quality_pulse_inds, ind)
        quality(ind) = true;
    end
end
end

function sigs = filter_signal(S, up)

%% Make output variable
sigs.fs = S.fs;
sigs.orig = S.v;

if up.options.do_filter
    
    %% Filter to remove high frequencies
    filt_characteristics = up.filtering.fiducial_points.elim_high_freqs;
    s_filt = elim_vhfs3(S, filt_characteristics);
    
    %% Filter to remove low frequencies
    if strcmp(up.analysis.no_of_pulses, 'multiple')
        filt_characteristics = up.filtering.fiducial_points.elim_low_freqs;
        s_filt = elim_vlfs(s_filt, filt_characteristics, up);
    else
        % remove baseline wander
        baseline = linspace(0, s_filt.v(end)-s_filt.v(1), length(s_filt.v)); baseline = baseline(:);
        s_filt.v = s_filt.v-baseline;
    end
    
    %% Output
    sigs.filt = s_filt.v;
    
end

end

function sigs = calculate_derivs(sigs, up)

%% Extract relevant signal (either filtered or not)
rel_sig.fs = sigs.fs;
if ~up.options.do_filter
    rel_sig.v = sigs.orig;
    s_g_filter_len = up.filtering.derivatives.s_g_filt_len_no_filtering;
else
    rel_sig.v = sigs.filt;
    s_g_filter_len = up.filtering.derivatives.s_g_filt_len_no_filtered;
end

%% Repeat single pulse wave to avoid edge effects
if strcmp(up.analysis.no_of_pulses, 'single')
    rel_sig.v = repmat(rel_sig.v, [up.paramSet.no_pulse_repeats,1]);
end

%% Calculate derivatives
dt = 1/rel_sig.fs;
sigs.first_d = savitzky_golay(rel_sig.v, 1, 9)./dt;
sigs.second_d = savitzky_golay(sigs.first_d, 1, 9)./dt;
sigs.third_d = savitzky_golay(sigs.second_d, 1, 9)./dt;

%% Retain original pulse wave
if strcmp(up.analysis.no_of_pulses, 'single')
    % select derivative values corresponding to the original pulse
    orig_len = length(rel_sig.v)/up.paramSet.no_pulse_repeats;
    orig_els = (floor(up.paramSet.no_pulse_repeats/2)*orig_len)+1 : ((floor(up.paramSet.no_pulse_repeats/2)+1)*orig_len);
    sigs.first_d = sigs.first_d(orig_els);
    sigs.second_d = sigs.second_d(orig_els);
    sigs.third_d = sigs.third_d(orig_els);
end

end

function pts = identify_fiducial_point_indices(sigs, pulses, up)

%% Identify fiducial point indices

% setup variables to store fiducial point indices
fid_pt_names = {'a', 'b', 'c', 'd', 'e', 'f', 's', 'dia', 'dic', 'p1pk', 'p2pk', 'p1in', 'p2in', 'ms', 'f1', 'f2', 'ms2'};
for fid_pt_no = 1 : length(fid_pt_names)
    eval(['pts.ind.' fid_pt_names{fid_pt_no} ' = nan(length(pulses.onsets)-1,1);'])
    eval(['pts.amp_norm.' fid_pt_names{fid_pt_no} ' = nan(length(pulses.onsets)-1,1);'])
    eval(['pts.amp.' fid_pt_names{fid_pt_no} ' = nan(length(pulses.onsets)-1,1);'])
    eval(['pts.t.' fid_pt_names{fid_pt_no} ' = nan(length(pulses.onsets)-1,1);'])
end

% setup buffer zones
buffer_tol.deb = 0.005; % in secs
buffer_tol.fin = 0.2; % in proportion of beat
buffer_p1 = [0.1, 0.18]; % in secs

% cycle through each pulse wave
for pulse_no = 1 : length(pulses.onsets)-1
    
    % extract data for this pulse wave
    curr_els = pulses.onsets(pulse_no):pulses.onsets(pulse_no+1);
    
    % (either filtered or not)
    if up.options.do_filter
        curr.sig = sigs.filt(curr_els);
    else
        curr.sig = sigs.orig(curr_els);
    end
    curr.derivs.first = sigs.first_d(curr_els);
    curr.derivs.second = sigs.second_d(curr_els);
    curr.derivs.third = sigs.third_d(curr_els);
    
    %% Pre-process
    
    % Make sure the pulse wave ends at the same amplitude as its onset.
    old.v = curr.sig; old.fs = sigs.fs;
    temp = eliminate_low_freq_from_single_beat(old, up);
    
    % Normalise pulse wave
    if up.options.normalise_pw
        
        % Normalise to occupy a range of 0 to 1
        temp.v = (temp.v-min(temp.v))/range(temp.v);
        
        % Ensure that signal commences at start of systole
        %     temp = align_pulse(temp, up);
        
    end
    
    curr.sig = temp.v;
    
    %% Identify fiducial points
    
    % find buffers
    initial_buffer = floor(buffer_tol.deb * sigs.fs);  % in samples
    end_buffer = length(curr.sig) - ceil(buffer_tol.fin * length(curr.sig));  % in samples
    
    % find f1 and f2
    temp_f1 = 1;
    temp_f2 = length(curr.sig);
    
    % find s
    temp_s = identify_s(curr);
    
    % find ms
    temp_ms = identify_ms(curr);
    
    % find a
    temp_a = identify_a(curr, initial_buffer);
    
    if isempty(temp_a)
        continue
    end
    
    % find b
    temp_b = identify_b(curr, temp_a, up);
    
    if isempty(temp_b)
        continue
    end
    
    % find p1
    p1_buffer = floor(buffer_p1 .* sigs.fs);  % in samples
    temp_p1 = identify_p1(curr, temp_b, temp_ms, p1_buffer, up.analysis.no_of_pulses, sigs.fs, up);
    
    % find e
    redo_log = 0;
    temp_e = identify_e(curr, temp_s, temp_ms, temp_b, end_buffer, redo_log);
    
    if isempty(temp_e)
        [temp_f, temp_c, temp_d, temp_dic, temp_dia, temp_p2] = deal(nan);
    else
        
        % find c
        temp_c = identify_c(curr, temp_b, temp_e);
        
        if isempty(temp_c)
            redo_log = 1;
            
            % re-find e
            old_temp_e = temp_e;
            temp_e = identify_e(curr, temp_s, temp_ms, end_buffer, redo_log);
            
            % re-find c
            temp_c = identify_c(curr, temp_b, temp_e);
        end
        
        % find f
        temp_f = identify_f(curr, temp_e, end_buffer);
        
        % find dic
        temp_dic = identify_dic(curr,temp_e);
        
        % find dia
        temp_dia = identify_dia(curr, temp_dic, temp_e, end_buffer, up);
        
        if isempty(temp_c)
            [temp_d, temp_p2] = deal(nan);
        else
            % find d
            temp_d = identify_d(curr, temp_c, temp_e);
            
            % find p2
            temp_p2 = identify_p2(curr, temp_d, temp_p1, temp_dic);
            
        end
        
    end
    
    % retain timings of original p1 and p2 estimates
    temp_p1in = temp_p1; temp_p2in = temp_p2;
    
    % make p1 or p2 coincident with the systolic peak
    [~, rel_el] = min(abs(temp_s-[temp_p1,temp_p2]));
    if rel_el == 1
        temp_p1 = temp_s;
    else
        temp_p2 = temp_s;
    end
    
    if ~isnan(temp_p2) && ~isnan(temp_p1)
        % make sure p2 is at a peak if necessary
        pks = find_pks_trs(curr.sig, 'pk');
        cutoff = mean([temp_p1, temp_p2]);
        possible_pks = find(pks > cutoff & pks < temp_e & curr.sig(pks) > curr.sig(temp_p2));
        if ~isempty(possible_pks)
            [~, temp_el] = max(curr.sig(pks(possible_pks)));
            temp_p2 = pks(possible_pks(temp_el));
        end
        % make sure p1 is at a peak if necessary
        pks = find_pks_trs(curr.sig, 'pk');
        cutoff = mean([temp_p1, temp_p2]);
        possible_pks = find(pks < cutoff & pks > temp_ms & curr.sig(pks) > curr.sig(temp_p1));
        if ~isempty(possible_pks)
            [~, temp_el] = max(curr.sig(pks(possible_pks)));
            temp_p1 = pks(possible_pks(temp_el));
        end
    end
    
    % store p1pk and p2pk
    temp_p1pk = temp_p1;
    temp_p2pk = temp_p2;
    
    % identify ms2
    temp_ms2 = identify_ms2(curr, temp_p1in, temp_p2pk, up);
    
    % calculate scale factor
    curr_range = range(curr.sig);
    scale_factor = 1/curr_range;
    
    % store points
    for fid_pt_no = 1 : length(fid_pt_names)
        eval(['curr_temp_el = temp_' fid_pt_names{fid_pt_no} ';']);
        if ~isnan(curr_temp_el)
            
            % - index of fiducial point
            sig_ind = curr_temp_el + curr_els(1)-1;
            eval(['pts.ind.' fid_pt_names{fid_pt_no} '(pulse_no) = sig_ind;'])
            
            % - amplitude of fiducial point
            if sum(strcmp(fid_pt_names{fid_pt_no}, {'a','b','c','d','e','f'}))
                amp = curr.derivs.second(curr_temp_el);
            elseif sum(strcmp(fid_pt_names{fid_pt_no}, {'s','dia','dic','p1pk','p1in','p2pk','p2in','f1','f2'}))
                amp = curr.sig(curr_temp_el);
            elseif sum(strcmp(fid_pt_names{fid_pt_no}, {'ms','ms2'}))
                amp = curr.derivs.first(curr_temp_el);
            end
            eval(['pts.amp.' fid_pt_names{fid_pt_no} '(pulse_no) = amp;'])
            amp_norm = amp*scale_factor;
            eval(['pts.amp_norm.' fid_pt_names{fid_pt_no} '(pulse_no) = amp_norm;'])
            
            % - timing of fiducial point
            t = (curr_temp_el-1)/sigs.fs;
            eval(['pts.t.' fid_pt_names{fid_pt_no} '(pulse_no) = t;'])
            clear amp_norm sig_ind amp t
        end
    end
    
    
    clear curr temp* curr_els empty_log
    
end
clear pulse_no


%% Ensure only pts are provided for those pulses with all pts available
pt_names = fieldnames(pts.ind);
include_pulse = true(size(pts.ind.dia));
for pt_name_no = 1 : length(pt_names)
    eval(['curr_pt_measures = pts.ind.' pt_names{pt_name_no} ';']);
    include_pulse(isnan(curr_pt_measures)) = false;    
end
for pt_name_no = 1 : length(pt_names)
    if strcmp(pt_names{pt_name_no}, 'f1') || strcmp(pt_names{pt_name_no}, 'f2')
        continue
    end
    eval(['pts.ind.' pt_names{pt_name_no} '(~include_pulse) = nan;']);
    eval(['pts.amp_norm.' pt_names{pt_name_no} '(~include_pulse) = nan;']);
    eval(['pts.amp.' pt_names{pt_name_no} '(~include_pulse) = nan;']);
    eval(['pts.t.' pt_names{pt_name_no} '(~include_pulse) = nan;']);
end

end

function pw_inds = calculate_pw_inds(fid_pts, sigs, pulses, ht, up)

%% Timings

% - Duration of pulse
pw_inds.T = fid_pts.t.f2-fid_pts.t.f1;

% - Instantaneous heart rate (bpm)
pw_inds.IHR = 60./pw_inds.T;

% - Time between systolic and diastolic peaks (secs)
pw_inds.delta_t = fid_pts.t.dia-fid_pts.t.s;

% - Crest time (secs)
pw_inds.CT = fid_pts.t.s-fid_pts.t.f1;

% - Crest time divided by height (s/m)
pw_inds.CT_div_ht = pw_inds.CT./ht;

% - Stiffness index (m/s)
pw_inds.SI = ht./pw_inds.delta_t;

% - Proportion of pulse wave duration which is spent on systolic upslope
pw_inds.prop_s = pw_inds.CT./pw_inds.T;

% - Duration of systole
pw_inds.t_sys = fid_pts.t.dic-fid_pts.t.f1;

% - Duration of diastole
pw_inds.t_dia = fid_pts.t.f2-fid_pts.t.dic;

% - Ratio of systolic to diastolic durations
pw_inds.t_ratio = pw_inds.t_sys./pw_inds.t_dia;

% - Proportion of pulse wave duration taken up by time between systolic and diastolic peaks
pw_inds.prop_delta_t = pw_inds.delta_t./pw_inds.T;

% - Time from p1 to diastolic peak
pw_inds.t_p1in_dia = fid_pts.t.dia-fid_pts.t.p1in;
pw_inds.t_p1pk_dia = fid_pts.t.dia-fid_pts.t.p1pk;

% - Time from p2 to diastolic peak
pw_inds.t_p2in_dia = fid_pts.t.dia-fid_pts.t.p2in;
pw_inds.t_p2pk_dia = fid_pts.t.dia-fid_pts.t.p2pk;

% - Time from 'b' to 'c'
pw_inds.t_b_c = fid_pts.t.c-fid_pts.t.b;

% - Time from 'b' to 'd'
pw_inds.t_b_d = fid_pts.t.d-fid_pts.t.b;

%% Amplitudes

% - Pulse amplitude
pw_inds.pulse_amp = fid_pts.amp.s-fid_pts.amp.f1;
pw_inds.pulse_amp_p1 = fid_pts.amp.p1in-fid_pts.amp.f1;
pw_inds.pulse_amp_p2 = fid_pts.amp.p2in-fid_pts.amp.f1;

% - Augmentation Pressure
pw_inds.AP = fid_pts.amp.p2pk-fid_pts.amp.p1in;

% - Agumentation Index
pw_inds.AI = 100*pw_inds.AP./pw_inds.pulse_amp;

% - Diastolic peak amplitude
pw_inds.dia_amp = fid_pts.amp.dia-fid_pts.amp.f1;

% - Reflection Index (calculated using systolic peak)
pw_inds.RI = pw_inds.dia_amp./pw_inds.pulse_amp;

% - Reflection Index (calculated using p1)
pw_inds.RI_p1 = pw_inds.dia_amp./pw_inds.pulse_amp_p1;

% - Reflection Index (calculated using p2)
pw_inds.RI_p2 = pw_inds.dia_amp./pw_inds.pulse_amp_p2;

% - Ratio of amplitudes of p2 and p1
pw_inds.ratio_p2_p1 = pw_inds.pulse_amp_p2./pw_inds.pulse_amp_p1;

%% Areas

% calculate area indices which require additional pulse wave analysis
for beat_no = 1:length(fid_pts.ind.f1)
    
    % skip if there isn't sufficient information for this pulse wave
    if isnan(fid_pts.ind.s(beat_no)), pw_inds.A1(beat_no) = nan; pw_inds.A2(beat_no) = nan; continue, end
    
    % find baseline of pulse wave (joining initial pulse onset to final pulse onset)
    if up.options.do_filter
        baseline = linspace(sigs.filt(fid_pts.ind.f1(beat_no)), sigs.filt(fid_pts.ind.f2(beat_no)), fid_pts.ind.f2(beat_no) - fid_pts.ind.f1(beat_no)+1); baseline = baseline(:);
    else
        baseline = linspace(sigs.orig(fid_pts.ind.f1(beat_no)), sigs.orig(fid_pts.ind.f2(beat_no)), fid_pts.ind.f2(beat_no) - fid_pts.ind.f1(beat_no)+1); baseline = baseline(:);
    end
    
    % - systolic area
    rel_pts = fid_pts.ind.f1(beat_no) : fid_pts.ind.dic(beat_no);
    baseline_pts = rel_pts - fid_pts.ind.f1(beat_no) + 1;
    if up.options.do_filter
        pw_inds.A1(beat_no) = sum(sigs.filt(rel_pts) - baseline(baseline_pts))/( sigs.fs*pw_inds.pulse_amp(beat_no));
    else
        pw_inds.A1(beat_no) = sum(sigs.orig(rel_pts) - baseline(baseline_pts))/( sigs.fs*pw_inds.pulse_amp(beat_no));
    end
    
    % - diastolic area
    rel_pts = fid_pts.ind.dic(beat_no) : fid_pts.ind.f2(beat_no);
    baseline_pts = rel_pts - fid_pts.ind.f1(beat_no) + 1;
    if up.options.do_filter
        pw_inds.A2(beat_no) = sum(sigs.filt(rel_pts) - baseline(baseline_pts))/(sigs.fs*pw_inds.pulse_amp(beat_no));
    else
        pw_inds.A2(beat_no) = sum(sigs.orig(rel_pts) - baseline(baseline_pts))/(sigs.fs*pw_inds.pulse_amp(beat_no));
    end

end
clear beat_no rel_beats

% - Ratio of diastolic to systolic area (called Inflection point area)
pw_inds.IPA = pw_inds.A2 ./ pw_inds.A1;

%% First Derivative

% - Maximum slope
pw_inds.ms = fid_pts.amp.ms;

% - Maximum slope divided by the pulse amplitude
pw_inds.ms_div_amp = fid_pts.amp.ms./pw_inds.pulse_amp;

%% Second Derivative

% - Amplitude of 'b' relative to 'a'
pw_inds.b_div_a = fid_pts.amp.b./fid_pts.amp.a;

% - Amplitude of 'c' relative to 'a'
pw_inds.c_div_a = fid_pts.amp.c./fid_pts.amp.a;

% - Amplitude of 'd' relative to 'a'
pw_inds.d_div_a = fid_pts.amp.d./fid_pts.amp.a;

% - Amplitude of 'e' relative to 'a'
pw_inds.e_div_a = fid_pts.amp.e./fid_pts.amp.a;

% - Amplitude of 'a' relative to pulse amplitude
pw_inds.a_div_amp = fid_pts.amp.a./pw_inds.pulse_amp;

% - Amplitude of 'b' relative to pulse amplitude
pw_inds.b_div_amp = fid_pts.amp.b./pw_inds.pulse_amp;

% - Amplitude of 'c' relative to pulse amplitude
pw_inds.c_div_amp = fid_pts.amp.c./pw_inds.pulse_amp;

% - Amplitude of 'd' relative to pulse amplitude
pw_inds.d_div_amp = fid_pts.amp.d./pw_inds.pulse_amp;

% - Amplitude of 'e' relative to pulse amplitude
pw_inds.e_div_amp = fid_pts.amp.e./pw_inds.pulse_amp;

% - Ageing index: original
pw_inds.AGI = pw_inds.b_div_a - pw_inds.c_div_a - pw_inds.d_div_a - pw_inds.e_div_a;

% - Ageing index: informal
pw_inds.AGI_inf = pw_inds.b_div_a - pw_inds.e_div_a;

% - Ageing index: modified
pw_inds.AGI_mod = pw_inds.b_div_a - pw_inds.c_div_a - pw_inds.d_div_a;

% Calculate SIs from second derivative slopes
pt1.t = fid_pts.t.b;
pt1.v = fid_pts.amp.b;
pt2.t = fid_pts.t.c;
pt2.v = fid_pts.amp.c;
pw_inds.slope_b_c = ((pt2.v - pt1.v)./fid_pts.amp.a)./(pt2.t-pt1.t);
pt2.t = fid_pts.t.d;
pt2.v = fid_pts.amp.d;
pw_inds.slope_b_d = ((pt2.v - pt1.v)./fid_pts.amp.a)./(pt2.t-pt1.t);

%% Indices calculated from multiple derivatives

% - Ratio of diastolic to systolic area (called Inflection point area) plus d-peak
pw_inds.IPAD = pw_inds.IPA + pw_inds.d_div_a;

% - Stiffness constant
pw_inds.k = fid_pts.amp.s ./ ((fid_pts.amp.s - fid_pts.amp.ms ) ./ pw_inds.pulse_amp );

%% Calculate median values of each pulse wave index (using only high quality pulse waves)

pw_ind_names = fieldnames(pw_inds);
for pw_ind_no = 1 : length(pw_ind_names)
    
    % extract data for this pulse wave index
    curr_pw_ind = pw_ind_names{pw_ind_no};
    eval(['raw_data = pw_inds.' curr_pw_ind ';']);
    
    % If there are multiple pulse waves then find the median value of this index
    if length(raw_data)==1
        mod_data.v = raw_data;
    else
        mod_data.raw = raw_data;
        mod_data.v = nanmedian(raw_data(pulses.quality(1:end-1)));
    end
    
    % store data
    eval(['pw_inds.' curr_pw_ind ' = mod_data;']);
    
    clear mod_data
end

end

function deriv = savitzky_golay(sig, deriv_no, win_size)

%% assign coefficients
% From: https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter#Tables_of_selected_convolution_coefficients
% which are calculated from: A., Gorry (1990). "General least-squares smoothing and differentiation by the convolution (Savitzky?Golay) method". Analytical Chemistry. 62 (6): 570?3. doi:10.1021/ac00205a007.

switch deriv_no
    case 0
        % - smoothing
        switch win_size
            case 5
                coeffs = [-3, 12, 17, 12, -3];
                norm_factor = 35;
            case 7
                coeffs = [-2, 3, 6, 7, 6, 3, -2];
                norm_factor = 21;
            case 9
                coeffs = [-21, 14, 39, 54, 59, 54, 39, 14, -21];
                norm_factor = 231;
            otherwise
                error('Can''t do this window size')
        end
    case 1
        % - first derivative
        switch win_size
            case 5
                coeffs = -2:2;
                norm_factor = 10;
            case 7
                coeffs = -3:3;
                norm_factor = 28;
            case 9
                coeffs = -4:4;
                norm_factor = 60;
            otherwise
                error('Can''t do this window size')
        end
        
    case 2
        % - second derivative
        switch win_size
            case 5
                coeffs = [2,-1,-2,-1,2];
                norm_factor = 7;
            case 7
                coeffs = [5,0,-3,-4,-3,0,5];
                norm_factor = 42;
            case 9
                coeffs = [28,7,-8,-17,-20,-17,-8,7,28];
                norm_factor = 462;
            otherwise
                error('Can''t do this window size')
        end
        
    case 3
        % - third derivative
        switch win_size
            case 5
                coeffs = [-1,2,0,-2,1];
                norm_factor = 2;
            case 7
                coeffs = [-1,1,1,0,-1,-1,1];
                norm_factor = 6;
            case 9
                coeffs = [-14,7,13,9,0,-9,-13,-7,14];
                norm_factor = 198;
            otherwise
                error('Can''t do this window size')
        end
        
    case 4
        % - fourth derivative
        switch win_size
            case 7
                coeffs = [3,-7,1,6,1,-7,3];
                norm_factor = 11;
            case 9 
                coeffs = [14,-21,-11,9,18,9,-11,-21,14];
                norm_factor = 143;
            otherwise
                error('Can''t do this window size')
        end
        
    otherwise
        error('Can''t do this order of derivative')        
end

if rem(deriv_no, 2) == 1
    coeffs = -1*coeffs;
end

A = [1,0];
filtered_sig = filter(coeffs, A, sig);
s=length(sig);
half_win_size = floor(win_size*0.5);
deriv=[filtered_sig(win_size)*ones(half_win_size,1);filtered_sig(win_size:s);filtered_sig(s)*ones(half_win_size,1)];
deriv = deriv/norm_factor;

end

function temp_s = identify_s(curr)

[~, temp_s] = max(curr.sig);

end

function temp_a = identify_a(curr, initial_buffer)

[~,filtered_ms] = max(curr.derivs.first);

pks = find_pks_trs(curr.derivs.second, 'pk');
rel_pks = pks(pks > initial_buffer & pks < filtered_ms); 
if isempty(rel_pks) && sum(pks<=initial_buffer) > 0   % Added in PulseAnalyse5
    rel_pks = pks(find(pks <= initial_buffer, 1, 'last'));
end
[~, temp_el] = max(curr.derivs.second(rel_pks));
temp_a = rel_pks(temp_el);

end

function temp_e = identify_e(curr, temp_s, temp_ms, temp_b, end_buffer, redo_log)

% Find local maxima in the second derivative
pks = find_pks_trs(curr.derivs.second, 'pk');
% Set an upper bound of 60 % of the PW duration
upper_bound = 0.6*length(curr.sig);   % const from: https://en.wikipedia.org/wiki/QT_interval#/media/File:QT_interval_corrected_for_heart_rate.png
% Set a lower bound of 'ms'
lower_bound = temp_ms;
% Identify the highest local maximum in the second derivative between these two bounds
rel_pks = pks(pks >= lower_bound & pks <= upper_bound);
[~, max_el] = max(curr.derivs.second(rel_pks));
% If this is the first local maximum in this search region ...
if max_el == 1
    % work out whether this has detected the "c" wave
    % - find out if there's an inflection point between "b" and this
    temp_trs = find_pks_trs(curr.derivs.third, 'tr');
    no_infl = sum(temp_trs > temp_b & temp_trs < rel_pks(max_el));
    % - if not then take the next peak
    if no_infl == 0
        % if there is 1 peak in this search region ...
        if length(rel_pks) < max_el+1   % Added in PulseAnalyse5
            % then take the next peak (slightly above the upper bound
            orig_el = find(pks >= lower_bound & pks <= upper_bound);
            rel_pks = pks(orig_el:orig_el+1);
        end
        rel_pk = rel_pks(max_el+1);
    else
        rel_pk = rel_pks(max_el);
    end
else
    rel_pk = rel_pks(max_el);
end
temp_e = rel_pk;

end

function temp_f = identify_f(curr, temp_e, end_buffer)

lower_bound = temp_e;
upper_bound = end_buffer;
trs = find_pks_trs(curr.derivs.second, 'tr');
possible_els = trs(trs >= lower_bound & trs <= upper_bound);

if isempty(possible_els)    % Added in PulseAnalyse5
    possible_els = trs(find(trs >=lower_bound, 1));
end

if isempty(possible_els)
    temp_f = nan;
else
    temp_f = possible_els(1);
end

end

function temp_b = identify_b(curr, temp_a, up)

% find b (PulseAnalyse6 onwards)

% Find local minima in second derivative
trs = find_pks_trs(curr.derivs.second, 'tr');
% define an upper bound as 25% of the duration of the signal
upper_bound = 0.25*length(curr.sig);
% find local minima between 'a' and this upper bound
temp = find(trs > temp_a & curr.derivs.second(trs) < 0 & trs < upper_bound);
% Identify the lowest of these local minima
[~, rel_el] = min(curr.derivs.second(trs(temp)));
temp = temp(rel_el);
temp_b = trs(temp); clear temp

if isempty(temp_b)
    if up.options.verbose, fprintf('\n - ''b'' outside range'); end
    
    % find local minima after 'a'
    temp = find(trs > temp_a & curr.derivs.second(trs) < 0, 1, 'first');
    temp_b = trs(temp); clear temp

end

end

function temp_d = identify_d(curr, temp_c, temp_e)

% Identify "d" as the lowest minimum of the second deriv between "c" and "e"
trs = find_pks_trs(curr.derivs.second, 'tr');
possible_trs = find(trs > temp_c & trs < temp_e);
if ~isempty(possible_trs)
    temp = trs(possible_trs);
    [~, temp_el] = min(curr.derivs.second(temp));
    temp_d = temp(temp_el); clear temp
else
    % unless there isn't a minimum, in which case it's an inflection, and
    % "d" is the same as "c"
    temp_d = temp_c;
end

end

function temp_c = identify_c(curr, temp_b, temp_e)

% Identify C as the highest local maximum on the second derivative between "b" and "e"
pks = find_pks_trs(curr.derivs.second, 'pk');
temp = find(pks > temp_b & pks < temp_e);
[~, rel_tr_el] = max(curr.derivs.second(pks(temp)));
temp_c = pks(temp(rel_tr_el)); clear temp rel_tr_el pks

% If there aren't any peaks that satisfy this criterion ...
if isempty(temp_c)
    % then identify C as the lowest local minimum on the third derivative
    % after "b" and before "e"
    trs = find_pks_trs(curr.derivs.third, 'tr');
    temp = find(trs > temp_b & trs < temp_e);
    [~, rel_tr_el] = min(curr.derivs.third(trs(temp)));
    if ~isempty(rel_tr_el)
        temp_c = trs(temp(rel_tr_el)); clear temp rel_tr_el trs
    end
end

end

function temp_dic = identify_dic(curr,temp_e)

temp_dic = temp_e;

end

function temp_dia = identify_dia(curr, temp_dic, temp_e, end_buffer, up)

% if there is a diastolic peak, then use that:
%  -  first peak in signal after "dic"
pks = find_pks_trs(curr.sig, 'pks');
temp_dia = pks(find(pks > temp_dic & pks < end_buffer, 1));

% if not, then ...
% % I tried (i) take first peak on first derivative after "e"
if isempty(temp_dia)
    pks = find_pks_trs(curr.derivs.first, 'pks');
    temp_dia = pks(find(pks > temp_e & pks < end_buffer, 1));
end
% % But the problem is that the first derivative isn't necessarily a peak at
% % the diastolic peak - it can be an inflection point. So:
% (ii) the soonest out of (a) first peak on first derivative after "e"
%                         (b) first min on third derivative after "e"
% if isempty(temp_dia)
%     pks = find_pks_trs(curr.derivs.first, 'pks');
%     temp_dia1 = pks(find(pks > temp_e, 1));
%     trs = find_pks_trs(curr.derivs.third, 'trs');
%     temp_dia2 = trs(find(trs > temp_e, 1));
%     temp_dia = min([temp_dia1, temp_dia2]);
% end

if isempty(temp_dia)
    if up.options.verbose, fprintf('\n - ''dia'' outside range'); end
    pks = find_pks_trs(curr.derivs.first, 'pks');
    temp_dia = pks(find(pks > temp_e, 1, 'first'));
end   

end

function temp_ms = identify_ms(curr)

% find max slope in DPPG
[~, temp_ms] = max(curr.derivs.first);

end

function temp_ms2 = identify_ms2(curr, temp_p1in, temp_p2_pk, up)

% Identify peaks in DPPG
d_pks = find_pks_trs(curr.derivs.first, 'pk');

% See whether any are between p1in and p2pk
rel_pks = find(d_pks > temp_p1in & d_pks < temp_p2_pk);
if ~isempty(rel_pks)
    [~, rel_pk_el] = max(d_pks(rel_pks));
    temp_ms2 = d_pks(rel_pks(rel_pk_el));
else
    temp_ms2 = round(mean([temp_p1in, temp_p2_pk]));
    if up.options.verbose, fprintf("\n Fix this"); end
end

end

function temp_p1 = identify_p1(curr, temp_b, temp_ms, buffer_p1, no_of_pulses, fs, up)

% find p1

% find local minima in the first derivative
fd_trs = find_pks_trs(curr.derivs.first, 'tr');
% find local maxima in the second derivative
sd_pks = find_pks_trs(curr.derivs.second, 'pk');

% Find the first local minimum in the first derivative after 0.1 s
current_buffer = buffer_p1(1);
temp = find(fd_trs > current_buffer,1);
% Find the second local minimum (or failing that, the first) in the first derivative after 'b'
temp2 = find(fd_trs > temp_b, 2);
if length(temp2) > 1
    temp2 = temp2(2);
end
% Take whichever comes first:
if temp2 < temp
    temp = temp2;
end
temp_p1 = fd_trs(temp);

% If this value for p1 is after the buffer of 0.18 s ...
if temp_p1 > buffer_p1(2)
    curr.derivs.fourth = savitzky_golay(curr.derivs.third, 1, 9);
    % Then find the first local minimum in the fourth derivative after 0.1 s
    fd_trs = find_pks_trs(curr.derivs.fourth, 'tr');
    temp = find(fd_trs > current_buffer,1);
    temp_p1 = fd_trs(temp); clear temp
end

% If this value for p1 is after the buffer of 0.18 s ...
if temp_p1 > buffer_p1(2)
    if up.options.verbose, fprintf('\n - P1 outside range'); end
    % Then find the last local minimum in the first derivative before 0.18 s
    temp_p1 = fd_trs(find(fd_trs <= current_buffer,1,'last'));
   
    % If this doesn't find temp_p1, then extend the buffer
    if isempty(temp_p1)
        temp_p1 = fd_trs(find(fd_trs > current_buffer,1,'first'));
    end
end



end

function temp_p2 = identify_p2(curr, temp_d, temp_p1, temp_dic)

% Find "p2" from the minimum value of the third derivative immediately before "d"
td_trs = find_pks_trs(curr.derivs.third, 'tr');
temp = find(td_trs < temp_d,1,'last'); clear d_pks
temp_p2 = td_trs(temp);

% unless c=d, in which case p2 will now be before p1, so take the minimum
% value of the third derivative immediately after "d"
if temp_p2 < temp_p1
    temp_p2 = td_trs(find(td_trs<temp_dic,1,'last'));
end

% check to see if there is a maximum in the signal between the current
% estimate of p2, and temp_dic. If so, take this instead
pks = find_pks_trs(curr.sig, 'pk');
temp = find(pks> temp_p2 & pks < temp_dic);
if length(temp) == 1
    temp_p2 = pks(temp);
elseif length(temp) == 2
    temp_p2 = pks(temp(2));
elseif length(temp) > 1
    fprintf('\nCheck this')
end
clear pks

end

function trs = find_pks_trs(sig,type)

if strcmp(type, 'tr')
    sig = -sig;
end

temp1 = find(sig(2:end-1) > sig(1:(end-2)) & sig(2:end-1) > sig(3:end) );

temp2 = find(sig(2:end-2) > sig(1:(end-3)) & sig(2:end-2)== sig(3:(end-1)) & sig(3:end-1) > sig(4:end) );

temp = unique([temp1; temp2]);

trs = temp+1;

end

function S_elim_lf = eliminate_low_freq_from_single_beat(sig, up)

old = sig.v;

% Correct for low frequency baseline drift in a single beat
correction_line = linspace(sig.v(1), sig.v(end), length(sig.v));
S_elim_lf.v = sig.v - correction_line' + sig.v(1);
S_elim_lf.fs = sig.fs;

end

function S_aligned = align_pulse(sig, up)

% Ensure that signal commences at start of systole

[~, align_el] = min(sig.v);
S_aligned.v = sig.v([align_el:end, 1:(align_el-1)]);
S_aligned.fs = sig.fs;

% add on one additional point so that you can define a second onset
S_aligned.v(end+1) = S_aligned.v(1);

end

function make_plots(sigs, pulse_no, fid_pts, up, plot_inds)
% make plot of individual beat if needed

%% - setup
paper_size = 1.5*[500,1050];
if ~up.options.plot_third_deriv, paper_size(2) = paper_size(2)*0.75; end
if up.options.plot_pw_only, paper_size(2) = paper_size(2)*0.37; end
figure('Position', [50,50, paper_size])
ftsize = 1.5*12; lwidth = 2;
sigs.t = [0:length(sigs.v)-1]/sigs.fs;

y_offset = 0.08;
if up.options.plot_third_deriv
    y_inc = 0.23;
    n_sub = 4;
elseif up.options.plot_pw_only
    y_inc = 0.8;
    n_sub = 1;
    y_offset = 0.17;
else
    y_inc = 0.31;
    n_sub = 3;
end

%% - plot sig
h_b(1) = subplot('Position', [0.21,y_offset+(n_sub-1)*y_inc,0.78,y_inc-0.01]);

% plot baseline curve
plot(sigs.t, sigs.v, 'b', 'LineWidth', lwidth); hold on,

if plot_inds && up.options.plot_areas
    
    % Areas: Systolic
    el = fid_pts.ind.dic(pulse_no) - fid_pts.ind.f1(pulse_no)+1;
    rel_els = 1:el-1;
    offset = 0.05*range(sigs.v);
    h = fill([sigs.t(rel_els), sigs.t(rel_els(end))], [sigs.v(rel_els); sigs.v(rel_els(1))], [1,0.8,0.8]);
    h.LineStyle = 'none';
    
    % Areas: Diastolic
    el2 = fid_pts.ind.f2(pulse_no) - fid_pts.ind.f1(pulse_no)+1;
    rel_els = el+1:el2;
    offset = 0.05*range(sigs.v);
    h = fill([sigs.t(rel_els), sigs.t(rel_els(1))], [sigs.v(rel_els); sigs.v(rel_els(end))], [0.8,0.8,1]);
    h.LineStyle = 'none';
    plot(sigs.t, sigs.v, 'b', 'LineWidth', lwidth)
    
end

% plot salient points
pt_names = {'dia', 'dic', 's'}; %, 'p1in', 'p2in'};
hspace0 = 0;
for pt_no = 1 : length(pt_names)
    eval(['curr_pt.el = fid_pts.ind.' pt_names{pt_no} '(pulse_no) - fid_pts.ind.f1(pulse_no)+1;']);
    if isnan(curr_pt.el)
        continue
    end
    
    % plot point
    curr_pt.v = sigs.v(curr_pt.el);
    curr_pt.t = sigs.t(curr_pt.el);
    plot(curr_pt.t, curr_pt.v, 'or')
    
    % annotate point
    curr_text = pt_names{pt_no};
    vspace0 = 0.12*range(sigs.v);
    switch curr_text
        case {'dia','p2'}
            hspace0 = 0.04;
            if (strcmp(curr_text, 'p2pk') || strcmp(curr_text, 'p2in')) && curr_pt.el < (fid_pts.s(pulse_no) - fid_pts.f1(pulse_no)+1)
                hspace0 = -0.04;
            elseif (strcmp(curr_text, 'p1pk') || strcmp(curr_text, 'p1in'))
                hspace0 = 0;
            end
        case 'dic'
            hspace0 = -0.01;
            vspace0 = -1*vspace0;
        case {'f1', 'p1','p1in'}
            hspace0 = -0.04;
        case 's'
            hspace0 = 0;
    end
    text(sigs.t(curr_pt.el)+hspace0, sigs.v(curr_pt.el)+vspace0 , curr_text,'FontSize', ftsize, 'Color', 'r', 'HorizontalAlignment', 'center');
    eval(['sigs.pts.' curr_text ' = curr_pt.el;']);
end

% set limits
curr_range = range(sigs.v);
if plot_inds
    factor1 = 0.4; factor2 = 0.2;
else
    factor1 = 0.15; factor2 = 0.15;
end
ylims = [min(sigs.v)-factor2*curr_range, max(sigs.v)+factor1*curr_range];

ylim(ylims)
curr_range = range(sigs.t);
xlims = [sigs.t(1)-0.05*curr_range, sigs.t(end)+0.1*curr_range];
xlim_labs = 0:0.25:0.25*ceil(sigs.t(end)*4);
xlim(xlims)

% set labels
ylab = ylabel({'Pulse', 'Wave'}, 'FontSize', ftsize, 'Rotation', 0);
set(ylab, 'Units', 'Normalized', 'Position', [-0.13, 0.5, 0]);
set(gca, 'FontSize', ftsize -2, 'XTick', xlim_labs, 'XGrid', 'on')
if ~up.options.normalise_pw
    ytick_vals = [round(min(sigs.v),3,'significant'), round(max(sigs.v),3,'significant')];
else
    ytick_vals = [];
end
set(gca, 'YTick', ytick_vals)
box off
if up.options.plot_pw_only
    xlabel('Time (s)', 'FontSize', ftsize)
else
    set(gca, 'XTickLabel', {})
end

% Plot indices
if plot_inds
    
    color = 0.4;
    
    % - delta T
    ht_annot = ylims(2)-0.11*range(ylims);
    plot(sigs.t(sigs.pts.s)*[1,1], [ht_annot, sigs.v(sigs.pts.s)], '--', 'color', color*ones(1,3))
    plot(sigs.t(sigs.pts.dia)*[1,1], [ht_annot, sigs.v(sigs.pts.dia)], '--', 'color', color*ones(1,3))
    normalised1  = coords_to_pos(sigs.t(sigs.pts.s), ht_annot);
    normalised2  = coords_to_pos(sigs.t(sigs.pts.dia), ht_annot);
    ah = annotation('doublearrow',[normalised1(1),normalised2(1)],[normalised1(2),normalised2(2)]);
    ah.Head1Length = 7; ah.Head2Length = 7;
    ah.Head1Width = 7; ah.Head2Width = 7;
    new_ht_annot= ht_annot +0.01*range(ylims);
    text(mean([sigs.t(sigs.pts.s),sigs.t(sigs.pts.dia)]), new_ht_annot, '\DeltaT','FontSize', ftsize, 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    
    % - Crest Time
    ht_annot = ylims(2)-0.11*range(ylims);
    plot(sigs.t(1)*[1,1], [ht_annot, sigs.v(1)], '--', 'color', color*ones(1,3))
    plot(sigs.t(sigs.pts.s)*[1,1], [ht_annot, sigs.v(sigs.pts.s)], '--', 'color', color*ones(1,3))
    normalised1  = coords_to_pos(sigs.t(1), ht_annot);
    normalised2  = coords_to_pos(sigs.t(sigs.pts.s), ht_annot);
    ah = annotation('doublearrow',[normalised1(1),normalised2(1)],[normalised1(2),normalised2(2)]);
    ah.Head1Length = 7; ah.Head2Length = 7;
    ah.Head1Width = 7; ah.Head2Width = 7;
    new_ht_annot= ht_annot +0.01*range(ylims);
    text(mean([sigs.t(1),sigs.t(sigs.pts.s)]), new_ht_annot, 'CT','FontSize', ftsize, 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    
    % - Reflection Index
    t_annot = sigs.t(end)+0.05;
    plot([sigs.t(end), t_annot], [1,1]*sigs.v(end), '--', 'color', color*ones(1,3))
    plot([sigs.t(sigs.pts.dia), t_annot], [1,1]*sigs.v(sigs.pts.dia), '--', 'color', color*ones(1,3))
    normalised1  = coords_to_pos(t_annot, sigs.v(end));
    normalised2  = coords_to_pos(t_annot, sigs.v(sigs.pts.dia));
    ah = annotation('doublearrow',[normalised1(1),normalised2(1)],[normalised1(2),normalised2(2)]);
    ah.Head1Length = 7; ah.Head2Length = 7;
    ah.Head1Width = 7; ah.Head2Width = 7;
    new_ht_annot= mean([sigs.v(sigs.pts.dia), sigs.v(end)]);
    text(sigs.t(end), new_ht_annot, 'RI','FontSize', ftsize, 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    
    % - Systolic time
    ht_annot = ylims(1)+0.02*range(ylims);
    plot(sigs.t(1)*[1,1], [ht_annot, sigs.v(1)], '--', 'color', color*ones(1,3))
    plot(sigs.t(sigs.pts.dic)*[1,1], [ht_annot, sigs.v(sigs.pts.dic)], '--', 'color', color*ones(1,3))
    normalised1  = coords_to_pos(sigs.t(1), ht_annot);
    normalised2  = coords_to_pos(sigs.t(sigs.pts.dic), ht_annot);
    ah = annotation('doublearrow',[normalised1(1),normalised2(1)],[normalised1(2),normalised2(2)]);
    ah.Head1Length = 7; ah.Head2Length = 7;
    ah.Head1Width = 7; ah.Head2Width = 7;
    new_ht_annot= ht_annot +0.01*range(ylims);
    text(mean([sigs.t(1),sigs.t(sigs.pts.dic)]), new_ht_annot, 'Systole','FontSize', ftsize, 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    
    % - Diastolic time
    ht_annot = ylims(1)+0.02*range(ylims);
    plot(sigs.t(end)*[1,1], [ht_annot, sigs.v(1)], '--', 'color', color*ones(1,3))
    plot(sigs.t(sigs.pts.dic)*[1,1], [ht_annot, sigs.v(sigs.pts.dic)], '--', 'color', color*ones(1,3))
    normalised1  = coords_to_pos(sigs.t(sigs.pts.dic), ht_annot);
    normalised2  = coords_to_pos(sigs.t(end), ht_annot);
    ah = annotation('doublearrow',[normalised1(1),normalised2(1)],[normalised1(2),normalised2(2)]);
    ah.Head1Length = 7; ah.Head2Length = 7;
    ah.Head1Width = 7; ah.Head2Width = 7;
    new_ht_annot= ht_annot +0.01*range(ylims);
    text(mean([sigs.t(sigs.pts.dic),sigs.t(end)]), new_ht_annot, 'Diastole','FontSize', ftsize, 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    
    if up.options.plot_areas
        % - Systolic area
        ht_annot = mean([sigs.v(1), sigs.v(sigs.pts.s)]);
        text(mean([sigs.t(sigs.pts.dic),sigs.t(1)]), ht_annot, 'A_s','FontSize', ftsize, 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        
        % - Diastolic area
        ht_annot = mean([sigs.v(1), sigs.v(sigs.pts.s)]);
        text(0.71*sigs.t(sigs.pts.s)+0.29*sigs.t(end), ht_annot, 'A_d','FontSize', ftsize, 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
    
end

if ~up.options.plot_pw_only
    
    %% - plot first derivative
    h_b(2) = subplot('Position', [0.21,y_offset+(n_sub-2)*y_inc,0.78,y_inc - 0.01]);
    
    % Plot x-axis
    plot([-10, 10], [0,0], '--k'); hold on
    
    % plot baseline curve
    plot(sigs.t, sigs.first_d, 'b', 'LineWidth', lwidth); hold on,
    
    % plot salient points
    pt_names = {'ms'}; %, 'ms2', 'dia'};
    curr_range = range(sigs.first_d);
    vspace0 = 0.1*curr_range;
    hspace0 = 0.05;
    for pt_no = 1 : length(pt_names)
        eval(['curr_fid_pt = fid_pts.ind.' pt_names{pt_no} '(pulse_no);'])
        if isnan(curr_fid_pt)
            continue
        end
        curr_pt.el = curr_fid_pt - fid_pts.ind.f1(pulse_no)+1;
        if isnan(curr_pt.el)
            continue
        end
        
        % plot point
        curr_pt.v = sigs.first_d(curr_pt.el);
        curr_pt.t = sigs.t(curr_pt.el);
        plot(curr_pt.t, curr_pt.v, 'or')
        
        % annotate point
        curr_text = pt_names{pt_no};
        text(sigs.t(curr_pt.el)+hspace0, sigs.first_d(curr_pt.el)+vspace0 , curr_text,'FontSize', ftsize, 'Color', 'r', 'HorizontalAlignment', 'center');
        eval(['sigs.pts.' curr_text ' = curr_pt.el;']);
        
    end
    
    % set limits
    curr_range = range(sigs.t);
    xlim(xlims)
    curr_range = range(sigs.first_d);
    ylim([min(sigs.first_d)-0.05*curr_range, max(sigs.first_d)+0.15*curr_range]); ylims = ylim;
    
    % set labels
    ylab = ylabel({'1st','derivative'}, 'FontSize', ftsize, 'Rotation', 0);
    set(ylab, 'Units', 'Normalized', 'Position', [-0.13, 0.5, 0]);
    set(gca, 'FontSize', ftsize -2, 'YTick', [])
    set(gca, 'XTick', xlim_labs, 'XTickLabel', {}, 'XGrid', 'on')
    if ~up.options.normalise_pw
        ytick_vals = [0, round(max(sigs.first_d),3,'significant')];
        set(gca, 'YTick', ytick_vals)
    else
        set(gca, 'YTick', 0);
    end
    box off
    
    % Plot indices
    if plot_inds
        
        % - ms
        ht_annot = mean([sigs.first_d(sigs.pts.ms), min(sigs.first_d)]);
        [~, temp] = min(sigs.first_d(sigs.pts.ms:sigs.pts.dic));
        el = temp + sigs.pts.ms -1;
        plot([sigs.t(sigs.pts.ms), sigs.t(el)], sigs.first_d(sigs.pts.ms)*[1,1], '--', 'color', color*ones(1,3))
        %     plot([sigs.t(sigs.pts.ms), sigs.t(el)], sigs.first_d(el)*[1,1], '--', 'color', color*ones(1,3))
        normalised1  = coords_to_pos(sigs.t(el), 0);
        normalised2  = coords_to_pos(sigs.t(el), sigs.first_d(sigs.pts.ms));
        ah = annotation('doublearrow',[normalised1(1),normalised2(1)],[normalised1(2),normalised2(2)]);
        curr_range = range(sigs.t);
        text(sigs.t(el)+0.06*curr_range, ht_annot, 'ms','FontSize', ftsize, 'Color', 'k', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
        
    end
    
    %% - plot Second derivative
    h_b(3) = subplot('Position', [0.21,y_offset+(n_sub-3)*y_inc,0.78,y_inc - 0.01]);
    
    % Plot x-axis
    plot([-10, 10], [0,0], '--k'); hold on
    
    % plot baseline curve
    curr_color = 0.4*0.5*[1,2,1];
    plot(sigs.t, sigs.second_d, 'b', 'LineWidth', lwidth); hold on,
    
    % plot salient points
    pt_names = {'a', 'b', 'c', 'd', 'e', 'f'};
    curr_range = range(sigs.second_d);
    vspace_const = 0.08;
    hspace0 = 0;
    for pt_no = 1 : length(pt_names)
        eval(['curr_pt.el = fid_pts.ind.' pt_names{pt_no} '(pulse_no) - fid_pts.ind.f1(pulse_no)+1;']);
        if isnan(curr_pt.el)
            continue
        end
        
        % plot point
        curr_pt.v = sigs.second_d(curr_pt.el);
        curr_pt.t = sigs.t(curr_pt.el);
        plot(curr_pt.t, curr_pt.v, 'or')
        
        % annotate point
        curr_text = pt_names{pt_no};
        switch curr_text
            case {'a','c', 'e'}
                vspace0 = (1.3*vspace_const)*curr_range;
            case {'b', 'd', 'f'}
                vspace0 = -1.3*vspace_const*curr_range;
        end
        text(sigs.t(curr_pt.el), sigs.second_d(curr_pt.el)+vspace0 , curr_text,'FontSize', ftsize, 'Color', 'r', 'HorizontalAlignment', 'center');
        eval(['sigs.pts.' curr_text ' = curr_pt.el;']);
    end
    
    % set limits
    curr_range = range(sigs.t);
    xlim(xlims)
    curr_range = range(sigs.second_d);
    ylim([min(sigs.second_d)-0.15*curr_range, max(sigs.second_d)+0.15*curr_range])
    
    % set labels
    ylab = ylabel({'2nd','derivative'}, 'FontSize', ftsize, 'Rotation', 0);
    set(ylab, 'Units', 'Normalized', 'Position', [-0.13, 0.5, 0]);
    set(gca, 'FontSize', ftsize -2, 'YTick', [], 'XTick', xlim_labs, 'XGrid', 'on')
    box off
    if ~up.options.plot_third_deriv
        xlabel('Time (s)', 'FontSize', ftsize)
    else
        set(gca, 'XTickLabel', {})
    end
    if ~up.options.normalise_pw
        ytick_vals = [round(min(sigs.second_d),3,'significant'), 0, round(max(sigs.second_d),3,'significant')];
        set(gca, 'YTick', ytick_vals)
    else
        set(gca, 'YTick', 0);
    end
    
    % Plot indices
    if plot_inds
        
        % - slope b-d
        plot(sigs.t([sigs.pts.b,sigs.pts.d]), sigs.second_d([sigs.pts.b, sigs.pts.d]), 'k', 'LineWidth', lwidth),
        curr_range = range(sigs.t);
        x = mean(sigs.t([sigs.pts.b,sigs.pts.d]));
        y = mean(sigs.second_d([sigs.pts.b,sigs.pts.b,sigs.pts.b,sigs.pts.d]));
        text(x+0.03*curr_range, y, 'slope_{b-d}','FontSize', ftsize, 'Color', 'k', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
        
    end
    
    %% - plot Third derivative
    if up.options.plot_third_deriv
        h_b(4) = subplot('Position', [0.21,y_offset+(n_sub-4)*y_inc,0.78,y_inc-0.01]);
        
        % Plot x-axis
        plot([-10, 10], [0,0], '--k'); hold on
        
        % plot baseline curve
        curr_color = 0.4*0.5*[1,2,1];
        plot(sigs.t, sigs.third_d, 'b', 'LineWidth', lwidth); hold on,
        
        % plot salient points
        pt_names = {'p1pk', 'p2pk'};
        curr_range = range(sigs.third_d);
        vspace_const = 0.08;
        hspace0 = 0;
        for pt_no = 1 : length(pt_names)
            eval(['curr_pt.el = fid_pts.ind.' pt_names{pt_no} '(pulse_no) - fid_pts.ind.f1(pulse_no)+1;']);
            if isnan(curr_pt.el)
                continue
            end
            
            % plot point
            curr_pt.v = sigs.third_d(curr_pt.el);
            curr_pt.t = sigs.t(curr_pt.el);
            plot(curr_pt.t, curr_pt.v, 'or')
            
            % annotate point
            curr_text = pt_names{pt_no};
            switch curr_text
                case {'p1pk'}
                    vspace0 = vspace_const*curr_range;
                case {'p2pk'}
                    vspace0 = -1*vspace_const*curr_range;
            end
            text(sigs.t(curr_pt.el), sigs.third_d(curr_pt.el)+vspace0 , curr_text,'FontSize', ftsize, 'Color', 'r', 'HorizontalAlignment', 'center');
            
        end
        
        % set limits
        curr_range = range(sigs.t);
        xlim(xlims)
        ylims = ylim; curr_range = range(sigs.third_d);
        ylim([min(sigs.third_d)-0.05*curr_range, max(sigs.third_d)+0.15*curr_range])
        
        % set labels
        xlabel('Time (s)', 'FontSize', ftsize)
        ylab = ylabel({'3rd','derivative'}, 'FontSize', ftsize, 'Rotation', 0);
        set(ylab, 'Units', 'Normalized', 'Position', [-0.13, 0.5, 0]);
        set(gca, 'FontSize', ftsize -2, 'YTick', [], 'XTick', xlim_labs, 'XGrid', 'on')
        box off
    end
    
end

linkaxes(h_b, 'x')
shg

if ~isempty(up.options.save_folder)
    savepath = [up.options.save_folder, up.options.save_file];
    if ~isempty(up.options.save_file)
        savepath = [savepath, '_'];
    end
    if ~plot_inds
        savepath = [savepath, 'FidPts'];
    else
        savepath = [savepath, 'PWInds'];
    end
    set(gcf,'color','w');
    set(gcf,'PaperUnits','centimeters');
    set(gcf,'PaperSize', [paper_size(1), paper_size(2)]./40);
    set(gcf,'PaperPosition',[0 0 paper_size(1) paper_size(2)]./40);
%     print(gcf,'-dpdf',savepath)
    print(gcf,'-depsc',savepath)
    print(gcf,'-dpng',savepath)
end
shg

end

function normalised  = coords_to_pos(x_coord, y_coord)

pos = get(gca, 'Position');
normalised(1) = (x_coord - min(xlim))/diff(xlim) * pos(3) + pos(1);
normalised(2) = (y_coord - min(ylim))/diff(ylim) * pos(4) + pos(2);
 
end

function subtracted = subtract_baseline(sig)

baseline = linspace(sig(1),sig(end),length(sig));

subtracted = sig(:) - baseline(:); 

end

function [norm, scale_factor] = normalise(sig)

norm = sig - min(sig); 
scale_factor = max(norm);
norm = norm / scale_factor;

end

function plot_results(sigs, pulses, fid_pts, pw_inds, up)

if ~up.options.do_plot || sum(isnan(fid_pts.ind.a)) == length(fid_pts.ind.a)
    return
end

% extract data for the first high quality pulse wave
pulse_no = 1;
if ~pulses.quality(pulse_no)
    temp = find(pulses.quality(1:end-1) & ~isnan(fid_pts.ind.a), 1);
    if ~isempty(temp)
        pulse_no = temp;
    end
end
curr_els = pulses.onsets(pulse_no):pulses.onsets(pulse_no+1);

% (either filtered or not)
if up.options.do_filter
    sigs.v = sigs.filt(curr_els);
else
    sigs.v = sigs.orig(curr_els);
end
sigs.first_d = sigs.first_d(curr_els);
sigs.second_d = sigs.second_d(curr_els);
sigs.third_d = sigs.third_d(curr_els);

% - subtract baseline from this pulse wave, and normalise (if specified)
if up.options.normalise_pw
    sigs.v = subtract_baseline(sigs.v);
    [sigs.v, scale_factor] = normalise(sigs.v);
    sigs.first_d = sigs.first_d./scale_factor;
    sigs.second_d = sigs.second_d./scale_factor;
    sigs.third_d = sigs.third_d./scale_factor;
end

% - Plot fiducial points
make_plots(sigs, pulse_no, fid_pts, up, 0)

% - Plot pulse wave indices
make_plots(sigs, pulse_no, fid_pts, up, 1)

end
