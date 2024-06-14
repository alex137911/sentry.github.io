% completeAlgorithm("C:\Users\acale\OneDrive\Documents\Waterloo BME\4B\BME 462\CapstoneProjectBase-physician_ui\data\mar13-2024\mar13\")

function completeAlgorithm(fpath) %fpath needs to have the backslash at the end
%loading and converting files
try
    eeg_file = dir(fullfile(strcat(fpath,'EEG'),'*.csv'));
    eeg_raw = readtable(fullfile(eeg_file.folder,eeg_file.name));
    eeg_raw = table2array(eeg_raw);
catch
    error('No EEG data is available')
end
try
    rw_file = dir(fullfile(strcat(fpath,'RightWrist'),'*.csv'));
    rw_raw = readtable(fullfile(rw_file.folder,rw_file.name));
    rw_raw = table2array(rw_raw);
    rw_recorded = true;
catch
    rw_recorded = false;
end
try
    lw_file = dir(fullfile(strcat(fpath,'LeftWrist'),'*.csv'));
    lw_raw = readtable(fullfile(lw_file.folder,lw_file.name));
    lw_raw = table2array(lw_raw);
    lw_recorded = true;
catch
    lw_recorded = false;
end
try
    ra_file = dir(fullfile(strcat(fpath,'RightAnkle'),'*.csv'));
    ra_raw = readtable(fullfile(ra_file.folder,ra_file.name));
    ra_raw = table2array(ra_raw);
    ra_recorded = true;
catch
    ra_recorded = false;
end
try
    la_file = dir(fullfile(strcat(fpath,'LeftAnkle'),'*.csv'));
    la_raw = readtable(fullfile(la_file.folder,la_file.name));
    la_raw = table2array(la_raw);
    la_recorded = true;
catch
    la_recorded = false;
end
try
    ml_model= load('trained_classifier.mat');
catch
    error('Classifier is unavailable');
end
if ~rw_recorded && ~lw_recorded && ~ra_recorded && ~la_recorded
    error('No IMU data is available')
end
if ~rw_recorded || ~lw_recorded
    error('Insufficient IMU data was recorded')
end

%recording parameters
eeg_sfreq = 256;
imu_sfreq = 100;

%filter eeg into frequency bands
eeg_filt = eeg_raw(:,2:end);
eeg_filt = lowpass(eeg_filt, 44, eeg_sfreq); %TODO: think about combining this into a bandpass filter
eeg_filt = highpass(eeg_filt, 0.5, eeg_sfreq);
delta = bandpass(eeg_filt, [0.5 4], eeg_sfreq);
theta = bandpass(eeg_filt, [4 8], eeg_sfreq);
alpha = bandpass(eeg_filt, [8 13], eeg_sfreq);
beta = bandpass(eeg_filt, [13 30], eeg_sfreq);
gamma = bandpass(eeg_filt, [30 44], eeg_sfreq);
%write to file
eeg_fnout = split(eeg_file.name,'.');
eeg_fnout = strcat(eeg_fnout{1},'_output.csv');
writematrix([eeg_raw(:,1) sum(delta,2) sum(theta,2) sum(alpha,2) sum(beta,2) sum(gamma,2)],fullfile(eeg_file.folder,eeg_fnout));
%continue processing eeg for detection
eeg_filt = eeg_filt(:,2);
delta = delta(:,2);
theta = theta(:,2);
alpha = alpha(:,2);
beta = beta(:,2);
gamma = gamma(:,2);
dar = delta./alpha;
dtr = delta./theta;
dtabr = (delta + theta)./(alpha+beta);

%extract features from MUSE data
freq_extractor = signalFrequencyFeatureExtractor("FrameSize", 7680, "FrameOverlapLength", 0, "SampleRate", eeg_sfreq, 'MeanFrequency', true, 'MedianFrequency', true, 'PeakLocation', true, 'Bandpower', true);
pow_extractor = signalFrequencyFeatureExtractor("FrameSize", 7680, "FrameOverlapLength", 0, "SampleRate", eeg_sfreq, "Bandpower", true);
delta_feats = extract(freq_extractor, delta);
theta_feats = extract(freq_extractor, theta);
alpha_feats = extract(freq_extractor, alpha);
beta_feats = extract(freq_extractor, beta);
gamma_feats = extract(freq_extractor, gamma);
dar_feats = extract(pow_extractor, dar);
dtr_feats = extract(pow_extractor, dtr);
dtabr_feats = extract(pow_extractor, dtabr);
global_feats = extract(pow_extractor, eeg_filt);

%TODO add a brief comment describing what this is
N = length(eeg_filt);
eeg_filt = reshape(eeg_filt(1:N-mod(N, 7680)), 7680, floor(N/7680));
delta = reshape(delta(1:N-mod(N, 7680)), 7680, floor(N/7680));
theta = reshape(theta(1:N-mod(N, 7680)), 7680, floor(N/7680));
alpha = reshape(alpha(1:N-mod(N, 7680)), 7680, floor(N/7680));
beta = reshape(beta(1:N-mod(N, 7680)), 7680, floor(N/7680));
gamma = reshape(gamma(1:N-mod(N, 7680)), 7680, floor(N/7680));
dar = reshape(dar(1:N-mod(N, 7680)), 7680, floor(N/7680));
dtr = reshape(dtr(1:N-mod(N, 7680)), 7680, floor(N/7680));

%TODO add a brief comment describing what this is
[psdx, freq] = pspectrum(delta, eeg_sfreq);                    
IntSpectrum = cumtrapz(freq, psdx); 
SEF_delta = zeros(1,length(psdx(1,:)));
for n = 1:length(psdx(1,:))
    [IntSpectrum_a, ia, ic] = unique(IntSpectrum(:,n));
    freq_t = freq(ia);
    SEF_delta(n) = interp1(IntSpectrum_a, freq_t, 0.95*IntSpectrum_a(end), 'linear');   
end
[psdx, freq] = pspectrum(theta, eeg_sfreq);                    
IntSpectrum = cumtrapz(freq, psdx); 
SEF_theta = zeros(1,length(psdx(1,:)));
for n = 1:length(psdx(1,:))
    [IntSpectrum_a, ia, ic] = unique(IntSpectrum(:,n));
    freq_t = freq(ia);
    SEF_theta(n) = interp1(IntSpectrum_a, freq_t, 0.95*IntSpectrum_a(end), 'linear');   
end
[psdx, freq] = pspectrum(alpha, eeg_sfreq);                    
IntSpectrum = cumtrapz(freq, psdx); 
SEF_alpha = zeros(1,length(psdx(1,:)));
for n = 1:length(psdx(1,:))
    [IntSpectrum_a, ia, ic] = unique(IntSpectrum(:,n));
    freq_t = freq(ia);
    SEF_alpha(n) = interp1(IntSpectrum_a, freq_t, 0.95*IntSpectrum_a(end), 'linear');   
end
[psdx, freq] = pspectrum(beta, eeg_sfreq);                    
IntSpectrum = cumtrapz(freq, psdx); 
SEF_beta = zeros(1,length(psdx(1,:)));
for n = 1:length(psdx(1,:))
    [IntSpectrum_a, ia, ic] = unique(IntSpectrum(:,n));
    freq_t = freq(ia);
    SEF_beta(n) = interp1(IntSpectrum_a, freq_t, 0.95*IntSpectrum_a(end), 'linear');   
end
[psdx, freq] = pspectrum(gamma, eeg_sfreq);                    
IntSpectrum = cumtrapz(freq, psdx); 
SEF_gamma = zeros(1,length(psdx(1,:)));
for n = 1:length(psdx(1,:))
    [IntSpectrum_a, ia, ic] = unique(IntSpectrum(:,n));
    freq_t = freq(ia);
    SEF_gamma(n) = interp1(IntSpectrum_a, freq_t, 0.95*IntSpectrum_a(end), 'linear');   
end

%TODO add a brief comment describing what this is
features = [delta_feats SEF_delta' theta_feats SEF_theta' alpha_feats SEF_alpha' beta_feats SEF_beta' gamma_feats SEF_gamma' dar_feats dtr_feats dtabr_feats global_feats];
features(end,:) = [];
features = table(features);

%predict sleep stages
[yfit,scores] = ml_model.trained_classifier.predictFcn(features);

%extratct periods of REM sleep
recording_start_time = eeg_raw(1,1);
rem_stamps = [];
rem_start_time = 0;
if yfit(1) =='REM'
    rem_stamps(end+1) = recording_start_time;
    rem_start_time = 1;
end
for i = 2:length(yfit)
    if yfit(i) == 'REM' && rem_start_time == 0
        rem_start_time = i;
        rem_stamps(end+1) = recording_start_time + i*30;
    end
    if yfit(i) ~= 'REM' && i>1 && yfit(i-1) == 'REM'
        rem_stamps(end+1) = recording_start_time + i*30;
        rem_start_time = 0;
    end
end
if isempty(rem_stamps)
    error('No REM sleep was detected in the recording')
end
if rem(length(rem_stamps),2) ~= 0
    error('Detected REM periods are missing an opening or closing timestamp')
end

%flag movement epoches of concern
if rw_recorded
    [rw_filt,rw_flagged] = AccelerationFlagging(rw_raw,rem_stamps,imu_sfreq);
    rw_fnout = split(rw_file.name,'.');
    rw_fnout = strcat(rw_fnout{1},'_output.csv');
    writematrix([rw_filt rw_flagged],fullfile(rw_file.folder,rw_fnout))
end
if lw_recorded
    [lw_filt,lw_flagged] = AccelerationFlagging(lw_raw,rem_stamps,imu_sfreq);
    lw_fnout = split(lw_file.name,'.');
    lw_fnout = strcat(lw_fnout{1},'_output.csv');
    writematrix([lw_filt lw_flagged],fullfile(lw_file.folder,lw_fnout))
end
if ra_recorded
    [ra_filt,ra_flagged] = AccelerationFlagging(ra_raw,rem_stamps,imu_sfreq);
    ra_fnout = split(ra_file.name,'.');
    ra_fnout = strcat(ra_fnout{1}, '_output.csv');
    writematrix([ra_filt ra_flagged],fullfile(ra_file.folder,ra_fnout))
end
if la_recorded
    [la_filt,la_flagged] = AccelerationFlagging(la_raw,rem_stamps,imu_sfreq);
    la_fnout = split(la_file.name,'.');
    la_fnout = strcat(la_fnout{1},'_output.csv');
    writematrix([la_filt la_flagged],fullfile(la_file.folder,la_fnout))
end



function [filtered,flagged] = AccelerationFlagging(raw,remStamps,sfreq)
    %SVM and correct for gravity
    filtered = [raw(:,1) abs(sqrt(raw(:,2).^2 + raw(:,3).^2 + raw(:,4).^2))];
    filtered(:,2) = filtered(:,2) - 1;

    %filter for frequencies of human movement
    [b,a] = butter(2,[0.5 20]/(sfreq/2),"bandpass");
    filtered(:,2) = filter(b,a,filtered(:,2));
    filtered(:,2) = abs(filtered(:,2));

    flagged = zeros(length(raw),1);
    for i = 1:length(remStamps)
        if rem(i,2) == 0
            continue
        end
        remAccel = filtered(filtered(:,1)>=remStamps(i),:);
        remAccel = remAccel(remAccel(:,1)<=remStamps(i+1),:);
        %TODO need to check that there is enough IMU data in the REM period
        sfreq = round(sfreq);
        idx1s = 1:sfreq:length(remAccel);
        backgroundAve = mean(remAccel(:,2));
        remf = remAccel;
        remf(remAccel(:,2)<0.1,2) = 0;
         for j = 1:length(idx1s)
            if j < length(idx1s)
                flagged(idx1s(j):idx1s(j+1)) = logical(remf(idx1s(j):idx1s(j+1),2) >= 2*backgroundAve);
            else
                flagged(idx1s(j):end) = logical(remf(idx1s(j):end,2) >= 2*backgroundAve);
            end
         end

        start_flagging = find(flagged(remAccel(:,1)) == 1);
        cur = 1;
        for j = 2:length(start_flagging)-1
            if start-flagging(j)-start_flagging(j-1) > 1
                curlen = j-1 - cur;
                if start_flagging(cur) - round(0.2*curlen) > 0
                    flagged(start_flagging(cur)-round(0.2*curlen):start_flagging(cur)) = 1;
                else
                    flagged(1:start_flagging(cur)) = 1;
                end
                if start_flagging(cur) + curlen + round(0.2*curlen) <= length(flagged)
                    flagged(start_flagging(cur):start_flagging(cur)+curlen+round(0.2*curlen)) = 1;
                else
                    flagged(startflagging(cur):end) = 1;
                end
                cur = j;
            end
        end
        for j = 1:length(start_flagging)-1
            if start_flagging(j+1) - start_flagging(j) > 50
                flagged(start_flagging(j):start_flagging(j+1)) = 1;
            end
        end
    end
end
end
