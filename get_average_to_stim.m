

% Experiment parameters:
fs = frequency_parameters.amplifier_sample_rate;            % recording sample rate
stim_freq = 20;                                             % frequency of stim (Hz)
chunk_size = 60;                                            % time chunk of read_Intan parsing (sec)

% Read raw .rhd file and store into base workspace:
read_Intan_RHD2000_file;  

% Check correct stimulation frequency:
fs_stim = frequency_parameters.board_dig_in_sample_rate;    % Digital IN sample rate

diffT = diff(board_dig_in_data);                            % Finds stim onsets and offsets (-1, 0, or 1)
exp_stim_freq = round(length(find(diffT==1))/chunk_size);   % Actual stim frequency (N of highs)           
on_idxs = find(diffT==1)+1;                                 % Corr. indices for recording data
off_idxs = find(diffT==-1)+1;                               % sometimes -1 is ON, and v/v...

t_n = min([length(on_idxs) length(off_idxs)]);              % sometimes n offsets ~= n onsets
n_samples_on = abs(min(off_idxs(1:t_n)-on_idxs(1:t_n)));    % N samples stim is ON


error = 0.0001;
t_step = diff(on_idxs) / fs_stim;                           % time steps b/w samples (sec)...
if any(t_step ~= (1/stim_freq))                             % ...should be equal to 1/stim_freq
    
    if any((t_step(t_step ~= (1/stim_freq)) - 1/stim_freq) > error)
        display('Stimulus frequency incorrect!')
        sprintf('Stim-freq is %d.\nFound different stim t-step:\n%s', stim_freq, num2str(unique(t_step(t_step ~= (1/stim_freq)))))
    end
    
end

n_samples_per_cycle = 1/stim_freq * fs;                 % N samples between each stim ON
n_samples_off = n_samples_per_cycle - n_samples_on;

n_pre = round(n_samples_off * 0.5);                                            % N samples to look before stim ON
n_post = round(n_samples_on + n_pre);                            % N samples to look after stim ON
win = n_pre + n_post + 1;                               % Window size of data to look at
dig_in = zeros(n_traces, win);
amp = zeros(n_traces, win);
figure()
for i=1:n_traces;
    
    % Filter the raw data:
    [b,a]=butter(4,[500/(fs/2)],'high');

    unfiltered_data = amplifier_data(i,:);
    filtered_data = filtfilt(b,a,unfiltered_data);
    
    tmp_dig_in = zeros(length(on_idxs)-2, win);
    tmp_amp = zeros(length(on_idxs)-2, win);
    % Grab chunk of data that corresponds to stim onset:
    for j=2:length(on_idxs)-1
        idx = on_idxs(j);                             % Indices of stim onset highs
        tmp_dig_in(j-1,:) = board_dig_in_data(idx-n_pre:idx+n_post);
        tmp_amp(j-1,:) = filtered_data(idx-n_pre:idx+n_post);
    end

    dig_in(i,:) = mean(tmp_dig_in, 1);
    amp(i,:) = mean(tmp_amp, 1);

%     figure()
%     subplot(2,1,1); 
%     plot(t_dig(idx-n_pre:idx+n_post), board_dig_in_data(idx-n_pre:idx+n_post));
%     title(sprintf('Stimulus onset, -%d and +%d samples', n_pre, n_post));
%     hold on;
%     subplot(2,1,2);
%     plot(t_amplifier(idx-n_pre:idx+n_post), filtered_data(idx-n_pre:idx+n_post))

    subplot(4,4,i)
    plot(t_amplifier(idx-n_pre:idx+n_post), filtered_data(idx-n_pre:idx+n_post))
    hold on

end

% Average across tetrodes:
tet = struct();
stim = struct();
last_idx = 1;
strt = 1;
for i=1:length(tetrodes)
    curr_tet = ['tetrode',num2str(i)];
    tet.(curr_tet) = mean(amp(strt:strt+3, :), 1);
    dig.(curr_tet) = mean(dig_in(strt:strt+3, :), 1);
    stim.on.(curr_tet) = find(mean(dig_in(strt:strt+3, :), 1), 1, 'first');
    stim.off.(curr_tet) = find(mean(dig_in(strt:strt+3, :), 1), 1, 'last');
    strt = strt+4;
end

% Plot tetrode averages, 0 is stim onset:
figure()
for i=1:length(tetrodes)
    
    % Plot filtered average:
    curr_tet = ['tetrode',num2str(i)];
    subplot(2,2,i)
    x = -1*n_pre:n_post; %t_amplifier(idx-n_pre:idx+n_post);
    plot(x, tet.(curr_tet))
    hold on
    plot(x, dig.(curr_tet))
    
    % Add vertical lines to show stim ON:
    ylim = get(gca,'ylim');  %Get x range 
    hold on
    plot([x(stim.on.(curr_tet)) x(stim.on.(curr_tet))], [ylim(1) ylim(2)],'k')
    
    hold on
    plot([x(stim.off.(curr_tet)) x(stim.off.(curr_tet))], [ylim(1) ylim(2)],'k')
    
end
