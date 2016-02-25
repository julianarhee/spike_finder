

[session_dir] = ...
    uigetdir('/media/juliana/TOSHIBA EXT/ephys', 'Select a session...');

x = inputdlg('How many tetrodes?','Tetrdes',1);
numTetrodes = str2num(x{:});
tetrodes = 1:1:numTetrodes;

if find(session_dir==' ')
    str = ['dir ''', session_dir ''''];
    f = {evalc(str)};
    f = regexp(f, '\n', 'split');
    files=[];
    for i=1:length(f{1})
        if ~isempty(strfind(f{1}{i}, '.rhd'))
            files = [files f{1}(i)];
        end
    end
else
    files = dir([session_dir '*.rhd']);
end

files = strtrim(files)';

D=struct();
for i=1:length(files)
    fparts = regexp(files{i}, '_', 'split');
    D(i).fname = fparts{1};
    %D(i).cond = fparts{2};
    D(i).current = str2num(fparts{2}(1:end-2));
    D(i).frequency = str2num(fparts{3}(1:end-2));
end

% GROUPING MAY DEPEND ON FREQUENCY OR CURRENT:
freqs = unique([D(:).frequency]);
currents = unique([D(:).current]);

groupby='frequency';

if strcmp(groupby, 'frequency')
    grouper = freqs;
    units = 'Hz';
elseif strcmp(groupby, 'current')
    grouper = currents;
    units = 'uA';
end

F = struct();
for i=1:length(grouper)
    curr_files = [];
    for j=1:length(D)
        if D(j).frequency == grouper(i)
            curr_files = [curr_files files(j)];
        end
    end
    F(i).fnames = curr_files;
end

path = session_dir;

for i=1:length(F)
    
    file = F(i).fnames(1);
    read_Intan_RHD2000_file_session(path, file)
    n_traces = size(amplifier_data,1);
    
    % do all the stuff in get_average_to_stim.m:
    fs = frequency_parameters.amplifier_sample_rate;            % recording sample rate
    stim_freq = freqs(i)                                             % frequency of stim (Hz)
    chunk_size = 60;     

%     % mk figure/output dir based on F groups (e.g., frequency):
%     outdir = fullfile(session_dir, [num2str(grouper(i)) units]);
%     
%     if ~isdir(outdir)
%         mkdir(outdir)
%     end
    
    % make the figures & save them to appropriate directory
    
    
    % SHOW AVG FOR EACH CHANNEL:
    
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
            sprintf('Stim-freq is %d.\nFound different stim t-step:\n%s',...
                stim_freq, num2str(unique(t_step(t_step ~= (1/stim_freq)))))
            stim_freq = round(1/mean(t_step));
            sprintf('Changing stim-freq to actual value: %d', stim_freq)
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
    for j=1:n_traces;

        % Filter the raw data:
        [b,a]=butter(4,[500/(fs/2)],'high');

        unfiltered_data = amplifier_data(j,:);
        filtered_data = filtfilt(b,a,unfiltered_data);

        tmp_dig_in = zeros(length(on_idxs)-2, win);
        tmp_amp = zeros(length(on_idxs)-2, win);
        % Grab chunk of data that corresponds to stim onset:
        for k=2:length(on_idxs)-1
            idx = on_idxs(k);                             % Indices of stim onset highs
            tmp_dig_in(k-1,:) = board_dig_in_data(idx-n_pre:idx+n_post);
            tmp_amp(k-1,:) = filtered_data(idx-n_pre:idx+n_post);
        end

        dig_in(j,:) = mean(tmp_dig_in, 1);
        amp(j,:) = mean(tmp_amp, 1);

    %     figure()
    %     subplot(2,1,1); 
    %     plot(t_dig(idx-n_pre:idx+n_post), board_dig_in_data(idx-n_pre:idx+n_post));
    %     title(sprintf('Stimulus onset, -%d and +%d samples', n_pre, n_post));
    %     hold on;
    %     subplot(2,1,2);
    %     plot(t_amplifier(idx-n_pre:idx+n_post), filtered_data(idx-n_pre:idx+n_post))

        subplot(4,4,j)
        plot(t_amplifier(idx-n_pre:idx+n_post), filtered_data(idx-n_pre:idx+n_post))
        hold on

    end
    
    % mk figure/output dir based on F groups (e.g., frequency):
    %outdir = fullfile(session_dir, [num2str(grouper(i)) units]);
    if strcmp(groupby, 'frequency') && grouper(i)~=stim_freq
        outdir = fullfile(session_dir, [num2str(stim_freq) units]);
    else
        outdir = fullfile(session_dir, [num2str(grouper(i)) units]);
    end
    
    if ~isdir(outdir)
        mkdir(outdir)
    end

    % save fig:
    figname = sprintf('%s_%s_each_channel', [num2str(D(i).current) 'uA'], [num2str(D(i).frequency) 'Hz']);
    outfile = strcat(outdir, '/', figname)
    savefig(fullfile(outdir, figname))
    
    
    % Average across tetrodes:
    tet = struct();
    stim = struct();
    last_idx = 1;
    strt = 1;
    for j=1:length(tetrodes)
        curr_tet = ['tetrode',num2str(j)];
        tet.(curr_tet) = mean(amp(strt:strt+3, :), 1);
        dig.(curr_tet) = mean(dig_in(strt:strt+3, :), 1);
        stim.on.(curr_tet) = find(mean(dig_in(strt:strt+3, :), 1), 1, 'first');
        stim.off.(curr_tet) = find(mean(dig_in(strt:strt+3, :), 1), 1, 'last');
        strt = strt+4;
    end

    % Plot tetrode averages, 0 is stim onset:
    figure()
    for j=1:length(tetrodes)

        % Plot filtered average:
        curr_tet = ['tetrode',num2str(j)];
        subplot(2,2,j)
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
    
    % save fig:
    figname = sprintf('%s_%s_each_tetrode', [num2str(D(i).current) 'uA'], [num2str(D(i).frequency) 'Hz']);
    savefig(fullfile(outdir, figname))
    
end