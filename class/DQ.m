classdef DQ
    %These code are developed by Dongqi
    
    properties
        Property1
    end
    
    methods(Static)
        
        function label(dataDir)
            % dataDir = "C:\Users\v-dongqihan\Downloads\autoCollectedSyllables\autoCollectedSyllables\";
            dataDir = "C:\Users\Zhehao\Dropbox (OIST)\My_Stimuli\Eleinff\";
            
            fg = figure ;
            
            listing = dir(dataDir);
            
            % -------------- FT options --------------
            wlen = 300;
            hop = 75;
            nfft = 300;
            
            fs = 32000;
            
            anal_win = blackmanharris(wlen, 'periodic');
            synth_win = hamming(wlen, 'periodic');
            
            % ------------------------------------
            
            curr_spectrogram = subplot(2, 2, 1);
            box on
            
            skip_10 = subplot(6, 2, 2);
            text(0, 0, 'Skip 10 samples', 'HorizontalAlignment', 'center')
            xlim([-5, 5])
            ylim([-5, 5])
            box on
            axis off
            
            not_a_syllable = subplot(6, 2, 2 + 2);
            text(0, 0, 'This is not a syllable', 'HorizontalAlignment', 'center')
            xlim([-5, 5])
            ylim([-5, 5])
            box on
            axis off
            
            back_to_previous = subplot(6, 2, 2+4);
            text(0, 0, 'previous label wrong, go back', 'HorizontalAlignment', 'center')
            xlim([-5, 5])
            ylim([-5, 5])
            box on
            axis off
            
            
            max_num_syllables  = 14;
            
            for i =  1:max_num_syllables
                genres(i) = subplot(4, floor(max_num_syllables / 2), max_num_syllables + i);
                box on
            end
            
            % --------------------------------------------
            
            for file_id = 3:length(listing)
                dataFile = listing(file_id).name;
                
                if ~strcmp(dataFile(max(1, end-3):end), '.mat' ) ||  strcmp(dataFile(1:min(7, length(dataFile))), 'labeled' )
                    continue
                end
                
                if exist(dataDir + "labeled_" + dataFile, 'file')
                    fprintf("%s already processed, skip\n", dataFile)
                    continue
                end
                
                dataPath = dataDir + dataFile;
                
                data = load(dataPath);
                
                labeled_data = data;
                
                labeled_data.wlen = wlen;
                labeled_data.hop = hop;
                labeled_data.nfft = nfft;
                labeled_data.fs = fs;
                
                for i =  1:max_num_syllables
                    labeled_data.label_example_syllable{i} = 0;
                    
                    % ------- in case you did it wrongly at, it can regret ------
                    prev_Ss_softplused{i} = [[nan, nan]; [nan, nan]];
                    prev_Ts{i} = [0, 1];
                    prev_Fs{i} = [0, 1];
                    
                    prev_prev_Ss_softplused{i} = [[nan, nan]; [nan, nan]];
                    prev_prev_Ts{i} = [0, 1];
                    prev_prev_Fs{i} = [0, 1];
                    
                end
                
                n_s = 0;
                
                % modified by Zhehao Cheng to fit sylinf/eleinf data format
                if isfield(data,'sylinf')
                    data.syllables = data.sylinf;
                elseif isfield(data,'eleinf')
                    data.syllables = data.eleinf;
                end
                
                
                while n_s < length(data.syllables)
                    
                    n_s = n_s + 1;
                    
                    syllable = data.syllables(n_s);
                    y = syllable.y / max(abs(syllable.y));
                    
                    [STFT, F, T] = DQ.stft(y, anal_win, hop, nfft, fs);
                    
                    S_softplused = bdsoftplus(abs(STFT));
                    labeled_data.syllables(n_s).S_softplused = S_softplused;
                    labeled_data.syllables(n_s).STFT = STFT;
                    labeled_data.syllables(n_s).T = T;
                    labeled_data.syllables(n_s).F = F;
                    
                    contourf(curr_spectrogram, T, F,  S_softplused, 20, 'LineStyle', 'none')
                    sound(syllable.y(1: min(fs, length(syllable.y))), fs) % play the sound for at longest 1 s
                    title(curr_spectrogram, sprintf("%s: %d / %d", dataFile, n_s, length(data.syllables)))
                    xlabel(curr_spectrogram, 'time')
                    ylabel(curr_spectrogram, 'frequency')
                    
                    previous_selected_axis = 'Hello, World';
                    
                    k = waitforbuttonpress;
                    if k == 0
                        if gca == not_a_syllable
                            labeled_data.syllables(n_s).label = -1;
                            continue
                        elseif gca == back_to_previous
                            if n_s <= 1
                                disp("Cannot go back because this is the first sample!")
                                n_s = 0;
                                continue
                            end
                            label = labeled_data.syllables(n_s - 1).label;
                            if label >= 0
                                contourf(genres(label), prev_prev_Ts{label}, prev_prev_Fs{label},  prev_prev_Ss_softplused{label}, 20, 'LineStyle', 'none')
                            end
                            n_s = n_s - 2;
                            continue
                        elseif gca == skip_10
                            for kk = 0 : 9
                                if n_s + kk <= length(data.syllables)
                                    labeled_data.syllables(n_s + kk).label = -1;
                                end
                            end
                            n_s = n_s + 9;
                            if n_s >= length(data.syllables)  % already end
                                save(dataDir + "labeled_" + dataFile, "labeled_data")
                            end
                            continue
                        else
                            for i = 1: max_num_syllables + 1
                                if i == max_num_syllables + 1
                                    disp("not correctly select a label!! Do again!")
                                    n_s = n_s - 1;
                                    continue
                                end
                                if genres(i) == gca
                                    labeled_data.syllables(n_s).label = i;
                                    break
                                end
                            end
                        end
                    end
                    
                    contourf(gca, T, F,  S_softplused, 20, 'LineStyle', 'none')
                    axis off
                    
                    if labeled_data.syllables(n_s).label > 0
                        
                        label = labeled_data.syllables(n_s).label;
                        labeled_data.label_example_syllable{label} = syllable;
                        
                        prev_prev_Ss_softplused{label} = prev_Ss_softplused{label};
                        prev_prev_Ts{label} = prev_Ts{label} ;
                        prev_prev_Fs{label} = prev_Fs{label} ;
                        
                        prev_Ss_softplused{label} = S_softplused;
                        prev_Ts{label} = T;
                        prev_Fs{label} = F;
                        
                    end
                    
                end
                
                fprintf("SAVING ------- %s \n", dataDir + "labeled_" + dataFile)
                save(dataDir + "labeled_" + dataFile, "labeled_data")
                
                fprintf("SAVED:  %s \n", dataDir + "labeled_" + dataFile)
                
            end
            
            
            
            function y = bdsoftplus(x)
                
                y = (1 - 2*(x < 0)).* log(1 + (1 - 2*(x < 0)) .* x);
                
            end
            
            
            
            function x = inversebdsoftplus(y)
                
                x = (1 - 2*(y < 0)).* exp(y) - (1 - 2*(y < 0));
                
            end
            
        end
        
        
        function labeled_data = simplelabel(fraginf)
            % dataDir = "C:\Users\v-dongqihan\Downloads\autoCollectedSyllables\autoCollectedSyllables\";
            dataDir = "C:\Users\Zhehao\Dropbox (OIST)\My_Stimuli\Eleinff\";
            
            fg = figure ;
            
            % -------------- FT options --------------
            wlen = 300;
            hop = 75;
            nfft = 300;
            
            fs = 32000;
            
            anal_win = blackmanharris(wlen, 'periodic');
            synth_win = hamming(wlen, 'periodic');
            
            % ------------------------------------
            
            curr_spectrogram = subplot(2, 10, 1);
            box on
            
            skip_10 = subplot(6, 10, 2);
            text(0, 0, 'Skip 10 ', 'HorizontalAlignment', 'center')
            xlim([-5, 5])
            ylim([-5, 5])
            box on
            axis off
            
            not_a_syllable = subplot(6, 10, 12);
            text(0, 0, 'Not a syllable', 'HorizontalAlignment', 'center')
            xlim([-5, 5])
            ylim([-5, 5])
            box on
            axis off
            
            back_to_previous = subplot(6, 10, 22);
            text(0, 0, 'Go back', 'HorizontalAlignment', 'center')
            xlim([-5, 5])
            ylim([-5, 5])
            box on
            axis off
            
            
            max_num_syllables  = 18;
            
            for i =  1:max_num_syllables
                genres(i) = subplot(2, 10, 2 + i);
                box on
            end
            
            % --------------------------------------------
            
            labeled_data.fraginf = fraginf;
            labeled_data.wlen = wlen;
            labeled_data.hop = hop;
            labeled_data.nfft = nfft;
            labeled_data.fs = fs;
            
            for i =  1:max_num_syllables
                labeled_data.label_example_syllable{i} = 0;
                
                % ------- in case you did it wrongly at, it can regret ------
                prev_Ss_softplused{i} = [[nan, nan]; [nan, nan]];
                prev_Ts{i} = [0, 1];
                prev_Fs{i} = [0, 1];
                
                prev_prev_Ss_softplused{i} = [[nan, nan]; [nan, nan]];
                prev_prev_Ts{i} = [0, 1];
                prev_prev_Fs{i} = [0, 1];
                
            end
            
            n_s = 0;
            
            
            
            while n_s < length(fraginf)
                
                n_s = n_s + 1;
                
                syllable = fraginf(n_s);
                y = syllable.y / max(abs(syllable.y));
                
                [STFT, F, T] = DQ.stft(y, anal_win, hop, nfft, fs);
                
                S_softplused = bdsoftplus(abs(STFT));
                labeled_fraginf(n_s).S_softplused = S_softplused;
                labeled_fraginf(n_s).STFT = STFT;
                labeled_fraginf(n_s).T = T;
                labeled_fraginf(n_s).F = F;
                
                contourf(curr_spectrogram, T, F,  S_softplused, 20, 'LineStyle', 'none')
                sound(syllable.y(1: min(fs, length(syllable.y))), fs) % play the sound for at longest 1 s
                title(curr_spectrogram, sprintf(" %d / %d", n_s, length(fraginf)))
                xlabel(curr_spectrogram, 'time')
                ylabel(curr_spectrogram, 'frequency')
                
                previous_selected_axis = 'Hello, World';
                
                k = waitforbuttonpress;
                if k == 0
                    if gca == not_a_syllable
                        labeled_data.fraginf(n_s).label = -1;
                        continue
                    elseif gca == back_to_previous
                        if n_s <= 1
                            disp("Cannot go back because this is the first sample!")
                            n_s = 0;
                            continue
                        end
                        label =  labeled_data.fraginf(n_s - 1).label;
                        if label >= 0
                            contourf(genres(label), prev_prev_Ts{label}, prev_prev_Fs{label},  prev_prev_Ss_softplused{label}, 20, 'LineStyle', 'none')
                        end
                        n_s = n_s - 2;
                        continue
                    elseif gca == skip_10
                        for kk = 0 : 9
                            if n_s + kk <= length(fraginf)
                                labeled_data.fraginf(n_s + kk).label = -1;
                            end
                        end
                        n_s = n_s + 9;
                        if n_s >= length(fraginf)  % already end
                            save(dataDir + "labeled_" + dataFile, "labeled_data")
                        end
                        continue
                    else
                        for i = 1: max_num_syllables + 1
                            if i == max_num_syllables + 1
                                disp("not correctly select a label!! Do again!")
                                n_s = n_s - 1;
                                continue
                            end
                            if genres(i) == gca
                                labeled_data.fraginf(n_s).label = i;
                                break
                            end
                        end
                    end
                end
                
                contourf(gca, T, F,  S_softplused, 20, 'LineStyle', 'none')
                axis off
                
                if  labeled_data.fraginf(n_s).label > 0
                    
                    label =  labeled_data.fraginf(n_s).label;
                    labeled_data.label_example_syllable{label} = syllable;
                    
                    prev_prev_Ss_softplused{label} = prev_Ss_softplused{label};
                    prev_prev_Ts{label} = prev_Ts{label} ;
                    prev_prev_Fs{label} = prev_Fs{label} ;
                    
                    prev_Ss_softplused{label} = S_softplused;
                    prev_Ts{label} = T;
                    prev_Fs{label} = F;
                    
                end
                
            end
            
            
            function y = bdsoftplus(x)
                
                y = (1 - 2*(x < 0)).* log(1 + (1 - 2*(x < 0)) .* x);
                
            end
            
            
            
            function x = inversebdsoftplus(y)
                
                x = (1 - 2*(y < 0)).* exp(y) - (1 - 2*(y < 0));
                
            end
            
        end
        
        
        function [x, t] = istft(STFT, awin, swin, hop, nfft, fs)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %          Inverse Short-Time Fourier Transform        %
            %               with MATLAB Implementation             %
            %                                                      %
            % Author: Ph.D. Eng. Hristo Zhivomirov        12/26/13 %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % function: [x, t] = istft(STFT, awin, swin, hop, nfft, fs)
            %
            % Input:
            % stft - STFT-matrix (only unique points, time
            %        across columns, frequency across rows)
            % awin - analysis window function
            % swin - synthesis window function
            % hop - hop size
            % nfft - number of FFT points
            % fs - sampling frequency, Hz
            %
            % Output:
            % x - signal in the time domain
            % t - time vector, s
            
            % signal length estimation and preallocation
            L = size(STFT, 2);          % determine the number of signal frames
            wlen = length(swin);        % determine the length of the synthesis window
            xlen = wlen + (L-1)*hop;    % estimate the length of the signal vector
            x = zeros(1, xlen);         % preallocate the signal vector
            
            % reconstruction of the whole spectrum
            if rem(nfft, 2)
                % odd nfft excludes Nyquist point
                X = [STFT; conj(flipud(STFT(2:end, :)))];
            else
                % even nfft includes Nyquist point
                X = [STFT; conj(flipud(STFT(2:end-1, :)))];
            end
            
            % columnwise IFFT on the STFT-matrix
            xw = real(ifft(X));
            xw = xw(1:wlen, :);
            
            % Weighted-OLA
            for l = 1:L
                x(1+(l-1)*hop : wlen+(l-1)*hop) = x(1+(l-1)*hop : wlen+(l-1)*hop) + ...
                    (xw(:, l).*swin)';
            end
            
            % scaling of the signal
            W0 = sum(awin.*swin);
            x = x.*hop/W0;
            
            % generation of the time vector
            t = (0:xlen-1)/fs;
            
        end
        
        function [STFT, f, t] = stft(x, win, hop, nfft, fs)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %              Short-Time Fourier Transform            %
            %               with MATLAB Implementation             %
            %                                                      %
            % Author: Ph.D. Eng. Hristo Zhivomirov        12/21/13 %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % function: [STFT, f, t] = stft(x, win, hop, nfft, fs)
            %
            % Input:
            % x - signal in the time domain
            % win - analysis window function
            % hop - hop size
            % nfft - number of FFT points
            % fs - sampling frequency, Hz
            %
            % Output:
            % STFT - STFT-matrix (only unique points, time
            %        across columns, frequency across rows)
            % f - frequency vector, Hz
            % t - time vector, s
            
            % representation of the signal as column-vector
            x = x(:);
            
            % determination of the signal length
            xlen = length(x);
            
            % determination of the window length
            wlen = length(win);
            
            % stft matrix size estimation and preallocation
            NUP = ceil((1+nfft)/2);     % calculate the number of unique fft points
            L = 1+fix((xlen-wlen)/hop); % calculate the number of signal frames
            STFT = zeros(NUP, L);       % preallocate the stft matrix
            
            % STFT (via time-localized FFT)
            for l = 0:L-1
                % windowing
                xw = x(1+l*hop : wlen+l*hop).*win;
                
                % FFT
                X = fft(xw, nfft);
                
                % update of the stft matrix
                STFT(:, 1+l) = X(1:NUP);
            end
            
            % calculation of the time and frequency vectors
            t = (wlen/2:hop:wlen/2+(L-1)*hop)/fs;
            f = (0:NUP-1)*fs/nfft;
            
        end
        

        function new_data = transform(data, fs, transform_type, transform_scale)
            
            
            switch transform_type
                case "speed_change"
                    % Simple make it plays faster. Duration and frequency change together
                    for i = 1 : size(data,2)  % i-th channel
                        new_data(:, i) = resample(data(:, i), 100, floor(transform_scale * 100));
                    end
                    
                case "frequency_change"
                    % Make the frequency change, while keep the duration unchanged
                    wlen = 400;
                    hop = 100;
                    nfft = 4000;
                    
                    anal_win = blackmanharris(wlen, 'periodic');
                    synth_win = hamming(wlen, 'periodic');
                    
                    n = 1 / transform_scale;
                    
                    for i = 1 : size(data,2)  % i-th channel
                        [STFT0{i}, F, T] = DQ.stft(data(:, i), anal_win, hop, nfft, floor(fs));  % F x T
                        
                        % Stretch the frequency distribution
                        nf = length(F);
                        if n <= 1
                            for j = 1 : length(T)
                                STFT{i}(:, j) = resample(STFT0{i}(1 : floor(nf * n), j), nf, floor(nf * n));
                            end
                        else
                            STFT{i} = STFT0{i} - STFT0{i};
                            for j = 1 : length(T)
                                STFT{i}(1 : floor(nf / n), j) = resample(STFT0{i}(:, j), floor(nf / n), nf);
                            end
                        end
                        
                        % inverse short-time Fourier transform
                        [new_data(:, i), ~] = DQ.istft(STFT{i}, anal_win, synth_win, floor(hop), nfft, fs);
                        
                    end
                    
                case "duration_change"
                    % Make the duration change, while keep the frequency unchanged
                    
                    new_data_0 = DQ.transform(data, fs, "frequency_change", transform_scale);
                    new_data = DQ.transform(new_data_0, fs, "speed_change", 1 / transform_scale);
                    
                    % normalize the power
                    new_data = mean(abs(data), 'all') / mean(abs(new_data), 'all')  * new_data;
                    
                case "add_noise"
                    % add Gaussian noise
                    new_data = data + transform_scale * mean(abs(data), 'all') * randn(size(data, 1), size(data, 2));
                    
                case "strength_change" % in db
                    % I modified here, Zhehao, 07082021
                    % !!!!!!!!!!!!!!!!!!!!!!!
                    
                    new_data = transform_scale ^ 0.5 * data;   % sound density is proportional to the square of wave amplitude
           
                    
                case "strength_change_db"  % very wierd
                    oldrms = rms(data);
                    
                    newrms = oldrms^transform_scale;
                    new_data = data*(newrms/oldrms);
                    
            
            end
            
            
            
            %% ---- Uncomment the below  plotting code for visualizing the spectrogram of old and new data ---------------
            wlen = 300;
            hop = 75;
            nfft = 300;
            anal_win = blackmanharris(wlen, 'periodic');
            synth_win = hamming(wlen, 'periodic');
            
            [STFT, F, T] = DQ.stft(data, anal_win, hop, nfft, fs);
            figure
            subplot(121)
            contourf(T, double(F), log(1+abs(STFT)), 30,'LineStyle','None')
            xlabel('Time(s)')
            ylabel('Frequency(Hz)')
            title("original spectogram")
            colorbar
            
            [STFT, F, T] = DQ.stft(new_data, anal_win, hop, nfft, fs);
            subplot(122)
            contourf(T, double(F), log(1+abs(STFT)), 30,'LineStyle','None')
            xlabel('Time(s)')
            ylabel('Frequency(Hz)')
            title("transformed spectogram ")
            colorbar
            
            
        end
        
        
    end
end

