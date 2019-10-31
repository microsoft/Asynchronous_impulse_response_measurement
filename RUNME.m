% Test script illustrating usage of asynchronous_IR_estimation.m
% 
% Reference:
% Gamper, H. (2017). "Clock drift estimation and compensation for 
% asynchronous impulse response measurements", in Proc. Workshop Hands-free 
% Speech Communication and Microphone Arrays (HSCMA), San Francisco, CA,
% USA.
% https://www.microsoft.com/en-us/research/publication/clock-drift-estimation-compensation-asynchronous-impulse-response-measurements/
%
% Created by Hannes Gamper, Microsoft Research.
% 1.0 - Oct 2019 Initial release
% 
% Please post comments to the Github page for this entry if you have any
% bugs or feature requests:
% https://github.com/microsoft/Asynchronous_impulse_response_measurement

%   Copyright 2019 Microsoft Corporation
%   
%   Permission is hereby granted, free of charge, to any person obtaining a 
%   copy of this software and associated documentation files (the "Software"), 
%   to deal in the Software without restriction, including without limitation 
%   the rights to use, copy, modify, merge, publish, distribute, sublicense, 
%   and/or sell copies of the Software, and to permit persons to whom the 
%   Software is furnished to do so, subject to the following conditions:
%   
%   The above copyright notice and this permission notice shall be included in 
%   all copies or substantial portions of the Software.
%   
%   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
%   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
%   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
%   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
%   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
%   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
%   DEALINGS IN THE SOFTWARE.

close all; clc; clearvars;

%% constants
fs = 8000;           % [Hz] sampling rate
signalType = 'tsp';  % 'tsp' (time-stretched pulse) | 'chirp' (exponential sine sweep)
FMIN = 0.01*fs;      % [Hz] chirp start frequency (ignored for 'tsp')
FMAX = 0.90*fs;      % [Hz] chirp stop frequency (ignored for 'tsp')
SNR = 30;            % [dB] signal-to-noise ratio
Nsamples = fs;       % length of test signal
reps = 4;            % number of test signal repetitions
sim_drift = -12.3456;% [samples/s] simulated clock drift rate

%% generate sweep and inverse sweep
[testseq, invsig] = generate_test_sequence(reps, Nsamples, 'tsp', fs, FMIN, FMAX);

%% simulate recording
% in actual experiment, sweepseq would be played back or recorded on DUT
% here, a noisy, reverberant, asynchronous recording of the test signal is simulated
fprintf('Simulating measurement...\n');
[rec, h] = simulate_asynchronous_measurement(testseq, Nsamples, fs, sim_drift, SNR);

%% analyse recording
% NOTE: a reference PLAYBACK clock and asynchronous RECORDING (on a DUT) are assumed
fprintf('Estimating IR and clock drift...\n');
[IRs, estimated_clock_drift_per_sample] = asynchronous_IR_estimation(rec, invsig, reps);
estimated_clock_drift = estimated_clock_drift_per_sample * fs;

%% plots
figure;
set(gcf, 'Unit', 'pixels', 'Position', [50, 50, 700, 700]);
Nchan = size(h,2);
lstyles = {'-', '--', '-'};
for ci = 1:Nchan
    for ii = 1:2
        if ii==1
            ir = h(:,ci);
        else
            ir = IRs(:,ci);
        end
        subplot(3, Nchan, ci);
        tvect = 0 : 1 / fs : (size(ir,1)-1)/fs;
        plot(1000*tvect, ir, 'LineStyle', lstyles{ii});
        ylim([-1.1,1.1]);
        grid on;
        hold on;
        title({
            sprintf('Channel %d', ci); 
            sprintf('simulated drift = %1.4f [samples/s]', sim_drift);
            sprintf('estimated drift = %1.4f [samples/s]', estimated_clock_drift)
            });
        xlabel('[ms]');
        ylabel('IR');
        if ii==1
            IR1 = fft(ir);
        else
            legend('reference', 'estimate', 'Location', 'SouthEast');
            
            IR = [fft(ir), IR1];
            IRerr = fft(ir) ./ IR1;
            fvect = linspace(0, fs, size(IRerr,1)+1);
            for spi = 1:2
                subplot(3, Nchan, ci + spi*2);
                ptmp = [IR, IRerr];
                if spi==1
                    for pii = 1:3
                        plot(fvect(1:end-1), 20*log10(abs(ptmp(:,pii))), 'LineStyle', lstyles{pii});
                        hold on;
                    end
                    ylabel('[dB]');
                    tstr = 'magnitude';
                    ylim([-33, 33]);
                else
                    for pii = 1:3
                        plot(fvect(1:end-1), angle(ptmp(:,pii)), 'LineStyle', lstyles{pii});
                        hold on;
                    end
                    ylabel('[rad]');
                    tstr = 'phase';
                    ylim([-1.1*pi,1.1*pi]);
                    set(gca, 'YTick', [-pi,0,pi]);
                    set(gca, 'YTickLabel', {'-\pi', 0, '\pi'});
                end
                xlim([0, fs/2]);
                xlabel('[Hz]');
                grid on;
                hold on;
                title(tstr);
                legend('reference', 'estimate', 'error', 'Location', 'SouthEast');
            end
        end
    end
end
