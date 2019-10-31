function [rec, h] = simulate_asynchronous_measurement(sig, N, fs, drift, SNR)
% SIMULATE_ASYNCHRONOUS_MEASUREMENT simulation of clock drift measurement
%   SIMULATE_ASYNCHRONOUS_MEASUREMENT(sig, N, fs, drift, SNR) simulates a
%   measurement of the test signal sig in a simulated reverberant
%   environment, subject to a clock drift rate drift [samples/s].
%
%   Inputs:
%      sig     ... audio signal
%      N       ... [samples] reverb tail length
%      fs      ... [Hz] sampling rate
%      drift   ... [samples/s] drift rate
%      SNR     ... [dB] signal-to-noise ratio (optional)
%
%   Outputs:
%      rec     ... simulated recording
%      h       ... reverb impulse response
%
% Created by Hannes Gamper, Microsoft Research.
% 1.0 - Oct 2019 Initial release

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

% simulate (room) impulse response
T60 = 0.2;                      % [s] T60 reverberation time
h = simple_reverb(N, T60 * fs);

% simulate recorded measurement
for ci = 1:size(h,2)
    c = filter(h(:,ci), 1, sig);
    if ci==1
        rec = zeros(size(c,1), size(h,2));
    end
    rec(:,ci) = c;
end

% resample to simulate clock drift
d_per_sample = drift / fs;
x = 1:size(rec,1);
xq = (0:1/(1 + d_per_sample):(length(x) - 1)/(1 + d_per_sample)) + 1;
rec = sinc_resampling(x, rec, xq);

if exist('SNR', 'var')
    eng = @(x) (mean(mean(abs(fft(x)),1)));
    n = randn(size(rec));
    n = n * (eng(rec) / eng(n)) * 10^(-SNR/20);
    assert(abs(SNR - 20*log10(eng(rec) / eng(n))) < 1);
    rec = rec + n;
end

end

%% simple reverb
function h = simple_reverb(N, T60_samples)
% SIMPLE_REVERB exponentially decaying random reverb
%   SIMPLE_REVERB(N, T60_samples) generates a reverb tail of length N
%   [samples] that decays by -60dB within T60_samples.

% constants
initial_delay = 5; % [samples] initial delay
Nchan = 2;

% inits
N = N - initial_delay;

ha = logspace(0, -3 * (N / T60_samples), N);
rng(8);
h = randn(N,Nchan);
h = h .* repmat(ha(:), 1, size(h,2));

% add initial delay
h = [zeros(initial_delay, Nchan); h];

% simple low-pass
h = filter([0.5,0.5], 1, h);
end
