function [sigseq, invsig] = generate_test_sequence(reps, N, sigtype, fs, fstart, fstop)
% GENERATE_SWEEP_SEQUENCE generates IR measurement test signals
%   GENERATE_SWEEP_SEQUENCE(reps, N, fs, FMIN, FMAX) generates a sequence
%   of sine sweps.
%
%   Inputs:
%   reps     ... number of sweep repetitions
%   N        ... [samples] length of individual sweeps
%   sigtype  ... 'tsp'   --> time-stretched pulse [1]
%                'chirp' --> exponential sine sweep [2]
%   fs       ... [Hz] sampling rate
%   fstart   ... [Hz] chirp start frequency
%   fstop    ... [Hz] chirp stop frequency
%
%   Outputs:
%   sigseq ... sequence of test signals
%   invsig ... inverse test signal
%
% References:
% [1] Farina, A. (2007) "Advancements in Impulse Response Measurements by 
% Sine Sweeps", in Proc. 122nd Conv. Audio Eng. Soc.
%
% [2] Meng, Q.; Sen, D.; Wang, S. & Hayes, L. (2008) "Impulse response 
% measurement with sine sweeps and amplitude modulation schemes", in Proc. 
% Int. Conf. Signal Process. and Commun. Syst.
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

% inits
T = N / fs;
t = linspace(0, T, N + 1);
t = t(1:end-1).';
logF = log(fstop/fstart);
K = T * 2 * pi * fstart / logF;
L = T / logF;

switch lower(sigtype)
    case 'chirp'
        % generate sweep
        sig = sin( K * (exp(t/L) - 1));
        
        % generate inverse sweep
        invsig = flipud(sig);
        w = K/L * exp(t/L);
        M = (2*pi*fstart) ./ w;
        G = sum(M.*(invsig.^2));
        invsig = invsig .* M;
        
        % normalise main peak to "1"
        invsig = invsig / G;
    case 'tsp'
        M = round(N / 4);
        k = 0:N/2;
        sig = exp(1i * 4 * M * pi * k.^2 / N^2).';
        sig = [sig; conj(flipud(sig(2:end-1)))];
        sig(N/2 + 1) = 1;
        sig = circshift(ifft(sig), -(N/2 - M));
        
        invsig = exp(-1i * 4 * M * pi * k.^2 / N^2).';
        invsig = [invsig; conj(flipud(invsig(2:end-1)))];
        invsig(N/2 + 1) = 1;
        invsig = circshift(ifft(invsig), (N/2 - M));
end

% generate sweep sequence
sigseq = repmat([sig; zeros(N,1)], reps, 1);

end
