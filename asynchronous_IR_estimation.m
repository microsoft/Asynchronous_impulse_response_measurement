function [IRs, d_per_sample] = asynchronous_IR_estimation(rec, invsig, reps, Npause)
%ASYNCHRONOUS_IR_ESTIMATION Matlab implementation of impulse response 
%estimation without clock synchronisation
% 
% IRstruct = ASYNCHRONOUS_IR_ESTIMATION(rec, invsweep, reps) returns 
% the impulse response and clock drift estimated from the recorded response 
% in rec.
%
% Inputs:
% rec          ... [L x M] recording samples, where M = # of input channels
% invsig       ... [N x 1] inverse test signal for IR estimation; assumed
%                  to be the same length as the test signal
% reps         ... number of test signal repetitions (minimum: 2)
% Npause       ... [samples] (optional) pause after each test signal;
%                  default: Npause = N
%
% Outputs:
% IRs          ... estimated impulse response(s)
% d_per_sample ... [samples^(-1)] estimated clock drift per sample
%
% Reference:
% Hannes Gamper (2017). "Clock drift estimation and compensation for 
% asynchronous impulse response measurements", in Proc. Workshop Hands-free 
% Speech Communication and Microphone Arrays (HSCMA), San Francisco, CA,
% USA.
% https://www.microsoft.com/en-us/research/wp-content/uploads/2017/03/Clock_drift_estimation_HSCMA_2017.pdf
%
% Created by Hannes Gamper, Microsoft Research.
% 1.0 - Oct 2019 Initial release
% 
% TODO:
% - implement proper gradient descent algorithm as alternative to fmincon

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

% constants
Ngrid = 201;              % number of grid search iterations (if fmincon is not available)
d_bounds = [-2,2];        % [samples] search bounds for drift estimation
d_min_per_sample = 10e-5; % [samples/sample] minimum drift considered non-zero

% check inputs
if size(rec,1)==1
    warning('Column vector expected as input signal.');
    rec = rec(:);
end
[L, Nchan] = size(rec);

if size(invsig,1)==1
    warning('Column vector expected as inverse test signal.');
    invsig = invsig(:);
end
assert(size(invsig,2)==1, 'S x 1 vector expected as inverse test signal.');
N = length(invsig);

if nargin < 4
    Npause = N;
end
NP = N + Npause;

if L < (NP * reps)
    warning('Input signal length should be >= (N + Npause) x reps. Padding input with zeros...');
    rec = [rec; zeros(2 * N * reps - L, Nchan)];
end

assert(reps > 1, 'At least 2 repetitions required for drift estimation.');

% chop signal into test sequences
filtbuf = buffering(rec, NP);
filtbuf = filtbuf(:,1:reps,:);

% filter with inverse sweep
filtbuf = filter(invsig, 1, filtbuf);

% to avoid getting stuck in a local minimum, estimate drift via
% cross-correlation across all microphones and repetitions
CC = fft(filtbuf(:,2:end,:)) .* conj(fft(filtbuf(:,1:end-1,:)));
cc = real(ifft(CC));
cc = fftshift(cc,1);
[~, maxind] = max(cc);
lags = (-(NP/2) : (NP/2 - 1));

% calculate starting point for minimisation
% for sub-sample resolution, consider using (parabolic) interpolation
d0 = lags(median(median(maxind)));

% estimate drift
d = min_search(@(d) (alignerr(filtbuf, d)), d0, d_bounds(1), d_bounds(2), Ngrid);

% check whether estimated drift is below minimum
if abs(d / NP) < d_min_per_sample
    warning('Drift %1.4f < %1.4f --> setting to zero.', (d/NP), d_min_per_sample);
    d = 0;
end

% estimate IRs
IRs = zeros(N, Nchan);
d_per_sample = d / NP; % clock drift per sample
d_per_sample_drift = d / (NP + d);

% resample raw rec to compensate for estimated drift
xrec = 1:size(rec,1);
xqrec = (0:1/(1 - d_per_sample_drift):(length(xrec) - 1)/(1 - d_per_sample_drift)) + 1;
rec_resampled = sinc_resampling(xrec, rec, xqrec);
recbuf = buffering(rec_resampled, NP);

% sanity check resampling & alignment
d2 = min_search(@(d) (alignerr(recbuf, d)), 0, d_bounds(1), d_bounds(2), Ngrid);
assert((d2 / NP) < d_min_per_sample, 'Resampling & alignment failed!');
    
for ci = 1:Nchan
    % average over runs
    buf = recbuf(:,:,ci);
    buf = mean(buf, 2);
    buf = filter(invsig, 1, buf);
    ir = buf(N+1:end);
    
    % store ir
    len = min(size(ir,1), N);
    IRs(1:len,ci) = ir(1:len);
end

end

%% buffer alignment error
function [err, y] = alignerr(x, d)

% inits
y = nan(size(x));
Nreps = size(x,2);

% delay individual runs
for n = 1:Nreps
    y(:,n,:) = fractional_delay(squeeze(x(:,n,:)), -d * n);
end

% compare to first run
refy = repmat(y(:,1,:),1,size(y,2));
err = nan_sum(median(nan_sum( (refy - y).^2, 1 ), 2), 3);
end

%% fractional delay via phase shift
function yd = fractional_delay(y, d)
% make sure length is even
ISODD = mod(size(y,1),2);
if ISODD
    y = [y; zeros(1, size(y,2))];
end

Y = fft(y);
Ym = abs(Y);
Ya = angle(Y);

a = linspace(0, -d * 2 * pi, size(y,1) + 1).';
a = a(1:end-1);

% apply phase shift
Ymd = Ym;
Yad = Ya + repmat(a, 1, size(Ya,2));
Yad(end/2 + 2 : end,:) = -flipud(Yad(2:end/2,:));

Yd = Ymd .* exp(1i * Yad);
yd = real(ifft(Yd));

if ISODD
    yd = yd(1:end-1,:);
end

end

%% buffering: replaces built-in Matlab buffer function
function buf = buffering(x, N)
[Nx, Nchan] = size(x);
Nb = ceil(Nx / N);
buf = [x; zeros(Nb * N - Nx, Nchan)];
buf = reshape(buf, N, [], Nchan);
end

%% nan_sum: replaces nansum in Matlab
function y = nan_sum(x, dim)
x(isnan(x)) = 0;
y = sum(x, dim);
end

%% search function minimum
function d = min_search(fun, x0, xl, xu, Ngrid)
try 
    % test whether fmincon is available on system
    fmincon(@sum, 0,[],[],[],[],[],[],[],optimset('Display', 'none', 'MaxFunEvals', 1));
    HAVE_FMINCON = true;
catch
    HAVE_FMINCON = false;
end

if HAVE_FMINCON
    % find function minimum using fmincon
    d = fmincon(fun, x0, [],[],[],[], x0 + xl, x0 + xu, [], optimset('Display', 'none'));
else
    % alternatively, use simple grid search
    warning('fmincon not available - using simple grid search instead...');
    xtest = x0 + linspace(xl, xu, Ngrid);
    xerrs = nan(1,length(xtest));
    for di = 1:length(xtest)
        xerrs(di) = fun(xtest(di));
    end
    [~,minind] = min(xerrs);
    d = xtest(minind);
end

end
