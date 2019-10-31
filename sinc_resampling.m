function yr = sinc_resampling(x, y, xq)
%SINC_RESAMPLING sinc interpolation
%   SINC_RESAMPLING(x, y, xq) interpolates the values y, sampled at x, at
%   the interpolated sample instances xq.
%
% NOTE: This is a brute-force implementation that can be quite slow for
% long input files.
%
% source: https://ccrma.stanford.edu/~jos/resample/resample.html
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

% constants
Nsinc = 512;  % [samples] length of sinc window

% ensure column vectors
x = x(:);
xq = xq(:);

% ensure x starts from "1"
x1 = x(1);
x = x - x1 + 1;
xq = xq - x1 + 1;

% inits
T = x(2) - x(1);
kwin = kaiser_win(2 * Nsinc + 1, 10).';
Nx = size(x, 1);
Nq = size(xq,1);
yr = zeros(Nq, size(y,2));

% limit xq to range of x (to avoid NaNs)
xq_lim = xq;
xq_lim(xq_lim < min(x)) = min(x);
xq_lim(xq_lim > max(x)) = max(x);

% find nearest sample instants
inds = interp1(x, 1:Nx, xq_lim, 'nearest');

for yi = 1:Nq
    ind = inds(yi);
    ptr1 = max(1, ind-Nsinc);
    ptr2 = min(Nx, ind+Nsinc);
    s = csin((xq(yi) - (ptr1:ptr2) * T)/T);
    
    % window sinc
    s = s .* kwin(max(1, Nsinc - ind + 2) : min(2*Nsinc+1, Nx - ind + Nsinc + 1));
    
    % perform interpolation
    yr(yi,:) = s * y(ptr1:ptr2,:);
end

end

%% sinc function
% substitute for built-in Matlab function
function s = csin(x)
pix = pi*x;
s = sin(pix) ./ pix;
s(x==0) = 1;
end
