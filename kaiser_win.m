function kwin = kaiser_win(N, bet)
%KAISER_WIN kaiser window
%   KAISER_WIN(N, bet) returns N x 1 Kaiser window with parameter bet.
%
% Reference: 
% https://ccrma.stanford.edu/~jos/sasp/Kaiser_Window.html
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

n = linspace(-(N-1)/2, (N-1)/2, N).';
I0 = @(bet) (besseli(0, bet));
kwin = I0(bet * sqrt(1 - (n / (N/2)).^2)) / I0(bet);
end

