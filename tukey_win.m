function twin = tukey_win(N, alph)
%TUKEY_WIN tukey window
%   TUKEY_WIN(N, alph) returns N x 1 Tukey window with parameter alph.
%
% Reference: 
% https://en.wikipedia.org/wiki/Window_function
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

ptr1 = floor(alph * (N-1) / 2);
ptr2 = ceil((N-1) * (1 - alph/2));

twin = ones(N,1);
n = 0:ptr1;
twin(1:ptr1+1) = 1/2 * (1 + cos(pi * (2*n / (alph*(N-1)) - 1)));
n = ptr2:N-1;
twin(ptr2+1:end) = 1/2 * (1 + cos(pi * (2*n / (alph*(N-1)) - 2/alph + 1)));
end

