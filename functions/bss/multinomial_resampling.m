% c = multinomial_resampling (w, m_out)
%     w = weights
%     c = counts (nb replications of each possible outcome)

% Copyright Notice
%
%    Copyright (C) 2016 CentraleSupelec
%
%    Author:  Julien Bect  <julien.bect@centralesupelec.fr>
%
%    Address questions, bug reports, feature requests, or any other
%    correspondance related to BSS to kriging-help@lists.sourceforge.net

% Copying Permission Statement
%
%    This file is part of BSS (https://github.com/jbect/bss).
%
%    BSS is free software; you can redistribute it and/or modify it under
%    the terms of the  GNU Lesser General Public License  as published by
%    the Free Software Foundation;  either version 2.1 of the License, or
%    (at your option) any later version.
%
%    BSS is distributed  in the hope that it will  be useful, but WITHOUT
%    ANY WARRANTY;  without even the implied  warranty of MERCHANTABILITY
%    or FITNESS  FOR A  PARTICULAR PURPOSE.  See the  GNU  Lesser General
%    Public License for more details.
%
%    You should  have received  a copy  of the  GNU Lesser General Public
%    License along with BSS;  if not, see <http://www.gnu.org/licenses/>.

function c = multinomial_resampling (w, m_out)

if ~ all (w >= 0)
    error ('Negative or NaN weights are not allowed.');
end

% Normalize weights
S = sum (w);
if S == 0
    error ('The sum of weights should be strictly positive.');
end
w = w / S;

% Input sample size
m_in = numel (w);

% Output sample size
if nargin < 2,
    % Default: same as input sample size
    m_out = m_in;
end

% Boundary of "cells" in [0; 1]
tR = cumsum(w(:));
tR(end) = 1;
tL = [0; tR(1:(end-1))];
num = 1:m_in;

% Remove cells of length zero
b = (tL == tR);
tR(b) = [];
%tL(b) = [];
num(b) = [];

u = sort (rand (m_out, 1), 'ascend');
c = zeros (size (w));

k1 = 1; j = 1;
while k1 <= m_out
    
    % Find in which cell the first "u" lies
    while u(k1) > tR(j),
        j = j + 1;
    end
    %j = find((us(k1) > tL) & (us(k1) <= tR), 1, 'first');
    
    % Find all the others u's in the same cell
    k2 = k1;
    while (k2 < m_out) && (u(k2 + 1) <= tR(j)),
        k2 = k2 + 1;
    end
    %k2 = find((us > tL(j)) & (us <= tR(j)), 1, 'last');
    
    % The corresponding particles are copies of particle j
    c(num(j)) = k2 - k1 + 1;
    
    % Next output particle ?
    k1 = k2 + 1;
end

end % function


%% some basic tests

%!% one input argument, not more, not less
%!error multinomial_resampling ();
%!test  multinomial_resampling (1);
%!test  multinomial_resampling (1, 1);
%!error multinomial_resampling (1, 1, 1);

%!% weights must be non-negative with a positive sum
%!error multinomial_resampling (-1);
%!error multinomial_resampling ([-1 2]);
%!error multinomial_resampling (0);
%!error multinomial_resampling ([0 0 0]);

%!% see if it doesn't crash for bigger sizes
%!test multinomial_resampling (ones(1, 50));
%!test multinomial_resampling (ones(1, 500));

%% test: special cases

%!% special case: only one particle (weights can be unnormalized)
%!test assert(isequal(multinomial_resampling (1), 1))
%!test assert(isequal(multinomial_resampling (2), 1))

%!% special case: two particles, the second with zero weight
%!test assert(isequal(multinomial_resampling ([1 0]), [2 0]))
%!test assert(isequal(multinomial_resampling ([2 0]), [2 0]))

%!% special case: three particles, the third one with zero weight
%!test
% c = multinomial_resampling ([1 1 0]);
% assert(isequal(c(3), 0));

%% test: valid output is for random input vectors

%!test  %% size of the ouput
%!  for n = 50:50:500,
%!     w = eps + rand(1, n);
%!     c = multinomial_resampling (w);
%!     assert(isequal(size(c), size(w)));
%!  end

%!test  %% sum of output values == nb particles
%!  for n = 50:50:500,
%!     w = eps + rand(1, n);
%!     c = multinomial_resampling (w);
%!     assert(isequal(sum(c), n));
%!  end

%!test  %% no NaNs of Infs
%!  for n = 50:50:500,
%!     w = eps + rand(1, n);
%!     c = multinomial_resampling (w);
%!     assert(~any(isnan(c)));
%!     assert(~any(isinf(c)));
%!  end
