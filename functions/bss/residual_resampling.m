% c = residual_resampling (w)
%     w = weights
%     c = counts (nb replications of each possible outcome)

% Copyright Notice
%
%    Copyright (C) 2016 CentraleSupelec
%
%    Authors:  Julien Bect  <julien.bect@centralesupelec.fr>
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

function N_out = residual_resampling (w, m_out)

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
    m_out = m_in;
end

% Expected number of replicates for each input particle
N_mean = w * m_out;

% Deterministic part of residual resampling
TOL = 10 * eps;
N_det = floor (N_mean + TOL);
m_det = sum (N_det);

% Unnormalized weights for the stochastic (multinomial) part
w_sto = max (0, N_mean - N_det);
if all (w_sto == 0)
    N_out = N_det;
    return;
end

% Stochastic part
m_sto = m_out - m_det;
N_sto = multinomial_resampling (w_sto, m_sto);

% Join determinitic and multinomial parts
N_out = N_det + N_sto;

end % function


%!% one input argument, not more, not less
%!error residual_resampling ();
%!test  residual_resampling (1);
%!test  residual_resampling (1, 0);
%!error residual_resampling (1, 0, 1);

%!% weights must be non-negative with a positive sum
%!error residual_resampling (-1);
%!error residual_resampling ([-1 2]);
%!error residual_resampling (0);
%!error residual_resampling ([0 0 0]);

%!% see if it doesn't crash for bigger sizes
%!test residual_resampling (ones(1, 50));
%!test residual_resampling (ones(1, 500));

%% test: special cases

%!% special case: only one particle (weights can be unnormalized)
%!test assert(isequal(residual_resampling (1), 1))
%!test assert(isequal(residual_resampling (2), 1))

%!% special case: two particles, the second with zero weight
%!test assert(isequal(residual_resampling ([1 0]), [2 0]))
%!test assert(isequal(residual_resampling ([2 0]), [2 0]))

%!% special case: three particles, the third one with zero weight
%!test
% c = residual_resampling ([1 1 0]);
% assert(isequal(c(3), 0));

%% test: valid output is for random input vectors

%!test  %% size of the ouput
%!  for n = 50:50:500,
%!     w = eps + rand(1, n);
%!     c = residual_resampling (w);
%!     assert(isequal(size(c), size(w)));
%!  end

%!test  %% sum of output values == nb particles
%!  for n = 50:50:500,
%!     w = eps + rand(1, n);
%!     c = residual_resampling (w);
%!     assert(isequal(sum(c), n));
%!  end

%!test  %% no NaNs of Infs
%!  for n = 50:50:500,
%!     w = eps + rand(1, n);
%!     c = residual_resampling (w);
%!     assert(~any(isnan(c)));
%!     assert(~any(isinf(c)));
%!  end
