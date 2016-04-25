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

function TN = TruncatedNormalDistribution (mu, sd, xL, xU)

if nargin == 0
    
    dim = 1;
    mu = 0;
    sd = 1;
    xL = -inf;
    xU = +inf;
    
else
    
    dim = size (mu, 2);
    
    if ~ isequal (size (mu), [1 dim])
        error ('Incorrect size for argument mu.');
    end
    
    if ~ isequal (size (sd), [1 dim])
        error('Incorrect size for argument sd.');
    end
    
    if ~ isequal (size (xL), [1 dim])
        error('Incorrect size for argument xL.');
    end
    
    if ~ isequal (size (xU), [1 dim])
        error('Incorrect size for argument xU.');
    end
    
end

if ~ all (xL < xU)
    error ('xL should be strictly smaller than xU');
end

a = (xL - mu) ./ sd;   pa = normpdf (a);
b = (xU - mu) ./ sd;   pb = normpdf (b);

normC = 1 ./ (normcdf (b) - normcdf (a));
alpha = (pb - pa) .* normC;

bpb = b .* pb;      apa = a .* pa;
bpb(pb == 0) = 0;   apa(pa == 0) = 0;

N_mean = mu - alpha .* sd;
N_var  = (sd .^ 2) .* (1 - alpha .^ 2 - normC .* (bpb - apa));
N_std  = sqrt (N_var);

if any (isnan (N_std)), keyboard; end

TN = struct(                                            ...
    'dim', dim, 'mu', mu, 'sd', sd, 'xL', xL, 'xU', xU, ...
    'mean', N_mean, 'var', N_var, 'std', N_std          );

TN = class (TN, 'TruncatedNormalDistribution');

end % function
