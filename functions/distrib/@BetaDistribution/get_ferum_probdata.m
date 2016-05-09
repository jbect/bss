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

function probdata = get_ferum_probdata (distrib)

dim = distrib.dim;

probdata.correlation = eye (dim);

probdata.marg = [     ...
    7 * ones(dim, 1)  ...  % beta distribution
    nan(dim, 2)       ...  % mean, std (not provided, since inputtype = 1)
    (mean (distrib))' ...  % starting point (use the mean)
    (distrib.alpha)'  ...  % p1 = alpha
    (distrib.beta)'   ...  % p2 = beta
    zeros(dim, 1)     ...  % p3 = 0 (left end-point)
    ones(dim, 1)      ...  % p4 = 1 (right end-point)
    ones(dim, 1)];         % inputtype = 1

end % function
