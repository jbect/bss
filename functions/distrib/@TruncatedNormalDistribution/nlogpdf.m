% NLOGPDF computes the opposite of the log-pdf of the distribution

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

function nlp = nlogpdf (distrib, x, exact)

% Points outside the support have infinite nlp
nlp = inf (size (x, 1), 1);          % Initialize with Infs everywhere
b = is_inside_support (distrib, x);  % Find out where we must really compute

% Compute up to an additive constant
t = bsxfun (@minus, x(b, :), distrib.mu);
t = bsxfun (@rdivide, t, distrib.sd);
nlp(b) = 0.5 * (sum (t .^ 2, 2));

% Compute exactly?
if (nargin < 3) || exact
    % log(2*pi)/2 = 0.91893853320467267 (double precision)
    nlp = nlp + distrib.dim * 0.91893853320467267 ...
        + sum (log (distrib.sigma)) - log (distrib.normC);
end

end % function
