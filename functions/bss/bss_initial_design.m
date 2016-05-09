% Initial design as proposed in our BSS journal paper

% Copyright Notice
%
%    Copyright (C) 2016 CentraleSupelec
%
%    Authors:  Julien Bect       <julien.bect@centralesupelec.fr>
%              Ling Li           <ling.li.supelec@gmail.com>
%              Emmanuel Vazquez  <emmanuel.vazquez@centralesupelec.fr>
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

function [x0, y0] = bss_initial_design (problem)

% Quantiles that define the bounds of the domain
epsilon = 1e-5;  proba = [epsilon; 1 - epsilon];

% Number of random LHS to draw
Q = 1e4;

% Size of the initial design
N0 = 5 * problem.dim;

% Construct a box
dom = quantile (problem.input_distrib, proba);
dom = stk_hrect (dom);
for j = 1:problem.dim
    dom.colnames{j} = sprintf ('x%d', j);
end

% Initial design
x0 = stk_sampling_maximinlhs (N0, [], dom, Q);
y0 = stk_feval (problem.cost_fun, x0);

end % function
