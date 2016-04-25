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

function [J, y_pred, stop_flag] = J_SUR_trick ...
    (model, x_obs, y_obs, x, w, threshold, SUR_params)

m = size (x, 1);      % number of integration points
n = size (x_obs, 1);  % number of observations

% Prepare a vector to store criterion values  (a vector of size m,
%  since the integration points also play the role of candidate points)
J = inf (m, 1);

% Compute predictions
y_pred = stk_predict (model, x_obs, y_obs, x);

% SAFER: fix the variance at observation points
b = ismember (x, x_obs, 'rows');
y_pred.var(b) = 0;

% Compute the current probability of misclassification
p_current = stk_pmisclass (threshold, y_pred);

if any (isnan (p_current))
    warning ('Hmmmm.... p_current contains some NaNs...');
end

% Value of the criterion if no new observation is added
if isempty (w)
    Gamma_current = (sum (p_current)) / m;
else
    Gamma_current = sum (w .* p_current);
end

% Stop if the current value of the integral criterion is below a threshold
if Gamma_current < 1e-150
    stop_flag = 1;
    return
end

% Compute the current *weighted* probability of misclassification
%  and sort the points according to this value, in decreasing order
if isempty (w)
    [wp_sort, is] = sort (p_current, 'descend');
else
    [wp_sort, is] = sort (w .* p_current, 'descend');
end

% Pruning: only keep points with a non-negligible
%  (current, weighted) probability of misclassification
cwp_sort = cumsum (wp_sort);
m0 = min (SUR_params.m0_max, ...
    find (cwp_sort >= SUR_params.frac * cwp_sort(end), 1, 'first'));
if (m0 > 10000)
    error (sprintf ('m0=%d: TOO LARGE :!!!', m0));
end
is = is(1:m0);
ns = length (is);

% Stop if there is no candidate point with a positive variance
if ns == 0,
    stop_flag = 2;
    return;
end

% Join observations and selected points
x0 = [x_obs; x(is, :)];

% Switch to a discrete index set
model0 = stk_model ('stk_discretecov', model, x0);

% Compute all the posterior quantities that we will need, once and for all
[ys_pred, ~, ~, Kpost_ys] = stk_predict (model0, (1:n)', y_obs, n + (1:ns)');

% Trick : make row deletions slightly faster, in the loop below
ys_pred.rownames = {};

% Loop over the set of candidate points
for k = 1:ns
    
    % set of points where the probability must be computed
    idx00 = 1:ns;
    idx00(k) = [];
    
    ys00_pred = ys_pred;
    ys00_pred(k, :) = [];
    
    Kpost_ys00_ynew = Kpost_ys(idx00, k);
    Kpost_ynew_ynew = Kpost_ys(k, k);
    
    % We use the current value as a proxy for the expected value for the points
    % that have not been selected
    p_expected = p_current;
    
    % Compute the pointwise criterion (integrand)
    p_expected(is(idx00)) = stk_pmisclass (threshold, ...
        ys00_pred, Kpost_ys00_ynew, Kpost_ynew_ynew);
    
    % Compute the SUR criterion for the k^th candidate point
    if isempty (w)
        J(is(k)) = (sum (p_expected)) / m;
    else
        J(is(k)) = sum (w .* p_expected);
    end
    
    if isnan (J(is(k))),
        warning ('Failed to compute a value of the SUR criterion (--> NaN).');
    end
end

stop_flag = 0;

end % function

%#ok<*WNTAG>
