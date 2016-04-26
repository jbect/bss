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

function [mcs, accept_rate, sig_RW] = bss_move ...
    (x0, obs, model, u, options, problem, sig_RW)

input_distrib = problem.input_distrib;
cost_fun = problem.cost_fun;  % Only used when use_gp = false

% At the first stage, sig_RW is empty: we use standard deviations
%  proportional to the standard deviations of the input distribution
if isempty (sig_RW)
    sig_RW = (std (problem.input_distrib)) * options.MH.init_std_ratio;
end

x1 = x0;  m = size (x1, 1);

if options.use_gp
    xi = obs.xi;  yi = obs.yi;           % Observations
    xi = {xi, stk_kreq_qr(model, xi)};   % Use an experimental feature of STK
    y1 = stk_predict (model, xi, yi, x1);
    y1 = [y1 stk_dataframe(sqrt(y1.var), {'std'})];
    g1 = g_fun (u, y1.mean, y1.std);
else
    y1 = cost_fun (x1);
    g1 = ones (m, 1);
end

accept_rate.expected = zeros (1, options.MH.nb_steps);
accept_rate.observed = zeros (1, options.MH.nb_steps);

for k = 1 : options.MH.nb_steps
    
    % Gaussian proposal distribution
    x2 = x1 + bsxfun (@times, randn (size (x1)), sig_RW);
    
    % Check domains
    if ~ all (is_inside_support (input_distrib, x1)),
        error ('This should never happen !');
    end
    b2 = is_inside_support (input_distrib, x2);
    
    if options.use_gp  % BSS
        
        % Predict the response on the proposed states
        y2 = stk_predict (model, xi, yi, x2);
        y2 = [y2 stk_dataframe(sqrt(y2.var), {'std'})]; %#ok<AGROW>
        
        % Compute the pdf g2 (pdf wrt the input distribution)
        g2 = zeros (m, 1);
        g2(b2) = g_fun (u, y2.mean(b2), y2.std(b2));
        
        % Ratio of target density values
        g_ratio = g2 ./ g1;
        
    else  % Subsim
        
        % Evaluate the response on the proposed states
        y2 = - Inf (m, 1);
        y2(b2) = cost_fun (x2(b2, :));
        
        % Compute the "relative pdf" g2 (pdf wrt the input distribution)
        g2 = zeros (m, 1);
        g2(b2) = double (y2(b2) > u);
        
        % Ratio of target density values
        g_ratio = g2;  % g1 = 1, always
    end
    
    % Twice the anti-logpdf of the input distribution
    Q2 = nlogpdf (input_distrib, x2, false);
    Q1 = nlogpdf (input_distrib, x1, false);
    
    % Acceptance ratio (symmetric RW)
    acceptance_ratio = min (1, (exp (Q1 - Q2)) .* g_ratio);
    
    % Accept / reject
    a = (rand (m, 1) < acceptance_ratio);
    
    % Copy accepted proposals
    x1(a, :) = x2(a, :);
    y1(a, :) = y2(a, :);
    if options.use_gp,
        g1(a) = g2(a);
    end
    
    accept_rate.expected(k) = (sum (acceptance_ratio)) / m;
    accept_rate.observed(k) = (sum (a)) / m;
    
    % Adaptation
    if options.MH.do_adapt
        direct = 2 * (accept_rate.expected(k) > options.MH.target_rate) - 1;
        sig_RW = sig_RW * (exp (options.MH.delta_init / k * direct));
    end
end

mcs = struct ('xt', x1, 'yt', y1, 'gx', g1, 'u_base', u);

end  % function
