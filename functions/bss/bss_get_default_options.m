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

function options = bss_get_default_options (options, problem)

% Default figure options
if ~ isfield (options, 'figs')
    options.figs = [];
end
options.figs = get_default_fig_options (options.figs);

% Bayes or not Bayes ?
if (~ isfield (options, 'use_gp')) || (isempty (options.use_gp))
    options.use_gp = true;
end

% COMMON: Monte Carlo sample size
if (~ isfield (options, 'samplesize')) || (isempty (options.samplesize))
    options.samplesize = 1000;
end

% COMMON: ratio between successive levels ("conditional probability")
if (~ isfield (options, 'p0')) || (isempty (options.p0))
    options.p0 = 0.10;
end

% COMMON: Metropolis-Hastings
if ~ isfield (options, 'MH')
    options.MH = [];
end
options.MH = get_default_MH_options (options.MH, problem);

% BAYES-ONLY: Size of the initial design
if (~ isfield (options, 'N0')) || (isempty (options.N0))
    options.N0 = 5 * problem.dim;
end

% BAYES-ONLY: model
if (~ isfield (options, 'model')) || (isempty (options.model))
    options.model = stk_model ('stk_materncov52_aniso', problem.dim);
    % note: we do not provide an automatic prior for the covariance parameters
end

% BAYES-ONLY: SUR parameters
if ~ isfield (options, 'SUR')
    options.SUR = [];
end
options.SUR = get_default_SUR_options (options.SUR);

% BAYES-ONLY: stopping threshold
if (~ isfield (options, 'stopping_thresh')) ...
        || (isempty (options.stopping_thresh))
    options.stopping_thresh = @(stage, final_stage, squaredCoV) ...
        bss_stopping_threshold (final_stage, squaredCoV, 0.05, 0.20);
end

% BAYES-ONLY: update covariance parameters inside SUR loop ?
if (~ isfield (options, 'update_inside_SUR_loop')) ...
        || (isempty (options.update_inside_SUR_loop))
    options.update_inside_SUR_loop = true;
end

end % function bss_get_default_options


function fo = get_default_fig_options (fo)

if (~ isfield (fo, 'disable_all')) || (isempty (fo.disable_all))
    fo.disable_all = false;
end

if (~ isfield (fo, 'dock')) || (isempty (fo.dock))
    fo.dock = false;  % Dockin does not work (yet) in Octave
end

if (~ isfield (fo, 'export')) || (isempty (fo.export))
    fo.export = false;
end

if (~ isfield (fo, 'range'))
    fo.range = [];
end

if (~ isfield (fo, 'save_prefix')) || (isempty (fo.save_prefix))
    fo.save_prefix = '';
end

end % function get_default_fig_options


function mh = get_default_MH_options (mh, problem)

% Proportionality factor for the initial instrumental distribution
%  ('initial' meaning at the first stage of the algorithm only; afterwards
%   the final sig_RW from the previous stage is used as the initial value)
if (~ isfield (mh, 'init_std_ratio')) || (isempty (mh.init_std_ratio))
    % Reminder: 2.38 / sqrt(d) is the optimal scaling for a Gaussian target
    %  in the case of Gaussian steps and large d according to RGG97
    mh.init_std_ratio = 2 / (sqrt (problem.dim));
end

% Number of steps for each particle (with or without adaptation)
if (~ isfield (mh, 'nb_steps')) || (isempty (mh.nb_steps))
    mh.nb_steps = 10;
end

% Adapt the standard deviation of MH steps ?
if (~ isfield (mh, 'do_adapt')) || (isempty (mh.do_adapt))
    mh.do_adapt = true;
end

% Adapt: Target acceptance rate
if (~ isfield (mh, 'target_rate')) || (isempty (mh.target_rate))
    mh.target_rate = 0.3;
end

% Adapt: Initial increment/decrement (on a log scale)
if (~ isfield (mh, 'delta_init')) || (isempty (mh.delta_init))
    mh.delta_init = log (2);
end

end % function get_default_MH_options


function sur = get_default_SUR_options (sur)

if (~ isfield (sur, 'm0_max')) || (isempty (sur.m0_max))
    sur.m0_max = 1000;
end

if (~ isfield (sur, 'frac')) || (isempty (sur.frac))
    sur.frac = 0.99;
end

if (~ isfield (sur, 'min_eval_per_stage')) ...
        || (isempty (sur.min_eval_per_stage))
    sur.min_eval_per_stage = 0;
end

if (~ isfield (sur, 'max_eval_per_stage')) ...
        || (isempty (sur.max_eval_per_stage))
    sur.max_eval_per_stage = 100;
end

end % function
