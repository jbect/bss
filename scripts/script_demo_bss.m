% script_bss_demo
%
% A short script that simply runs BSS with default parameters on one example.
%
% This should be fine for showing demos.

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


%% Example

CASE_NUM = 3; % cantilever beam example
problem = select_example (CASE_NUM);


%% BSS options

options = bss_get_default_options ([], problem);
%options.figs.dock = true;   % Matlab-only

% Override default settings to force more evaluations on this simple case
options.stopping_thresh = @(stage, final_stage, squaredCoV) ...
    bss_stopping_threshold (final_stage, squaredCoV, 0.005, 0.02);

warning('OFF', 'STK:stk_predict:NegativeVariancesSetToZero');


%% Run BSS

tic;

stage_data = bss (problem, options);

fprintf ('Pf_estimate = %.3e\n', stage_data(end).Pf);
fprintf ('Pf_true = %.3e\n ', problem.Pf_true);

t0 = toc;  fprintf ('elapsed time = %.1f seconds\n', t0);


%% Display (part of) stage_data

fprintf ('\n\nStage summary\n');

disp (stk_dataframe ([[stage_data.n_evals]; [stage_data.u_base]; ...
    [stage_data.u_target]]', {'n_evals', 'u_base', 'u_target'}))
