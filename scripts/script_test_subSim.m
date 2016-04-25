% Try the SUBSET SIMULATION algorithm on all test cases and compute the ratio
% between the estimated probability of failure and the reference value (for
% visual inspection).
%
% Note: "subset simulation" is obtained as a special case of BSS (with
% options.use_gp set to false), so this is actually a test of the bss() function
% in disguise.

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

%%
% Common options for all runs

fig_opts = struct('disable_all', true);

options = struct(                       ...
    'p0',                     0.1,      ...  % ratio between stages
    'samplesize',             64000,    ...  % MC samplesize
    'figs',                   fig_opts, ...  % see script_set_fig_opts.m
    'use_gp',                 false     );   % no GP -> plain subsim


%%
% Run all examples

result = struct();

LIST_CASES = [1:6 101];

for k = 1:length(LIST_CASES),
    
    tic;
    
    CASE_NUM = LIST_CASES(k);
    problem = select_example (CASE_NUM);
    stage_data = bss(problem, options);
    
    result(k).t        = toc;
    result(k).pf_estim = stage_data(end).Pf;
    result(k).pf_ref   = problem.Pf_true;
    result(k).problem  = problem;
end


%%
% Display summary

fprintf('\n\nSUMMARY\n=======\n\n');

for k = 1:length(LIST_CASES)
    
    fprintf('case #%03d:  ', LIST_CASES(k));
    fprintf('DIM=%d  ', result(k).problem.dim);
    fprintf('pf_estim=%.2e  ', result(k).pf_estim);
    fprintf('pf_ref=%.2e  ', result(k).pf_ref);
    fprintf('ratio=%.2e  ', result(k).pf_estim / result(k).pf_ref);
    fprintf('\n');
    
end
