% Try the BAYESIAN SUBSET SIMULATION algorithm on all test cases and compute the
% ratio between the estimated probability of failure and the reference value
% (for visual inspection).

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

clear all; close all; clc

s = RandStream ('mt19937ar', 'Seed', 0);
RandStream.setGlobalStream (s);


%%
% Run all examples

result = struct();

LIST_CASES = [1:6 101];

for k = 1:length(LIST_CASES),
    
    tic;
    
    CASE_NUM = LIST_CASES(k);
    problem = select_example (CASE_NUM);
    
    options = bss_get_default_options ([], problem);
    options.figs.disable_all = true;
    
    stage_data = bss (problem, options);
    
    result(k).t          = toc;
    result(k).pf_estim   = stage_data(end).Pf;
    result(k).pf_ref     = problem.Pf_true;
    result(k).problem    = problem;
    result(k).stage_data = stage_data;
end


%%
% Display summary

fprintf('\n\nSUMMARY\n=======\n\n');

for k = 1:length(LIST_CASES),
    
    fprintf('case #%03d:  ', LIST_CASES(k));
    fprintf('DIM=%d  ', result(k).problem.dim);
    fprintf('pf_estim=%.2e  ', result(k).pf_estim);
    fprintf('pf_ref=%.2e  ', result(k).pf_ref);
    fprintf('ratio=%.2e  ', result(k).pf_estim / result(k).pf_ref);
    fprintf('\n');
    
end
