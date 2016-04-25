% bss_subsim_stage ...

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

function [mcs, stage_data, final_stage] = bss_subsim_stage ...
    (stage, mcs, stage_data, problem, options)

xt = mcs.xt;
yt = mcs.yt;
final_stage = false;  % is it the final stage ?

if stage == 1,
    Pf_base = 1.00;
    sig_RW = [];
else
    Pf_base = stage_data(stage - 1).Pf;
    sig_RW = stage_data(stage - 1).sig_RW;
end

% visu BEFORE choosing a target threshold
bss_figure(21, options.figs, xt, [], [], [stage_data(1:stage-1).u_target], ...
    [], problem, stage, []); drawnow;

% UPDATE TARGET THRESHOLD
u_target = quantile(yt, 1 - options.p0);
p1 = mean(yt > problem.u_final);
if (p1 > options.p0 ^ 2) && (p1 < options.p0)
    u_target = quantile(yt, 1 - sqrt(p1));
end

if u_target >= problem.u_final,
    u_target = problem.u_final;
    final_stage = true;
end

g_next = (yt > u_target);
estim_ratio = mean(g_next);
Pf_target = Pf_base * estim_ratio;

% visu AFTER choosing a target threshold
bss_figure(21, options.figs, xt, [], [], [stage_data(1:stage-1).u_target], ...
    u_target, problem, stage, Pf_target); drawnow; ax = axis();

fprintf('>>>>> estim_ratio = %.5f\n', estim_ratio);
fprintf('>>>>> Pf_target = %.5f\n', Pf_target);

% save stage results
stage_data(stage).estim_ratio = estim_ratio;
stage_data(stage).Pf = Pf_target;
stage_data(stage).u_target = u_target;

% don't produce samples for the next stage if we're at the final stage
if final_stage, return; end

% Reweight/resample
nb_replicates = residual_resampling (g_next);
xt = replicate_particles (xt, nb_replicates);

% visu AFTER resampling (killing !)
bss_figure(21, options.figs, xt, [], [], [stage_data(1:stage-1).u_target], ...
    u_target, problem, stage, Pf_target, ax); drawnow;

[mcs, accept_rate, sig_RW] = bss_move ...
    (xt, [], [], u_target, options, problem, sig_RW);

stage_data(stage).accept_rate = accept_rate;
stage_data(stage).sig_RW = sig_RW;

fprintf('>>>>> expected acceptance rates:  ');
fprintf('%.2f%%  ', accept_rate.expected * 100); fprintf('\n');
fprintf('>>>>> observed acceptance rates:  ');
fprintf('%.2f%%  ', accept_rate.observed * 100); fprintf('\n');

end % function
