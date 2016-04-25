% BSS... FIXME: missing documentation
%
% Fields of the 'options' structure:
%
%   .use_gp                    true: BSS, false: plain subset simulation
%
%   .figs                      Figure-related options
%                              (range, disable_all, export, dock, save_prefix)
%
%   .p0                        Conditional proba between successive stages
%   .samplesize                Monte Carlo sample size (a.k.a. population size)
%
%   .MH.init_std_ratio
%   .MH.nb_steps
%   .MH.do_adapt
%   .MH.target_rate
%   .MH.delta_init
%
%   .model                     Random process model
%   .N0                        Size of the initial DoE
%   .update_inside_SUR_loop    Update covariance parameters inside SUR loop ?
%
%   .SUR.m0_max                Maximal value for m0 in the "m0 trick"
%   .SUR.frac                  Paremeter 'rho' in the "m0 trick " (BSS paper)
%   .SUR.max_eval_per_stage    Maximal number of evaluation per stage

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

function stage_data = bss (problem, options, x0, y0)

options = bss_get_default_options (options, problem);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THRESHOLDS... %
%%%%%%%%%%%%%%%%%

% u_final:    final target threshold
%
% u_base:     base threshold for the current stage
%             (the one that has been used to define the current sampling
%             density; equal to -Inf at the first stage)
%
% u_target:   target threshold, will be the base threshold for the next stage
%             (u_target = u_final at the last stage, < u_final otherwise)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZATION %
%%%%%%%%%%%%%%%%%%

% Generate a sample from the input distribution for a start
mcs.xt = random (problem.input_distrib, options.samplesize);
mcs.yt = feval (problem.cost_fun, mcs.xt);    % "cheating" (visu only)
mcs.gx = ones (options.samplesize, 1);        % relative weights
mcs.u_base = -Inf;                            % base threshold

if options.use_gp,  % Specific inits for bss
    
    if nargin < 3,
        [x0, y0] = bss_initial_design (problem);
    end
    
    % Initial design
    obs = struct ('xi', x0, 'yi', y0);
    
    % Estimate covariance parameters on the initial dataset
    model = options.model;
    model.param = stk_param_estim (options.model, x0, y0);
    
elseif nargin > 2,
    
    warning ('Input arguments x0 and y0 are ignored in SubSim mode.');
    
end

stage = 0;  final_stage = false;  stage_data = struct ();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN LOOP %
%%%%%%%%%%%%%

while ~ final_stage,
    
    % a new stage begins !
    stage = stage + 1;
    stage_data(stage).u_base = mcs.u_base;
    stage_data(stage).u_target = [];
    
    fprintf('\n\n####################\n');
    fprintf('##    stage %2d    ##\n', stage);
    fprintf('####################\n\n');
    fprintf('base threshold = %.3f\n\n', mcs.u_base);
    
    if options.use_gp,
        
        [obs, mcs, model, stage_data, final_stage] = bss_stage ...
            (stage, obs, mcs, model, stage_data, problem, options);
        
    else
        
        [mcs, stage_data, final_stage] = bss_subsim_stage ...
            (stage, mcs, stage_data, problem, options);
        
    end
    
end

end % function bss %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bss_stage %
%%%%%%%%%%%%%

function [obs, mcs, model, stage_data, final_stage] = bss_stage ...
    (stage, obs, mcs, model, stage_data, problem, options)

xt = mcs.xt;  yt = mcs.yt;  gx = mcs.gx;
nt = size (xt, 1);  assert (isequal (size (yt), [nt, 1]));

xi = obs.xi;  yi = obs.yi;
ni = size (xi, 1);  assert (isequal (size (yi), [ni, 1]));

final_stage = false;  % is it the final stage ?

u_target = [];

if stage == 1,
    Pf_base = 1.00;
    sig_RW = [];
    squaredCoV_prev = 0;
else
    Pf_base = stage_data(stage - 1).Pf;
    sig_RW = stage_data(stage - 1).sig_RW;
    squaredCoV_prev = stage_data(stage - 1).squaredCoV;
end

if ~ options.update_inside_SUR_loop,
    % update parameters once, before entering the SUR loop
    model.param = stk_param_estim (options.model, xi, yi);
    % TODO: try several starting points ?
end

stage_history.u_target            = [];
stage_history.Pf_target_estim     = [];
stage_history.u_final             = [];
stage_history.Pf_final_estim      = [];
stage_history.stopping_crit       = [];
stage_history.stopping_thresh     = [];
stage_history.squaredCoV          = [];
stage_history.integratedPMisclass = [];
stage_history.averagePMisclass    = [];

% STEPWISE UNCERTAINTY REDUCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('[[stage %d, stepwise uncertainty reduction]]\n', stage);

% Indices of all the points (in xt) selected during this stage
idx_added = [];

k = 0;

while k < options.SUR.max_eval_per_stage,
    
    % PREDICTIONS ON THE MONTE CARLO SAMPLING POINTS %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    yp_sample = stk_predict(model, xi, yi, xt);
    yp_sample = [yp_sample ...
        stk_dataframe(sqrt(yp_sample.var), {'std'})];
    
    % UPDATE TARGET THRESHOLD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ~final_stage,
        
        % At the final stage, u_target = u_final.
        % Otherwise, u_target is a moving target !
        [u_target, g_next] = bss_select_next_threshold ...
            (yp_sample, gx, options.p0, u_target, problem.u_final);
        final_stage = (u_target == problem.u_final);
        
    end
    
    stage_history.u_target(end+1) = u_target;
    stage_history.u_final(end+1)  = problem.u_final;
    
    if final_stage
        fprintf('>>> final stage !  ==>  ');
        fprintf('u_target is now locked to u_final = %.3f.\n', u_target);
    else
        fprintf('>>> u_target = %.3f', u_target);
        fprintf('  <  u_final = %.3f\n', problem.u_final);
        assert(u_target < problem.u_final);
    end
    
    % THRESHOLD EXCEEDANCE PROBABILITIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Current estimate of the proba of exceeding u_target
    %  (u_target is the level we currently trying to reach)
    [q_excess, p_excess] = stk_distrib_normal_cdf ...
        (u_target, yp_sample.mean, yp_sample.std);
    p_misclass = min (p_excess, q_excess);
    Pf_target_estim = mean(p_excess ./ gx) * Pf_base;
    fprintf ('>>> Pf_target_estim = %.3e\n', Pf_target_estim);
    
    % Current estimate of the proba of exceeding u_final
    %  (u_final is the final target level)
    [~, ggg_final] = stk_distrib_normal_cdf ...
        (problem.u_final, yp_sample.mean, yp_sample.std);
    Pf_final_estim = (mean (ggg_final ./ gx)) * Pf_base;
    fprintf ('>>> (Pf_final_estim = %.3e)\n', Pf_final_estim);
    
    % ERROR ESTIMATES FOR THE TARGET THRESHOLD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % We use a recurrence relation on the (squared) coefficient of variation.
    
    % Squared coefficient of variation for the ratio
    %   alphaBSS_t / alphaBSS_{t-1}
    squaredCoV_incr = (Pf_base / Pf_target_estim)^2 ...
        * var(g_next ./ gx) / length(gx);
    
    % Squared coefficient of variation for alphaBSS_t
    squaredCoV = squaredCoV_incr + (1 + squaredCoV_incr) * squaredCoV_prev;
    fprintf ('>>> Estimated CoV = %.2f%%\n', sqrt (squaredCoV) * 100);
    
    integratedPMisclass = Pf_base * (mean (p_misclass ./ gx));
    relativeIntegratedPMisclass = integratedPMisclass / Pf_target_estim;
    fprintf ('>>> Relative integrated PMisclass (RIPM) = %.2f%%\n', ...
        relativeIntegratedPMisclass * 100);
    
    % STOPPING CRITERION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % A measure of the relative error induced by the use of a GP model
    %   about the probability at the target level (u_target)
    stopping_crit = relativeIntegratedPMisclass;
    
    % Compute the stopping threshold for the current stage
    stopping_thresh = options.stopping_thresh ...
        (stage, final_stage, squaredCoV);
    
    % COMPUTE THE SAMPLING CRITERION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [J, yp, stop_flag] = J_SUR_trick ...
        (model, xi, yi, xt, 1 ./ gx, u_target, options.SUR);
    
    % note: we could save some computation time by moving this AFTER the block
    % "DECIDE IF WE KEEP ADDING POINTS" but we would lose some figures
    
    % SAVE SOME INTRA-STAGE HISTORY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    stage_history.Pf_target_estim(k + 1)     = Pf_target_estim;
    stage_history.Pf_final_estim(k + 1)      = Pf_final_estim;
    stage_history.squaredCoV(k + 1)          = squaredCoV;
    stage_history.integratedPMisclass(k + 1) = integratedPMisclass;
    stage_history.averagePMisclass(k + 1)    = mean(p_misclass);
    stage_history.stopping_crit(k + 1)       = stopping_crit;
    stage_history.stopping_thresh(k + 1)     = stopping_thresh;
    
    % FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ~ options.figs.disable_all,
        
        % FIXME: Make it possible to disable figures selectively
        
        if k == 0,
            
            % previous design points + new sample points (without threshold)
            bss_figure(21, options.figs, xt, xi, [], ...
                [stage_data(1:stage-1).u_target], [], ...
                problem, stage, []);
            drawnow;
            
            % previous design points + new sample points (with a threshold)
            bss_figure(21, options.figs, xt, xi, [], ...
                [stage_data(1:stage-1).u_target], u_target, ...
                problem, stage, Pf_target_estim);
            drawnow;
            
        end
        
        %%% Intra-stage history
        % bss_figure (26, options.figs, stage_history);  % probabilities
        % bss_figure (25, options.figs, stage_history);  % thresholds
        % bss_figure (29, options.figs, stage_history);  % squaredCoV
        % bss_figure (50, options.figs, stage_history);  % pmisclass
        % bss_figure (27, options.figs, stage_history);  % stopping criterion
        % drawnow;
        
        %%% Specialized 2D plots
        % bss_figure (23, options.figs, xt, J);           % sampling crit
        % bss_figure (24, options.figs, xt, yp);          % predictor
        % bss_figure (28, options.figs, xt, p_misclass);  % pmisclass
        % drawnow;
        
    end
    
    % DECIDE IF WE KEEP ADDING POINTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Forced stopping if we cannot compute the SUR criterion
    fprintf ('>>> SUR stop_flag = %d  ', stop_flag);
    b0 = (stop_flag > 0);
    if b0
        fprintf (' [MUST STOP]\n');
    else
        fprintf (' [CAN CONTINUE]\n');
    end
    
    % First stopping condition: RIPM must be low enough
    fprintf ('>>> stopping_crit(RIPM)=%.2f%%  ', stopping_crit * 100);
    b1 = (stopping_crit <= stopping_thresh);
    if b1
        fprintf ('<=  stopping_thresh=%.2f%%  [STOP ALLOWED]\n', stopping_thresh * 100);
    else
        fprintf ('>  stopping_thresh=%.2f%%  [MUST CONTINUE]\n', stopping_thresh * 100);
    end
    
    % Second stopping condition: k >= k_min
    k_min = options.SUR.min_eval_per_stage;
    b2 = (k >= k_min);
    if b2
        fprintf ('>>> k=%d >= k_min=%d  [STOP ALLOWED]\n', k, k_min);
    else
        fprintf ('>>> k=%d < k_min=%d [MUST CONTINUE]\n', k, k_min);
    end
    
    if b0 || (b1 && b2)
        fprintf ('>>>   --> all stopping condition are satisfied,\n');
        fprintf ('>>>   --> stop adding evaluation points.\n');
        fprintf ('SUR phase completed ');
        fprintf ('(%d points have been evaluated)\n', k);
        break;
    else
        fprintf ('>>>   --> at least one stopping condition is not satisfied,\n');
        fprintf ('>>>   --> pick a new evaluation point\n');
    end
    
    % PICK A NEW EVALUATION POINT & EVALUATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [ignore_JJ_min, I_new] = min (J); %#ok<ASGLU>
    
    if ismember (xt(I_new, :), xi, 'rows'),
        warning('This point has already been evaluated.'); %#ok<WNTAG>
        keyboard; break;
    end
    
    % Indices of all the points (in xt) selected during this stage
    idx_added = [idx_added; I_new];
    
    xi = [xi; xt(I_new, :)];
    yi = [yi; yt(I_new, :)];
    
    % youplaboum, we have one more evaluation in store
    nb_evals = size(xi, 1);
    assert(isequal(size(xi), [nb_evals problem.dim]));
    assert(isequal(size(yi), [nb_evals 1]));
    
    % FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ~ options.figs.disable_all,
        
        %%% Predicted f versus true f
        % bss_figure (14, options.figs, ...
        %    yi, yt, yp, [stage_data(1:(stage-1)).u_target],  ...
        %    u_target);
        
        %%% PExcess versus true f
        % bss_figure (15, options.figs, ...
        %    yt, p_excess, idx_added, u_target);
        
        bss_figure(21, options.figs, xt, xi, idx_added, ...
            [stage_data(1:stage-1).u_target], u_target, problem, stage, ...
            Pf_target_estim);
        
        drawnow;
        
    end
    
    % UPDATE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if options.update_inside_SUR_loop,
        % update parameters once, before entering the SUR loop
        model.param = stk_param_estim (options.model, xi, yi);
    end
    
    k = k + 1;
    
    fprintf ('>>> nb_evals = %d ', nb_evals);
    fprintf ('(stage=%d, k=%d)\n\n', stage, k);
    
end % while (k <= options.SUR.max_eval_per_stage)

% FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~ options.figs.disable_all,
    
    bss_figure (21, options.figs, xt, xi, [], ...
        [stage_data(1:stage-1).u_target], u_target, problem, stage, ...
        Pf_target_estim); drawnow;
    
    bss_figure (220, options.figs, xt, xi, idx_added, ...
        [stage_data(1:stage-1).u_target], u_target, problem, stage, ...
        Pf_target_estim); drawnow;
    
    bss_figure (221, options.figs, xt, xi, idx_added, ...
        [stage_data(1:stage-1).u_target], u_target, problem, stage, ...
        Pf_target_estim); drawnow;
    
end

% Update the threshold %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf ('\n[[stage %d, updating threshold]]\n', stage);

% target threshold
if ~ final_stage, % update if we're not at the final stage
    fprintf (['>>> The current threshold (stage %d) has been ' ...
        'updated.\n'], stage);
    fprintf (['>>> === SUR was ran using the temporary threshold' ...
        '%.3f.\n'], u_target);
    
    [u_target, g_next] = bss_select_next_threshold ...
        (yp_sample, gx, options.p0, u_target, problem.u_final);
    
    final_stage = (u_target == problem.u_final);
    
    fprintf (['>>> === The next sample will be simulated using the ' ...
        'updated threshold u_target = %.3f.\n'], u_target);
end

stage_data(stage).u_target = u_target;

% update Pf estimate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n[[stage %d, updating Pf estimate...]]\n', stage);

% the estimate the product of all stagewise estimated ratios
if final_stage,
    [~, ggg] = stk_distrib_normal_cdf (problem.u_final, yp_sample.mean, yp_sample.std);
else
    ggg = g_next;
end
estim_ratio = mean(ggg ./ gx);
Pf_target = Pf_base * estim_ratio;

fprintf('>>> estim_ratio = %.5f\n', estim_ratio);
fprintf('>>> Pf_target = %.5f\n', Pf_target);

% save stage results
stage_data(stage).estim_ratio = estim_ratio;
stage_data(stage).Pf = Pf_target;
stage_data(stage).squaredCoV = squaredCoV;
stage_data(stage).n_evals = k;

% don't produce samples for the next stage if we're at the final stage
if final_stage, return; end

% SMC: generate a new MC sample %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n[[stage %d, simulating a new sample]]\n', stage);

g_next = g_fun(u_target, yp_sample.mean, yp_sample.std);

% Reweight/resample
nb_replicates = residual_resampling (g_next ./ gx);
xt_1 = replicate_particles (xt, nb_replicates);

obs = struct('xi', xi, 'yi', yi);

[mcs, accept_rate, sig_RW] = bss_move ...
    (xt_1, obs, model, u_target, options, problem, sig_RW);

fprintf('>>> expected acceptance rates:  ');
fprintf('%.2f%%  ', accept_rate.expected * 100); fprintf('\n');
fprintf('>>> observed acceptance rates:  ');
fprintf('%.2f%%  ', accept_rate.observed * 100); fprintf('\n');

stage_data(stage).accept_rate = accept_rate;
stage_data(stage).sig_RW = sig_RW;

mcs.yt = feval(problem.cost_fun, mcs.xt); % CHEATING

end % function

%#ok<*AGROW>
