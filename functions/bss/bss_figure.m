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

function bss_figure(fignum, fig_opts, varargin)

%#ok<*DEFNU>

if fig_opts.disable_all, return; end

if ~ isfield (fig_opts, 'color')
    fig_opts.color = true;
else
    assert (islogical (fig_opts.color));
end

dock = fig_opts.dock;
export = fig_opts.export;

if fig_opts.color
    EVAL_NEW_STYLE = {'k*', 'MarkerfaceColor', 'r', 'MarkerSize', 8};
    EVAL_OLD_STYLE = {'ko', 'MarkerfaceColor', 'g', 'MarkerSize', 6};
    CONTOUR_TL = {'Color', 'k', 'LineWidth', 2};
    CONTOUR_FL = {'Color', 'r', 'LineStyle', '-.', 'LineWidth', 1};
    MCSAMPLE = {'b.'};
else
    EVAL_NEW_STYLE = {'k^', 'MarkerfaceColor', 'k', 'MarkerSize', 8};
    EVAL_OLD_STYLE = {'ko', 'MarkerfaceColor', 0.7 * [1 1 1], 'MarkerSize', 6};
    CONTOUR_TL = {'Color', 'k', 'LineWidth', 2};
    CONTOUR_FL = {'Color', 'k', 'LineStyle', '-.', 'LineWidth', 1};
    MCSAMPLE = {'k.'};
end

LABEL_STYLE = {'FontWeight', 'bold', 'FontSize', 10};
AXIS_COL = 0.5 * [1 1 1];

eval(sprintf('figure%2d(varargin{:})', fignum));

% various helper functions%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function bss_raise_fig(fignum, flag_dock)
    
    flag_dock = dock && ((nargin == 1) || flag_dock);
    raise_fig(fignum, flag_dock);
    cla; hold on;
    
    end

    function h = bss_subplot(varargin)
    
    h = subplot(varargin{:});
    cla; hold on;
    set(h, 'Units', 'normalized');
    
    end

% fig_ax_ay %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [ax, ay] = fig_ax_ay(xt, xi, input_distrib)
    
    if (isfield (fig_opts, 'axis')) && (~ isempty (fig_opts.axis))
        ax = fig_opts.axis(1:2);
        ay = fig_opts.axis(3:4);
        return
    end
    
    sigma = std(input_distrib);
    
    ax = [min(xt(:, 1)) max(xt(:, 1))];
    ay = [min(xt(:, 2)) max(xt(:, 2))];
    
    %         if ~isempty(xi),
    %             ax(1) = min(ax(1), min(xi(:, 1)));
    %             ax(2) = max(ax(2), max(xi(:, 1)));
    %             ay(1) = min(ay(1), min(xi(:, 2)));
    %             ay(2) = max(ay(2), max(xi(:, 2)));
    %         end
    
    % enlarge to avoid points on the boundary
    dx = ax(2) - ax(1); ax = ax + 0.02 * dx * [-1 1];
    dy = ay(2) - ay(1); ay = ay + 0.02 * dy * [-1 1];
    
    % turn ellipses into circles !
    dx = (ax(2) - ax(1)) / sigma(1);
    dy = (ay(2) - ay(1)) / sigma(2);
    r = dy / dx;
    if r > 1,
        ddx = (r - 1) * dx;
        ax = ax + 0.5 * ddx * sigma(1) * [-1 1];
    else
        ddy = (1/r - 1) * dy;
        ay = ay + 0.5 * ddy * sigma(2) * [-1 1];
    end
    
    amin = min(input_distrib);
    ax(1) = max(ax(1), amin(1));
    ay(1) = max(ay(1), amin(2));
    
    amax = max(input_distrib);
    ax(2) = min(ax(2), amax(1));
    ay(2) = min(ay(2), amax(2));
    
    end % function fig_ax_ay


    function [xg1, xg2, zg] = compute_xg_zg ...
        (ax, ay, input_distrib, cost_fun)
    
    GRID_NX = 200;
    
    xg1 = linspace (ax(1), ax(2), GRID_NX);
    xg2 = linspace (ay(1), ay(2), GRID_NX);
    [X1, X2] = meshgrid (xg1, xg2);
    x = [X1(:) X2(:)];
    
    b = is_inside_support (input_distrib, x);
    
    zg = nan (GRID_NX^2, 1);
    if iscell (cost_fun)
        xi = cost_fun{1};
        yi = cost_fun{2};
        model = cost_fun{3};
        xt = x(b, :);
        zpred = stk_predict (model, xi, yi, xt);
        zg(b) = zpred.mean;
    else
        zg(b) = cost_fun (x(b, :));
    end
    zg = reshape (zg, GRID_NX, GRID_NX);
    
    end % function compute_xg_zg


    function plot_background_f (xg1, xg2, z, range)
    
    pcolor (xg1, xg2, z);
    shading interp;
    
    if (nargin >= 5) && (~ isempty (range)),
        caxis (range);
    end
    
    colormap (flipud (bone (1024)));
    colorbar;
    
    end % function plot_background_f


% FIGURE 14 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function figure14(varargin) %#ok<DEFNU>
    
    bss_raise_fig(14);
    figure_yp_yt(varargin{:});
    
    end % function figure14

    function figure_yp_yt (yi, yt, yp, u_target_serie, u_target)
    
    ytt = [yi; yt];
    ym = [yi; yp.mean];
    
    plot (ytt, ym, 'b.' );
    
    ymin = min( min(ytt), min(ym) );
    ymax = max( max(ytt), max(ym) );
    plot( [ymin ymax], [ymin ymax], 'r--' );
    
    if nargin > 2
        y1 = min([u_target_serie u_target]);
        y2 = max(max(ytt), max(ym));
        delta = y2 - y1; ax = [y1 - delta, y2 + 0.1*delta];
        axis([ax ax]);
        plot([ax(1) u_target], u_target*[1 1], 'r--' );
        plot(u_target*[1 1], [ax(1) u_target], 'r--' );
        if ~isempty(u_target_serie),
            for j = 1:length(u_target_serie),
                plot([ax(1) u_target_serie(j)], ...
                    u_target_serie(j) * [1 1], 'g--' );
                plot(u_target_serie(j) * [1 1], ...
                    [ax(1) u_target_serie(j)], 'g--' );
            end
        end
    end
    
    xlabel('true f (unknown)', LABEL_STYLE{:});
    ylabel('predicted f', LABEL_STYLE{:});
    hold off;
    
    end % end figure_yp_yt

% FIGURE 15 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function figure15(varargin) %#ok<DEFNU>
    
    bss_raise_fig(15);
    figure_pn_yt(varargin{:});
    
    end % function figure15

    function figure_pn_yt(yt, p_excess, idx_new, u_target)
    
    % plot all sample points
    plot(yt, p_excess, 'b.');
    % threshold
    plot(u_target*[1, 1], ylim, 'r--', 'LineWidth', 2);
    % points that were added during this stage
    plot(yt(idx_new), p_excess(idx_new), EVAL_NEW_STYLE{:});
    
    xlabel('true f (unknown)', LABEL_STYLE{:});
    ylabel('P( f > u_{target} )', LABEL_STYLE{:});
    hold off;
    
    end % function figure_pn_yt

% FIGURE 17 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     function figure17(xt, xi, input_distrib) %#ok<DEFNU>
%
%         if size(xt, 2) ~= 2, return; end
%         bss_raise_fig(17);
%
%         plot(xt(:, 1), xt(:, 2), 'b.'); hold on;
%         plot(xi(:, 1), xi(:, 2), EVAL_NEW_STYLE{:});
%
%         Q = 200; angle = linspace(0, 2*pi, Q)';
%         mu = mean(input_distrib);
%         sigma = std(input_distrib);
%         for k = 1:30,
%             R = sqrt(chi2inv(1 - 10^(-k), 2));
%             xx = repmat(mu, Q, 1) + R * [cos(angle) sin(angle)] * diag(sigma);
%             plot(xx(:, 1), xx(:, 2), '--', 'Color', AXIS_COL);
%         end
%
%         [ax, ay] = fig_ax_ay(xt, input_distrib);
%
%         plot(ax, mu(2)*[1 1], '--', 'Color', AXIS_COL);
%         plot(mu(1)*[1 1], ay, '--', 'Color', AXIS_COL);
%         axis([ax ay]);
%
%         xlabel('x_1', LABEL_STYLE{:});
%         ylabel('x_2', LABEL_STYLE{:});
%         hold off;
%
%     end % function figure17

% Common subfunctions for figure 21, 220, 221, 30 %%%%%%%%%%%%%%%%%%%%%%%%%

    function trace_contour_target_final(xg1, xg2, z, u_target, u_final)
    
    if ~ isequal (u_target, u_final),
        
        % Intermediate stage / target level
        [C, h] = contour (xg1, xg2, z, ...
            [u_target u_target], CONTOUR_TL{:});
        %clabel(C, h);
        
        % Intermediate stage / final level
        [C, h] = contour (xg1, xg2, z, ...
            [u_final u_final], CONTOUR_FL{:});
        %clabel(C, h);
    else
        
        % Final stage / final level
        [C, h] = contour (xg1, xg2, z, ...
            [u_final u_final], CONTOUR_TL{:});
        %clabel(C, h);
        
    end
    
    end % function trace_contour_target_final


    function plot_design_points(xt, xi, idx_new)
    
    xn = double (xt(idx_new, :));
    xi = setdiff (double (xi), xn, 'rows');
    
    plot(xi(:, 1), xi(:, 2), EVAL_OLD_STYLE{:});
    plot(xn(:, 1), xn(:, 2), EVAL_NEW_STYLE{:});
    
    end % function plot_design_points


    function plot_contour_targets_prev (xg1, xg2, z, u_target_prev)
    
    if isscalar (u_target_prev)
        % Workaround to draw a single contour line
        u_target_prev = [u_target_prev u_target_prev];
    end
    
    [C, h] = contour(xg1, xg2, z, ...
        u_target_prev, 'Color', 'k', 'LineStyle', '--');
    % clabel(C, h);
    
    end


    function xylabel ()
    
    if export % prepare for psfrags
        xlabel('x1', LABEL_STYLE{:});
        ylabel('x2', LABEL_STYLE{:});
    else
        xlabel('x_1', LABEL_STYLE{:});
        ylabel('x_2', LABEL_STYLE{:});
    end
    
    end % function xylabel


    function draw_input_distrib_contours (input_distrib)
    
    % FIXME: ceci ne marche que pour les gaussiennes !!!
    
    Q = 200;
    angle = linspace (0, 2 * pi, Q)';
    mu = mean (input_distrib);
    sigma = std (input_distrib);
    for k = 1:30,
        proba = 10^(-k);
        R = sqrt (chi2inv (1 - proba, 2));
        xx = repmat (mu, Q, 1) ...
            + R * [cos(angle) sin(angle)] * (diag (sigma));
        plot (xx(:, 1), xx(:, 2), ':', 'Color', AXIS_COL);
    end
    
    end % function draw_input_distrib_contours


    function eos_title ...
        (fignum, stage, xi, Pf_target_estim, u_target, u_final)
    
    nb_evals = size (xi, 1);
    
    if export,  % do not display title
        
        prefix = sprintf ('%s-fig%d-%d', ...
            fig_opts.save_prefix, fignum, stage);
        
        ppm = get (gcf, 'PaperPositionMode');
        set (gcf, 'PaperPositionMode', 'auto');
        print (gcf, '-depsc', [prefix '.eps']);
        set (gcf, 'PaperPositionMode', ppm);
        
        % Instead of a title, write some info in a txt file
        F = fopen ([prefix '.txt'], 'wt');
        fprintf (F, 'End of stage #%d\n', stage);
        fprintf (F, 'total nb evals = %d\n', nb_evals);
        fprintf (F, 'prob estim = %.2e', Pf_target_estim);
        fprintf (F, 'u_target = %.2e\n', u_target);
        fprintf (F, 'u_final = %.2e\n', u_final);
        fclose (F);
        
    else  % display title
        
        h = title (sprintf (['End of stage #%d  /  %d evals  /  ' ...
            'prob estim = %.2e'], stage, nb_evals, Pf_target_estim));
        set (h, 'FontSize', 12, 'FontWeight', 'b', 'Color', 'r');
        
    end
    
    end % function eos_title


% FIGURE 21 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function figure21(varargin)
    
    xt = varargin{1};
    if size(xt, 2) ~= 2, return; end
    
    bss_raise_fig(21);
    figure_x1_x2_A(varargin{:});
    
    end % function figure21

    function figure_x1_x2_A ...
        (xt, xi, idx_new, u_target_prev, u_target, ...
        problem, stage, Pf_target_estim)
    
    cost_fun      = problem.cost_fun;
    input_distrib = problem.input_distrib;
    u_final       = problem.u_final;
    
    draw_input_distrib_contours (input_distrib);
    
    [ax, ay] = fig_ax_ay(xt, xi, input_distrib);
    
    [xg1, xg2, z] = compute_xg_zg (ax, ay, input_distrib, cost_fun);
    
    % Color background
    if fig_opts.color
        plot_background_f (xg1, xg2, z, fig_opts.range);
    end
    
    %[C, h] = contourf(xg1, xg2, z, 50);
    %colormap(autumn(1024)); colorbar; shading interp;
    
    plot(xt(:, 1), xt(:, 2), 'b.');
    
    plot_contour_targets_prev(xg1, xg2, z, u_target_prev)
    trace_contour_target_final(xg1, xg2, z, u_target, u_final);
    
    mu = mean(input_distrib);
    plot(ax, mu(2)*[1 1], '--', 'Color', AXIS_COL);
    plot(mu(1)*[1 1], ay, '--', 'Color', AXIS_COL);
    
    plot_design_points(xt, xi, idx_new);
    
    axis([ax ay]);  xylabel ();
    
    box off;  box on;  % ???
    
    str = sprintf('stage #%d', stage);
    if ~isempty(xi)
        str = sprintf('%s  /  %d evals', str, size(xi, 1));
    end
    if ~isempty(Pf_target_estim)
        str = sprintf('%s  /  prob estim = %.2e', str, Pf_target_estim);
    end
    title(str, 'FontSize', 12, 'FontWeight', 'b', 'Color', 'b');
    
    %         persistent count
    %         if export
    %             if isempty(count), count = 0; else count = count + 1; end
    %             export_fig(sprintf('%s-fig21-%d.jpg', ...
    %                 fig_opts.save_prefix, count), '-jpeg', '-m3', '-q99');
    %         end
    
    end % function figure_x1_x2_A

% FIGURE 220 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function figure220 (varargin)
    
    xt = varargin{1};
    if size(xt, 2) ~= 2, return; end
    
    bss_raise_fig (220);
    figure_x1_x2_B (varargin{:});
    
    end % function figure220

    function figure_x1_x2_B ...
        (xt, xi, idx_new, u_target_prev, u_target, ...
        problem, stage, Pf_target_estim)
    
    cost_fun      = problem.cost_fun;
    input_distrib = problem.input_distrib;
    u_final       = problem.u_final;
    
    [ax, ay] = fig_ax_ay(xt, xi, input_distrib);
    
    [xg1, xg2, z] = compute_xg_zg (ax, ay, input_distrib, cost_fun);
    
    % Color background (no background in 'bw' mode)
    if fig_opts.color
        plot_background_f (xg1, xg2, z, fig_opts.range);
    end
    
    plot_contour_targets_prev (xg1, xg2, z, u_target_prev)
    trace_contour_target_final (xg1, xg2, z, u_target, u_final);
    plot_design_points (xt, xi, idx_new);
    axis ([ax ay]);  xylabel ();
    
    eos_title (220, stage, xi, Pf_target_estim, u_target, u_final);
    
    end % function figure_x1_x2_B


% FIGURE 221 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function figure221 (varargin)
    
    xt = varargin{1};
    if size(xt, 2) ~= 2, return; end
    
    bss_raise_fig (221);
    figure_x1_x2_C (varargin{:});
    
    end % function figure221

    function figure_x1_x2_C ...
        (xt, xi, idx_new, u_target_prev, u_target, ...
        problem, stage, Pf_target_estim)
    
    cost_fun      = problem.cost_fun;
    input_distrib = problem.input_distrib;
    u_final       = problem.u_final;
    
    [ax, ay] = fig_ax_ay (xt, xi, input_distrib);
    
    [xg1, xg2, z] = compute_xg_zg (ax, ay, input_distrib, cost_fun);
    
    % Color background (no background in 'bw' mode)
    if fig_opts.color
        plot_background_f (xg1, xg2, z, fig_opts.range);
    end
    
    draw_input_distrib_contours (input_distrib);
    
    % Plot Monte Carlo sample
    plot (xt(:, 1), xt(:, 2), MCSAMPLE{:});
    
    trace_contour_target_final (xg1, xg2, z, u_target, u_final);
    axis([ax ay]);  xylabel ();
    
    eos_title (221, stage, xi, Pf_target_estim, u_target, u_final);
    
    end % function figure_x1_x2_C


% COLORFUL SCATTER-PLOTS ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function scatterc (x, z)
    x = double (x);  assert (size (x, 2) == 2);
    scatter (x(:, 1), x(:, 2), [], z, '.');
    end


% FIGURE 23 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function figure23 (x, Jx)
    if size(x, 2) ~= 2,  return;  end  % 2D only
    bss_raise_fig (23);
    
    b = ~ (isinf (Jx) | isnan (Jx));
    if sum(b) > 0,  scatterc (x(b, :), Jx(b));  end
    
    title('sampling criterion');
    end


% FIGURE 24 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function figure24 (x, zpred)
    if size(x, 2) ~= 2,  return;  end   % 2D only
    bss_raise_fig (24);
    
    scatterc (x, zpred.mean);
    
    title('kriging mean');
    end


% FIGURE 28 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function figure28 (x, pmisclass)
    if size(x, 2) ~= 2,  return;  end    % 2D only
    bss_raise_fig (28);
    
    plot (x(:, 1), x(:, 2), 'k.');  % Do we need this one?
    p = max(1e-3, pmisclass);
    scatterc (x, log10 (2 * p));
    
    title('log10(2 * pmisclass)');
    end


% FIGURE 25 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function figure25(stage_history)
    
    bss_raise_fig(25); hold off; cla;
    
    plot(stage_history.u_target,  'ro-'); hold on;
    plot(stage_history.u_final,   'k--');
    
    all_u = [stage_history.u_target(:); stage_history.u_final(:)];
    
    u_min = min(all_u); u_max = max(all_u); delta_u = u_max - u_min;
    if delta_u > 0,
        axis([xlim [u_min u_max] + delta_u * 0.05 * [-1 1]]);
    end
    
    legend('target', 'target2', 'final');
    title('thresholds');
    
    end % function figure25


% FIGURE 26 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function figure26(stage_history)
    
    bss_raise_fig(26); hold off; cla;
    
    semilogy(stage_history.Pf_target_estim,  'ro-'); hold on;
    semilogy(stage_history.Pf_final_estim,   'k--');
    
    legend('target', 'target2', 'final');
    title('probabilities');
    
    end % function figure26


% FIGURE 27 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function figure27(stage_history)
    
    bss_raise_fig(27); hold off; cla;
    
    semilogy(stage_history.stopping_crit, 'ro-'); hold on;
    semilogy(stage_history.stopping_thresh, 'k--');
    
    legend('criterion', 'threshold', 'Location', 'SouthWest');
    title('stopping criterion');
    
    end % function figure27


% FIGURE 29 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function figure29(stage_history)
    
    bss_raise_fig(29); hold off; cla;
    
    gp_err_rel = stage_history.integratedPMisclass ...
        ./ stage_history.Pf_target_estim;
    
    semilogy(sqrt(stage_history.squaredCoV), 'bo--');
    hold on; semilogy(gp_err_rel, 'ro--');
    legend('CoV', 'iPMisclass');
    set(gcf, 'Name', 'squaredCoV');
    
    end % function figure29


% FIGURE 50 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function figure50(stage_history)
    
    bss_raise_fig(50); hold off; cla;
    
    semilogy(stage_history.averagePMisclass, 'ko--');
    
    set(gcf, 'Name', 'averagePMisclass');
    
    end


% FIGURE 30 = 21 + 22 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     function figure30(varargin)
%
%         xt = varargin{2};
%         if size(xt, 2) ~= 2, return; end
%
%         bss_raise_fig(30, false);
%         set(gcf, 'PaperPositionMode', 'auto');
%
%         % resize window to get a decent aspect ratio
%         ss = get(0, 'ScreenSize'); w = ss(3); h = ss(4); ww = 950; hh = 400;
%         set(gcf, 'Position', [floor((w - ww)/2) floor((h - hh)/2) ww hh]);
%
%         h = bss_subplot(1, 2, 1);
%         set(h, 'OuterPosition', [0 0 0.5 1]);
%         figure_x1_x2_A (varargin{:});
%
%         h = bss_subplot(1, 2, 2);
%         set(h, 'OuterPosition', [0.5 0 0.5 1]);
%         figure_x1_x2_B (varargin{:});
%         colorbar('off');
%
%         persistent count
%         if isempty(count), count = 0; else count = count + 1; end
%         print(30, sprintf('fig30-%d.eps', count), '-depsc2', '-painters');
%         %export_fig(sprintf('fig30-%d.png', count), '-png', '-m4');
%
%     end % function figure30

end % function
