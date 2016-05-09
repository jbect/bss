% Initial design as proposed in our BSS journal paper

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

% Copying Permission Statement  (this file)
%
%    To the extent possible under law, Julien Bect have waived  all copy-
%    right  and related  or neighboring rights to select_example.m.  This
%    work is published from France.
%
%    License: CC0  <http://creativecommons.org/publicdomain/zero/1.0/>

function [problem, range] = select_example (CASE_NUM)

switch CASE_NUM
    
    case {1, 2} %--- Dynamic response of a non-linear oscillator --------------
        
        DIM = 6;
        
        switch CASE_NUM
            
            case 1 % easy case
                
                % Input distribution: IDD Gaussian with
                mu    = [1.0  1.0  0.1  0.5  0.6  1.0 ];
                sigma = [0.05 0.10 0.01 0.05 0.10 0.20];
                
                % Reference probability
                pf_ref = 9.09 * 1e-6;
                
            case 2 % harder case
                
                % Input distribution: IDD Gaussian with
                % (note: the only difference is in the 5th input)
                mu    = [1.0  1.0  0.1  0.5  0.45  1.0 ];
                sigma = [0.05 0.10 0.01 0.05 0.075 0.20];
                
                % Reference probability
                pf_ref = 1.564 * 1e-8;
                
        end
        
        % Limit-state function
        LSS_expr = '- f_dynaresp(x1, x2, x3, x4, x5, x6)';
        
        input_distrib = TruncatedNormalDistribution ...
            (mu, sigma, zeros(1, DIM), inf(1, DIM));
        
        problem = struct(                   ...
            'name',          'dynaresp',    ...
            'dim',           DIM,           ...
            'cost_fun',      @f_dynaresp,   ...
            'input_distrib', input_distrib, ...
            'u_final',       0,             ...
            'Pf_true',       pf_ref         );
        
        
    case 3 % --- cantilever beam ----------------------------------------------
        
        DIM = 2;
        
        MU_1 = 1e-3;  % en Mpa  (1000 Pa)
        MU_2 = 300;   % em mm   [[thicker than the beam in R&S 93 !!!]]
        
        % Input distribution: IDD Gaussian with
        mu    = [MU_1 MU_2];
        sigma = mu .* [0.2 0.10];
        
        % Reference probability
        pf_ref = 3.9406e-06;
        
        input_distrib = TruncatedNormalDistribution ...
            (mu, sigma, zeros(1, DIM), inf(1, DIM));
        
        problem = struct(                   ...
            'name',          'cantilbeam',  ...
            'dim',           2,             ...
            'cost_fun',      @f_cantilbeam, ...
            'input_distrib', input_distrib, ...
            'u_final',       0,             ...
            'Pf_true',       pf_ref         );
        
        % Limit-state function
        LSS_expr = '- f_cantilbeam(x1, x2)';
        
        % fixed range for colorbar
        range = [0 100];
        
        
    case 4 %--- four branch series system -------------------------------------
        
        % Input distribution: IDD Gaussian with
        mu    = [0 0];
        sigma = [1 1];
        
        % Reference probability
        pf_ref = 5.4950e-9;
        
        problem = struct(                                    ...
            'name',          'fourbranch',                   ...
            'dim',            2,                             ...
            'cost_fun',       @f_fourbranch,                 ...
            'input_distrib',  NormalDistribution(mu, sigma), ...
            'u_final',        4,                             ...
            'Pf_true',        pf_ref                         );
        
        % Limit-state function
        LSS_expr = sprintf ('%.20e - f_fourbranch ([x1(:) x2(:)])', ...
            problem.u_final);
                        
        
    case 6 %--- 1D toy example ------------------------------------------------
        
        DIM = 1;
        
        pf_ref = 1e-5;
        
        % Input distribution: IDD Gaussian with
        mu    = 0;
        sigma = 1;
        
        problem = struct(                              ...
            'name',           'gauss1db',              ...
            'dim',            DIM,                     ...
            'cost_fun',       @(x)(x.^2),              ...
            'input_distrib',  NormalDistribution(mu, sigma),           ...
            'u_final',        norminv(1 - pf_ref/2)^2, ...
            'Pf_true',        pf_ref);
        
        % Limit-state function
        LSS_expr = sprintf('%.20e - x1(:) .^ 2', problem.u_final);
        
        
    case 101 %--- Waarts E1 ----------------------------------------------------
        
        % LSS = x1 - x2 (linear problem !)
        
        DIM = 2;
        
        pf_ref = normcdf(0, 5, sqrt(2));
        
        % valeur exacte du beta:     3.535533905932737
        % valeur donnee par Waarts:  3.54
        
        % Input distribution: IDD Gaussian with
        mu    = [7.0 2.0];
        sigma = [1.0 1.0];
        
        problem = struct(                                    ...
            'name',           'WaartsE1',                    ...
            'dim',            DIM,                           ...
            'cost_fun',       @(x)(x(:, 2) - x(:, 1)),       ...
            'input_distrib',  NormalDistribution(mu, sigma), ...
            'u_final',        0,                             ...
            'Pf_true',        pf_ref                         );
        
        % Limit-state function
        LSS_expr = 'x1(:) - x2(:)';
        
        
    case {201, 202, 203} %--- Test scalar distrib ------------------------------
        
        % Test scalar distributions, with x -> x as a cost function.
        % Teh threshold is chosen such that pf_ref = 1e-12;
        
        pf_ref = 1e-12;
        
        switch CASE_NUM
            case 201
                P = NormalDistribution (10, 2);
            case 202
                P = TruncatedNormalDistribution (10, 2, 9, inf);
            case 203
                P = BetaDistribution (2, 5);
        end
        
        name = sprintf ('2xx-%s', class (P));
        
        problem = struct (                              ...
            'name',           name,                     ...
            'dim',            1,                        ...
            'cost_fun',       @(x) x,                   ...
            'input_distrib',  P,                        ...
            'u_final',        quantile (P, 1 - pf_ref), ...
            'Pf_true',        pf_ref                    );
        
        % Limit-state function
        LSS_expr = sprintf ('%.20e - x1(:)', problem.u_final);
        
    otherwise
        
        error('Are you drunk ?');
        
end

if ~ exist ('range', 'var')
    range = [];
end

problem.LSS_expr = LSS_expr;

end
