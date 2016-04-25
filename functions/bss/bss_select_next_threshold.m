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

function [u_target, g_next] = bss_select_next_threshold ...
    (yp_sample, gx, p0, u_target, u_final)

% update target threshold
[u_target, g_next] = bss_threshold_dichotomy_ ...
    (yp_sample.mean, yp_sample.std, gx, p0, u_target, u_final);

fprintf ('>>> updated target threshold: u_target = %.3f\n', u_target);

% variance monitoring
cv1 = std (g_next ./ gx);
cv2 = std ((yp_sample > u_target) ./ gx);

fprintf ('>>> cv1=%f, cv2=%f\n', cv1, cv2);

end



function [uu_median, g_next] = bss_threshold_dichotomy_...
    (ya, ystd, g_curr, p0, propu, u_final)

if isempty(propu), % initial guess
    propu = quantile(ya, 1 - p0);
end

theta = 1e-9;

% test with u_final, first
g_next = g_fun(u_final, ya, ystd);
condi_prob = mean(g_next ./ g_curr);
if condi_prob > p0,
    uu_median = u_final;
    return;
end

% change p0 if we're close to the final threshold
if condi_prob > p0^2,
    p0 = sqrt(condi_prob);
end

g_next = g_fun(propu, ya, ystd);
condi_prob = mean(g_next ./ g_curr);

if condi_prob == p0
    uu_median = propu;
    return
end

if condi_prob > p0
    prop_u = propu;
    while condi_prob > p0
        prop_u = prop_u + 1;
        g_next = g_fun(prop_u, ya, ystd);
        condi_prob = mean(g_next ./ g_curr);
    end
    prop_l = propu;
elseif condi_prob < p0
    prop_l = propu;
    while condi_prob < p0
        prop_l = prop_l - 1;
        g_next = g_fun(prop_l, ya, ystd);
        condi_prob = mean(g_next ./ g_curr);
    end
    prop_u = propu;
end

uu_median = prop_l + (prop_u - prop_l)/2;
g_next = g_fun(uu_median, ya, ystd);
condi_uu = mean(g_next ./ g_curr);

count_iter = 0;

while (prop_u - prop_l) > theta %abs( condi_uu - p0 ) > theta
    if condi_uu == p0
        break;
    elseif condi_uu > p0
        prop_l = uu_median;
    else %if condi_uu < p0
        prop_u = uu_median;
    end
    uu_median = (prop_l + prop_u)/2;
    g_next = g_fun(uu_median, ya, ystd);
    condi_uu = mean(g_next ./ g_curr);
    
    count_iter = count_iter + 1;
    if count_iter == 500, keyboard, end
end

end % function
