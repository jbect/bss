% Copyright Notice
%
%    Copyright (C) 2016 CentraleSupelec
%
%    Author:  Julien Bect  <julien.bect@centralesupelec.fr>
%             Ling Li      <ling.li.supelec@gmail.com>
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
%    right and related or neighboring rights to f_fourbranch.m. This work
%    is published from France.
%
%    License: CC0  <http://creativecommons.org/publicdomain/zero/1.0/>

function z = f_fourbranch (x)

n = size (x, 1);
assert (isequal (size (x), [n 2]));

x = double (x);

z = - min([                                            ...
    3 + 0.1*(x(:,1)-x(:,2)).^2-(x(:,1)+x(:,2))/sqrt(2), ...
    3 + 0.1*(x(:,1)-x(:,2)).^2+(x(:,1)+x(:,2))/sqrt(2), ...
    (x(:,1)-x(:,2)) + 6/sqrt(2),                        ...
    (x(:,2)-x(:,1)) + 6/sqrt(2) ], [], 2);
