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

BSS_ROOT = fileparts (mfilename ('fullpath'));

addpath (fullfile (BSS_ROOT, 'scripts'));
addpath (fullfile (BSS_ROOT, 'examples'));

FUNC_DIR = fullfile (BSS_ROOT, 'functions');
addpath (FUNC_DIR);
addpath (fullfile (FUNC_DIR, 'distrib'));
addpath (fullfile (FUNC_DIR, 'bss'));

cd (fullfile (FUNC_DIR, 'bss'));
mex replicate_particles.c
cd (BSS_ROOT);

% STK
addpath (fullfile (BSS_ROOT, 'stk'));
stk_init;
