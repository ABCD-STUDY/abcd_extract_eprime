function val = mmil_getfield(s,tag,defval)
%function val = mmil_getfield(s,tag,[defval])
%
% Purpose: get field value from struct
%   if not field or isempty, use defval
%
% Required Input:
%   s: struct
%   tag: field name
%
% Optional Input:
%   defval: default value
%
% Created:  04/18/11 by Don Hagler
% Last Mod: 10/19/21 by Don Hagler
%

if ~mmil_check_nargs(nargin,2), return; end;
if ~exist('defval','var'), defval = []; end;

if isfield(s,tag)
  val = getfield(s,tag);
  if isempty(val), val = defval; end
else
  val = defval;
end

