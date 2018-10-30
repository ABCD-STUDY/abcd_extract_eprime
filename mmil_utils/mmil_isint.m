function res=mmil_isint(v)
%function res=mmil_isint(v)
%
% Purpose: tests whether input value is integer
%   (equal to value rounded to integer)
%
% Created:  04/01/09 by Don Hagler
% Last Mod: 02/16/11 by Don Hagler
%
%  copied from Ziad Saad's isint
%

if any(~isfinite(v))
  res = 0;
  return;
end;

fv = fix(v);
df = v - fv;
if (~df)
  res = 1;
else
  res = 0;
end;

