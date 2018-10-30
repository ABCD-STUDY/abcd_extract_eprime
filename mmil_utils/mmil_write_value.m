function mmil_write_value(fid,value)
%function mmil_write_value(fid,value)
%
% Purpose: writes a value to open file specified by fid
%   handles multiple data types, including vectors and cell arrays
%
% Created:  02/25/11 by Don Hagler
% Last Mod: 10/07/13 by Don Hagler
%

if isfloat(value)
  if length(value)~=1, fprintf(fid,'['); end;
  for i=1:length(value)
    if isinf(value(i))
      if sign(value(i))<0
        fprintf(fid,'-Inf');
      else
        fprintf(fid,'Inf');
      end;
    elseif mmil_isint(value(i))
      fprintf(fid,'%d ',value(i));
    else
      fprintf(fid,'%0.8f ',value(i));
    end;
  end;
  if length(value)~=1, fprintf(fid,']'); end;
elseif isnumeric(value) | islogical(value)
  if length(value)==1
    fprintf(fid,'%d',value);
  else
    fprintf(fid,'[%s]',sprintf('%d ',value));
  end;
elseif ischar(value)
  fprintf(fid,'''%s''',value);
elseif isempty(value)
  fprintf(fid,'[]');
elseif iscell(value)
  fprintf(fid,'{');
  for k=1:length(value)
    mmil_write_value(fid,value{k});
    if k<length(value), fprintf(fid,' '); end;
  end;
  fprintf(fid,'}');
else
  error('unrecognized type');
end;

