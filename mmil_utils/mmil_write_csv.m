function mmil_write_csv(fname,data,varargin)
%function mmil_write_csv(fname,data,[options])
%
% Required Input:
%  fname: output file name
%  data: 2D data matrix of numbers or cell array of numbers and/or strings
%
% Optional Parameters:
%   'row_labels': cell array of row labels
%     {default: []}
%   'col_labels': cell array of column labels
%     {default: []}
%   'separator': separator between columns
%     {default: ','}
%   'firstcol_label': label for first column
%     {default: []}
%
% Created:  05/23/09 by Don Hagler
% Last Mod: 12/12/13 by Don Hagler
%   

if ~mmil_check_nargs(nargin,2), return; end;

parms = mmil_args2parms(varargin,{...
  'row_labels',[],[],...
  'col_labels',[],[],...
  'separator',',',[],...
  'firstcol_label',[],[],...
});

[nrows,ncols] = size(data);

if ~isempty(parms.row_labels) & length(parms.row_labels) ~= nrows
  error('number of row labels (%d) does not match first dim of data (%d)',...
    length(parms.row_labels),nrows);
end;
if ~isempty(parms.col_labels) & length(parms.col_labels) ~= ncols
  error('number of col labels (%d) does not match second dim of data (%d)',...
    length(parms.col_labels),ncols);
end;

data_is_cellarray = iscell(data);

fid = fopen(fname,'w');
if ~isempty(parms.col_labels)
  if ~isempty(parms.row_labels)
    fprintf(fid,'%s%s',parms.firstcol_label,parms.separator);
  end;
  for j=1:ncols
    if j>1, fprintf(fid,'%s',parms.separator); end
    fprintf(fid,'"%s"',parms.col_labels{j});
  end;
  fprintf(fid,'\n');
end;
for i = 1:nrows
  if ~isempty(parms.row_labels)
    fprintf(fid,'%s%s',parms.row_labels{i},parms.separator);
  end;
  for j = 1:ncols
    if j>1, fprintf(fid,'%s',parms.separator); end
    if data_is_cellarray, tvar = data{i,j};
    else  tvar = data(i,j);
    end
    if isempty(tvar)
    elseif ischar(tvar)
      fprintf(fid,'"%s"',tvar);
    else
      try
        mmil_write_value(fid,tvar);
      catch me
        error('element {%d,%d} of input data matrix is invalid:\n%s',...
          i,j,me.message);
      end;
    end;
  end
  fprintf(fid,'\n');
end
fclose(fid);

