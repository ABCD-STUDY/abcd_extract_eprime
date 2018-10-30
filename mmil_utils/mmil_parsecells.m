function celldat = mmil_parsecells(celldat,rows,cols)
%function celldat = mmil_parsecells(celldat,rows,cols)
%
% Purpose: parse cell matrix of numbers and strings
%   e.g., the output of mmil_readtext
%
% Required Input:
%   celldat: cell matrix of unknown types
%
% Optional Input:
%   rows: vector of row numbers to be parsed
%     {default = 'all'}
%   cols: vector of column numbers to be parsed
%     {default = 'all'}
%
% Output: cell matrix of numbers, strings, cell arrays, numeric arrays
%
% Created:  04/01/09 by Jason Sherfey
% Prev Mod: 10/07/13 by Don Hagler
% Last Mod: 01/30/18 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2, rows = 'all'; end
if nargin < 3, cols = 'all'; end

if ~isnumeric(rows) || max(rows) > size(celldat,1), rows = 'all'; end
if ~isnumeric(cols) || max(cols) > size(celldat,2), cols = 'all'; end  
 
if strcmp(rows,'all'), rows = 1:size(celldat,1); end
if strcmp(cols,'all'), cols = 1:size(celldat,2); end

celldat = celldat(rows,cols);

nr = size(celldat,1);
nc = size(celldat,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:nr
  for j = 1:nc
    dat = celldat{i,j};
    try if isnan(dat), continue; end; end
    if isnumeric(dat), continue; end
    if ischar(dat)
      tmp = dat;
      try
        if any(regexp(dat,'[\[\]]'))
          tmp = eval(dat);
        end
      end
      if isnumeric(tmp),    celldat = updatecell(celldat,i,j,tmp); continue; end
      if strcmp(dat,'pwd'), celldat = updatecell(celldat,i,j,pwd); continue; end
      [dat status] = iscellarr(dat);
      if status, celldat = updatecell(celldat,i,j,dat); continue; end
      [dat status] = isnumarr(dat);
      if status, celldat = updatecell(celldat,i,j,dat); continue; end
      % do not convert to numeric if negative sign is not at beginning
      if ~isempty(regexp(dat,'-')) && isempty(regexp(strtrim(dat),'^-')), continue; end;
      % convert to numeric if all characters are numbers, spaces, brackets, periods, or negative signs
      if length(regexp(dat,'[\d\s\[\].-]'))==length(dat) && ~isempty(str2num(dat))
        celldat = updatecell(celldat,i,j,str2num(dat));
        continue;
      end      
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dat status] = iscellarr(dat)
status = 0;
if isempty(regexp(dat,'^{.*}','match')), return; end
val = regexp(dat,'[^,{}'']+','match');
keep = {};
for n = 1:length(val)
  [tmpval s] = isnumarr(val{n});
  if s==1, keep{end+1} = tmpval; continue; end
  if s==2, keep(end+1:end+length(tmpval)) = tmpval; continue; end
  if isempty(regexp(val{n},'\w')), continue; end
  keep(end+1) = regexp(val{n},'[^'']+','match');
end
dat = keep;
status = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dat status] = isnumarr(dat)
status = 0;
brk = regexp(dat,'[\[^\]$]');
len = length(brk);
if len < 2 || mod(len,2)~=0, return; end
if ~(brk(1)==1 && brk(end)==length(dat)), return; end
if len > 2
  tmp = dat;
  dat = {};
  for n = 1:2:len
    dat{end+1} = str2num(tmp(brk(n):brk(n+1)));
  end
else
  dat = str2num(dat);
end
if isnumeric(dat), status = 1; end
if iscell(dat),    status = 2; end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function data = updatecell(data,i,j,value)
data{i,j} = value;
