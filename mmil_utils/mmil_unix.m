function [status,result,errmsg]=mmil_unix(cmd,batch_limit)
%function [status,result,errmsg]=mmil_unix(cmd,[batch_limit])
%
% Purpose: execute a unix command
%   moves MATLAB paths to bottom of LD_LIBRARY_PATH stack
%   to simulate non-MATLAB shell environment.
%
% Required Input:
%   cmd : UNIX command-line to execute
%
% Optional Input:
%   batch_limit: number of lines to execute at a time
%     {default = 250}
%
% Output:
%   status: 0 if successful, 1 if error
%   result: stdout (or stdout +  stderr if nargout==2)
%   errmsg: stderr
%
% See also: unix, jsystem
%
% Created:  08/01/07 by Ben Cipollini (as smartunix)
% Prev Mod: 08/05/11 by Don Hagler
% Prev Mod: 07/15/19 by Feng Xue
% Prev Mod: 07/02/20 by Don Hagler
% Last Mod: 07/06/20 by Don Hagler
%

status = 0; result = []; errmsg = [];

if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('batch_limit','var'), batch_limit = 250; end;
if isempty(cmd), error('cmd is empty'); end;
shell = '/bin/tcsh -c';

% remove MATLAB directories from LD_LIBRARY_PATH
path = getenv('LD_LIBRARY_PATH');
paths = regexp(path,':','split');
newpaths = paths(cellfun(@isempty,regexp(paths,'matlab|/mcr/')));
if isempty(newpaths)
  newpath = '';
else
  newpath = [sprintf('%s:',newpaths{1:end-1}) newpaths{end}];
end;
setenv('LD_LIBRARY_PATH',newpath);

% split command into multiple lines
cmd_list = regexp(cmd,'\n','split');
ncmds = length(cmd_list);
if ncmds<=1
  [status,result,errmsg] = jsystem(cmd,shell);
else
  nbatches = ceil(ncmds/batch_limit);
  for b=1:nbatches
    j = 1 + batch_limit*(b-1);
    k = min(j + batch_limit - 1,ncmds);
    tmp_cmds = cmd_list(j:k);
    tmp_cmd = sprintf('%s\n',tmp_cmds{:});
    [status,tmp_result,tmp_errmsg] = jsystem(tmp_cmd,shell);
    result = [result tmp_result];
    errmsg = [errmsg tmp_errmsg];
    if status, break; end;
  end;
end;

% restore LD_LIBRARY_PATH
setenv('LD_LIBRARY_PATH',path);

if nargout==2
  result = [result errmsg];
end;

