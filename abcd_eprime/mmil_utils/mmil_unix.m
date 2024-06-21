function [status,result,errmsg]=mmil_unix(cmd,batch_limit,jflag)
%function [status,result,errmsg]=mmil_unix(cmd,[batch_limit],[jflag])
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
%   jflag: [0|1] if 1, use jsystem; if 0, use unix
%     {default = 1}
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
% Prev Mod: 04/26/21 by Don Hagler
% Last Mod: 07/27/23 by Don Hagler
%

status = 0; result = []; errmsg = [];

if ~mmil_check_nargs(nargin,1), return; end;
if ~exist('batch_limit','var') || isempty(batch_limit), batch_limit = 250; end;
if ~exist('jflag','var') || isempty(jflag), jflag = 1; end;

if isempty(cmd), error('cmd is empty'); end;
%shell = '/bin/tcsh -c';
shell = 'tcsh -c';

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
  if jflag
    [status,result,errmsg] = jsystem(cmd,shell);
  else
    [status,result] = unix(cmd);
  end    
else
  nbatches = ceil(ncmds/batch_limit);
  for b=1:nbatches
    j = 1 + batch_limit*(b-1);
    k = min(j + batch_limit - 1,ncmds);
    tmp_cmds = cmd_list(j:k);
    tmp_cmd = sprintf('%s\n',tmp_cmds{:});
    if jflag
      [status,tmp_result,tmp_errmsg] = jsystem(tmp_cmd,shell);
    else
      [status,tmp_result] = unix(tmp_cmd);
      tmp_errmsg = [];
    end
    if isempty(result)
      result = tmp_result;
    elseif ~isempty(tmp_result)
      result = sprintf('%s\n%s',result,tmp_result);
    end
    if isempty(errmsg)
      errmsg = tmp_errmsg;
    elseif ~isempty(tmp_errmsg)
      errmsg = sprintf('%s\n%s',errmsg,tmp_errmsg);
    end
    if status, break; end;
  end;
end;

% restore LD_LIBRARY_PATH
setenv('LD_LIBRARY_PATH',path);

if nargout==2
  result = [result errmsg];
end;

