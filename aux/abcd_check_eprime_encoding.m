function fname = abcd_check_eprime_encoding(fname,outdir,outstem,forceflag)
%function fname = abcd_check_eprime_encoding(fname,[outdir],[outstem],[forceflag])
%
% Purpose: check whether eprime file is plain text or unicode UTF-16
%   and if so, convert to ASCII
%
% Required parameters:
%   fname: full path name of eprime output file
%
% Optional parameters:
%   outdir: output directory
%     used if conversion to ASCII is required
%     {default = pwd}
%   outdir: output directory
%     used if conversion to ASCII is required
%     {default = pwd}
%   outstem: output file stem
%     if empty, will use filestem of fname
%     {default = []}
%   forceflag: [0|1] overwrite existing output
%     {default = 0}
%
% Created:  12/21/16 by Don Hagler
% Prev Mod: 01/06/17 by Don Hagler
% Prev Mod: 04/08/19 by Feng Xue
% Prev Mod: 07/02/20 by Don Hagler
% Last Mod: 08/28/20 by Don Hagler
% Last Mod: 05/24/24 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;

if ~exist('outdir','var') || isempty(outdir), outdir = pwd; end;
if ~exist('outstem','var'), outstem = []; end;

if ~exist(fname,'file')
  error('file %s not found',fname);
end;

[fpath,fstem,fext] = fileparts(fname);
if isempty(outstem)
  outstem = fstem;
end
outstem = abcd_clean_fstem(outstem);

cmd = sprintf('file "%s"',fname);
[status,result] = mmil_unix(cmd);
if status
  error('cmd %s failed:\n%s',cmd,result);
end;

if ~isempty(regexp(result,'UTF-16'))
  fname_out = sprintf('%s/%s_ASCII%s',outdir,outstem,fext);
  if ~exist(fname_out,'file') || forceflag
    mmil_mkdir(outdir);
    cmd = sprintf('iconv -f UTF-16 -t ASCII --output="%s" "%s"',...
      fname_out,fname);
    [status,result] = mmil_unix(cmd);
    if status
      error('cmd %s failed:\n%s',cmd,result);
    end;
  end;
  fname = fname_out;
end;

