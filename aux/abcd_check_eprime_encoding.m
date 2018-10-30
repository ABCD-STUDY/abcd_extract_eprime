function fname = abcd_check_eprime_encoding(fname,outdir,forceflag)
%function fname = abcd_check_eprime_encoding(fname,[outdir],[forceflag])
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
%   forceflag: [0|1] overwrite existing output
%     {default = 0}
%
% Created:  12/21/16 by Don Hagler
% Last Mod: 01/06/17 by Don Hagler
%

if ~mmil_check_nargs(nargin,1), return; end;

if ~exist('outdir','var') || isempty(outdir), outdir = pwd; end;

if ~exist(fname,'file')
  error('file %s not found',fname);
end;

[fpath,fstem,fext] = fileparts(fname);

cmd = sprintf('file %s',fname);
[s,r] = unix(cmd);
if s
  error('cmd %s failed:\n%s',cmd,r);
end;

if ~isempty(regexp(r,'UTF-16'))
  fname_out = sprintf('%s/%s_ASCII%s',outdir,fstem,fext);
  if ~exist(fname_out,'file') || forceflag
    mmil_mkdir(outdir);
    cmd = sprintf('iconv -f UTF-16 -t ASCII --output=%s %s',...
      fname_out,fname);
    [s,r] = unix(cmd);
    if s
      error('cmd %s failed:\n%s',cmd,r);
    end;
  end;
  fname = fname_out;
end;

