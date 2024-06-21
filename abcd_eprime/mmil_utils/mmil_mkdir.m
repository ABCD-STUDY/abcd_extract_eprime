function mmil_mkdir(outdir)
%function mmil_mkdir(outdir)
%
% NOTE: gives error if unable to create outdir
%
% Created:  04/04/11 by Don Hagler
% Last Mod: 04/04/11 by Don Hagler
%

[succ,msg] = mkdir(outdir);
if ~succ, error('failed to create output directory %s:\n%s',outdir,msg); end;
