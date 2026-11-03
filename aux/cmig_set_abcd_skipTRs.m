function skipTRs = cmig_set_abcd_skipTRs(scan_info,varargin)
%function skipTRs = cmig_set_abcd_skipTRs(scan_info,[options])
%
% Purpose: set skipTRs based on scanner manufacturer and software version
%          from BIDS scan JSON sidecar fields
%
% Required Input:
%   scan_info: struct loaded from BIDS scan JSON sidecar
%     relevant fields: Manufacturer, SoftwareVersions
%
% Optional Input ('key',value):
%   'skipTRs_GE_DV25': skipTRs for GE scanners with DV25 or Orchestra software
%     {default = 4}
%   'skipTRs_GE': skipTRs for GE scanners with other software versions
%     {default = 15}
%   'skipTRs_Philips': skipTRs for Philips scanners
%     {default = 8}
%   'skipTRs_Siemens': skipTRs for Siemens scanners
%     {default = 8}
%   'fname_events': BIDS _events.tsv, used only when Manufacturer is missing
%     or unrecognized. Its first row is a dummy event whose duration is the
%     time the skipped TRs occupy, so skipTRs is that duration over the TR
%     {default = []}
%
% Output:
%   skipTRs: number of TRs to skip at the start of each run
%
% Created:  05/11/2026 by Don Hagler
% Prev Mod: 05/11/2026 by Don Hagler
% Prev Mod: 09/03/2026 by Don Hagler
% Last Mod: 10/03/2026 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Based on:
%   abcd_set_skipTRs
%     Created:  11/23/2020 by Don Hagler
%     Last Mod: 06/20/2025 by Don Hagler

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

skipTRs = 0;
if ~mmil_check_nargs(nargin,1), return; end

parms = mmil_args2parms(varargin,{...
  'skipTRs_GE_DV25',4,[],...
  'skipTRs_GE',15,[],...
  'skipTRs_Philips',8,[],...
  'skipTRs_Siemens',8,[],...
  'fname_events',[],[],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

manufacturer = '';
if isfield(scan_info,'Manufacturer'), manufacturer = scan_info.Manufacturer; end
sw_ver = '';
if isfield(scan_info,'SoftwareVersions'), sw_ver = scan_info.SoftwareVersions; end

if ~isempty(regexpi(manufacturer,'GE'))
  if ~isempty(regexp(sw_ver,'DV25|Orchestra','once'))
    skipTRs = parms.skipTRs_GE_DV25;
  else
    skipTRs = parms.skipTRs_GE;
  end
elseif ~isempty(regexpi(manufacturer,'Philips'))
  skipTRs = parms.skipTRs_Philips;
elseif ~isempty(regexpi(manufacturer,'Siemens'))
  skipTRs = parms.skipTRs_Siemens;
else
  % no usable Manufacturer, as with the release's mmps_mproc sidecars. The
  %   dummy event carries the same information: its duration is the time the
  %   skipped TRs occupy, and duration over TR reproduces every branch above
  skipTRs = skipTRs_from_dummy(parms.fname_events,scan_info);
  if isempty(skipTRs)
    fprintf('%s: WARNING: unsupported manufacturer "%s" and no usable dummy event; setting skipTRs = 0\n',...
      mfilename,manufacturer);
    skipTRs = 0;
  else
    fprintf('%s: manufacturer "%s" not usable; skipTRs = %d from the dummy event\n',...
      mfilename,manufacturer,skipTRs);
  end
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function skipTRs = skipTRs_from_dummy(fname_events,scan_info)
  % duration of the events file's dummy row, divided by the TR. Found by its
  %   trial_type: in old-style files it is the first row, in new-style files
  %   scanner_trigger can come first. Empty if that cannot be established, so
  %   the caller can fall back to 0 with a warning rather than a wrong number
  skipTRs = [];
  if isempty(fname_events) || ~exist(fname_events,'file'), return; end
  if ~isfield(scan_info,'RepetitionTime') || isempty(scan_info.RepetitionTime), return; end
  TR = scan_info.RepetitionTime;
  if ~isnumeric(TR) || TR <= 0, return; end
  try
    tbl = mmil_readtext(fname_events,'\t');
  catch
    return;
  end
  if size(tbl,1) < 2 || size(tbl,2) < 4, return; end
  cols = tbl(1,:);
  i_dur  = find(strcmp(cols,'duration'),1);
  i_type = find(strcmp(cols,'trial_type'),1);
  if isempty(i_dur) || isempty(i_type), return; end
  i_row = [];
  for k = 2:size(tbl,1)
    trial_type = regexprep(num2str_if_needed(tbl{k,i_type}),'"','');
    if strcmpi(strtrim(trial_type),'dummy'), i_row = k; break; end
  end
  if isempty(i_row), return; end
  dur = tbl{i_row,i_dur};
  if ~isnumeric(dur), dur = str2double(dur); end
  if ~isfinite(dur) || dur < 0, return; end
  nTRs = dur / TR;
  % a non-integer means the assumption does not hold for this file
  if abs(nTRs - round(nTRs)) > 1e-6, return; end
  skipTRs = round(nTRs);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function str = num2str_if_needed(val)
  if ischar(val), str = val; else str = num2str(val); end
return;
