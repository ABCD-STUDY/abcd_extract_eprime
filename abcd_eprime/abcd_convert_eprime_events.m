function [eprime_nruns,eprime_runs,errcode] = abcd_convert_eprime_events(fname_eprime,varargin)
%function [eprime_nruns,eprime_runs,errcode] = abcd_convert_eprime_events(fname_eprime,[options])
%
% Required Input:
%   fname_eprime: eprime file name
%
% Optional Input (key,'value'):
%   'outdir': output directory
%     {default = pwd}
%   'outstem': output file stem
%     if not supplied, will use file stem from fname_eprime
%       for intermediate timing files
%     if supplied, will use this for final output file name
%       and intermediate files
%     {default = []}
%   'TR': repetition time (seconds)
%     {default = 0.8}
%   'numTRs': number of TRs
%     {default = 500}
%   'skipTRs': number of TRs at beginning of each run to ignore
%     {default = 0}
%   'minfrac': minimum fraction of a TR to register an event
%     {default = 0.5}
%   'run': starting run number (e.g., if there were preceding scans)
%     {default = 1}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%
%  Output: 
%    eprime_nruns: number of valid runs in the E-prime file  
%    errcode: [0|1] whether the file was successfully processed
%
% Created:  04/16/24 by Don Hagler
% Prev Mod: 05/01/24 by Don Hagler
% Prev Mod: 05/21/24 by Don Hagler
% Prev Mod: 05/22/24 by Don Hagler
% Prev Mod: 05/24/24 by Don Hagler
% Last Mod: 06/13/24 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eprime_nruns = 0; errcode = 0;
if ~mmil_check_nargs(nargin,1), return; end

% check input
parms = check_input(fname_eprime,varargin);

% extract events
[eprime_nruns,eprime_runs,errcode] = extract_events(parms);
if errcode, return; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fname_eprime,options)
  parms = mmil_args2parms(options,{...
    'fname_eprime',fname_eprime,[],...
    ...
    'outdir',pwd,[],...
    'outstem',[],[],...
    'TR',0.8,[],...
    'numTRs',500,[],...
    'skipTRs',0,[0,1000],...
    'minfrac',0.5,[0,1],...
    'forceflag',false,[false true],...
    ...
    'run',1,[1,100],...
    'taskname',[],[],...
    'verbose',true,[false true],...
  });
  
  % check eprime file
  if ~exist(parms.fname_eprime,'file'), error('file %s not found',parms.fname_eprime); end

  % check outstem
  if isempty(parms.outstem)
    parms.empty_outstem_flag = 1;
    [~,parms.outstem] = fileparts(parms.fname_eprime);
  else
    parms.empty_outstem_flag = 0;
  end
  % remove problematic characters
  parms.outstem = abcd_clean_fstem(parms.outstem);
  
  % check taskname
  if isempty(parms.taskname)
    %% todo: read taskname from eprime file
    error('input taskname required for now');    
  end

  % create outdir
  mmil_mkdir(parms.outdir);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [eprime_nruns,eprime_runs,errcode] = extract_events(parms)
  eprime_nruns = 0;
  errcode = 0;
  % set parameters for extracting from eprime files
  tp = [];
  tp.outstem = parms.outstem;
  tp.outdir = parms.outdir;
  tp.numTRs = parms.numTRs;
  tp.TR = parms.TR;
  tp.minfrac = parms.minfrac;
  tp.nskipTRs = 0; % stimuli start after skipped TRs
  tp.verbose = parms.verbose;
  tp.forceflag = parms.forceflag;
  args = mmil_parms2args(tp);
  % create stimulus timing files
  switch upper(parms.taskname)
    case 'MID'
      [eprime_nruns,eprime_runs,errcode] = abcd_extract_eprime_mid(parms.fname_eprime,args{:});
    case 'SST'
      [eprime_nruns,eprime_runs,errcode] = abcd_extract_eprime_sst(parms.fname_eprime,args{:});
    case 'NBACK'
      [eprime_nruns,eprime_runs,errcode] = abcd_extract_eprime_nback(parms.fname_eprime,args{:});
    otherwise
      error('extract_events not implemented for %s',parms.taskname);
  end
  if errcode, return; end

  if parms.verbose
    fprintf('%s: exporting events...\n',mfilename);
  end
  for r=1:eprime_nruns
    run = parms.run - 1 + r;
    eprime_run = r;
    export_events(run,eprime_run,parms);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function export_events(run,eprime_run,parms)
  if ~mmil_check_nargs(nargin,3), return; end

  % run: run number (e.g., 1 or 2) within session
  if ~exist('run','var') || isempty(run), run = 1; end
  % eprime_run: run number (e.g., 1 or 2) within eprime file
  if ~exist('eprime_run','var') || isempty(eprime_run), eprime_run = 1; end

  % set time for dummy trial
  skip_time = parms.skipTRs * parms.TR;

  % get condition names
  switch parms.taskname
    case 'MID'
      stim_labels = abcd_set_contrasts_mid();
    case 'SST'
      stim_labels = abcd_set_contrasts_sst();
    case 'nBack'
      stim_labels = abcd_set_contrasts_nback();
  end

  if parms.empty_outstem_flag
    fname_out = sprintf('%s/%s_scan%d_events.tsv',...
      parms.outdir,parms.taskname,run);
  else
    fname_out = sprintf('%s/%s_%s_scan%d_events.tsv',...
      parms.outdir,parms.outstem,parms.taskname,run);
  end
  if ~exist(fname_out,'file') || parms.forceflag
    if parms.verbose
      fprintf('%s: exporting events to %s...\n',mfilename,fname_out);
    end
    % find fsl txt files for each condition
    onsets = [];
    durations = [];
    resp_times = [];
    trial_types = {};
    for c=1:length(stim_labels)
      condname = stim_labels{c};
      flist = dir(sprintf('%s/%s_scan%d_%s_fsl.txt',parms.outdir,parms.outstem,eprime_run,condname));
      if isempty(flist)
        if parms.verbose
          fprintf('%s: WARNING: stim file does not exist for run %d %s\n',...
            mfilename,run,condname);
        end
        continue;
      end
      % read file with event timing
      fname_in = sprintf('%s/%s',parms.outdir,flist(1).name);
      data_in = load(fname_in);
      if isempty(data_in)
        if parms.verbose
          fprintf('%s: WARNING: stim file %s is empty\n',...
            mfilename,fname_in);
        end
        continue;
      end
      % find matching rt file
      fname_rt = regexprep(fname_in,'_fsl\.txt','_rt.txt');
      if ~exist(fname_rt,'file')
        if parms.verbose
          fprintf('%s: WARNING: response time file %s not found\n',...
            mfilename,fname_rt);
        end
        data_rt = nan([size(data_in,1),1]);
      else
        data_rt = load(fname_rt);
      end
      tmp_onsets = data_in(:,1) + skip_time;
      tmp_durations = data_in(:,2);
      tmp_resp_times = data_rt;
      tmp_trial_types = repmat({condname},[length(tmp_onsets),1]);
      onsets = cat(1,onsets,tmp_onsets);
      durations = cat(1,durations,tmp_durations);
      resp_times = cat(1,resp_times,tmp_resp_times);
      trial_types = cat(1,trial_types,tmp_trial_types);
    end
    % append dummy event
    onsets = cat(1,onsets,0);
    durations = cat(1,durations,skip_time);
    resp_times = cat(1,resp_times,NaN);
    trial_types = cat(1,trial_types,{'dummy'});
    % sort by onset
    [onsets,ind_sort] = sort(onsets);
    durations = durations(ind_sort);
    resp_times = resp_times(ind_sort);
    trial_types = trial_types(ind_sort);
    % create output tsv file
    fid = fopen(fname_out,'wt');
    if fid<0, error('failed to open %s for writing',fname_out); end
    fprintf(fid,'onset\tduration\tresponse_time\ttrial_type\n');
    for i=1:length(onsets)
      if isnan(resp_times(i))
        fprintf(fid,'%0.3f\t%0.3f\t"n/a"\t"%s"\n',...
          onsets(i),durations(i),trial_types{i});
      else
        fprintf(fid,'%0.3f\t%0.3f\t%0.3f\t"%s"\n',...
          onsets(i),durations(i),resp_times(i),trial_types{i});
      end
    end
    fclose(fid);
  end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


