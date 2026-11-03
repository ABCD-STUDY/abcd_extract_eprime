function [eprime_nruns,eprime_runs,errcode,behav,errmsg] = abcd_extract_eprime_sst(fname,varargin)
%function [eprime_nruns,eprime_runs,errcode,behav,errmsg] = abcd_extract_eprime_sst(fname,[options])
%
% Purpose: extract condition time courses
%   from eprime data files for SST task
%
% Required input:
%   fname: name of input eprime file
%
% Optional input:
%   'outdir': output directory
%     {default = pwd}
%   'outstem': output file stem
%     if empty, will use filestem of fname
%     {default = []}
%   'TR': scan repetition time (s)
%     {default = 0.8}
%   'numTRs': number of repetitions
%     {default = 500}
%   'minfrac': minimum fraction of a TR to register an event
%     {default = 0.5}
%   'nskipTRs': number of dummy TRs before the first retained frame; sets
%     start_time after a single trigger (the timing files are relative to
%     start_time and are not trimmed)
%     {default = []: 0 for a behavioral-only experiment (detected from
%      ExperimentName); else computed from Manufacturer/SoftwareVersion
%      via cmig_set_abcd_skipTRs() if given; else 0 with a warning}
%   'Manufacturer': DICOM/BIDS Manufacturer field, used to compute
%     nskipTRs when nskipTRs is [] (ignored otherwise)
%     {default = []}
%   'SoftwareVersion': DICOM SoftwareVersion(s) field, used to compute
%     nskipTRs when nskipTRs is [] (ignored otherwise)
%     {default = []}
%   'switch_thresh': percent correct below which to swap button assignements
%     {default = 0.3}
%   'min_valid_trials': minimum number of trials with both a response and a
%     correct-response value before a button-switch decision is attempted;
%     below this, the accuracy estimate is too noisy to act on and no
%     switch is attempted
%     {default = 10}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%   'timing_files_flag': [0|1] generate timing files
%     {default = 1}
%   'ste_recode_flag': [] auto-detect from SST experiment version (true
%     if the first digit is >=2, else false); or explicit true/false to
%     override the auto-detected value unconditionally (e.g. to compare
%     both interpretations for the same session)
%     {default = []}
%   'verbose': [0|1] display messages
%     {default = 1}
%
% Output: 
%    eprime_nruns: number of valid runs in the e-prime in the file  
%    errcode: [0|1] whether the file was successfully processed
%    behav: behavioral data
%    errmsg: string describing error if errcode=1
%
% Created:  10/07/16 by Don Hagler
% Prev Mod: 02/12/19 by Dani Cornejo
% Prev Mod: 03/12/19 by Feng Xue
% Prev Mod: 08/05/19 by Octavio Ruiz
% Prev Mod: 05/03/21 by Don Hagler
% Prev Mod: 08/31/21 by Don Hagler
% Prev Mod: 01/12/22 by Octavio Ruiz
% Prev Mod: 02/10/22 by Octavio Ruiz
% Prev Mod: 04/22/22 by Don Hagler
% Prev Mod: 06/14/24 by Don Hagler
% Prev Mod: 06/25/24 by Don Hagler
% Prev Mod: 08/19/25 by Bader Chaarani
% Prev Mod: 08/21/25 by Don Hagler
% Prev Mod: 09/15/25 by Don Hagler
% Prev Mod: 06/29/26 by Don Hagler
% Prev Mod: 07/06/26 by Don Hagler
% Prev Mod: 07/08/26 by Don Hagler
% Prev Mod: 07/23/26 by Don Hagler
% Prev Mod: 07/24/26 by Don Hagler
% Prev Mod: 07/27/26 by Don Hagler
% Prev Mod: 07/28/26 by Don Hagler
% Prev Mod: 07/29/26 by Don Hagler
% Prev Mod: 07/30/26 by Don Hagler
% Prev Mod: 07/31/26 by Don Hagler
% Prev Mod: 10/05/26 by Don Hagler
% Last Mod: 10/21/26 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE: Based on make-SST-refs.py
%         provided by Michael Michael Riedel, Ph.D. from FIU
%       Using style from abcd_extract_eprime_mid.m
%         provided by Donald J Hagler, Ph.D.
%       Includes code for glitch and violator variables by Octavio Ruiz
%         adapted from python code provided by Sage Hahn from UVM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize outputs
eprime_nruns = []; eprime_runs = []; errcode = 0; behav = []; errmsg = [];

% check arguments 
if ~mmil_check_nargs(nargin,1), return; end;

% check input parameters
parms = check_input(fname,varargin); 

% create output directory
mmil_mkdir(parms.outdir);

% create structure with info for each event
[event_info,start_time,errcode,errmsg,parms] = get_event_info(parms);
if errcode, return; end;

% auto-detect ste_recode_flag from SST experiment version if not
% explicitly specified (true if the first digit is >=2, else false)
if isempty(parms.ste_recode_flag)
  version_str = mmil_getfield(event_info(1),'version','');
  version_parts = regexp(version_str,'\.','split');
  if ~isempty(version_parts) && ~isnan(str2double(version_parts{1})) && ...
      str2double(version_parts{1})>=2
    parms.ste_recode_flag = true;
  else
    parms.ste_recode_flag = false;
  end
end

% check for no responses
[nresp, errcode, errmsg] = check_no_responses(event_info,parms);
if errcode, return; end;
        
% switch buttons if necessary 
%   file _switched.csv will be created 
[event_info,parms.switch_flag,errcode,errmsg] = sst_switch(event_info,parms); 
if errcode, return; end;

% replace and/or add events names that were not coded
%   correctly by E-prime
%   e.g., late or early responses
event_info = revise_sst_types(event_info,parms);

% get behav data and write to csv file 
[behav,eprime_runs,event_info,errcode,errmsg] = get_behavioral_data(event_info,parms);
if errcode, return; end;

parms.eprime_runs = eprime_runs; 
eprime_nruns = length(eprime_runs);
parms.eprime_nruns = eprime_nruns; 
if parms.verbose, fprintf('%s: number of e-prime runs = %d \n',mfilename,eprime_nruns); end

% when STE trials are kept as a nuisance condition rather than recoded
% as Go, give them a real timing condition instead of no representation
% at all (see revise_sst_types's ste_recode_flag=false branch). Appended
% here, after get_behavioral_data() has already run, so it does not also
% pick up early_stop/STE and duplicate/collide with the dedicated
% STE_counts_total/STE_percent_total fields get_behavioral_data() already
% computes directly from ssd_dur/ssd_rt.
if ~parms.ste_recode_flag
  parms.condnames{end+1} = 'early_stop';
  parms.typenames{end+1} = 'STE';
end
% zero-SSD stop trials with an SSD-window response (ZeroSSDResp), neither go
%   nor stop, are a nuisance condition like early_stop; appended here for the
%   same reason (counted by get_behavioral_data() as ZeroSSDResp_counts)
parms.condnames{end+1} = 'zero_ssd_response';
parms.typenames{end+1} = 'ZeroSSDResp';
parms.nconds = length(parms.condnames);
parms.ntypes = length(parms.typenames);

if parms.timing_files_flag && eprime_nruns
  % create files for each condition
  for i=1:parms.nconds
    type  = parms.typenames{i};
    cond  = parms.condnames{i};

    % find all instances of current type of event
    ind_type = find(strcmp(type,{event_info.type}));
    % raw run, same source as behavioral counting (not a timestamp guess)
    run = [event_info(ind_type).run];

    if ismember(cond,{'early_stop','zero_ssd_response'})
      % STE trials kept as a nuisance condition (ste_recode_flag=false),
      % and ZeroSSDResp trials; 'early_stop' would otherwise regex-extract
      % eventname='stop' below and collide with correct_stop/incorrect_stop's
      % own logic
      parms.mindur = 0.08;
      % get onset
      onset = [event_info(ind_type).('stop_onset_time')];
      % set offset from the full event span (see revise_sst_types)
      offset = onset + [event_info(ind_type).('ste_dur')];
      % get resp_times: response occurred during the SSD window
      resp_times = [event_info(ind_type).('ssd_rt')];
    else
      pat = '\w_(?<event>\w+)';
      k = regexp(cond,pat,'names');
      eventname = k.event;

      switch eventname
        case 'go'
          parms.mindur = 0.08;
          % get onset
          onset = {event_info(ind_type).('go_onset_time')};
          % get offset
          offset = {event_info(ind_type).('go_offset_time')};
          % get resp_times
          resp_times = {event_info(ind_type).('go_rt_adjusted')};
          % check offsets match onsets
          [onset,offset,resp_times,run] = check_offsets(onset,offset,resp_times,run,eventname,parms);
        case 'stop'
          parms.mindur = 0;
          % get onset
          onset = [event_info(ind_type).('stop_start_time')];
          % set offset
          offset = onset+20;
          % get resp_times
          resp_times = [event_info(ind_type).('stop_rt_adjusted')];
      end
    end
    % assign each event to its own (raw run) run's timing files
    [ind_start,event_ref_time] = assign_run(run,start_time);
    % create files for each scan
    write_files(cond,onset,offset,resp_times,ind_start,event_ref_time,parms);
  end 
end; %timing_files_flag

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fname,options)
  parms = mmil_args2parms(options,{...
    'fname',fname,[],...
    ...
    'outdir',pwd,[],...
    'outstem',[],[],...
    'TR',0.8,[],...
    'numTRs',500,[],...
    'minfrac',0.5,[0,1],...
    'nskipTRs',[],[],...
    'Manufacturer',[],[],...
    'SoftwareVersion',[],[],...
    'switch_thresh',0.3,[0,0.5],...
    'min_valid_trials',10,[],...
    'forceflag',false,[false true],...
    'timing_files_flag',true,[false true],...
    'verbose',true,[false true],...
    'ste_recode_flag',[],[],...
    ...
    'eprime_runs',1:2,[],...
    'eprime_nruns',0,0:2,...
    'expected_numtrials',182,[],...
    'expected_start_delay_GE',12,[],...
    'expected_start_delay_Siemens',6.4,[],...
    'expected_init_delay_GE',0,[],...
    'expected_init_delay_Siemens',0,[],...
    'start_delay_Siemens',6.4,[],...
    'max_start_diff',0.5,[],...
    'max_init_diff',5.0,[],...
    ...
    'colnames',  {'NARGUID','SessionDate','SessionTime','ExperimentName','ExperimentVersion',...
                  'Subject','Procedure[Block]','Block','Trial',...
                  'GetReady.RTTime','Wait4Scanner.RTTime','Wait4ScannerB.RTTime','BeginFix.StartTime','BeginFix.OnsetTime',...
                  'Go.RT', 'Fix.RT', 'StopSignal.RT','Go.Duration', 'SSDDur',...
                  'StopSignal.Duration','SSD.RT', 'Go.OnsetTime', 'Go.OffsetTime',  'Go.ACC' , 'Go.RESP',...
                  'Fix.RESP', 'SSD.OnsetTime','StopSignal.StartTime','SSD.OffsetTime','StopSignal.ACC', 'Go.CRESP',...
                  'StopSignal.OnsetTime','StopSignal.OffsetTime','StopSignal.RESP', 'SSD.ACC','TrialCode', 'SSD.RESP', 'Stimulus' },[],...
    'fieldnames',{'narguid','date','time','experiment','version',...
                  'subject','procedure_block','block','trial',...
                  'getready_rttime','wait1_rttime','wait2_rttime','beginfix_st','beginfix_onset',...
                  'go_rt','fix_rt','stop_rt','go_dur','ssd_dur',...
                  'stop_dur','ssd_rt','go_onset_time','go_offset_time' , 'go_acc','go_resp',...
                  'fix_resp', 'ssd_onset_time','stop_start_time','ssd_offset_time','stop_acc','go_cresp',...
                  'stop_onset_time','stop_offset_time','stop_resp', 'ssd_acc', 'type', 'ssd_resp', 'stimulus' },[],...
    'typenames',{'CorrectGo', 'CorrectLateGo', 'IncorrectGo', 'IncorrectLateGo','NoRespGo',...
                 'CorrectStop', 'IncorrectStop'},[],...
    'condnames',{'correct_go','correctlate_go', 'incorrect_go', 'incorrectlate_go' ...
                  'noresp_go', 'correct_stop', 'incorrect_stop'},[],...
  });

  % check conditions
  parms.nconds = length(parms.condnames);
  parms.ntypes = length(parms.typenames);
  parms.mindur = parms.minfrac*parms.TR;
  if parms.nconds ~= parms.ntypes
    error('condnames and typenames length mismatch');
  end;
  % check input file
  if ~exist(parms.fname,'file')
    error('file %s not found',parms.fname);
  end;
  [fdir,fstem,fext] = fileparts(parms.fname);
  % remove problematic characters
  if isempty(parms.outstem)
    parms.outstem = abcd_clean_fstem(fstem);
  end;
  % compute time at the start/end of each TR. Direct multiplication
  % (0:numTRs-1)*TR, not linspace(TR,numTRs*TR,numTRs) -- linspace's
  % formula introduces a ~1 ULP rounding difference per grid point
  % relative to k*TR (confirmed: at k=302, TR=0.8, linspace gives
  % 241.59999999999999, direct multiplication gives 241.60000000000002),
  % which occasionally flips the strict '<'/'>' TR-mask comparison right
  % at a trial-onset tie exactly on a TR boundary -- common since MID/
  % SST/NBACK trial timing is often designed to land there. Python's
  % grid (np.arange(numTRs)*TR + TR) already used direct multiplication;
  % matching it here removes the discrepancy.
  parms.TR_onset  = (0:parms.numTRs-1) * parms.TR;
  parms.TR_offset = parms.TR_onset + parms.TR;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [event_info,start_time,errcode,errmsg,parms] = get_event_info(parms)

  event_info=[]; start_time=[]; errcode = 0; errmsg = [];
  
  try 
    % write event info to file
    fname_csv = abcd_check_eprime_sprdsh(parms.fname, parms.colnames,...
                  parms.fieldnames, parms.outdir, parms.outstem, parms.forceflag, parms.verbose);
    event_info = mmil_csv2struct(fname_csv);
  catch me
    if parms.verbose
      fprintf('%s: ERROR: failed to read or interpret e-prime file %s:\n%s\n',...
        mfilename,parms.fname,me.message);
    end
    errcode = 1;
    errmsg = 'failed to read or interpret e-prime file';
    return;
  end
  
  % check experiment
  experiment = mmil_getfield(event_info(1),'experiment');
  if isempty(regexpi(experiment,'sst'))
    if parms.verbose
      fprintf('%s: ERROR: wrong experiment name in e-prime file %s: %s\n',...
        mfilename,parms.fname,experiment);
    end
    errcode = 1;
    errmsg = 'wrong experiment name';
    return;
  else
    if parms.verbose, fprintf('%s: experiment name: %s\n',mfilename,experiment); end
  end;

  % check for behavioral or GE-specific experiment
  behav_flag = 0; ge_flag = 0;
  if ~isempty(regexpi(experiment,'behavioral'))
    behav_flag = 1;
    expected_start_delay = [0,0];
    expected_init_delay = [0,0];
  elseif ~isempty(regexp(experiment,'_GE'))
    ge_flag = 1;
    expected_start_delay = parms.expected_start_delay_GE*[1,1];
    expected_init_delay = parms.expected_init_delay_GE*[1,1];
  else
    expected_start_delay = parms.expected_start_delay_Siemens*[1,1];
    expected_init_delay = parms.expected_init_delay_Siemens*[1,1];
  end

  % resolve nskipTRs=[]: 0 for a behavioral-only experiment; else computed
  % from Manufacturer/SoftwareVersion if given; else 0 with a warning.
  % Missing imaging data is not special-cased -- it isn't behavioral, so it
  % takes this same path. Resolved here (not in the caller) so it's
  % available below for the single-trigger start_time formula.
  if isempty(parms.nskipTRs)
    if behav_flag
      parms.nskipTRs = 0;
    elseif ~isempty(parms.Manufacturer) || ~isempty(parms.SoftwareVersion)
      % cmig_set_abcd_skipTRs() expects char fields -- an unsupplied
      % Manufacturer/SoftwareVersion is [] (empty double) here, not ''
      manufacturer = parms.Manufacturer; if isempty(manufacturer), manufacturer = ''; end
      software_version = parms.SoftwareVersion; if isempty(software_version), software_version = ''; end
      scan_info = struct('Manufacturer',manufacturer,...
                          'SoftwareVersions',software_version);
      parms.nskipTRs = cmig_set_abcd_skipTRs(scan_info);
    else
      parms.nskipTRs = 0;
      if parms.verbose
        fprintf('%s: WARNING: nskipTRs not supplied and no Manufacturer/SoftwareVersion given -- defaulting nskipTRs=0\n',mfilename);
      end
    end
  end

  % get trigger and start times
  trig_time = []; start_time = [];
  if isfield(event_info,'getready_rttime') % GE
    % get initial trigger times
    ind_trig = find(~cellfun(@isempty,{event_info.getready_rttime}));
    ind_trig_first = ind_trig(find([100,diff(ind_trig)]>1));
    trig_time = [event_info(ind_trig_first).getready_rttime];
    % exclude zeros
    trig_time = nonzeros(trig_time)';
    % NOTE: for GE, 16 triggers are sent/recorded
    %       set start_time to last trigger
    ind_trig_last = ind_trig(find([diff(ind_trig),100]>1));
    start_time = [event_info(ind_trig_last).getready_rttime];
    % exclude zeros
    start_time = nonzeros(start_time)';
  elseif isfield(event_info,'wait1_rttime') || isfield(event_info,'wait2_rttime') % Siemens/Philips
    % get initial trigger times
    trig_time = [];
    if isfield(event_info,'wait1_rttime')
      ind_trig = find(~cellfun(@isempty,{event_info.wait1_rttime}));
      if ~isempty(ind_trig)
        trig_time = event_info(ind_trig(1)).wait1_rttime;
      end
    end
    if isfield(event_info,'wait2_rttime')
      ind_trig = find(~cellfun(@isempty,{event_info.wait2_rttime}));
      if ~isempty(ind_trig)
        trig_time = [trig_time,event_info(ind_trig(1)).wait2_rttime];
      end
    end
    % exclude zeros
    trig_time = nonzeros(trig_time)';
    % set start time relative to trigger time    
    if behav_flag
      % NOTE: for behavioral only, only one trigger is sent/recorded
      %       set start_time to trigger
      start_time = trig_time;
    else
      % NOTE: for single-trigger files, use nskipTRs*TR instead of the
      %       fixed start_delay_Siemens constant (numerically the same
      %       for genuine Siemens/Philips hardware)
      start_time = trig_time + 1000*parms.nskipTRs*parms.TR;
    end
  end

  % get init_time
  ind_init = find(~cellfun(@isempty,{event_info.beginfix_onset}));
  init_time = [event_info(ind_init).beginfix_onset];
  % exclude zeros
  init_time = nonzeros(init_time)';

  % remove extra triggers
  if length(trig_time) > length(init_time)
    trig_time = trig_time(1:length(init_time));
  end  
  % remove extra starts
  if length(start_time) > length(init_time)
    start_time = start_time(1:length(init_time));
  end

  if ge_flag
    % check whether start delay diff is an integer multiple of TR
    %   indicating missed triggers
    if length(start_time)==1
      ex_start_delay = expected_start_delay(1);
    else
      ex_start_delay = expected_start_delay;
    end
    start_delay_diff = (start_time - trig_time) - 1000*ex_start_delay;
    start_delay_diff_TR = start_delay_diff / (1000 * parms.TR);
    start_delay_diff_TR_int = round(start_delay_diff_TR);
    idx_miss = find(ismember(start_delay_diff_TR_int,[1:3]) & abs(start_delay_diff_TR - start_delay_diff_TR_int) <= 0.01);
    if ~isempty(idx_miss)
      if parms.verbose
        for i=1:length(idx_miss)
          fprintf('%s: WARNING: adjusting start_time for run %d by %0.2f seconds due to indications of missed triggers\n',...
            mfilename,idx_miss(i),start_delay_diff(idx_miss(i))/1000);
        end
      end
      % adjust start_time by subtracting start_delay_diff
      start_time(idx_miss) = start_time(idx_miss) - start_delay_diff(idx_miss);
      %% NOTE: start_time represents the start of the non-dummy imaging volume
      %%       init_time is when the stimulus starts
    end
  end

  % check for disparity between trig_time and init_time
  abcd_check_eprime_timing(trig_time,start_time,init_time,...
    expected_start_delay,expected_init_delay,...
    parms.outdir,parms.outstem,...
    parms.max_start_diff,parms.max_init_diff,...
    parms.verbose);

  % remove non-events
  ind_events = find(~cellfun(@isempty,{event_info.type})); 
  event_info = event_info(ind_events); 

  % check go_resp
  go_resp = replace_cells({event_info.go_resp});
  fix_resp = replace_cells({event_info.fix_resp});
  stop_resp = replace_cells({event_info.stop_resp});
  ssd_resp = replace_cells({event_info.ssd_resp});
  go_cresp = replace_cells({event_info.go_cresp});

  if any(cellfun(@isstr,go_resp)) || any(cellfun(@isstr,fix_resp)) ||...
     any(cellfun(@isstr,stop_resp)) || any(cellfun(@isstr,ssd_resp)) ||...
     any(cellfun(@isstr,go_cresp))
    if ~any(cellfun(@isstr,go_cresp))
      if parms.verbose
        fprintf('%s: ERROR: string go_resp values without string go_cresp values in e-prime file %s\n',...
          mfilename,parms.fname);
      end
      errcode = 1;
      errmsg = 'string go_resp values without string go_cresp values';
      return;
    end;
    % get response numbers and names for unique, non-emptpy correct responses
    uniq_go_cresp = unique(go_cresp(find(~cellfun(@isempty,go_cresp))));
    num_cresp = length(uniq_go_cresp);
    uniq_go_cresp_nums = [];
    uniq_go_cresp_names = [];
    for i=1:length(uniq_go_cresp)
      %   allow either '1{LEFTARROW}' or '1,{LEFTARROW}'
      k = regexp(uniq_go_cresp{i},'(?<num>\d+),?{(?<name>\w+)}','names');
      if ~isempty(k)
        uniq_go_cresp_nums(i) = str2num(k.num);
        uniq_go_cresp_names{i} = k.name;
      else
        % check for single character r, g, b, or y
        k = regexp(uniq_go_cresp{i},'^(?<char>[rgby])$','names');
        if ~isempty(k)
          uniq_go_cresp_nums(i) = i;
          uniq_go_cresp_names{i} = k.char;
        else
          % check for 'LEFTARROW' or 'RIGHTARROW'
          if all(ismember(uniq_go_cresp,{'LEFTARROW','RIGHTARROW'}))
            uniq_go_cresp_nums(i) = find(strcmp(uniq_go_cresp{i},{'LEFTARROW','RIGHTARROW'}));
            uniq_go_cresp_names{i} = uniq_go_cresp{i};
          else
            if parms.verbose
              fprintf('%s: ERROR: string go_cresp with unexpected pattern (%s) in e-prime file %s\n',...
                mfilename,uniq_go_cresp{i},parms.fname);
            end
            errcode = 1;
            errmsg = sprintf('string go_cresp with unexpected pattern (%s)',uniq_go_cresp{i});
            return;
          end
        end
      end        
    end
    if parms.verbose
      fprintf('%s: WARNING: replacing go_resp strings with numeric\n', mfilename);
    end;
    for i=1:length(event_info)
      % assign numbers to each go_resp/fix_resp/stop_resp/ssd_resp/go_cresp
      event_info(i).go_resp = resolve_resp_value(go_resp{i},uniq_go_cresp_names,uniq_go_cresp_nums);
      event_info(i).fix_resp = resolve_resp_value(fix_resp{i},uniq_go_cresp_names,uniq_go_cresp_nums);
      event_info(i).stop_resp = resolve_resp_value(stop_resp{i},uniq_go_cresp_names,uniq_go_cresp_nums);
      event_info(i).ssd_resp = resolve_resp_value(ssd_resp{i},uniq_go_cresp_names,uniq_go_cresp_nums);
      event_info(i).go_cresp = resolve_resp_value(go_cresp{i},uniq_go_cresp_names,uniq_go_cresp_nums);
    end
  end
  
  % reset variables indicating accurate responses
  event_info = sst_acc(event_info);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vals = replace_cells(vals)
  if any(cellfun(@iscell,vals))
    ind_cell = find(cellfun(@iscell,vals));
    vals(ind_cell) = cellfun(@(x) x{1},vals(ind_cell),'UniformOutput',false);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = resolve_resp_value(val,uniq_go_cresp_names,uniq_go_cresp_nums)
  if ~isstr(val), return; end
  % try extracting an embedded digit directly, via the SAME regex used to
  % build uniq_go_cresp_nums/names, applied to the value itself -- handles
  % a self-referential value already in 'N,{NAME}'/'N{NAME}' format (e.g.
  % go_cresp when Go.RESP is plain numeric but Go.CRESP is 'N,{NAME}');
  % such a value is never bracket-unwrapped upstream, so it never matches
  % a bare name via strcmp below
  k = regexp(val,'(?<num>\d+),?{(?<name>\w+)}','names');
  if ~isempty(k)
    val = str2num(k.num);
    return;
  end
  % fall back to name lookup (bare {NAME} unwrapped upstream, r/g/b/y, or
  % bare LEFTARROW/RIGHTARROW)
  k = find(strcmp(val,uniq_go_cresp_names));
  val = uniq_go_cresp_nums(k);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check for no responses
function [nresp, errcode, errmsg] = check_no_responses(event_info,parms)
  nresp = 0; errcode = 0; errmsg = [];
  % check go_resp, fix_resp, and stop_resp
  go_resp = {event_info.go_resp};
  fix_resp = {event_info.fix_resp};
  stop_resp = {event_info.stop_resp};
  ssd_resp = {event_info.ssd_resp};
  % count number of responses (any type)
  nresp = nnz(~cellfun(@isempty,go_resp)) +...
          nnz(~cellfun(@isempty,fix_resp)) +...
          nnz(~cellfun(@isempty,stop_resp)) +...
          nnz(~cellfun(@isempty,ssd_resp));
  % if no responses, set errcode and errmsg
  if nresp == 0
    if parms.verbose, fprintf('%s: ERROR: no responses in %s\n',mfilename,parms.fname); end
    errcode = 1;
    errmsg = 'no responses';
    return;
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function event_info = sst_acc(event_info)
  for i=1:length(event_info)
    if event_info(i).go_resp == event_info(i).go_cresp
      event_info(i).go_acc=1;
    elseif event_info(i).go_resp ~= event_info(i).go_cresp
      event_info(i).go_acc=0;
    end
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [event_info,switch_flag,errcode,errmsg] = sst_switch(event_info,parms)
  [switch_flag,acc1,errcode,errmsg] = sst_switch_flag(event_info,parms);
  if errcode, return; end;

  if switch_flag
    codes = derive_switch_codes({event_info.go_cresp});
    if isempty(codes)
      if parms.verbose
        fprintf('%s: go_cresp does not yield a unique response-code pair, skipping switch\n',mfilename);
      end
      switch_flag = 0;
    else
      event_info = sst_switch_event(event_info,codes);
      [~,acc2,~] = sst_switch_flag(event_info,parms);
      if acc2 <= acc1
        event_info = sst_switch_event(event_info,codes);
        switch_flag = 0;
      elseif ~isequal(codes,[1 2]) && parms.verbose
        fprintf('%s: switched using derived codes %g<->%g\n',mfilename,codes(1),codes(2));
      end
    end
  end

  % write to csv
  if switch_flag 
    fname_csv_out  = sprintf('%s/%s_events_switched.csv',...
        parms.outdir,parms.outstem); 
    if ~exist(fname_csv_out,'file') || parms.forceflag
      mmil_struct2csv(event_info,fname_csv_out)
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function event_info = sst_switch_event(event_info,codes)
  for i=1:length(event_info)
    event_info(i).go_resp = switch_val(event_info(i).go_resp,codes);
    event_info(i).fix_resp = switch_val(event_info(i).fix_resp,codes);
    event_info(i).stop_resp = switch_val(event_info(i).stop_resp,codes);
    event_info(i).ssd_resp = switch_val(event_info(i).ssd_resp,codes);
  end
  event_info = sst_acc(event_info);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = switch_val(val,codes)
  if val==codes(1)
    val=codes(2);
  elseif val==codes(2)
    val=codes(1);
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function codes = derive_switch_codes(cresp_cell)
  % Derive the two response codes to flip from the correct-response
  % field's own unique values (device/site wiring varies, e.g. some
  % Philips sites use 3/4 instead of 1/2), rather than assuming 1/2.
  % Returns [code_a,code_b] when exactly 2 unique non-empty values exist;
  % [] when indeterminate (0, 1, or more than 2 values).
  vals = cresp_cell(~cellfun(@isempty,cresp_cell));
  uniq = unique(cell2mat(vals));
  if numel(uniq)==2
    codes = uniq;
  else
    codes = [];
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [switch_flag,accuracy,errcode,errmsg] = sst_switch_flag(event_info,parms)
  switch_flag = 0;
  errcode = 0; errmsg = [];
  accuracy = NaN;

  % go_acc may retain a stale raw Go.ACC value for trials where go_resp is
  % empty (sst_acc only overwrites go_acc when go_resp and go_cresp are
  % both non-empty), so count trials with a real response/correct-response
  % pair directly rather than trusting go_acc's denominator.
  go_resp = {event_info.go_resp};
  go_cresp = {event_info.go_cresp};
  ind_valid = find(~cellfun(@isempty,go_resp) & ~cellfun(@isempty,go_cresp));
  n_valid = length(ind_valid);
  if n_valid < parms.min_valid_trials
    if parms.verbose, fprintf('%s: only %d valid trials (< %d), skipping switch\n',...
      mfilename,n_valid,parms.min_valid_trials); end
    return;
  end

  correct = 0;
  for i=ind_valid
    if isequal(go_resp{i},go_cresp{i})
      correct = correct+1;
    end
  end
  correct_total = n_valid;
  accuracy = 100*correct/correct_total;
  if parms.verbose, fprintf('%s: accuracy = %0.1f%%\n',mfilename,accuracy); end
  
  if accuracy < 100*parms.switch_thresh
    if parms.verbose, fprintf('%s: accuracy < %0.1f%%, switching button responses\n',...
      mfilename,100*parms.switch_thresh); end
    switch_flag = 1;
  end;  
return 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function event_info = revise_sst_types(event_info,parms)
  % correct_go
  go_cresp = {event_info.go_cresp};
  go_resp = {event_info.go_resp};  
  correct_go_ind = [];
  for i=1:size(go_cresp,2)
    if go_resp{i} == go_cresp{i}
      correct_go_ind = [correct_go_ind i]; 
      event_info(i).type = 'CorrectGo';
    end 
  end
  % correct_go_late 
  fix_resp = {event_info.fix_resp}; 
  correct_go_late_ind = [];
  for i=1:size(go_cresp,2)  
    if fix_resp{i} == go_cresp{i} 
      if go_resp{i} == go_cresp{i}
      else 
        correct_go_late_ind = [correct_go_late_ind i]; 
        event_info(i).type = 'CorrectLateGo';
     end 
    end 
  end
  % no_resp_go 
  go_event_ind = find(~cellfun(@isempty,go_cresp)); 
  noresp_ind = find(cellfun(@isempty,go_resp)); 
  nofix_ind = find(cellfun(@isempty,fix_resp));
  noresp_nofix_ind = [intersect(noresp_ind,nofix_ind)];  
  no_resp_go_ind = [intersect(noresp_nofix_ind,go_event_ind)];   
  for i=no_resp_go_ind
    event_info(i).type = 'NoRespGo';
  end
  % incorrect_go
  incorrect_go_ind = [];
  incorrect_go_late_ind = [];
  for i=1:size(go_cresp,2)  
    if go_resp{i} ~= go_cresp{i}
      incorrect_go_ind = [incorrect_go_ind i]; 
      event_info(i).type = 'IncorrectGo'; %legit 1~=2
    elseif go_resp{i} == go_cresp{i}
    else 
      if fix_resp{i} ~= go_cresp{i}
        incorrect_go_late_ind = [incorrect_go_late_ind i]; 
        event_info(i).type = 'IncorrectLateGo';
      end 
    end 
  end
  % correct & incorrect _stop
  ssd_rt = {event_info.ssd_rt};
  ssd_dur = {event_info.ssd_dur};
  stop_resp = {event_info.stop_resp};  
  fix_resp = {event_info.fix_resp};
  ssd_resp = {event_info.ssd_resp};
  % note how many STE trials will be recoded as Go trials
  if parms.verbose
    idx_val = find(~cellfun(@isempty, ssd_rt));
    idx_nonzero = cell2mat(ssd_rt(idx_val)) > 0;
    idx_ste = idx_val(idx_nonzero);
    if ~isempty(idx_ste)
      if parms.ste_recode_flag
        fprintf('%s: recoding %d STE trials as Go trials\n',mfilename,length(idx_ste));
      else
        fprintf('%s: labeling %d STE trials as early_stop (nuisance)\n',mfilename,length(idx_ste));
      end
    end
  end
  % ensure ste_dur exists as a field on every event, even when this file
  % has zero STE or ZeroSSDResp trials -- MATLAB's dynamic field access below
  % (event_info(ind_type).('ste_dur')) errors if the field was never
  % assigned anywhere in the struct array, not just for the empty indices
  [event_info.ste_dur] = deal([]);
  for i=1:size(ssd_dur,2)
    if ssd_dur{i} > 0  % non-zero SSD stop trial
      if ssd_rt{i} > 0
        % Stop Too Early (STE): participant responded during the SSD window
        if parms.ste_recode_flag
          % legacy behavior: recode STE trial as Go trial
          event_info(i).go_resp = event_info(i).ssd_resp;
          event_info(i).go_rt = event_info(i).ssd_rt;
          event_info(i).go_onset_time = event_info(i).ssd_onset_time;
          %% NOTE: missing SSD.OffsetTime, so it is unclear what is correct value for this
          event_info(i).go_offset_time = event_info(i).stop_offset_time;
          %% NOTE: unclear what should be used for Go duration for these trials
          event_info(i).go_dur = event_info(i).stop_dur + event_info(i).ssd_dur;
          event_info(i).stop_dur = [];
          % determine Go.CRESP from the stimulus image
          stimulus_path = event_info(i).stimulus;
          if contains(stimulus_path, 'LeftArrow')
            event_info(i).go_cresp = 1;
          elseif contains(stimulus_path, 'RightArrow')
            event_info(i).go_cresp = 2;
          end
          % recalculate Go.ACC to determine if the Go trial is correct or incorrect
          if event_info(i).go_resp == event_info(i).go_cresp
            event_info(i).go_acc = 1; % Correct Go
          else
            event_info(i).go_acc = 0; % Incorrect Go
          end
          % recode the trial type to reflect the new Go status
          if event_info(i).go_acc == 1
             event_info(i).type = 'CorrectGo';
          else
             event_info(i).type = 'IncorrectGo';
          end
        else
          % default: label as nuisance type; exclude from both go and stop counts
          % preserve the full STE event span (SSD window + stop signal)
          % before stop_dur is cleared, for early_stop's own timing files
          event_info(i).ste_dur = event_info(i).stop_dur + event_info(i).ssd_dur;
          event_info(i).stop_dur = [];
          event_info(i).type = 'STE';
        end

      elseif stop_resp{i} > 0
        event_info(i).type = 'IncorrectStop';
      elseif fix_resp{i} > 0
        event_info(i).type = 'IncorrectStop';
      else
        event_info(i).type = 'CorrectStop';
      end
    elseif isequal(ssd_dur{i}, 0)
      if ~isempty(ssd_resp{i})
        % zero-SSD trial with an SSD-window response (ZeroSSDResp): neither
        % go nor stop; stop_dur cleared to leave it out of the stop trials,
        % its event span kept in ste_dur for its own timing files
        event_info(i).ste_dur = event_info(i).stop_dur + event_info(i).ssd_dur;
        event_info(i).stop_dur = [];
        event_info(i).type = 'ZeroSSDResp';
      elseif stop_resp{i} > 0
        event_info(i).type = 'IncorrectStop';
      elseif fix_resp{i} > 0
        event_info(i).type = 'IncorrectStop';
      else
        event_info(i).type = 'CorrectStop';
      end
    end
  end
  
  % write output file
  fname_csv_out = sprintf('%s/%s_events_revised.csv',parms.outdir,parms.outstem);
  if ~exist(fname_csv_out,'file') || parms.forceflag
    mmil_struct2csv(event_info,fname_csv_out);
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [ind_start,event_ref_time] = assign_run(run,start_time)
  % uses each event's own raw run number directly, instead of picking
  % the nearest-preceding start_time (which could silently move an event
  % into a different run's timing files)
  if ~isempty(run)
    ind_start = run;
    event_ref_time = start_time(run);
  else
    ind_start = [];
    event_ref_time = [];
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [onset,offset,resp_times,run] = check_offsets(onset,offset,resp_times,run,eventname,parms)
  ind_empty = find(cellfun(@isempty,onset) | cellfun(@isempty,offset));
  if ~isempty(ind_empty)
    if parms.verbose
      fprintf('%s: WARNING: %s event has %d onsets and %d offsets\n',...
        mfilename,eventname,nnz(~cellfun(@isempty,onset)),nnz(~cellfun(@isempty,offset)));
      % NOTE: this may indicate file truncation
    end
    ind_keep = setdiff([1:length(onset)],ind_empty);
    onset = onset(ind_keep);
    offset = offset(ind_keep);
    resp_times = resp_times(ind_keep);
    run = run(ind_keep);
  end;
  onset = cell2mat(onset);
  offset = cell2mat(offset);
  if iscell(resp_times), resp_times = cell2mat(resp_times); end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function write_files(eventname,onset,offset,resp_times,ind_start,event_ref_time,parms)
 nscans = parms.eprime_nruns;
 scans = parms.eprime_runs;  
  
  for s=scans
    ind_scan = find(ind_start == s);
    fname_1D = sprintf('%s/%s_scan%d_%s.1D',...
      parms.outdir,parms.outstem,s,eventname);
    fname_txt = sprintf('%s/%s_scan%d_%s.txt',...
      parms.outdir,parms.outstem,s,eventname);
    fname_block = sprintf('%s/%s_scan%d_%s_block.txt',...
      parms.outdir,parms.outstem,s,eventname);  
    fname_block_dur = sprintf('%s/%s_scan%d_%s_block_dur.txt',...
      parms.outdir,parms.outstem,s,eventname);  
    fname_fsl = sprintf('%s/%s_scan%d_%s_fsl.txt',...
      parms.outdir,parms.outstem,s,eventname);
    fname_rt = sprintf('%s/%s_scan%d_%s_rt.txt',...
      parms.outdir,parms.outstem,s,eventname);

    if ~exist(fname_1D,'file') || ~exist(fname_txt,'file') || ...
       ~exist(fname_block,'file') || ~exist(fname_fsl,'file') ||...
       ~exist(fname_rt,'file') || parms.forceflag
      if ~isempty(onset(ind_scan))
        % subtract start time and convert to seconds
        rel_onset  = (onset(ind_scan) - event_ref_time(ind_scan))/1000;
        rel_offset = (offset(ind_scan) - event_ref_time(ind_scan))/1000;
        % get resp_times for this scan
        rel_resp_times = resp_times(ind_scan)/1000;
        % loop over events
        TR_mask = zeros(parms.numTRs,1);
        for k=1:length(rel_onset)
          ind_rep = find(rel_onset(k)<(parms.TR_offset-parms.mindur) &...
                         rel_offset(k)>(parms.TR_onset+parms.mindur));
          TR_mask(ind_rep) = 1;
        end;
        % no TRs removed: rel_onset is from start_time, already the first
        %   retained frame (nskipTRs x TR after a single trigger), so the
        %   grid starts there; trimming nskipTRs more (as from 07/31/26)
        %   put the .1D and .txt nskipTRs TRs early
        % calculate event times from TR_onset (0-based from start_time)
        event_times = (find(TR_mask)-1)*parms.TR;
        % convert time units to multiples of TRs
        event_TRs = event_times/parms.TR;
        % write to 1D file
        fid = fopen(fname_1D,'wt');
        if fid<0, error('failed to open %s for writing',fname_1D); end;
        for k=1:length(TR_mask)
          fprintf(fid,'%d\n',TR_mask(k));
        end;
        fclose(fid);
        % write to txt file with stim times in seconds
        fid = fopen(fname_txt,'wt');
        if fid<0, error('failed to open %s for writing',fname_txt); end;
        fprintf(fid,'%s\n',sprintf('%0.2f ',event_TRs*parms.TR));
        fclose(fid);
        % write to block txt file
        fid = fopen(fname_block,'wt');  
        if fid<0, error('failed to open %s for writing',fname_block); end;    
        fprintf(fid,'%s\n',sprintf('%0.2f ',rel_onset));
        fclose(fid);
        fid = fopen(fname_block_dur,'wt');  
        if fid<0, error('failed to open %s for writing',fname_block_dur); end;
        block_len = mean(rel_offset- rel_onset); 
        fprintf(fid,'%s\n',sprintf('%0.2f ',block_len));
        fclose(fid);
        % write to fsl file
        fid = fopen(fname_fsl,'wt');
        if fid<0, error('failed to open %s for writing',fname_fsl); end;
        for i=1:length(rel_onset)
          fprintf(fid, '%0.2f %0.2f 1 \n',rel_onset(i), rel_offset(i)-rel_onset(i));
        end
        fclose(fid);
        % write to rt file
        fid = fopen(fname_rt,'wt');
        if fid<0, error('failed to open %s for writing',fname_rt); end;
        for i=1:length(rel_resp_times)
          fprintf(fid, '%0.3f\n',rel_resp_times(i));
        end
        fclose(fid);
      else
        TR_mask = zeros(parms.numTRs,1);
        % write to 0s 1D file
        fid = fopen(fname_1D,'wt');
        if fid<0, error('failed to open %s for writing',fname_1D); end;
        for k=1:length(TR_mask)
          fprintf(fid,'%d\n',TR_mask(k));
        end;
        fclose(fid);
        % write * txt file
        fid = fopen(fname_txt,'wt');
        if fid<0, error('failed to open %s for writing',fname_txt); end;
        fprintf(fid,'*\n');
        fclose(fid);   
        % write * txt file block 
        fid = fopen(fname_block,'wt');
        if fid<0, error('failed to open %s for writing',fname_block); end;
        fprintf(fid,'*\n');
        fclose(fid);  
        fid = fopen(fname_block_dur,'wt');
        if fid<0, error('failed to open %s for writing',fname_block); end;
        fprintf(fid,'0\n');
        fclose(fid); 
        % write empty fsl file
        fid = fopen(fname_fsl,'wt');
        if fid<0, error('failed to open %s for writing',fname_fsl); end;
        fclose(fid);    
        % write empty rt file
        fid = fopen(fname_rt,'wt');
        if fid<0, error('failed to open %s for writing',fname_rt); end;
        fclose(fid);    
      end
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [behav,runs_ok,event_info,errcode,errmsg] = get_behavioral_data(event_info,parms) 
  errcode = 0; errmsg = [];

  behav = []; 
  behav.('SubjID') = []; behav.('VisitID') = []; 
  behav.version = mmil_getfield(event_info(1),'version','UNKNOWN');
  if parms.verbose
    fprintf('%s: experiment version = %s\n',mfilename,behav.version);
  end
  behav.switch_flag = parms.switch_flag;
  behav.perform_flag = 1;

  blocks = [event_info.block]; 
  trials = [event_info.trial]; 
  % determine version of task
  n_blocks = size(unique(blocks),2);
  n_trials = size(unique(trials),2);
  if n_trials > 2, 
    blocks = [event_info.block]; 
  else 
    blocks = [event_info.trial];
  end
  
  % check if runs are truncated 
  runs = blocks; 
  nruns = length(unique(blocks)); 
  runs_ok = []; 
  for i=1:nruns
    run_len = length(find(runs==i));     
    if run_len == parms.expected_numtrials
      runs_ok = [runs_ok i];  
    else 
      if parms.verbose
        fprintf('%s: WARNING: run %d is short: %d trials found (expected %d)\n',...
          mfilename,i,run_len,parms.expected_numtrials);
      end
    end 
  end
  nruns = length(runs_ok); 
  
  all_types = {event_info.type};
  if nruns == 1
    ind_ok = find(runs==runs_ok);
    all_types = all_types(ind_ok);
    blocks = blocks(ind_ok); 
    event_info = event_info(ind_ok);
    runs = runs(ind_ok);
  elseif nruns == 0
   if parms.verbose, fprintf('%s: ERROR: no valid e-prime runs in %s\n',mfilename,parms.fname); end
   errcode = 1;
   errmsg = 'no valid e-prime runs';
   return;
  end
  behav.nruns = nruns;

  % store raw run assignment on event_info so timing-file writing below
  % can reuse this same value directly
  for k=1:length(event_info)
    event_info(k).run = runs(k);
  end

  % event durations
  go_dur = {event_info.go_dur}; 
  time_go_dur = [event_info.go_dur];
  time_go_dur = time_go_dur(1); 
  stop_dur = {event_info.stop_dur}; 
  go_event_ind = find(~cellfun(@isempty,go_dur));
  stop_event_ind = find(~cellfun(@isempty,stop_dur));
  % other trials, neither go nor stop: STE (when not recoded as Go) and
  %   ZeroSSDResp
  other_event_ind = find(ismember({event_info.type},{'STE','ZeroSSDResp'}));
  % total numbers of trials    
  behav.('Go_total') = size(go_event_ind,2);
  behav.('Stop_total') = size(stop_event_ind,2);
  behav.('Total') = size(go_event_ind,2)+size(stop_event_ind,2)+...
    size(other_event_ind,2);
  
  % number of trials per run
  for j=1:2 %runs 
    % go
    type_counts_run_go = sprintf('Go_total_run_%d',j); 
    behav.(type_counts_run_go) = 'NaN';  
    % stop
    type_counts_run_stop = sprintf('Stop_total_run_%d',j);
    behav.(type_counts_run_stop) = 'NaN';  
    % all
    behav.(sprintf('Total_run_%d',j)) = 'NaN';
  end
  for j=runs_ok %runs 
    % go
    type_counts_run_go = sprintf('Go_total_run_%d',j);
    ind_type_run = [intersect(find(blocks==j),go_event_ind)];   
    behav.(type_counts_run_go) = size(ind_type_run,2); 
    n_run = size(ind_type_run,2);
    % stop
    type_counts_run_stop = sprintf('Stop_total_run_%d',j);
    ind_type_run = [intersect(find(blocks==j),stop_event_ind)];   
    behav.(type_counts_run_stop) = size(ind_type_run,2); 
    n_run = n_run + size(ind_type_run,2);
    % all: go, stop and other
    ind_type_run = intersect(find(blocks==j),other_event_ind);
    behav.(sprintf('Total_run_%d',j)) = n_run + size(ind_type_run,2);
  end
  
  % results for each condition
  for i=1:parms.nconds
    % total counts per type
    type  = parms.typenames{i};         
    type_counts_total = sprintf('%s_counts_total',type);
    ind_type = find(strcmp(type,all_types)); 
    behav.(type_counts_total) = size(ind_type,2); 
    % total counts per type per run  
    for j=1:2 % runs        
      type_counts_run = sprintf('%s_counts_run_%d',type,j);   
      behav.(type_counts_run) = 'NaN'; 
    end 
    for j=runs_ok       
      type_counts_run = sprintf('%s_counts_run_%d',type,j); 
      ind_type_run = [intersect(find(blocks==j),ind_type)];  
      behav.(type_counts_run) = size(ind_type_run,2); 
    end 
    
    % selecting go or stop event 
    cond  = parms.condnames{i};  
    pat = '\w_(?<event>\w+)'; 
    k = regexp(cond,pat,'names');
    eventname = k.event;    
    switch eventname  
      case 'go'
        event_ind = go_event_ind; 
        event_rt = {event_info.('go_rt')}; 
        if strcmp(cond,'correctlate_go')
           event_rt_fix = {event_info.('fix_rt')};
           event_rt = cellfun(@(x) x+time_go_dur, event_rt_fix,'UniformOutput',false);
        end 
        if strcmp(cond,'incorrectlate_go')
           event_rt_fix = {event_info.('fix_rt')};
           event_rt = cellfun(@(x) x+time_go_dur, event_rt_fix,'UniformOutput',false);
        end 
      case 'stop'
        event_ind = stop_event_ind;
        event_rt = {event_info.('stop_rt')};
    end 
    % total percent per type 
    type_percent_total = sprintf('%s_percent_total',type);  
    behav.(type_percent_total) = size(ind_type,2)./size(event_ind,2);  
    % percent per type per run  
    for j=1:2 % runs 
      type_counts_run = sprintf('%s_percent_run_%d',type,j); 
      behav.(type_counts_run) = 'NaN';  
    end
    for j=runs_ok 
      type_counts_run = sprintf('%s_percent_run_%d',type,j); 
      ind_type_run = [intersect(find(blocks==j),ind_type)];
      ind_type_percent_run = [intersect(find(blocks==j),event_ind)]; 
      behav.(type_counts_run) = size(ind_type_run,2)./size(ind_type_percent_run,2);  
    end
    % total rt per type 
    type_rt_total = sprintf('%s_rt_total',type); 
    all_rt = cell2mat(event_rt(ind_type)); 
    behav.(type_rt_total) = mean(all_rt(all_rt~=0));  
    type_rt_std_total = sprintf('%s_rt_std_total',type);
    behav.(type_rt_std_total) = std(all_rt(all_rt~=0)); 
    
    % total rt per type per run 
    for j=1:2 % runs 
      type_rt_runs = sprintf('%s_rt_run_%d',type,j);
      behav.(type_rt_runs) = 'NaN';  
      type_rt_std_runs = sprintf('%s_rt_std_run_%d',type,j);
      behav.(type_rt_std_runs) = 'NaN';
    end    
    for j=runs_ok % runs 
      type_rt_runs = sprintf('%s_rt_run_%d',type,j);
      ind_type_run = [intersect(find(blocks==j),ind_type)];
      all_rt = cell2mat(event_rt(ind_type_run));  
      behav.(type_rt_runs) = mean(all_rt(all_rt~=0)); 
      type_rt_std_runs = sprintf('%s_rt_std_run_%d',type,j);
      behav.(type_rt_std_runs) = std(all_rt(all_rt~=0));
    end  
  end

  % ssd and ssrt
  ssd_stop_inds = cellfun(@(t)(~isempty(strfind(t,'Stop'))), {event_info.type});
  behav.('SSD_mean_total') = mean([event_info(ssd_stop_inds).ssd_dur]);
  behav.('SSRT_mean_total') = behav.CorrectGo_rt_total-behav.SSD_mean_total;

  % STE (stop-too-early) trial count: ssd_dur>0 and ssd_rt>0, computed
  % directly from the raw per-trial fields rather than type=='STE', since
  % ssd_dur/ssd_rt are untouched by revise_sst_types() regardless of
  % ste_recode_flag (only stop_dur/go_*/type are recoded there) -- so
  % this count is correct whether or not STE trials were recoded as Go
  ssd_dur_all = {event_info.ssd_dur};
  ssd_rt_all  = {event_info.ssd_rt};
  ste_ind = find(cellfun(@(d) ~isempty(d) && d>0, ssd_dur_all) & ...
                 cellfun(@(r) ~isempty(r) && r>0, ssd_rt_all));
  behav.('STE_counts_total') = size(ste_ind,2);
  behav.('STE_percent_total') = size(ste_ind,2)./size(stop_event_ind,2);
  for j=1:2 % runs
    behav.(sprintf('STE_counts_run_%d',j)) = 'NaN';
    behav.(sprintf('STE_percent_run_%d',j)) = 'NaN';
  end
  for j=runs_ok
    ind_run = intersect(find(blocks==j),ste_ind);
    stop_run = intersect(find(blocks==j),stop_event_ind);
    behav.(sprintf('STE_counts_run_%d',j)) = size(ind_run,2);
    behav.(sprintf('STE_percent_run_%d',j)) = size(ind_run,2)./size(stop_run,2);
  end

  % ZeroSSDResp trial count: zero-SSD trials with an SSD-window response
  zr_ind = find(strcmp({event_info.type},'ZeroSSDResp'));
  behav.('ZeroSSDResp_counts_total') = size(zr_ind,2);
  for j=1:2 % runs
    behav.(sprintf('ZeroSSDResp_counts_run_%d',j)) = 'NaN';
  end
  for j=runs_ok
    ind_run = intersect(find(blocks==j),zr_ind);
    behav.(sprintf('ZeroSSDResp_counts_run_%d',j)) = size(ind_run,2);
  end

  % rm meaningless fields
  rmfields_names = {'NoRespGo_rt_total','NoRespGo_rt_run_1','NoRespGo_rt_run_2', ...
                    'NoRespGo_rt_std_total','NoRespGo_rt_std_run_1','NoRespGo_rt_std_run_2', ...
                    'CorrectStop_rt_total','CorrectStop_rt_run_1', ...
                    'CorrectStop_rt_std_total','CorrectStop_rt_std_run_1', ...
                    'CorrectStop_rt_run_2',...
                    'CorrectStop_rt_std_run_2'};
                  
  behav = rmfield(behav,rmfields_names); 
  
  %behav flag 
  if behav.CorrectGo_percent_total < 0.6   
    if parms.verbose, fprintf('%s: CorrectGo_percent_total < 60%% \n',mfilename); end
    behav.perform_flag = 0; 
  end 
  if behav.IncorrectGo_percent_total > 0.3
    if parms.verbose, fprintf('%s: IncorrectGo_percent_total > 30%% \n',mfilename); end
    behav.perform_flag = 0; 
  end 
  late_go = behav.CorrectLateGo_percent_total+behav.IncorrectLateGo_percent_total;
  if late_go > 0.3 
    if parms.verbose, fprintf('%s: LateGo > 30%% \n',mfilename); end
    behav.perform_flag = 0; 
  end
  if behav.NoRespGo_percent_total > 0.3
    if parms.verbose, fprintf('%s: NoRespGo_percent_total > 30%% \n',mfilename); end
    behav.perform_flag = 0; 
  end 
  if behav.CorrectStop_percent_total < 0.2 || behav.CorrectStop_percent_total > 0.8
    if parms.verbose, fprintf('%s: CorrectStop_percent_total < 20%% or > 80%% \n',mfilename); end
    behav.perform_flag = 0; 
  end 
  if behav.Go_total < 290
    if parms.verbose, fprintf('%s: Go_total trials less than 290 \n',mfilename); end
    behav.perform_flag = 0; 
  end
  if behav.Stop_total < 40
    if parms.verbose, fprintf('%s: Stop_total trials less than 40 \n',mfilename); end
    behav.perform_flag = 0; 
  end
  if behav.CorrectGo_rt_total < behav.IncorrectStop_rt_total
    if parms.verbose, fprintf('%s: CorrectGo_rt_total < IncorrectStop_rt_total \n',mfilename);end
    behav.perform_flag = 0;
  end;
  if parms.verbose, fprintf('%s: perform_flag = %.0f\n',mfilename,behav.perform_flag); end

  %% check for problematic conditions
  
  behav.SSRT_integrated_total = NaN;
  behav.violatorflag = NaN;
  behav.glitchflag = NaN;
  behav.zeroSSDcount = NaN;
  
  tfmri_sst_all_beh_total_issrt = NaN;
  tfmri_sst_beh_glitchflag = NaN;
  tfmri_sst_beh_glitchcnt = NaN;
  tfmri_sst_beh_0SSDcount = NaN;
  tfmri_sst_beh_0SSD_flag = NaN;  % called "tfmri_sst_beh_0SSD>20flag" in Sage's script
  tfmri_sst_beh_violatorflag = NaN;
  
  % set_correct_go
  manynans = num2cell(NaN(1,length(event_info)));
  [event_info.('correct_go_response')] = manynans{:};

  inds = cellfun(@(x,y)(~isempty(x) & (x==y)),...
                  {event_info.go_resp}, {event_info.go_cresp}, 'UniformOutput',0);
  inds = find(cellfun(@(x)(isequal(x,1)), inds));
  [event_info(inds).correct_go_response] = deal(1.0);
 
  inds = cellfun(@(x,y,z)( isempty(x) & (y==z)),...
                  {event_info.go_resp}, {event_info.go_cresp}, {event_info.fix_resp},'UniformOutput',0);
  inds = find(cellfun(@(x)(isequal(x,1)), inds));
  [event_info(inds).correct_go_response] = deal(1.0);

  inds = cellfun(@(x,y,z)(~isempty(x) & ~isequal(x,y) & ~isempty(strfind(z,'Go'))),...
                  {event_info.go_resp}, {event_info.go_cresp}, {event_info.type}, 'UniformOutput',0);
  inds = find(cellfun(@(x)(isequal(x,1)), inds));
  [event_info(inds).correct_go_response] = deal(0.0);
  
  inds = cellfun(@(x,y,z,r)( isempty(x) & ~isequal(y,z) & ~isempty(strfind(r,'Go'))),...
                  {event_info.go_resp}, {event_info.go_cresp}, {event_info.fix_resp}, {event_info.type},...
                  'UniformOutput',0);
  inds = find(cellfun(@(x)(isequal(x,1)), inds));
  [event_info(inds).correct_go_response] = deal(0.0);

  inds = cellfun(@(x,y,z)( isempty(x) & isempty(y) & ~isempty(strfind(z,'Go'))),...
                  {event_info.go_resp}, {event_info.fix_resp}, {event_info.type}, 'UniformOutput',0);
  inds = find(cellfun(@(x)(isequal(x,1)), inds));
  [event_info(inds).correct_go_response] = deal('omission');

  % set_correct_stop
  manynans = num2cell(NaN(1,length(event_info)));
  [event_info.('correct_stop')] = manynans{:};

  inds = cellfun(@(x,y,z,r)( isempty(x) & isempty(y) & isempty(z) & ~isempty(strfind(r,'Stop')) ),...
            {event_info.stop_resp}, {event_info.fix_resp}, {event_info.ssd_resp}, {event_info.type}, 'UniformOutput',0);
  inds = find(cellfun(@(x)(isequal(x,1)), inds));
  [event_info(inds).correct_stop] = deal(1.0);

  inds = cellfun(@(x,y,z,r)((~isempty(x) | ~isempty(y) | ~isempty(z)) & ~isempty(strfind(r,'Stop')) ),...
            {event_info.stop_resp}, {event_info.fix_resp}, {event_info.ssd_resp}, {event_info.type}, 'UniformOutput',0);
  inds = find(cellfun(@(x)(isequal(x,1)), inds));
  [event_info(inds).correct_stop] = deal(0.0);

  % set_correct_go_rt
  inds = cellfun(@(x,y,z)(~isempty(x) & isempty(y) & ~isempty(strfind(z,'Go')) ),...
                  {event_info.fix_resp}, {event_info.go_resp}, {event_info.type}, 'UniformOutput',0 );
  inds = find(cellfun(@(x)(isequal(x,1)), inds));
  [event_info.go_rt_adjusted] = event_info.go_rt;
  for i = 1:length(inds)
    j = inds(i);
    event_info(j).go_rt_adjusted = event_info(j).go_dur + event_info(j).fix_rt;
  end

  % set correct_stop_rt
  [event_info.('stop_rt_adjusted')] = event_info.stop_rt;

  % set to stop signal duration + fix rt, as was answered during fix
  inds = cellfun(@(x,y,z,r)(~isempty(x) & isempty(y) & ~isempty(strfind(z,'Stop')) & isequal(r,0) ),...
            {event_info.fix_resp}, {event_info.stop_resp}, {event_info.type}, {event_info.correct_stop}, 'UniformOutput',0);
  inds = find(cellfun(@(x)(isequal(x,1)), inds));
  for i = 1:length(inds)
    j = inds(i);
    event_info(j).stop_rt_adjusted = event_info(j).stop_dur + event_info(j).fix_rt;
  end

  % set updated for response during stop
  inds = cellfun(@(x,y,z)((~isempty(x) | ~isempty(y)) & ~isempty(z)),...
            {event_info.stop_resp}, {event_info.fix_resp}, {event_info.stop_rt_adjusted}, 'UniformOutput',0);
  inds = find(cellfun(@(x)(isequal(x,1)), inds));
  for i = 1:length(inds)
    j = inds(i);
    event_info(j).stop_rt_adjusted = event_info(j).stop_rt_adjusted + event_info(j).ssd_dur;
  end

  % null go_rt_adjusted/stop_rt_adjusted to NaN for true no-response
  % trials (NoRespGo, CorrectStop) -- neither is touched by any override
  % above, so they'd otherwise stay at their raw (0) value
  inds = cellfun(@(x,y,z)( isempty(x) & isempty(y) & ~isempty(strfind(z,'Go')) ),...
            {event_info.go_resp}, {event_info.fix_resp}, {event_info.type}, 'UniformOutput',0);
  inds = find(cellfun(@(x)(isequal(x,1)), inds));
  [event_info(inds).go_rt_adjusted] = deal(NaN);

  inds = cellfun(@(x,y,z,r)( isempty(x) & isempty(y) & isempty(z) & ~isempty(strfind(r,'Stop')) ),...
            {event_info.stop_resp}, {event_info.fix_resp}, {event_info.ssd_resp}, {event_info.type}, 'UniformOutput',0);
  inds = find(cellfun(@(x)(isequal(x,1)), inds));
  [event_info(inds).stop_rt_adjusted] = deal(NaN);

  % check_omissions
  inds = cellfun(@(x,y,z)( isempty(x) & isempty(y) & ~isempty(strfind(z,'Go')) ),...
            {event_info.go_resp}, {event_info.fix_resp}, {event_info.type}, 'UniformOutput',0);
  inds = find(cellfun(@(x)(isequal(x,1)), inds));
  n_missing = length(inds);

  inds = cellfun(@(z)( ~isempty(strfind(z,'Go')) ), {event_info.type}, 'UniformOutput',0);
  inds = find(cellfun(@(x)(isequal(x,1)), inds));
  n_go_trials = length(inds);

  if n_missing == n_go_trials
    if parms.verbose
      fprintf('%s: WARNING: all trials in this file are ommissions\n',mfilename);
    end
  end

  % separate just go trials
  inds = cellfun( @(x)( ~isempty(strfind(x,'Go')) ), {event_info.type}, 'UniformOutput',0 );
  inds = find(cellfun(@(x)(isequal(x,1)), inds));
  go_trials = event_info(inds);

  % set omissions to max go.rt - if primary rt is_null means omission
  inds = cellfun( @(x)( strcmp(x,'omission') ), {go_trials.correct_go_response}, 'UniformOutput',0 );
  inds = find(cellfun(@(x)(isequal(x,1)), inds));
  go_rt_max = max( [go_trials.go_rt_adjusted] );
  for i = 1:length(inds)
    j = inds(i);
    go_trials(j).go_rt_adjusted = go_rt_max;
  end
  % sort go trials
  sorted_go = sort([go_trials.go_rt_adjusted]);

  % get stop trials
  inds = cellfun( @(x)(~isempty(strfind(x,'Stop'))), {event_info.type}, 'UniformOutput',0 );
  inds = find(cellfun(@(x)(isequal(x,1)), inds));
  stop_trials = event_info(inds);

  % calculate prob stop failure
  prob_stop_failure = 1 - ...
    sum( cellfun(@(x)(isequal(x,1)),{stop_trials.correct_stop})) / length( stop_trials );

  % calc go.rt index
  index = prob_stop_failure * length(sorted_go);

  % if prob_stop_failure is 1, use the max go.rt
  if ceil(index) == length(sorted_go)
    index = length(sorted_go);
  else
    index = [floor(index) ceil(index)]+1;
  end

  % calc SSRT
  mean_ssd = mean( [stop_trials.ssd_dur] );
  try
    ssrt = mean(sorted_go(index)) - mean_ssd;
  catch
      fprintf('index = %f\n', index);
      fprintf('len sorted go = %f\n', length(sorted_go));
      fprintf('prob stop fail = %f\n', prob_stop_failure);
      ssrt = NaN;
  end
  tfmri_sst_all_beh_total_issrt = ssrt;
  if parms.verbose, fprintf('%s: tfmri_sst_all_beh_total_issrt = %.4f\n',...
      mfilename,tfmri_sst_all_beh_total_issrt); end

  % number of glitched trials (stop and ZeroSSDResp trials, excluding STE
  %   go/nuisance rows)
  glitch_trials = event_info(cellfun(@(x)(~isempty(strfind(x,'Stop')) | ...
    strcmp(x,'ZeroSSDResp')), {event_info.type}));
  inds = arrayfun(@(x,y)( (x < 50) & (x > 0) & (y <= 50)), [glitch_trials.ssd_rt], [glitch_trials.ssd_dur]);
  tfmri_sst_beh_glitchcnt = sum(inds);
  if parms.verbose, fprintf('%s: tfmri_sst_beh_glitchcnt = %.0f\n',...
      mfilename,tfmri_sst_beh_glitchcnt); end
  
  tfmri_sst_beh_glitchflag = (tfmri_sst_beh_glitchcnt > 0);
  if parms.verbose, fprintf('%s: tfmri_sst_beh_glitchflag = %.0f\n',...
      mfilename,tfmri_sst_beh_glitchflag); end

  %  def get_0SSD_cnt(subj)  and   get_0SSD_flag(subj):
  inds = arrayfun(@(x)( x == 0 ), [event_info.ssd_dur]);
  tfmri_sst_beh_0SSDcount = sum(inds);
  if parms.verbose, fprintf('%s: tfmri_sst_beh_0SSDcount = %.0f\n',...
      mfilename,tfmri_sst_beh_0SSDcount); end
  
  tfmri_sst_beh_0SSD_flag = (tfmri_sst_beh_0SSDcount > 20);  % Caled tfmri_sst_beh_0SSD>20flag in Sage's script
  if parms.verbose, fprintf('%s: tfmri_sst_beh_0SSD_flag = %.0f\n',...
      mfilename,tfmri_sst_beh_0SSD_flag); end

  inds = cellfun(@(x,y)( (isequal(x,0) | isequal(x,1))  & ~isempty(strfind(y,'Go')) ),...
            {event_info.correct_go_response}, {event_info.type}, 'UniformOutput',0);
  inds = find(cellfun(@(x)(isequal(x,1)), inds));
  go_rt = mean( [event_info(inds).go_rt_adjusted] );

  
  inds = cellfun(@(x)( isequal(x,0) ), {event_info.correct_stop}, 'UniformOutput',0);
  inds = find(cellfun(@(x)(isequal(x,1)), inds));
  stop_rt = mean( [event_info(inds).stop_rt_adjusted] );
  
  if parms.verbose
    fprintf('%s:   go_rt = %.4f\n', mfilename,go_rt);
    fprintf('%s: stop_rt = %.4f\n', mfilename,stop_rt);
  end

  tfmri_sst_beh_violatorflag = (stop_rt > go_rt);
  if parms.verbose, fprintf('%s: tfmri_sst_beh_violatorflag = %.0f\n',...
      mfilename,tfmri_sst_beh_violatorflag); end

  behav.SSRT_integrated_total = tfmri_sst_all_beh_total_issrt;
  behav.violatorflag = tfmri_sst_beh_violatorflag;
  behav.glitchflag = tfmri_sst_beh_glitchflag;
  behav.zeroSSDcount  = tfmri_sst_beh_0SSDcount;

  % write to csv
  fname_csv_out = sprintf('%s/%s_behavioral.csv',parms.outdir,parms.outstem); 
  if ~exist(fname_csv_out,'file') || parms.forceflag
    mmil_struct2csv(behav,fname_csv_out)
  end;
return; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function get_behavioral_data_empty(parms) 

  behav = []; 
  behav.('SubjID') = []; behav.('VisitID') = []; 
  behav.version = [];
  behav.switch_flag = 0; 
  behav.perform_flag = 0;
  behav.nruns = 0; 
  
  % total numbers of trials    
  behav.('Go_total') = 'NaN'; 
  behav.('Stop_total') = 'NaN'; 
  behav.('Total') = 'NaN'; 
  % number of trials per run    
  for j=1:2 %runs 
    % go
    type_counts_run_go = sprintf('Go_total_run_%d',j);
    behav.(type_counts_run_go) = 'NaN'; 
    % stop
    type_counts_run_stop = sprintf('Stop_total_run_%d',j);  
    behav.(type_counts_run_stop) = 'NaN';  
    % all
    behav.(sprintf('Total_run_%d',j)) = 'NaN';
  end
  
  % results for each condition
  for i=1:parms.nconds
    % total counts per type
    type  = parms.typenames{i};         
    type_counts_total = sprintf('%s_counts_total',type);  
    behav.(type_counts_total) = 'NaN'; 
    % total counts per type per run  
    for j=1:2 % runs        
      type_counts_run = sprintf('%s_counts_run_%d',type,j);  
      behav.(type_counts_run) = 'NaN'; 
    end 

    % total percent per type 
    type_percent_total = sprintf('%s_percent_total',type);  
    behav.(type_percent_total) = 'NaN'; 
    % percent per type per run  
    for j=1:2 % runs 
      type_counts_run = sprintf('%s_percent_run_%d',type,j); 
      behav.(type_counts_run) = 'NaN'; 
    end
    % total rt per type 
    type_rt_total = sprintf('%s_rt_total',type); 
    behav.(type_rt_total) = 'NaN'; 
    type_rt_std_total = sprintf('%s_rt_std_total',type);
    behav.(type_rt_std_total) = 'NaN'; 
    
    % total rt per type per run 
    for j=1:2 % runs 
      type_rt_runs = sprintf('%s_rt_run_%d',type,j);  
      behav.(type_rt_runs) = 'NaN'; 
      type_rt_std_runs = sprintf('%s_rt_std_run_%d',type,j);
      behav.(type_rt_std_runs) = 'NaN'; 
    end     
  end 
  % ssd and ssrt 
  behav.('SSD_mean_total')= 'NaN'; 
  behav.('SSRT_mean_total') = 'NaN';  
  % rm meaningless fields 
  rmfields_names = {'NoRespGo_rt_total','NoRespGo_rt_run_1','NoRespGo_rt_run_2', ...
                    'NoRespGo_rt_std_total','NoRespGo_rt_std_run_1','NoRespGo_rt_std_run_2', ...
                    'CorrectStop_rt_total','CorrectStop_rt_run_1', ...
                    'CorrectStop_rt_std_total','CorrectStop_rt_std_run_1', ...
                    'CorrectStop_rt_run_2',...
                    'CorrectStop_rt_std_run_2'};
                  
  behav = rmfield(behav,rmfields_names); 
  % write to csv
  fname_csv_out = sprintf('%s/%s_behavioral.csv',parms.outdir,parms.outstem); 
  if ~exist(fname_csv_out,'file') || parms.forceflag
    mmil_struct2csv(behav,fname_csv_out)
  end;
  
return; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
