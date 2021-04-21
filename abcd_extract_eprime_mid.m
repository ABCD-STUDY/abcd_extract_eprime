function [eprime_nruns,errcode,behav,errmsg] = abcd_extract_eprime_mid(fname,varargin)
%function [eprime_nruns,errcode,behav,errmsg] = abcd_extract_eprime_mid(fname,[options])
%
% Purpose: extract condition time courses
%   from eprime data files for MID task
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
%   'nskipTRs': number of TRs to remove from beginning of time course
%     {default = 0}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%   'timing_files_flag': [0|1] generate timing files
%     {default = 1}
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
% Prev Mod: 01/23/19 by Dani Cornejo
% Prev Mod: 08/05/19 by Octavio Ruiz
% Prev Mod: 11/03/20 by Don Hagler
% Prev Mod: 03/31/21 by Emma Pearson, Anthony Juliano, and Bader Chaarani at UVM
% Prev Mod: 04/06/21 by Don Hagler
% Last Mod: 04/16/21 by Don Hagler

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE: based on abcd_mid_extract.m
%       provided by Mary Soules (mfield@med.umich.edu) 10/03/16

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize outputs
eprime_nruns = []; errcode = 0; behav = []; errmsg = [];

% check arguments 
if ~mmil_check_nargs(nargin,1), return; end;

% check input parameters
parms = check_input(fname,varargin); 

% create output directory
mmil_mkdir(parms.outdir);

% create struct with info for each event
[event_info,start_time,all_types,errcode,errmsg] = get_event_info(parms);  
if errcode, return; end;
        
% get behavioral data and write it to a csv
%   also we will return this table to the calling program 
[behav,eprime_runs,errcode,errmsg] = get_behav_feedback(event_info,start_time,parms,all_types);
if errcode, return; end;

parms.eprime_runs = eprime_runs; 
eprime_nruns = length(eprime_runs);
parms.eprime_nruns = eprime_nruns;
if parms.verbose, fprintf('%s: number of e-prime runs = %d \n',mfilename,eprime_nruns); end

if parms.timing_files_flag && eprime_nruns  
  % separate file for each combination of condnames and stimnames (events)
  for i=1:parms.nconds
    type = parms.typenames{i};
    cond = parms.condnames{i};
    ind_type = find(strcmp(type,all_types));
    for j=1:parms.nstims
      stim = parms.stimnames{j};
      eventname = sprintf('%s_%s',cond,stim); 
      % get times of stim onset and offset
      switch stim
        case 'pos_feedback'
          % find events with positive feedback (acc = 1)
          acc = [event_info(ind_type).acc];
          ind_keep = find(acc==1);
          onset = {event_info(ind_type(ind_keep)).feedback_onset};
          offset = {event_info(ind_type(ind_keep)).feedback_offset};
          [onset,offset] = check_offsets(onset,offset,eventname,parms);
          % find most recent start time for each event
          [ind_start,event_ref_time] = set_ref(onset,start_time);
          % create files for each scan
          write_files(eventname,onset,offset,ind_start,event_ref_time,parms)
        case 'neg_feedback'
          % find events with negative feedback (acc = 0)
          acc = [event_info(ind_type).acc];
          ind_keep = find(acc==0);
          onset = {event_info(ind_type(ind_keep)).feedback_onset};
          offset = {event_info(ind_type(ind_keep)).feedback_offset};
          [onset,offset] = check_offsets(onset,offset,eventname,parms);
          % find most recent start time for each event
          [ind_start,event_ref_time] = set_ref(onset,start_time);
          % create files for each scan
          write_files(eventname,onset,offset,ind_start,event_ref_time,parms)    
        otherwise
          ind_keep = [1:length(ind_type)];
          onset = {event_info(ind_type(ind_keep)).cue_onset}; 
          offset = {event_info(ind_type(ind_keep)).probe_offset}; 
          [onset,offset] = check_offsets(onset,offset,eventname,parms);
          % find most recent start time for each event
          [ind_start,event_ref_time] = set_ref(onset,start_time); 
          % create files for each scan
          write_files(eventname,onset,offset,ind_start,event_ref_time,parms)    
      end 
    end
  end
end

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
    'nskipTRs',0,[],...
    'forceflag',false,[false true],...
    'timing_files_flag',true,[false true],...
    'verbose',true,[false true],...
    'eprime_runs',1:2,[],...
    'eprime_nruns',0,0:2,...
    'min_run_resps',20,[],...
    ...
    'colnames',  {'NARGUID','SessionDate','SessionTime','ExperimentName','MIDVERSION',...
                  'Block','SubTrial','Condition','PrepTime.OnsetTime','PrepTime.OffsetTime',...
                  'Cue.OnsetTime','Anticipation.OnsetTime','Probe.OnsetTime','Feedback.OnsetTime',...
                  'Cue.OffsetTime','Anticipation.OffsetTime','Probe.OffsetTime','Feedback.OffsetTime','ResponseCheck',...
                  'TextDisplay1.RTTime','Feedback.RTTime','Probe.RESP','TextDisplay1.RESP','Feedback.RESP',...
                  'prbacc','prbrt','RunMoney',... 
                  'Anticipation.Duration','Cue.Duration','Probe.Duration','Probe.OffsetTime'},[],...
    'fieldnames',{'narguid','date','time','experiment','version',...
                  'run','trial','type','prep_onset','prep_offset',...
                  'cue_onset','antic_onset','probe_onset','feedback_onset',...
                  'cue_offset','antic_offset','probe_offset','feedback_offset','response_check',...
                  'text_display_rttime','feedback_rttime','probe_resp','text_display_resp','feedback_resp',...
                  'acc','rt','money',...
                  'antic_dur','cue_dur','probe_dur','probe_offset'},[],...
    'typenames',{'SmallReward','LgReward','SmallPun','LgPun','Triangle'},[],...
    'condnames',{'small_reward','large_reward','small_loss','large_loss','neutral'},[],...
    'stimnames',{'antic','pos_feedback','neg_feedback'},[],...
  });

  parms.nconds = length(parms.condnames);
  parms.ntypes = length(parms.typenames);
  parms.nstims = length(parms.stimnames);
  parms.mindur = parms.minfrac*parms.TR;
  if parms.nconds ~= parms.ntypes
    error('condnames and typenames length mismatch');
  end;
  if ~exist(parms.fname,'file')
    error('file %s not found',parms.fname);
  end;
  [fdir,fstem,fext] = fileparts(parms.fname);
  % remove problematic characters
  if isempty(parms.outstem)
    parms.outstem = abcd_clean_fstem(fstem);
  end;
  % calculate time at the end of each TR
  parms.TR_offset = linspace(parms.TR,parms.numTRs * parms.TR,parms.numTRs);
  parms.TR_onset = parms.TR_offset - parms.TR;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ind_start,event_ref_time] = set_ref(onset,start_time)
  if ~isempty(onset)
    d = bsxfun(@minus,onset,start_time');
    d(d<0) = Inf;
    [d,ind_start] = min(d,[],1);
    event_ref_time  = start_time(ind_start);
  else
    ind_start = start_time;
    event_ref_time = 0;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [onset,offset] = check_offsets(onset,offset,eventname,parms)
  ind_empty = find(cellfun(@isempty,onset) | cellfun(@isempty,offset));
  if ~isempty(ind_empty)
    if parms.verbose
      fprintf('%s: WARNING: %s event has %d onsets but %d offsets\n',...
        mfilename,eventname,nnz(~cellfun(@isempty,onset)),nnz(~cellfun(@isempty,offset)));
    end
    ind_keep = setdiff([1:length(onset)],ind_empty);
    onset = onset(ind_keep);
    offset = offset(ind_keep);
  end;
  onset = cell2mat(onset);
  offset = cell2mat(offset);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function write_files(eventname,onset,offset,ind_start,event_ref_time,parms)
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
    if ~exist(fname_1D,'file') || ~exist(fname_txt,'file') || ...
       ~exist(fname_block,'file') || ~exist(fname_fsl,'file') || parms.forceflag
      if ~isempty(onset(ind_scan))
        % subtract start time and convert to seconds
        rel_onset  = (onset(ind_scan) - event_ref_time(ind_scan))/1000;
        rel_offset = (offset(ind_scan) - event_ref_time(ind_scan))/1000;
        % loop over events
        TR_mask = zeros(parms.numTRs,1);
        for k=1:length(rel_onset)
          ind_rep = find(rel_onset(k)<(parms.TR_offset-parms.mindur) &...
                         rel_offset(k)>(parms.TR_onset+parms.mindur));
          TR_mask(ind_rep) = 1;
        end;
        % calculate event times from TR_onset
        event_times = parms.TR_onset(find(TR_mask));
        % convert time units to multiples of TRs
        event_TRs = event_times/parms.TR;
        % remove  nskipTRs
        if parms.nskipTRs>0
          TR_mask = TR_mask(parms.nskipTRs+1:end);
          event_times = event_times(parms.nskipTRs+1:end) - parms.nskipTRs*parms.TR;
          event_TRs = event_times/parms.TR; 
        end;
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
        block_len = mean(rel_offset - rel_onset); 
        fprintf(fid,'%s\n',sprintf('%0.2f ',block_len));
        fclose(fid);
        % write to fsl file
        fid = fopen(fname_fsl,'wt');
        if fid<0, error('failed to open %s for writing',fname_fsl); end;
        for i=1:length(rel_onset)
          fprintf(fid, '%0.2f %0.2f 1 \n',rel_onset(i), rel_offset(i)-rel_onset(i));
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
      end
    end;
  end;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [event_info,start_time,all_types,errcode,errmsg] = get_event_info(parms)
  
  event_info=[]; start_time=[]; all_types=[];
  errcode = 0; errmsg = [];

  try 
    % write event info to file
    fname_csv = abcd_check_eprime_sprdsh(parms.fname, parms.colnames,...
                  parms.fieldnames, parms.outdir, parms.forceflag, parms.verbose);
    event_info = mmil_csv2struct(fname_csv);
  catch me
    if parms.verbose, fprintf('%s: ERROR: failed to read e-prime file %s:\n%s\n',...
      mfilename,parms.fname,me.message); end
    errcode = 1;
    errmsg = 'failed to read e-prime file';
    return; 
  end
  
  % check experiment
  experiment = mmil_getfield(event_info(1),'experiment');
  if isempty(regexpi(experiment,'mid'))
    if parms.verbose, fprintf('%s: ERROR: wrong experiment name in e-prime file %s: %s\n',...
      mfilename,parms.fname,experiment); end
    errcode = 1;
    errmsg = 'wrong experiment name';
    return;
  else
    if parms.verbose, fprintf('%s: experiment name: %s\n',mfilename,experiment); end
  end;

  % get start times
  ind_start = find(~cellfun(@isempty,{event_info.prep_onset}));
  start_time = [event_info(ind_start).prep_onset];
  % remove non-events
  all_types = {event_info.type};
  ind_events = find(~cellfun(@isempty,all_types));
  event_info = event_info(ind_events);
  all_types = {event_info.type}; 
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  function get_behav_empty(parms) 
%  
%    behav = [];
%    behav.('SubjID') = []; behav.('VisitID') = []; 
%    behav.switch_flag = 0; 
%    behav.perform_flag = 0;
%    behav.nruns = 0; 
%    
%    counts_total = sprintf('total_numtrials'); 
%    behav.(counts_total) = NaN;
%    for j=1:2 % runs
%      counts_total_run = sprintf('total_numtrials_run_%d',j);
%      behav.(counts_total_run) = NaN; 
%    end 
%    
%    for i=1:parms.nconds
%      type = parms.typenames{i}; 
%      cond = parms.condnames{i};
%  
%      type_counts_total = sprintf('%s_total_bothruns_numtrials',cond); 
%      behav.(type_counts_total) = NaN; 
%      type_counts_total = sprintf('%s_total_bothruns_mean_rt',cond); 
%      behav.(type_counts_total) = NaN;
%      type_counts_total = sprintf('%s_total_bothruns_std_rt',cond); 
%      behav.(type_counts_total) = NaN;
%      
%      for j=1:2 % runs
%        type_counts_run = sprintf('%s_total_run_%d_numtrials',cond,j);   
%        behav.(type_counts_run) = NaN;
%        type_counts_run = sprintf('%s_total_run_%d_mean_rt',cond,j);
%        behav.(type_counts_run) = NaN;
%        type_counts_run = sprintf('%s_total_run_%d_std_rt',cond,j);
%        behav.(type_counts_run) = NaN;  
%      end 
%  
%      for j=1:parms.nstims
%        stim = parms.stimnames{j}; 
%        eventname = sprintf('%s_%s',cond,stim);
%  
%        switch stim
%          case 'pos_feedback'
%            
%            type_counts_stim = sprintf('%s_bothruns_numtrials',eventname);  
%            behav.(type_counts_stim) = NaN;  
%            type_acc_stim = sprintf('%s_bothruns_rate',eventname);
%            behav.(type_acc_stim) = NaN; 
%            type_rt_stim = sprintf('%s_bothruns_mean_rt',eventname);
%            behav.(type_rt_stim) = NaN; 
%            type_rt_stim = sprintf('%s_bothruns_std_rt',eventname);
%            behav.(type_rt_stim) = NaN; 
%  
%            for j=1:2 % runs
%              type_counts_run = sprintf('%s_run_%d_numtrials',eventname,j);  
%              behav.(type_counts_run) = NaN; 
%              type_acc_run = sprintf('%s_run_%d_rate',eventname,j);
%              behav.(type_acc_run) = NaN; 
%              type_rt_run = sprintf('%s_run_%d_mean_rt',eventname,j);    
%              behav.(type_rt_run) = NaN; 
%              type_rt_run = sprintf('%s_run_%d_std_rt',eventname,j);
%              behav.(type_rt_run) = NaN; 
%            end
%  
%          case 'neg_feedback'
%            
%            type_counts_stim = sprintf('%s_bothruns_numtrials',eventname);         
%            behav.(type_counts_stim) = NaN; 
%            type_acc_stim = sprintf('%s_bothruns_rate',eventname);
%            behav.(type_acc_stim) = NaN; 
%            type_rt_stim = sprintf('%s_bothruns_mean_rt',eventname);
%            behav.(type_rt_stim) = NaN;  
%            type_rt_stim = sprintf('%s_bothruns_std_rt',eventname);
%            behav.(type_rt_stim) = NaN; 
%          
%            for j=1:2 % runs
%              type_counts_run = sprintf('%s_run_%d_numtrials',eventname,j);   
%              behav.(type_counts_run) = NaN; 
%              type_acc_run = sprintf('%s_run_%d_rate',eventname,j);
%              behav.(type_acc_run) = NaN;        
%              type_rt_run = sprintf('%s_run_%d_mean_rt',eventname,j);   
%              behav.(type_rt_run) = NaN; 
%              type_rt_run = sprintf('%s_run_%d_std_rt',eventname,j);
%              behav.(type_rt_run) = NaN;         
%            end
%  
%        end
%      end 
%    end
%    
%    % collapse across small or large 
%    collapse_cond = {'reward', 'loss'}; 
%    
%    for i=1:length(collapse_cond)
%      cond = collapse_cond{i}; 
%      type_counts_total = sprintf('hilo%s_total_bothruns_numtrials',cond); 
%      behav.(type_counts_total) = NaN; 
%      type_rt_total = sprintf('hilo%s_total_bothruns_mean_rt',cond); 
%      behav.(type_rt_total) = NaN; 
%      type_rt_total = sprintf('hilo%s_total_bothruns_std_rt',cond); 
%      behav.(type_rt_total) = NaN; 
%      
%      for j=1:2 % runs
%       type_counts_run = sprintf('hilo%s_total_run_%d_numtrials',cond,j);    
%       behav.(type_counts_run) = NaN;
%       type_rt_run = sprintf('hilo%s_total_run_%d_mean_rt',cond,j);
%       behav.(type_rt_run) = NaN;
%       type_rt_run = sprintf('hilo%s_total_run_%d_std_rt',cond,j);
%       behav.(type_rt_run) = NaN;
%      end 
%     
%      for j=1:parms.nstims
%        stim = parms.stimnames{j}; 
%        eventname = sprintf('%s_%s',cond,stim);   
%      
%        switch stim
%          case 'pos_feedback'
%      
%            type_counts_stim = sprintf('hilo%s_bothruns_numtrials',eventname);
%            behav.(type_counts_stim) = NaN; 
%            type_acc_stim = sprintf('hilo%s_bothruns_rate',eventname);
%            behav.(type_acc_stim) = NaN; 
%            type_rt_stim = sprintf('hilo%s_bothruns_mean_rt',eventname);
%            behav.(type_rt_stim) = NaN; 
%            type_rt_stim = sprintf('hilo%s_bothruns_std_rt',eventname);
%            behav.(type_rt_stim) = NaN;  
%            
%            for j=1:2 % runs
%              type_counts_run = sprintf('hilo%s_run_%d_numtrials',eventname,j);   
%              behav.(type_counts_run) = NaN; 
%              type_acc_run = sprintf('hilo%s_run_%d_rate',eventname,j);
%              behav.(type_acc_run) = NaN; 
%              type_rt_run = sprintf('hilo%s_run_%d_mean_rt',eventname,j);  
%              behav.(type_rt_run) = NaN;   
%              type_rt_run = sprintf('hilo%s_run_%d_std_rt',eventname,j);
%              behav.(type_rt_run) = NaN;
%            end 
%      
%          case 'neg_feedback'
%      
%            type_counts_stim = sprintf('hilo%s_bothruns_numtrials',eventname);
%            behav.(type_counts_stim) = NaN; 
%            type_acc_stim = sprintf('hilo%s_bothruns_rate',eventname);
%            behav.(type_acc_stim) = NaN;  
%            type_rt_stim = sprintf('hilo%s_bothruns_mean_rt',eventname);  
%            behav.(type_rt_stim) = NaN; 
%            type_rt_stim = sprintf('hilo%s_bothruns_std_rt',eventname);
%            behav.(type_rt_stim) = NaN; 
%            
%            for j=1:2 % runs
%              type_counts_run = sprintf('hilo%s_run_%d_numtrials',eventname,j);     
%              behav.(type_counts_run) = NaN; 
%              type_acc_run = sprintf('hilo%s_run_%d_rate',eventname,j);
%              behav.(type_acc_run) = NaN;        
%              type_rt_run = sprintf('hilo%s_run_%d_mean_rt',eventname,j);  
%              behav.(type_rt_run) = NaN;    
%              type_rt_run = sprintf('hilo%s_run_%d_std_rt',eventname,j);
%              behav.(type_rt_run) = NaN;
%            end
%            
%        end
%      end
%    end
%    
%    % money 
%    behav.('earnings_total') = NaN; 
%    for j=1:2 % runs
%      type_counts_run = sprintf('earnings_run_%d',j);  
%      behav.(type_counts_run) = NaN;   
%    end 
%    
%    % write to csv
%    fname_csv_out = sprintf('%s/%s_behavioral.csv',parms.outdir,parms.outstem); 
%    if ~exist(fname_csv_out,'file') || parms.forceflag
%      mmil_struct2csv(behav,fname_csv_out);
%    end;
%  return; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [behav,runs_ok,errcode,errmsg] = get_behav_feedback(event_info,start_time,parms,all_types)
  errcode = 0; errmsg = [];
  
  behav = [];
  behav.('SubjID') = []; behav.('VisitID') = []; 
  behav.switch_flag = 0; 
  behav.feedback_flag = 1;
  behav.perform_flag = 1;

  % check for variable in 2017 version
  fbflag = isfield(event_info, 'feedback_resp');
  
  % check if second runs is truncated 
  runs = [event_info.run];  
  nruns = length(unique(runs)); 

  runs_ok = []; 
  for i=1:nruns
    run_len = length(find(runs==i));   
    if run_len == 50 % require expected run length
     runs_ok = [runs_ok i];  
    else 
     if parms.verbose, fprintf('%s: run %d is short: %d trials found (<50 trials) \n',mfilename,i,run_len); end
    end 
  end
  nruns = length(runs_ok); 

  if nruns == 1
    new_info = find(runs==runs_ok);
    all_types = all_types(new_info);  
    event_info = event_info(new_info);
    runs = runs(new_info); 
  elseif nruns == 0
    if parms.verbose, fprintf('%s: ERROR: no valid e-prime runs in %s\n',mfilename,parms.fname); end
    errcode = 1;
    errmsg = 'no valid e-prime runs';
    return;
  end  
  behav.nruns = nruns; 
  
  counts_total = sprintf('total_numtrials'); 
  behav.(counts_total) = length(all_types);
  for j=1:2 % runs
    counts_total_run = sprintf('total_numtrials_run_%d',j);
    behav.(counts_total_run) = NaN; 
  end 
  for j=runs_ok % runs
    ind_type_run = find(runs==j); 
    counts_total_run = sprintf('total_numtrials_run_%d',j);
    behav.(counts_total_run) = length(ind_type_run); 
  end 
  
  for i=1:parms.nconds
    type = parms.typenames{i}; 
    cond = parms.condnames{i};
    ind_type = find(strcmp(type,all_types)); 

    type_counts_total = sprintf('%s_total_bothruns_numtrials',cond); 
    behav.(type_counts_total) = length(ind_type); 
    type_counts_total = sprintf('%s_total_bothruns_mean_rt',cond);
    if fbflag
      rt_keep = get_rt(event_info(ind_type));
      behav.(type_counts_total) = mean(rt_keep);
    else
      behav.(type_counts_total) = NaN;
    end;
    type_counts_total = sprintf('%s_total_bothruns_std_rt',cond); 
    if fbflag
      behav.(type_counts_total) = std(rt_keep);
    else
      behav.(type_counts_total) = NaN;
    end;
    
    for j=1:2 % runs
      type_counts_run = sprintf('%s_total_run_%d_numtrials',cond,j);   
      behav.(type_counts_run) = NaN;
      type_counts_run = sprintf('%s_total_run_%d_mean_rt',cond,j);
      behav.(type_counts_run) = NaN;
      type_counts_run = sprintf('%s_total_run_%d_std_rt',cond,j);
      behav.(type_counts_run) = NaN;  
    end 
    for j=runs_ok % runs
      type_counts_run = sprintf('%s_total_run_%d_numtrials',cond,j);   
      ind_type_run = [intersect(find(runs==j),ind_type)];    
      behav.(type_counts_run) = length(ind_type_run); 
      type_counts_run = sprintf('%s_total_run_%d_mean_rt',cond,j);
      if fbflag
        rt_keep = get_rt(event_info(ind_type_run));
        behav.(type_counts_run) = mean(rt_keep);
      else
        behav.(type_counts_run) = NaN;
      end;
      type_counts_run = sprintf('%s_total_run_%d_std_rt',cond,j);
      if fbflag
        behav.(type_counts_run) = std(rt_keep);
      else
        behav.(type_counts_run) = NaN;
      end;
    end 

    acc = [event_info(ind_type).acc];

    for j=1:parms.nstims
      stim = parms.stimnames{j}; 
      eventname = sprintf('%s_%s',cond,stim);

      acc = [event_info(ind_type).acc]; 

      switch stim
        case 'pos_feedback'

          ind_keep = find(acc==1); 
          type_counts_stim = sprintf('%s_bothruns_numtrials',eventname);  
          if length(ind_keep) < 4
            if parms.verbose, fprintf('%s: %s < 4 \n',mfilename,type_counts_stim); end
            behav.feedback_flag = 0;
          end 
          
          behav.(type_counts_stim) = length(ind_keep); 
          type_acc_stim = sprintf('%s_bothruns_rate',eventname);
          behav.(type_acc_stim) = length(ind_keep)./length(ind_type); 
          type_rt_stim = sprintf('%s_bothruns_mean_rt',eventname);
          rt_keep = get_rt(event_info(ind_type(ind_keep)));
          behav.(type_rt_stim) = mean(rt_keep);  
          type_rt_stim = sprintf('%s_bothruns_std_rt',eventname);
          behav.(type_rt_stim) = std(rt_keep);

          for j=1:2 % runs
            type_counts_run = sprintf('%s_run_%d_numtrials',eventname,j);  
            behav.(type_counts_run) = NaN; 
            type_acc_run = sprintf('%s_run_%d_rate',eventname,j);
            behav.(type_acc_run) = NaN; 
            type_rt_run = sprintf('%s_run_%d_mean_rt',eventname,j);    
            behav.(type_rt_run) = NaN; 
            type_rt_run = sprintf('%s_run_%d_std_rt',eventname,j);
            behav.(type_rt_run) = NaN; 
          end
          for j=runs_ok
            type_counts_run = sprintf('%s_run_%d_numtrials',eventname,j); 
            ind_type_run = [intersect(find(runs==j),ind_type)];  
            ind_keep_run = [intersect(find(runs==j),ind_type(ind_keep))];    
            behav.(type_counts_run) = length(ind_keep_run); 
            type_acc_run = sprintf('%s_run_%d_rate',eventname,j);
            behav.(type_acc_run) = length(ind_keep_run)./length(ind_type_run); 
            type_rt_run = sprintf('%s_run_%d_mean_rt',eventname,j);
            rt_keep = get_rt(event_info(ind_keep_run));
            behav.(type_rt_run) = mean(rt_keep);
            type_rt_run = sprintf('%s_run_%d_std_rt',eventname,j);
            behav.(type_rt_run) = std(rt_keep);
          end

        case 'neg_feedback'
          
          ind_keep = find(acc==0); 
          type_counts_stim = sprintf('%s_bothruns_numtrials',eventname);  
          if length(ind_keep) < 4
            if parms.verbose, fprintf('%s: %s < 4 \n',mfilename,type_counts_stim); end
            behav.feedback_flag = 0;
          end 
          
          behav.(type_counts_stim) = length(ind_keep); 
          type_acc_stim = sprintf('%s_bothruns_rate',eventname);
          behav.(type_acc_stim) = length(ind_keep)./length(ind_type); 
          type_rt_stim = sprintf('%s_bothruns_mean_rt',eventname);
          if fbflag
            rt_keep = get_rt(event_info(ind_type(ind_keep)));
            behav.(type_rt_stim) = mean(rt_keep);
          else
            behav.(type_rt_stim) = NaN;
          end
          type_rt_stim = sprintf('%s_bothruns_std_rt',eventname);
          if fbflag
            behav.(type_rt_stim) = std(rt_keep);
          else
            behav.(type_rt_stim) = NaN;
          end;
 
          for j=1:2 % runs
            type_counts_run = sprintf('%s_run_%d_numtrials',eventname,j);  
            behav.(type_counts_run) = NaN; 
            type_acc_run = sprintf('%s_run_%d_rate',eventname,j);
            behav.(type_acc_run) = NaN;
            type_rt_run = sprintf('%s_run_%d_mean_rt',eventname,j);    
            if fbflag
              rt_keep = get_rt(event_info(ind_type_run));
              behav.(type_rt_run) = mean(rt_keep);
            else
              behav.(type_rt_run) = NaN;
            end
            type_rt_run = sprintf('%s_run_%d_std_rt',eventname,j);
            if fbflag
              behav.(type_rt_run) = std(rt_keep);
            else
              behav.(type_rt_run) = NaN;
            end
          end
          for j=runs_ok
            type_counts_run = sprintf('%s_run_%d_numtrials',eventname,j); 
            ind_type_run = [intersect(find(runs==j),ind_type)];  
            ind_keep_run = [intersect(find(runs==j),ind_type(ind_keep))];    
            behav.(type_counts_run) = length(ind_keep_run); 
            type_acc_run = sprintf('%s_run_%d_rate',eventname,j);
            behav.(type_acc_run) = length(ind_keep_run)./length(ind_type_run);
            type_rt_run = sprintf('%s_run_%d_mean_rt',eventname,j);
            if fbflag
              rt_keep = get_rt(event_info(ind_keep_run));
              behav.(type_rt_run) = mean(rt_keep);
            else
              behav.(type_rt_run) = NaN;
            end;
            type_rt_run = sprintf('%s_run_%d_std_rt',eventname,j);
            if fbflag
              behav.(type_rt_run) = std(rt_keep);
            else
              behav.(type_rt_run) = NaN;
            end
          end
          
      end
    end 
  end
  
  % collapse across small or large 
  collapse_cond = {'reward', 'loss'};
  collapse_ind_type.reward = [];
  collapse_ind_type.loss = [];
  
  for i=1:parms.nconds
    type = parms.typenames{i}; 
    cond = parms.condnames{i};  
    ind_type = find(strcmp(type,all_types)); 
    if findstr('reward', cond) > 0
      collapse_ind_type.reward = [collapse_ind_type.reward ind_type]; 
    elseif findstr('loss', cond) > 0
      collapse_ind_type.loss = [collapse_ind_type.loss ind_type]; 
    end
  end  
  
  for i=1:length(collapse_cond)
    cond = collapse_cond{i}; 
    type_counts_total = sprintf('hilo%s_total_bothruns_numtrials',cond); 
    behav.(type_counts_total) = size(collapse_ind_type.(cond),2); 
    type_rt_total = sprintf('hilo%s_total_bothruns_mean_rt',cond);
    if fbflag
      rt_keep = get_rt(event_info(ind_type));
      behav.(type_rt_total) = mean(rt_keep);
    else
      behav.(type_rt_total) = NaN;
    end
    type_rt_total = sprintf('hilo%s_total_bothruns_std_rt',cond); 
    if fbflag
      behav.(type_rt_total) = std(rt_keep);
    else
      behav.(type_rt_total) = NaN;
    end
    
    for j=1:2 % runs
      type_counts_run = sprintf('hilo%s_total_run_%d_numtrials',cond,j);    
      behav.(type_counts_run) = NaN;
      type_rt_run = sprintf('hilo%s_total_run_%d_mean_rt',cond,j);
      behav.(type_rt_run) = NaN;
      type_rt_run = sprintf('hilo%s_total_run_%d_std_rt',cond,j);
      behav.(type_rt_run) = NaN;
    end 
    for j=runs_ok
      type_counts_run = sprintf('hilo%s_total_run_%d_numtrials',cond,j); 
      ind_type_run = [intersect(find(runs==j),collapse_ind_type.(cond))];    
      behav.(type_counts_run) = length(ind_type_run); 
      type_rt_run = sprintf('hilo%s_total_run_%d_mean_rt',cond,j);      
      if fbflag
        rt_keep = get_rt(event_info(ind_type_run));
        behav.(type_rt_run) = mean(rt_keep);
      else
        behav.(type_rt_run) = NaN;
      end
      type_rt_run = sprintf('hilo%s_total_run_%d_std_rt',cond,j);
      if fbflag
        behav.(type_rt_run) = std(rt_keep);
      else
        behav.(type_rt_run) = NaN;
      end
    end 
    
    for j=1:parms.nstims
      stim = parms.stimnames{j}; 
      eventname = sprintf('%s_%s',cond,stim);   
    
      acc = [event_info(collapse_ind_type.(cond)).acc];  
    
      switch stim
        case 'pos_feedback'
    
          ind_keep = find(acc==1); 
          ind_type = collapse_ind_type.(cond); 
          type_counts_stim = sprintf('hilo%s_bothruns_numtrials',eventname);
          behav.(type_counts_stim) = length(ind_keep); 
          type_acc_stim = sprintf('hilo%s_bothruns_rate',eventname);
          behav.(type_acc_stim) = length(ind_keep)./length(ind_type); 
          type_rt_stim = sprintf('hilo%s_bothruns_mean_rt',eventname);
          rt_keep = get_rt(event_info(ind_type(ind_keep)));
          behav.(type_rt_stim) = mean(rt_keep); 
          type_rt_stim = sprintf('hilo%s_bothruns_std_rt',eventname);
          behav.(type_rt_stim) = std(rt_keep); 
          
          for j=1:2 % runs
            type_counts_run = sprintf('hilo%s_run_%d_numtrials',eventname,j);   
            behav.(type_counts_run) = NaN; 
            type_acc_run = sprintf('hilo%s_run_%d_rate',eventname,j);
            behav.(type_acc_run) = NaN; 
            type_rt_run = sprintf('hilo%s_run_%d_mean_rt',eventname,j);  
            behav.(type_rt_run) = NaN;   
            type_rt_run = sprintf('hilo%s_run_%d_std_rt',eventname,j);
            behav.(type_rt_run) = NaN;
          end 
          for j=runs_ok
            type_counts_run = sprintf('hilo%s_run_%d_numtrials',eventname,j);  
            ind_type_run = [intersect(find(runs==j),ind_type)];  
            ind_keep_run = [intersect(find(runs==j),ind_type(ind_keep))];    
            behav.(type_counts_run) = length(ind_keep_run); 
            type_acc_run = sprintf('hilo%s_run_%d_rate',eventname,j);
            behav.(type_acc_run) = length(ind_keep_run)./length(ind_type_run); 
            type_rt_run = sprintf('hilo%s_run_%d_mean_rt',eventname,j);
            rt_keep = get_rt(event_info(ind_keep_run));
            behav.(type_rt_run) = mean(rt_keep);    
            type_rt_run = sprintf('hilo%s_run_%d_std_rt',eventname,j);
            behav.(type_rt_run) = std(rt_keep);
          end 

        case 'neg_feedback'
    
          ind_keep = find(acc==0); 
          ind_type = collapse_ind_type.(cond); 
          type_counts_stim = sprintf('hilo%s_bothruns_numtrials',eventname);
          behav.(type_counts_stim) = length(ind_keep); 
          type_acc_stim = sprintf('hilo%s_bothruns_rate',eventname);
          behav.(type_acc_stim) = length(ind_keep)./length(ind_type); 
          type_rt_stim = sprintf('hilo%s_bothruns_mean_rt',eventname);
          if fbflag
            rt_keep = get_rt(event_info(ind_type(ind_keep)));
            behav.(type_rt_stim) = mean(rt_keep);
          else
            behav.(type_rt_stim) = NaN;
          end
          type_rt_stim = sprintf('hilo%s_bothruns_std_rt',eventname);
          if fbflag
            behav.(type_rt_stim) = std(rt_keep);
          else
            behav.(type_rt_stim) = NaN;
          end
          
          for j=1:2 % runs
            type_counts_run = sprintf('hilo%s_run_%d_numtrials',eventname,j);     
            behav.(type_counts_run) = NaN; 
            type_acc_run = sprintf('hilo%s_run_%d_rate',eventname,j);
            behav.(type_acc_run) = NaN;        
            type_rt_run = sprintf('hilo%s_run_%d_mean_rt',eventname,j);  
            behav.(type_rt_run) = NaN;    
            type_rt_run = sprintf('hilo%s_run_%d_std_rt',eventname,j);
            behav.(type_rt_run) = NaN;
          end
          for j=runs_ok
            type_counts_run = sprintf('hilo%s_run_%d_numtrials',eventname,j);  
            ind_type_run = [intersect(find(runs==j),ind_type)];  
            ind_keep_run = [intersect(find(runs==j),ind_type(ind_keep))];    
            behav.(type_counts_run) = length(ind_keep_run); 
            type_acc_run = sprintf('hilo%s_run_%d_rate',eventname,j);
            behav.(type_acc_run) = length(ind_keep_run)./length(ind_type_run);
            type_rt_run = sprintf('hilo%s_run_%d_mean_rt',eventname,j);
            if fbflag
              rt_keep = get_rt(event_info(ind_keep_run));
              behav.(type_rt_run) = mean(rt_keep);
            else
              behav.(type_rt_run) = NaN;
            end
            type_rt_run = sprintf('hilo%s_run_%d_std_rt',eventname,j);
            if fbflag
              behav.(type_rt_run) = std(rt_keep);
            else
              behav.(type_rt_run) = NaN;
            end
          end
      end        
    end
  end
  
  % money 
  money = [event_info.money]; 
  behav.('earnings_total') = sum(money); 
  for j=1:2 % runs
    type_counts_run = sprintf('earnings_run_%d',j);  
    behav.(type_counts_run) = NaN;   
  end   
  for j=runs_ok
    type_counts_run = sprintf('earnings_run_%d',j);  
    behav.(type_counts_run) = sum(money(find(runs==j)));   
  end   

  % check responses
  [trial_resp_array, valid_trial_resp_array] = check_responses(event_info);
  if trial_resp_array(1) < parms.min_run_resps || trial_resp_array(2) < parms.min_run_resps
    behav.perform_flag = 0;
  else
    behav.perform_flag = 1;
  end
  for j=1:2
    trialnumtrials = sprintf('trialresp_run_%d_numtrials',j);
    behav.(trialnumtrials) = trial_resp_array(j);
    validtrialnumtrials = sprintf('valid_trialresp_run_%d_numtrials',j);
    behav.(validtrialnumtrials) = valid_trial_resp_array(j);
  end
  trialtotalnumtrials = sprintf('trialresp_bothruns_numtrials');
  behav.(trialtotalnumtrials) = trial_resp_array(3);
  validtrialtotalnumtrials = sprintf('validtrialresp_bothruns_numtrials',j);
  behav.(validtrialtotalnumtrials) = valid_trial_resp_array(3);
  
  % write to csv
  fname_csv_out = sprintf('%s/%s_behavioral.csv',parms.outdir,parms.outstem); 
  if ~exist(fname_csv_out,'file') || parms.forceflag
    mmil_struct2csv(behav,fname_csv_out)
  end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check if response is premature
function premature = get_premature(event_info)
  premature = 0;
  response_check = mmil_getfield(event_info,'response_check');
  if strcmp(event_info.response_check,'You pressed too soon!'), premature = 1; end
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check button press
function [trial_resp_array, valid_trial_resp_array] = check_responses(event_info)
  runs = [event_info.trial];  
  nruns = length(unique(runs));
  trial_resp_array = []; valid_trial_resp_array = [];
  nresp_run1 = 0; nresp_valid_run1 = 0;
  nresp_run2 = 0; nresp_valid_run2 = 0;
  for k=1:length(event_info)
    premature = get_premature(event_info(k));
    response = ~isempty(mmil_getfield(event_info(k),'probe_resp')) ||...
               ~isempty(mmil_getfield(event_info(k),'text_display_resp')) ||...
               ~isempty(mmil_getfield(event_info(k),'feedback_resp'));
    if event_info(k).run == 1
      if response || premature, nresp_run1 = nresp_run1 + 1; end;
      if response && ~premature, nresp_valid_run1 = nresp_valid_run1 + 1; end;
    elseif event_info(k).run == 2
      if response || premature, nresp_run2 = nresp_run2 + 1; end;
      if response && ~premature, nresp_valid_run2 = nresp_valid_run2 + 1; end;
    end;
  end
  trial_resp_array = [nresp_run1,nresp_run2,nresp_run1 + nresp_run2];
  valid_trial_resp_array = [nresp_valid_run1,nresp_valid_run2,nresp_valid_run1 + nresp_valid_run2];
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get valid response times
%  function vals = get_rt(event_info)
%    vals = [];
%    if ~isfield(event_info,'rt_any')
%      vals = {event_info.rt};
%    else
%      vals = {event_info.rt_any};
%      for k=1:length(event_info)
%        premature = get_premature(event_info(k));
%        if premature, vals{k} = []; end;
%      end;
%    end;
%    if iscell(vals)
%      i_valid = find(~cellfun(@isempty,vals));
%      vals = cell2mat(vals(i_valid));
%    end;
%    vals = nonzeros(vals(~isnan(vals)));
%  return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get valid response times
function vals = get_rt(event_info)
  vals = [];
  if ~isfield(event_info, 'feedback_resp')
    vals = {event_info.rt};
  else
    vals = [];
    for k=1:length(event_info)
      premature = get_premature(event_info(k));
      if premature
        vals{k} = [];
      elseif event_info(k).acc
        vals{k} = event_info(k).rt;
      elseif event_info(k).text_display_rttime ~= 0
        vals{k} = event_info(k).text_display_rttime - event_info(k).probe_onset;
      elseif event_info(k).feedback_rttime ~= 0
        vals{k} = event_info(k).feedback_rttime - event_info(k).probe_onset;
      end
    end;
  end;
  if iscell(vals)
    i_valid = find(~cellfun(@isempty,vals));
    vals = cell2mat(vals(i_valid));
  end;
  vals = nonzeros(vals(~isnan(vals)));
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
