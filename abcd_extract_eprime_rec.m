function [errcode,behav,errmsg] = abcd_extract_eprime_rec(fnamewm,fnamerec,varargin)
%function [errcode,behav,errmsg] = abcd_extract_eprime_rec(fnamewm,fnamerec,[options])
%
% Purpose: generate behav data for n-back REC files 
%
% Required input:
%   fnamewm: name of input eprime file WM (scanner run). 
%   fnamerec: name of input eprime file REC (recall run). 
%
% Optional input:
%   'outdir': output directory
%     {default = pwd}
%   'outstem': output file stem
%     if empty, will use filestem of fname
%     {default = []}
%   'switch_thresh': percent correct below which to swap button assignements
%     {default = 0.3}
%   'forceflag': [0|1] overwrite existing output
%     {default = 0}
%   'verbose': [0|1] display messages
%     {default = 1}
%
% Output: 
%    errcode: [0|1] whether the file was successfully processed
%    behav: behavioral data
%    errmsg: string describing error if errcode=1
%
% Created : 11/16/17 by Dani Cornejo 
% Prev Mod: 01/23/19 by Dani Cornejo
% Prev Mod: 05/21/20 by Octavio Ruiz
% Prev Mod: 08/28/20 by Don Hagler
% Last Mod: 11/03/20 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this section:
% Based on abcd_extract_eprime_nback.m 
% Created : 01/06/17 by Jose Teruel 
% Prev Mod: 09/13/17 by Don Hagler
% Last Mod: 11/08/17 by Dani Cornejo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize outputs
errcode = 0; behav = []; errmsg = [];

% check arguments 
if ~mmil_check_nargs(nargin,2), return; end;

% check input parameters
parms = check_input(fnamewm,fnamerec,varargin);

% create output directory
mmil_mkdir(parms.outdir);

if ~parms.errcode_nback 
  % create struct with info for each event in WM
  [event_info_wm,event_info_proc_wm,errcode,errmsg] = get_event_info_wm(parms);
  if errcode, return; end;
 
  % create struct with info for each event REC
  [event_info_rec,event_info_proc_rec,errcode,errmsg] = get_event_info_rec(parms); 
  if errcode, return; end;
  
  % switch buttons if necesary for REC 
  [event_info_rec,~,parms.switch_flag,errcode,errmsg] = rec_switch(event_info_rec,event_info_proc_rec,parms); 
  if errcode, return; end;
  
  % get behav data and write it to a csv 
  behav = get_behavioral_data_rec(event_info_wm,event_info_rec,parms); 
else
  % get behav data and write it to a csv
  get_behavioral_data_rec_empty(parms);
  if parms.verbose, fprintf('%s: ERROR: empty behavioral file created due to error with n-back file %s\n',...
    mfilename,parms.fnamewm); end
  errcode = 1;
  errmsg = 'error with n-back file';
end % if errcode_nback 

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parms = check_input(fnamewm,fnamerec,options)
  parms = mmil_args2parms(options,{...
    'fnamewm',fnamewm,[],...
    'fnamerec',fnamerec,[],...
    ...
    'outdir',pwd,[],...
    'outstem',[],[],...
    'switch_thresh',0.3,[0,0.5],...
    'forceflag',false,[false true],...
    'verbose',true,[false true],...
    'errcode_nback',0,0:1,...
    ...
    'colnameswm', {'NARGUID','SessionDate','SessionTime','ExperimentName','ExperimentVersion',...
                  'Procedure[Block]','BlockType','StimType','Stim.OnsetTime','Stim.OffsetTime','Stim.ACC',...
                  'Cue2Back.OnsetTime','Cue2Back.OffsetTime','CueTarget.OnsetTime','CueTarget.OffsetTime',...
                  'CueFix.OnsetTime','CueFix.OffsetTime','CueFix.StartTime',...
                  'Fix15sec.OnsetTime', 'Fix15sec.OffsetTime','TargetType','Stim.RT',...
                  'CorrectResponse','Stim.RESP','Stimulus'},[],...
    'fieldnameswm', {'narguid','date','time','experiment','version',...
                  'procedure_block','block_type','stim_type','stim_onset','stim_offset','stim_acc',...
                  'cue2back_onset','cue2back_offset','cue0back_onset','cue0back_offset',...
                  'cuefix_onset','cuefix_offset','cuefix_start',...
                  'fixation_onset', 'fixation_offset','target_type','stim_rt',...
                  'correct_response','stim_resp','stim'},[],...
    'colnamesrec', {'NARGUID','SessionDate','SessionTime','ExperimentName','ExperimentVersion',...
                  'Block','Stim.ACC[Block]','Stim.CRESP[Block]','Stim.RESP[Block]',...
                  'Stim.RT[Block]','StimType[Block]','Stimulus[Block]'},[],...
    'fieldnamesrec', {'narguid','date','time','experiment','version',...
                  'block','stim_acc_block','stim_cresp_block','stim_resp_block',...
                  'stim_rt_block','stim_type_block','stim_block'},[],...    
    'recnames', {'newneutface','newplace','newposface','newnegface','oldposface',...
                 'oldplace','oldnegface','oldneutface'},[],...              
    'typenames',{'2-Back','0-Back'},[],...
    'condnames',{'2_back','0_back'},[],...
    'stimnames',{'posface','neutface','negface','place'},[],...
    'extra_stimnames',{'target','lure','nonlure'},[],...
    'cues',{'Cue0BackPROC','Cue2BackPROC'},[]...
    'cuenames',{'cue0back','cue2back'},[],...
    'procedures',{'TRSyncPROC', 'TRSyncPROCR2'},[]...
  });
  parms.ncues = length(parms.cues);
  parms.nprocedures = length(parms.procedures);  
  parms.nconds = length(parms.condnames);
  parms.ntypes = length(parms.typenames);
  parms.nstims = length(parms.stimnames);
  parms.nextrastims = length(parms.extra_stimnames);
  parms.nrec = length(parms.recnames); 
  if parms.nconds ~= parms.ntypes
    error('condnames and typenames length mismatch');
  end;
  if ~exist(parms.fnamewm,'file')
    error('%s: file %s not found',mfilename,parms.fnamewm);
  end;
  if ~exist(parms.fnamerec,'file')
    error('%s: file %s not found',mfilename,parms.fnamerec);
  end;
  [~,fstem,~] = fileparts(parms.fnamerec);
  % remove problematic characters
  if isempty(parms.outstem)
    parms.outstem = abcd_clean_fstem(fstem);
  end;
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

function [event_info_wm,event_info_proc_wm,errcode,errmsg] = get_event_info_wm(parms)
  
  event_info_wm=[]; event_info_proc_wm = [];
  errcode = 0; errmsg = [];

  try
    % write event info to file
    fnamewm_csv = abcd_check_eprime_sprdsh(parms.fnamewm, parms.colnameswm,...
           parms.fieldnameswm, parms.outdir, parms.forceflag, parms.verbose);
    event_info_wm = mmil_csv2struct(fnamewm_csv);
  catch me
    if parms.verbose, fprintf('%s: ERROR: failed to read e-prime file %s:\n%s\n',...
      mfilename,parms.fnamewm,me.message); end
    errcode = 1;
    errmsg = 'failed to read e-prime file';
    return;
  end
  
  % remove non-events
  all_types = {event_info_wm.block_type}; 
  ind_events = find(~cellfun(@isempty,all_types));
  event_info_proc_wm = event_info_wm;
  event_info_wm = event_info_wm(ind_events);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [event_info_rec,event_info_proc_rec,errcode,errmsg] = get_event_info_rec(parms)
  
  event_info_rec = []; event_info_proc_rec = []; 
  errcode = 0; errmsg = [];
  
  try 
    % write event info to file
    fnamerec_csv = abcd_check_eprime_sprdsh(parms.fnamerec, parms.colnamesrec,...
            parms.fieldnamesrec, parms.outdir, parms.forceflag, parms.verbose);
    event_info_rec = mmil_csv2struct(fnamerec_csv); 
  catch me
    if parms.verbose, fprintf('%s: ERROR: failed to read e-prime file %s:\n%s\n',...
      mfilename,parms.fnamerec,me.message); end
    errcode = 1;
    errmsg = 'failed to read e-prime file';
    return;
  end
  
  % check experiment
  experiment = mmil_getfield(event_info_rec(1),'experiment');
  if isempty(regexpi(experiment,'rec'))
    if parms.verbose, fprintf('%s: ERROR: wrong experiment name in e-prime file %s: %s\n',...
      mfilename,parms.fnamerec,experiment); end
    errcode = 1;
    errmsg = 'wrong experiment name';
    return;
  else
    if parms.verbose, fprintf('%s: experiment name: %s\n',mfilename,experiment); end
  end;
  
  % remove non-events
  all_types = {event_info_rec.stim_type_block}; 
  ind_events = find(~cellfun(@isempty,all_types));
  event_info_proc_rec = event_info_rec;
  event_info_rec = event_info_rec(ind_events);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [event_info,event_info_proc,switch_flag,errcode,errmsg] = rec_switch(event_info,event_info_proc,parms)

  [event_info, event_info_proc] = rec_switch_resp_format(event_info,event_info_proc,parms); 

  [switch_flag,acc1,errcode,errmsg] = rec_switch_flag(event_info,parms);  
  if errcode, return; end;
   
  if switch_flag 
    [event_info, event_info_proc] = rec_switch_event(event_info,event_info_proc);
    [~,acc2,~] = rec_switch_flag(event_info,parms);
    if acc2 < acc1
      [event_info, event_info_proc] = rec_switch_event(event_info,event_info_proc);
       switch_flag = 0; 
    end 
  end 
  
  if switch_flag
    % write to csv
    fname_csv_out  = sprintf('%s/%s_events_switched.csv',...
      parms.outdir,parms.outstem); 
    if ~exist(fname_csv_out,'file') || parms.forceflag
      mmil_struct2csv(event_info_proc,fname_csv_out);
    end
  end 
  
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [event_info, event_info_proc] = rec_switch_event(event_info,event_info_proc)

  for i=1:length({event_info.stim_resp_block})
    if event_info(i).stim_resp_block == 1
      event_info(i).stim_resp_block = 2;
    elseif event_info(i).stim_resp_block == 2
      event_info(i).stim_resp_block = 1;
    end
    if event_info(i).stim_resp_block == event_info(i).stim_cresp_block
      event_info(i).stim_acc_block = 1;
    else 
      event_info(i).stim_acc_block = 0;
    end  
  end
  for i=1:length({event_info_proc.stim_resp_block})
    if event_info_proc(i).stim_resp_block == 1
      event_info_proc(i).stim_resp_block = 2; 
    elseif event_info_proc(i).stim_resp_block == 2
      event_info_proc(i).stim_resp_block = 1;
    end
    if (event_info_proc(i).stim_acc_block >= 0)
      if event_info_proc(i).stim_resp_block == event_info_proc(i).stim_cresp_block
         event_info_proc(i).stim_acc_block = 1;
      else 
         event_info_proc(i).stim_acc_block = 0;
      end 
    end
  end % for proc
  
return; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [switch_flag,accuracy,errcode,errmsg] = rec_switch_flag(event_info,parms) 
  
  errcode = 0; errmsg = [];
  switch_flag = 0; 
  correct_resp = {event_info.stim_cresp_block};
  resp = {event_info.stim_resp_block}; 
  correct_total = size(cell2mat({event_info.stim_cresp_block}),2);   
  correct = 0;  
  
  for i=1:correct_total
    if correct_resp{i}==resp{i}
      correct = correct+1;
    end;
  end;
  accuracy = 100*correct/correct_total;
  if parms.verbose, fprintf('%s: accuracy = %0.1f%%\n',mfilename,accuracy); end

  if accuracy == 0
    if parms.verbose, fprintf('%s: ERROR: accuracy equals zero for %s\n',mfilename,parms.fnamerec); end
    errcode = 1;
    errmsg = 'accuracy equals zero';
    return; 
  elseif accuracy < 100*parms.switch_thresh
    if parms.verbose, fprintf('%s: accuracy < %0.1f%%, switching button responses\n',...
                              mfilename,100*parms.switch_thresh); end
    switch_flag = 1;
  end;
  
return; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [event_info, event_info_proc] = rec_switch_resp_format(event_info,event_info_proc,parms) 

  correct_resp = {event_info.stim_cresp_block};
  correct_resp_proc = {event_info_proc.stim_cresp_block};
  resp = {event_info.stim_resp_block}; 
  resp_proc = {event_info_proc.stim_resp_block};
  
  if strcmp(correct_resp{1},'RIGHTARROW') || strcmp(correct_resp{1},'LEFTARROW')
    if parms.verbose, fprintf('%s: switching response format \n',mfilename); end
    for i=1:length(correct_resp)
      if strcmp(correct_resp{i},'RIGHTARROW') 
        event_info(i).stim_cresp_block = 1;
      elseif strcmp(correct_resp{i},'LEFTARROW') 
        event_info(i).stim_cresp_block = 2;
      end
    end
    for i=1:length(resp)
      if strcmp(resp{i},'RIGHTARROW') 
        event_info(i).stim_resp_block = 1;
      elseif strcmp(resp{i},'LEFTARROW') 
        event_info(i).stim_resp_block = 2;
      end
    end
    for i=1:length(correct_resp_proc)
      if strcmp(correct_resp_proc{i},'RIGHTARROW') 
        event_info_proc(i).stim_cresp_block = 1;
      elseif strcmp(correct_resp_proc{i},'LEFTARROW') 
        event_info_proc(i).stim_cresp_block = 2;
      end
    end
    for i=1:length(resp_proc)
      if strcmp(resp_proc{i},'RIGHTARROW') 
        event_info_proc(i).stim_resp_block = 1;
      elseif strcmp(resp_proc{i},'LEFTARROW') 
        event_info_proc(i).stim_resp_block = 2;
      end
    end
  end %if 'ARROW'
return; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function behav = get_behavioral_data_rec(event_info_wm,event_info_rec,parms) 
  
  behav = []; 
  behav.('SubjID') = []; behav.('VisitID') = []; 
  behav.switch_flag = parms.switch_flag;
  
  pictures = []; 
  
  all_types_wm = {event_info_wm.block_type}; 
  all_stims_wm = {event_info_wm.stim_type}; 
  all_targets_wm = {event_info_wm.target_type};
  all_pics_wm = {event_info_wm.stim}; 
    
  for i=1:parms.nconds
    type = parms.typenames{i}; 
    cond = parms.condnames{i}; 
    ind_type = find(strcmp(type,all_types_wm));
    
    for j=1:parms.nstims
      stim = parms.stimnames{j}; 
      ind_stim = find(strcmp(stim,lower(all_stims_wm)));
    
      for k=1:parms.nextrastims
        extra = parms.extra_stimnames{k};
        ind_extra = find(strcmp(extra,lower(all_targets_wm))); 
        eventname = sprintf('block_%s_%s_%s',cond,stim,extra); 
        ind = intersect(intersect(ind_type,ind_stim),ind_extra);
        pics = unique(all_pics_wm(ind));  
      
        if strcmp(extra,'nonlure')
          pics_1 = unique(all_pics_wm(ind));
          event_lure = sprintf('block_%s_%s_%s',cond,stim,'lure');
          pics_2 = pictures.(event_lure);
          pics = setdiff(pics_1,pics_2);
          pictures.(eventname)= pics; 
        else
          pictures.(eventname)= pics;   
        end    
      end 
    end 
  end %ncond 
   
  all_rec_types = {event_info_rec.stim_type_block};
  all_stim_acc = {event_info_rec.stim_acc_block};
  all_rec_pics = {event_info_rec.stim_block};
  
  % exclude events with empty acc
  i_valid = find(~cellfun(@isempty,all_stim_acc));
  all_rec_types = all_rec_types(i_valid);
  all_stim_acc = cell2mat(all_stim_acc(i_valid));
  all_rec_pics = all_rec_pics(i_valid);

  for l=1:parms.nrec
    rectype = parms.recnames{l};
    ind_type_rec = find(strcmp(rectype,all_rec_types)); 
    eventname = sprintf('%s_acc',rectype);
    event_acc = all_stim_acc(ind_type_rec);
    %behav.(eventname) =  event_acc;
    eventname = sprintf('%s_hr',rectype);
    behav.(eventname) =  mean(event_acc);
    eventname = sprintf('%s_fa',rectype);
    behav.(eventname) =  1-mean(event_acc);
    eventname = sprintf('%s_pic',rectype);
    event_pics = all_rec_pics(ind_type_rec);  
    %behav.(eventname) =  event_pics;
    
    if strfind(rectype,'old')
      for i=1:parms.nconds
        cond = parms.condnames{i}; 
        event_pics_cond = []; event_acc_cond = []; 
        for j=1:length(event_pics)
          found_stim = find_field(pictures,event_pics(j),cond); 
          if found_stim 
            event_pics_cond = [event_pics_cond event_pics(j)]; 
            event_acc_cond = [event_acc_cond event_acc(j)]; 
          end 
        end
        eventname = sprintf('%s_%s_acc',rectype,cond);  
        %behav.(eventname) = event_acc_cond; 
        eventname = sprintf('%s_%s_pic',rectype,cond);  
        %behav.(eventname) = event_pics_cond; 
        eventname = sprintf('%s_%s_hr',rectype,cond);
        behav.(eventname) =  mean(event_acc_cond);
        eventname = sprintf('%s_%s_fa',rectype,cond);
        behav.(eventname) =  1-mean(event_acc_cond);
      end 
    end %if_old
  end
  
  for i=1:parms.nstims
    stims = parms.stimnames{i};
    eventname = sprintf('%s_pr',stims); 
    hr_name = sprintf('old%s_hr',stims); %oldposface_hr
    fa_name = sprintf('new%s_fa',stims);
    pr = behav.(hr_name) - behav.(fa_name); 
    behav.(eventname) = pr; 
    eventname = sprintf('%s_br',stims); 
    br = (behav.(fa_name)/(1-(behav.(hr_name)-behav.(fa_name))))-0.5; 
    behav.(eventname) = br;
    
    eventname = sprintf('%s_dprime',stims); 
    if behav.(hr_name) == 1
      hrd = 0.9948; 
    elseif behav.(hr_name) == 0
      hrd = 0.0052;
    else 
      hrd = behav.(hr_name);
    end 
    if behav.(fa_name) == 1
      fad = 0.9948; 
    elseif behav.(fa_name) == 0
      fad = 0.0052;
    else 
      fad = behav.(fa_name);
    end
    dprime = norminv(hrd) - norminv(fad);       
    %Dprime = z(Happy_HRd) - z(Happy_FAd)
    behav.(eventname) = dprime;      
  end 
   
  % write to csv
  fname_csv_out = sprintf('%s/%s_behavioral.csv',parms.outdir,parms.outstem); 
  if ~exist(fname_csv_out,'file') || parms.forceflag
    mmil_struct2csv(behav,fname_csv_out)
  end;
  
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function get_behavioral_data_rec_empty(parms) 
  
  behav = []; 
  behav.('SubjID') = []; behav.('VisitID') = []; 
  behav.switch_flag = 0;
  
  for l=1:parms.nrec
    rectype = parms.recnames{l}; 
    eventname = sprintf('%s_hr',rectype);
    behav.(eventname) =  'NaN'; 
    eventname = sprintf('%s_fa',rectype);
    behav.(eventname) =  'NaN'; 
    
    if strfind(rectype,'old')
      for i=1:parms.nconds
        cond = parms.condnames{i}; 
        eventname = sprintf('%s_%s_hr',rectype,cond);
        behav.(eventname) =  'NaN'; 
        eventname = sprintf('%s_%s_fa',rectype,cond);
        behav.(eventname) =  'NaN'; 
      end 
    end %if_old
  end
  
  for i=1:parms.nstims
    stims = parms.stimnames{i};
    eventname = sprintf('%s_pr',stims);
    behav.(eventname) = 'NaN'; 
    eventname = sprintf('%s_br',stims); 
    behav.(eventname) = 'NaN'; 
    eventname = sprintf('%s_dprime',stims); 
    behav.(eventname) = 'NaN'; 
  end 
   
  % write to csv
  fname_csv_out = sprintf('%s/%s_behavioral.csv',parms.outdir,parms.outstem); 
  if ~exist(fname_csv_out,'file') || parms.forceflag
    mmil_struct2csv(behav,fname_csv_out)
  end;
  
return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function found_stim = find_field(pictures,image,looking)
  fields = fieldnames(pictures); 
  found_field = []; 
  for p=1:numel(fields)
    pic_field = [pictures.(fields{p})];
    if ismember(image,pic_field)
     found_field = [found_field fields{p}]; 
    end 
  end
  found_stim = strfind(found_field,looking); 
return; 
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  
