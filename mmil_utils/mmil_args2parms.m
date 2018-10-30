function parms = mmil_args2parms(options, options_schema, strict)
%function parms = mmil_args2parms(options, options_schema, [strict])
%
% Purpose:
%   Take a set of argument-value pairs and a set of argument-value
%   defaults, and create a 'parms' structure where all arguments are fields
%   of the parms object.
%
% Input Arguments:
%   options         : the 'optional arguments' from a function
%   options_schema  : cell array containing 3 values per known 'option':
%                     - option name
%                     - default value
%                     - allowed values
%                         vector of true/false values
%                         vector of min/max values
%                         vector of allowed values (more than 2 elements)
%                         cell array of allowed values
%                         empty to specify no limitations.
%  strict           : whether to fail if options not specified in the
%                       options_schema are found.
%   {default = true}
%
% Output Arguments:
%   parms           : structure containing named options as fields
%
% NOTE: do not pass empty cell {} for default value
%
% See also: mmil_parms2args
%
% Created:  07/15/07 by Ben Cipollini 
% Prev Mod: 10/09/13 by Don Hagler
% Last Mod: 06/01/17 by Don Hagler
%

if ~mmil_check_nargs(nargin,2), return; end;

if (0 ~= mod(length(options),2))     %   Validate that the # of args is even
  error('list of arguments must be even (must have name/value pair)');
end;
if (0 ~= mod(length(options_schema),3))    %   Validate the options_schema info is right
  error('programming error: list of default arguments must be divisible by three (must have name/default/allowed triplets)');
end;
if (~exist('strict', 'var'))
  strict = true;
end;
parms   = [];

%   Rip all cell arguments into non-cell arguments.
for index = 1:length(options)
  if (iscell(options{index}) && ~iscell(options{index}{1}))
    options{index} = { options{index} }; 
  end;
end;

valid_fields    = options_schema(1:3:end);
input_fields    = options(1:2:end);
if ~isempty(valid_fields)
  unknown_fields  = setdiff(input_fields, valid_fields);
else
  unknown_fields = [];
end;

%   Validate that there are no extraneous params sent in
if (strict && ~isempty(unknown_fields))
  error('the following unrecogized options were passed in: %s', sprintf('%s ',unknown_fields{:}));
end;

if (~isempty( options ))
  parms=struct(options{:});
end;

%   This allows 'pass-through' of parameters; 
%
%   remove any fields that are unknown
%
%   unless no schema is defined.
%
if (~strict && ~isempty(parms) && ~isempty(options_schema))
  parms = rmfield(parms,unknown_fields);
end;

%   Check arg values and set defaults
for f=1:3:length(options_schema)
  %   The value has been set explicitly by the caller;
  %   Validate the input parameters by the 'range' field
  if  (isfield(parms,options_schema{f}))
    param_name		= options_schema{f};
	  param_value		= getfield(parms,param_name);
    param_range   = options_schema{f+2};

    % no value was specified,
    if isempty(param_value)
      parms = setfield(parms, options_schema{f}, options_schema{f+1});
    % no range was specified,
    elseif isempty(param_range)
      continue;
   	%	param range is a cell array of strings; make sure the current value is within that range.
    elseif iscell(param_range)
      num_flag = 0;
      char_flag = 0;
      for i=1:length(param_range)
        if isnumeric(param_range{i}), num_flag=1; end;
        if ischar(param_range{i}), char_flag=1; end;
      end;
      if num_flag && char_flag
        error('type of parameter range (cell array of numbers and strings) specified for parameter ''%s'' is currently unsupported', options_schema{f});
      elseif char_flag
        if iscell(param_value)
          for i=1:length(param_value)
            if ~ischar(param_value{i})
              error('parameter ''%s'' must be string or cell array of strings', options_schema{f});
            end;
          end;
        elseif ~ischar(param_value)
          error('parameter ''%s'' must be string or cell array of strings', options_schema{f});
        else
          param_value = {param_value};
        end;
        if length(find(ismember(param_value,param_range))) ~= length(param_value)
          error('parameter ''%s'' value must be one of the following: { %s}', ...
  			    param_name, sprintf('''%s'' ',param_range{:}));
        end;
      elseif num_flag
        param_range = cell2mat(param_range);
        if ~isnumeric(param_value)
          error('parameter ''%s'' must be numeric', options_schema{f});
        end;
        if length(find(ismember(param_value,param_range))) ~= length(param_value)
          error('parameter ''%s'' value must be one of the following: { %s}', ...
  			    param_name,sprintf('%d ',param_range));
        end;
      else
        error('type of parameter range specified for parameter ''%s'' is currently unsupported', options_schema{f});
      end;
    %  param range is logical and has two elements (i.e. true/false)
    elseif islogical(param_range) && length(param_range)==2
      if ~ismember(param_value,param_range)
        error('parameter %s value must be true (1) or false (0)', options_schema{f});
      end;
    %  param range is numeric and has two elements (i.e. min and max)
    elseif isnumeric(param_range) && length(param_range)==2
      %  param range is numeric or logical, and within a specified range,
      if ~isempty(find(param_value < param_range(1))) || ...
         ~isempty(find(param_value > param_range(2)))
        if int64(param_range(1))==param_range(1) && int64(param_range(2))==param_range(2)
          error('parameter %s value must be between %d and %d',...
            options_schema{f},param_range(1),param_range(2));
        else
          error('parameter %s value must be between %0.4f and %0.4f',...
            options_schema{f},param_range(1),param_range(2));
        end;
      end;
    %  param range is numeric and has more than two elements (allowed values)
    elseif isnumeric(param_range)
      if ~isnumeric(param_value)
        error('parameter ''%s'' must be numeric', options_schema{f});
      end;
      if length(find(ismember(param_value,param_range))) ~= length(param_value)
        error('parameter ''%s'' value must be one of the following: [ %s]', ...
  			  param_name,sprintf('%d ',param_range));
      end;
    %  param range is of a type we currently don't support.
    else
      error('type of parameter range specified for parameter ''%s'' is currently unsupported', options_schema{f});
    end;
  % field not found, so set the default value.
  else
    parms = setfield(parms, options_schema{f}, options_schema{f+1});
  end;
end;

