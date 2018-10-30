#!/usr/bin/env python3

import sys, io, os
import datetime, dateutil.parser

# Import pandas, avoiding warning
import warnings
warnings.simplefilter(action='ignore', category=UserWarning)
import pandas as pd
pd.set_option('display.width', 1024)
pd.set_option('max_colwidth',   256)

import glob, json

from difflib import SequenceMatcher
import copy


#---------------------------------------------------------------------------------------------------
Info_Prot = {
    'fname':       None,
    'encoding':    None,
    'sep_tab':     None,
    'quoted_rows': None,
    'hoffs':       None,
    'n_rows':      None,
    'n_cols':      None,
    'exp_in_file': None,
    'exper':       None,
    'exp_datime':  None,
    'delay' :      float('nan'),
    'exp_t0':      '0000-00-00 00:00:00',
    'tdiff' :      float('nan'),
    'fout_ok' :       None,
    'fout_msg':       '',
    'diagnos':       -256,
    'pGUIDmatch':     0.0,
    'naming_ok':       0,
    'fname_exp_match': 0,
    'ok' : None,
    'msg': ''
}

# Task keys and corresponding file-name identifiers
tsk_fid_dict = {'MID':   ['MID'],
                'nBack': ['NBACK','WM','REC'],
                'SST':   ['SST'] }

encd_optns_list = ['utf-8','utf-16']

# Acceptable extreme values
rows_n_min      =  2
rows_head_n_max =  5
cols_n_min      = 10

Verbose = False
#---------------------------------------------------------------------------------------------------


# =============================================================================================================================================
# ---------------------------------------------------------------------------------------------------------------
def program_description():
    print()
    print('Read text file(s) containing a E-Prime spreadsheet(s), detect encoding and format, interpret content,')
    print("check if experiment matches file name, extract experiment's date and time and other information.")
    print('This program can also locate and check a set of files under a given directory or path.')
    print('                                                      Octavio Ruiz.  2017jun05-nov01, 2018feb19-jun25')
    print('Usage:')
    print('  ./eprime_sprdsht_get.py                          Print this help')
    print('  ./eprime_sprdsht_get.py file Summary             Read file, print summary of file encoding and contents')
    print('  ./eprime_sprdsht_get.py file Info                Read file, print diagnostics and file information')
    print('  ./eprime_sprdsht_get.py file InfoCheckName       Read file, print diagnostics and file information,')
    print('                                                   check that file name matches experiment in spreadsheet')
    print('  ./eprime_sprdsht_get.py file ExportFile outfile  Read file, and write a new file containing the interpreted EPrime data in our standard format')
    print('                                                   (tab-separated, no comments before header line, extension = .txt)')
    print('                                                   outfile should be entered with no extension')
    print()
    print('  ./eprime_sprdsht_get.py dir                                                           Read files in directory and print a table with diagnostics, experiment date, time, and filename')
    print('  ./eprime_sprdsht_get.py dir PickFile                                                  Idem, suppress listing, pick file with best format and diagnostics, print diagnostics and file information')
    print('  ./eprime_sprdsht_get.py dir PickFile "YYYYmmdd HHMMSS"                                Idem, and calculate difference (in minutes) relative to given time')
    print('  ./eprime_sprdsht_get.py dir PickFile "YYYYmmdd HHMMSS" pGUID                          Idem, order file names by specified subject')
    print('  ./eprime_sprdsht_get.py dir PickFile "YYYYmmdd HHMMSS" pGUID Task Info                Idem, filter for files containing specified task, pick file with smallest pGUID and time differences')
    print('                                                                                        Suppress table listing, print diagnostics and file information')
    print('  ./eprime_sprdsht_get.py dir PickFile "YYYYmmdd HHMMSS" pGUID Task ExportFile outfile  Idem, write picked file using our standard format')
    print('                                                                                        (tab-separated, no comments before header line, extension = .txt)')
    print('                                                                                        outfile should be entered with no extension')
    print('Where')
    print('  Task is  "" or one of', tsk_fid_dict.keys() )
    print()
    print('Info, when argument is a file, prints:')
    print('    diagnos ,  pGUIDmatch ,  naming_ok ,  exp_t0 ,  task ,  full_file_fname')
    print()
    print('Info, when arguments are  dir PickFile "YYYYmmdd HHMMSS" ,  prints:')
    print('    diagnos ,  pGUIDmatch ,  naming_ok ,  exp_t0 ,  task ,  time_diff(minutes) ,  full_file_fname')
    print()
    print('diagnos is a sum of values from the following code:')
    print('  diag = 0  File not found')
    print('  diag > 0  File was found and read succesfully. Encoding and format are:')
    print('         1   =>   Encoding = utf-8')
    print('         2   =>   Encoding = utf-16')
    print('         4   =>   Separator = Tab, instead of comma')
    print('         8   =>   Quoted rows')
    print('        16   =>   One or more rows before column-names header')
    print('        64   =>   Experiment in spreadsheet matches file name (returned only if option = "FileNameCheck")')
    print('  diag < 0  File was found and read, but it is not acceptable:')
    print('   -1..-16   =>   Format as described for diag > 0')
    print('       -32   =>   Unable to recognize experiment in spreadsheet')
    print('       -64   =>   Experiment in spreadsheet does not match file name, or practice experiment (returned only if option = "FileNameCheck")')
    print('      -128   =>   Start_time_info not found or unable to extract')
    print('      -256   =>   Un-diagnosed error')
    print('diagnos is returned to the shell as exit-status code, that can be checked with "echo $?"')
    print()
    print('naming_ok = 1  =>  file name starts with "NDAR_INV"')
    print()
    print('Spreadsheets with practice experiments are considered not valid.')
    print()
    print('Examples:')
    print('  ./eprime_sprdsht_get.py   .  PickFile  "20170520 164900"  ""  ""  Info')
    print('  ./eprime_sprdsht_get.py   .  PickFile  "20170520 164900"  NDAR_INVJPLWZ1Z0  ""   Info')
    print('  ./eprime_sprdsht_get.py   .  PickFile  "20170520 164900"  NDAR_INVJPLWZ1Z0  MID  Info')
    print('  ./eprime_sprdsht_get.py   /space/syn07/1/data/ABCD/DAL_ABCD_QC/aux_incoming/CUB/NDAR_INVJPLZW1Z0/baseline_year_1_arm_1/sst-exported_NDAR_INVJPLZW1Z0_baseline_year_1_arm_1_SessionC')
    print('  ./eprime_sprdsht_get.py  "/space/syn07/1/data/ABCD/DAL_ABCD_QC/aux_incoming/*/*INVB9CDPZUA*/*/*exported*"')
    print()
# ---------------------------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------
def File_Diagnostics( EPrime_Info ):
     diagnos  = 0
     diagnos += encd_optns_list.index( EPrime_Info['encoding'] ) + 1  if  EPrime_Info['encoding'] in encd_optns_list  else  0
     diagnos += EPrime_Info['sep_tab']     *  4
     diagnos += EPrime_Info['quoted_rows'] *  8
     diagnos += 16 if (EPrime_Info['hoffs'] > 0) else 0

     EPrime_Info['diagnos'] = diagnos

     return EPrime_Info
# ---------------------------------------------------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------
def command_line_get_variables():
    cmnd_syntx_ok = True
    fname = ''
    optn  = ''
    fname_out = ''
    subj = ''
    ref_time  = ''
    task  = '' 
    optn2 = ''

    if 2 <= len(sys.argv) and len(sys.argv) <= 8:
        fname = sys.argv[1]

        if len(sys.argv) == 2:
            optn = 'ListFiles'

        elif len(sys.argv) >= 3:
            optn = sys.argv[2]

            if optn in ['Summary', 'Info', 'InfoCheckName']:
                pass

            elif optn == 'ExportFile' and len(sys.argv) == 4:
                fname_out = sys.argv[3]

            elif optn == 'PickFile':
                if len(sys.argv) == 3:
                    pass

                elif len(sys.argv) == 4:
                    ref_time = sys.argv[3]

                elif len(sys.argv) == 5:
                    ref_time = sys.argv[3]
                    subj = sys.argv[4]
                    pass

                elif len(sys.argv) == 6:
                    cmnd_syntx_ok = False
                    program_description()

                elif len(sys.argv) >= 7:
                    ref_time = sys.argv[3]
                    subj = sys.argv[4]
                    task = sys.argv[5]

                    if len(sys.argv) == 7 and sys.argv[6] in ['Info', 'InfoCheckName']:
                        optn2 = sys.argv[6]

                    elif len(sys.argv) == 8 and sys.argv[6] == 'ExportFile':
                        optn2     = sys.argv[6]
                        fname_out = sys.argv[7]

                    else:
                        cmnd_syntx_ok = False
                        program_description()
                else:
                    cmnd_syntx_ok = False
                    program_description()
            else:
                cmnd_syntx_ok = False
                program_description()
    else:
        cmnd_syntx_ok = False
        program_description()

    return cmnd_syntx_ok, fname, optn, fname_out, ref_time, subj, task, optn2
# ---------------------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------------------
def file_read( fname ):
    # Read file into a text buffer.  Determine encoding (UTF-16 or -8), number of rows,
    # if rows are encased by quotation marks, and what are the field separators
    buf = None
    encoding  = None
    txt_lngth = None
    num_lines = None
    sep_tab   = None
    itab      = None
    quoted_rows = None
    iquot       = None
    line_typic_seps_num = None

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    txt = ''
    ienc = 0
    txt_lines = []

    for i, encd in enumerate(encd_optns_list):
        try:
            ff = open(fname, encoding=encd)
            txt = ff.read()
            ff.close()
            
            txt_lngth = len(txt)
            txt_lines = txt.split('\n')
            num_lines = len(txt_lines)
            if num_lines >= rows_n_min:   # and got to this point:
                encoding = encd
                ienc = i+1
                break
        except:
            # File reading failed => not the right encoding, try another one
            pass

    if not encoding:
        return buf, encoding, txt_lngth, num_lines, sep_tab, itab, quoted_rows, iquot, line_typic_seps_num
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Tab- or comma-separated values?
    if txt.count('\t') > cols_n_min:
        sep = '\t'
        sep_tab = True
        itab = 1
    else:
        sep = ','
        sep_tab = False
        itab = 2
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    #            Are rows encased by quotation marks? - If so, remove quotes
    # check a "typical" line, after min. number of initial rows
    i_typ = int( (rows_n_min + num_lines) / 2 )
    line = txt_lines[i_typ]
    if line.startswith('"') and line.endswith('"'):
        quoted_rows = True
        iquot = 1
    else:
        quoted_rows = False
        iquot = 0

    if quoted_rows:
        txt = ''
        for line in txt_lines:
            txt = txt + line.strip('"') + '\n'
        txt_lines = txt.split('\n')
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # line_typic_seps_num = txt_lines[i_typ].count( sep )

    # Many times, when an EPrime file is uploaded at the site, its encoding is Microsoft, not Unix
    # When I read those files, there are blank lines between each spreadsheet line.
    # Here I check for this case
    # print()
    # print('i_typ=', i_typ )
    # print( [i  for i in range(i_typ-2, i_typ+3)] )
    # print( [txt_lines[i].count( sep )  for i in range(i_typ-2, i_typ+3)] )
    # print( pd.Series( [txt_lines[i].count( sep )  for i in range(i_typ-2, i_typ+3)] ) )
    # print( pd.Series( [txt_lines[i].count( sep )  for i in range(i_typ-2, i_typ+3)] ).min() )
    # print()
    line_typic_seps_num  = pd.Series( [txt_lines[i].count( sep )  for i in range(i_typ-2, i_typ+3)] ).min()

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Store "de-quoted" text as a readable memory object, that can be read into a pd.DataFrame
    buf = io.StringIO( txt )
    buf.flush()
    buf.seek(0)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    return buf, encoding, txt_lngth, num_lines, sep_tab, itab, quoted_rows, iquot, line_typic_seps_num
# ---------------------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------------------
def sprdsht_read( sprdsh, sep_tab, num_lines, line_typic_seps_num ):
    # Some files have comments in the first line or lines, instead of the column names.
    # Read buffer using increasing row offsets until starting a set of column names consistent with an EPrime file:
    # first row defining  number of columns > cols_n_min  and  number of columns ~ line_typic_seps_num
    data  = pd.DataFrame()
    hoffs = None
    ok  = True
    msg = ''

    for hoffs in range(0, rows_head_n_max):
        data = pd.DataFrame()
        sprdsh.seek(0)
        try:
            if sep_tab:
                data = pd.read_csv(sprdsh, header=hoffs, sep='\t')
            else:
                data = pd.read_csv(sprdsh, header=hoffs )
        except:
            pass

        if (len(data.columns) > cols_n_min) and (len(data.columns) >= (0.9 * line_typic_seps_num)):
            # If row contains column names, there must be none or very few empty (unnamed) cells
            if ['Unnamed' in s for s in data.columns].count(True) < 5:
                break

    # Many times, when an EPrime file is uploaded at the site, its encoding is Microsoft, not Unix
    # When I read those files, there are blank lines between each spreadsheet line.
    if line_typic_seps_num > 0:

        if (len(data) > 0.8 * num_lines) and (len(data.columns) >= (0.8 * line_typic_seps_num)):
            ok = True
        else:
            ok = False
            msg = "Error: unable to interpret file (1). "

    else:
        if (len(data) > 0.4 * num_lines) and len(data.columns) > cols_n_min:
            ok = True
        else:
            ok = False
            msg = "Error: unable to interpret file (2). "

    return data, hoffs, ok, msg
# ---------------------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------------------
def ExperimentCheck( data, fname ):
    exp_in_file = None
    exper = 'Unknown'
    fname_exp_match = 0
    ok  = True
    msg = ''
    exp_diagnos = 0
    
    if data.shape[0] >= rows_n_min and data.shape[1] >= cols_n_min:
        pass
    else:
        ok = False
        msg = 'Unable to interpret file as an EPrime spreadsheet. '
        exp_diagnos = -32
        return exp_in_file, exper, fname_exp_match, ok, msg, exp_diagnos


    # Find experiment reported in spreadsheet
    if 'ExperimentName' in data.columns[0]:   # It can be 'ExperimentName' or '0-ExperimentName' or ...?...
        ExpNameCol = data.columns[0]
        exp_in_file = data[ExpNameCol].unique()[0]

        if Verbose:
            print( 'exp_in_file:', exp_in_file )

        if isinstance(exp_in_file, str):   # Some files have a number instead of an experim.description

            if 'practice' in exp_in_file.lower():
                exper = 'Practice'
                ok = False
                msg = 'Practice experiment. '
                exp_diagnos = -64

            else:
                exp_fnameID_list = []
                for tsk in tsk_fid_dict.keys():
                    if tsk.lower() in exp_in_file.lower():
                        exper = tsk
                        exp_fnameID_list = tsk_fid_dict[tsk]
                        break

                if Verbose:
                    print('exp_fnameID_list =', exp_fnameID_list )

                if len(fname):
                    # Check if experiment matches file name identifier
                    for exp_fn in exp_fnameID_list:
                        if fname.lower().rfind(exp_fn.lower()) > (len(fname)-13):   # (look for task id near the end of file name)
                            fname_exp_match = 1
                            break
                else:
                    pass

                if len(fname):
                    if fname_exp_match:
                        exp_diagnos =  64
                # else:
                #   Not an error; but we can check later fname_exp_match if requested.

        else:
            ok = False
            msg = 'Invalid experiment information in spreadsheet. '
            exp_diagnos = -32

    else:
        ok = False
        msg = 'Unable to recognize experiment in spreadsheet. '
        exp_diagnos = -32
    
    return exp_in_file, exper, fname_exp_match, ok, msg, exp_diagnos
# ---------------------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------------------
#                                                     Main function
def EPrime_Info_and_Data_get( fname ):
    Info = copy.deepcopy( Info_Prot )
    data = pd.DataFrame()
    exp_diagnos = 0

    # # Does file name contain official pGUID format?
    # Info['naming_ok'] = 0
    # f_nopath_name = os.path.basename( fname )
    # if f_nopath_name.startswith('NDAR_INV') and f_nopath_name.count('_') == 2:
    #     Info['naming_ok'] = 1

    # Does file name complies with ABCD format: official pGUID format + _taskname ?
    Info['naming_ok'] = 0
    f_nopath_name = os.path.basename( fname )
    if f_nopath_name.startswith('NDAR_INV') and f_nopath_name.count('_') == 2:
        pGUID_short = f_nopath_name.replace('NDAR_INV','')
        if len(pGUID_short.split('_')[0]) == 8:
            Info['naming_ok'] = 1

    # ---------------------------------------------------------------------------------------------------------------
    #                                 Read file and determine encoding and format

    sprdsh, encoding, txt_lngth, lines_num, sep_tab, itab, quoted_rows, iquot, line_typic_seps_num  =  file_read( fname )
    
    Info['fname']       = fname
    Info['encoding']    = encoding
    Info['sep_tab']     = sep_tab
    Info['quoted_rows'] = quoted_rows

    if Verbose:
        print('encoding =', encoding, ',   txt_lngth =', txt_lngth, ',   lines_num =', lines_num,
           ',   sep_tab =', sep_tab, ',   quoted_rows =', quoted_rows, ',   line_typic_seps_num =', line_typic_seps_num )

    if not encoding:
        Info['ok']  = False
        Info['msg'] = 'Error: unable to read file'
        Info['diagnos'] = 0
        if Verbose:
            print( Info['msg'], '\n')
        return Info, data
    # ---------------------------------------------------------------------------------------------------------------

    # ---------------------------------------------------------------------------------------------------------------
    #                                  Read spreadsheet into a pandas DataFrame

    data, hoffs, ok, msg  =  sprdsht_read( sprdsh, sep_tab, lines_num, line_typic_seps_num )

    Info['hoffs']  = hoffs
    Info['n_rows'] = data.shape[0]
    Info['n_cols'] = data.shape[1]
    Info['ok']  = ok
    Info['msg'] = msg

    if ok:
        # Further checks: discard plain-text files that are not really spreadsheets
        if Info['n_rows'] < rows_n_min or Info['n_cols'] < cols_n_min:
            Info['ok'] = False
            Info['msg'] += 'Invalid num. of rows or cols. '

    if not Info['ok']:
        if Verbose:
            print('sprdsh_ok =', ok )
            print('Info =', Info, '\n' )

        Info = File_Diagnostics( Info )
        Info['diagnos'] = -Info['diagnos']
        return Info, data


    if Verbose:
        print('Spreadsheet:  shape =', data.shape, ',  hoffs =', hoffs )
        print('First two and last two rows:')
        print( data.iloc[[0,1,-2,-1]] )
        print()
    # ---------------------------------------------------------------------------------------------------------------


    # ---------------------------------------------------------------------------------------------------------------
    #             Check that file contains a known task.

    exp_in_file, exper, fname_exp_match, ok, msg, exp_diagnos  =  ExperimentCheck( data, fname )

    Info['exp_in_file']     = exp_in_file
    Info['exper']           = exper
    Info['fname_exp_match'] = fname_exp_match
    Info['ok']  = ok
    Info['msg'] = Info['msg'] + msg

    if Verbose:
        print('experim_check_ok =', Info['ok'] )

    if not Info['ok']:
        if Verbose:
            print('Info =', Info, '\n' )

        Info = File_Diagnostics( Info )

        # Add exp_diagnos to File_Diagnostics.
        # If exp_diagnos < 0, make the file-encoding number negative
        if exp_diagnos > 0:
            Info['diagnos'] = exp_diagnos + Info['diagnos']
        else:
            Info['diagnos'] = exp_diagnos - Info['diagnos']
        return Info, data


    if Verbose:
        print('Task reported in file:', exper, '. ', msg)
        print()
    # ---------------------------------------------------------------------------------------------------------------


    # ---------------------------------------------------------------------------------------------------------------
    #                               Extract time information and task starting time

    exp_date_str = data.iloc[0]['SessionDate']
    exp_time_str = data.iloc[0]['SessionTime']
    delay = float('nan')
    try:
        exp_datime = dateutil.parser.parse( exp_date_str + ' ' + exp_time_str )
    except:
        exp_diagnos += -128
        Info['ok']  = False
        Info['msg'] += 'Unable to extract or interpret date or time. '

    if Verbose:
        print( Info, '\n')

    if Info['ok']:
        Info['exp_datime'] = exp_datime

        # if '_REC' in fname or '_WM' in fname:
        #     # File contains an nBack recall experiment, performed outside the scanner,
        #     # therefore there is no waiting for a sync.pulse from the scanner to begin the experiment
        #     Info['exp_t0'] = Info['exp_datime']
        #     ok = True

        if  Info['exper'] == 'nBack' and ('Rec' in Info['exp_in_file'] or 'REC' in Info['exp_in_file']):
            # File contains an nBack recall experiment, performed outside the scanner,
            # therefore there is no waiting for a sync.pulse from the scanner to begin the experiment
            Info['delay']  = float('nan')
            Info['exp_t0'] = Info['exp_datime']
            ok = True

        else:
            if 'GetReady.RTTime' in data.columns:
                delay = data.iloc[0]['GetReady.RTTime']
                if Verbose:
                    print( data.iloc[[0,1,-2,-1]][['SessionDate', 'SessionTime', 'GetReady.RTTime']] )

            elif 'Wait4Scanner.RTTime' in data.columns:
                delay = data.iloc[0]['Wait4Scanner.RTTime']
                if Verbose:
                    print( data.iloc[[0,1,-2,-1]][['SessionDate', 'SessionTime', 'Wait4Scanner.RTTime']] )
            else:
                delay = float('nan')
                if Verbose:
                    print( data.iloc[[0,1,-2,-1]][['SessionDate', 'SessionTime']] )
                    print('Unable to find GetReady.RTTime or Wait4Scanner.RTTime in data.columns')

            if delay == delay:   # not nan
                try:
                    if isinstance( delay, str ):
                        delay = int( delay )
                    delay = int( delay / 1000 )   # seconds
                    Info['delay']  = delay
                    Info['exp_t0'] = Info['exp_datime'] + datetime.timedelta(seconds=delay)
                    ok = True
                except:
                    Info['delay']  = float('nan')
                    Info['exp_t0'] = '0000-00-00 00:00:00'
                    ok = False
                    exp_diagnos += -128
                    msg = 'Unable to assess starting time. '
            else:
                Info['delay']  = float('nan')
                Info['exp_t0'] = '0000-00-00 00:00:00'
                ok = False
                exp_diagnos += -128
                msg = 'Unable to extract internal delay. '

        Info['ok']  = ok
        Info['msg'] = Info['msg'] + msg

        if Verbose:
            print( Info, '\n')

    else:
        if Verbose:
            print( Info['msg'], '\n')
            print('fname:', fname )
            print('exp_date_str:', exp_date_str )
            print('exp_time_str:', exp_time_str )
            print()
    # ---------------------------------------------------------------------------------------------------------------


    # ---------------------------------------------------------------------------------------------------------------
    #                                          Encode diagnostics

    # Summarize file-encoding in variable Info['diagnos']
    Info = File_Diagnostics( Info )

    # If exp_diagnos < 0, combine the negative summary of the file-encoding number with exp_diagnos
    if exp_diagnos < 0:
        Info['diagnos'] = exp_diagnos - Info['diagnos']
    # ---------------------------------------------------------------------------------------------------------------

    return Info, data
# ---------------------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------------------
def ExportFile( data, fname_out ):
    ok = False
    msg = ''

    fxpos = fname_out.rfind('.')
    if fxpos > 0:
        if fxpos > (len(fname_out) - 5):   # The rightmost '.' is near the end of filename, and thus consistent with an extension
            fname_bas = fname_out[0:fxpos]
        else:
            fname_bas = fname_out
    else:
        if fxpos < 0:
            fname_bas = fname_out
        else:
            msg = 'Error: invalid output file name'
            return ok, msg

    if len(fname_bas) > 0:
        fname_out = fname_bas + '.txt'

    if Verbose:
        print( '\nWriting tab-separated file:', fname_out, '...' )

    try:
        data.to_csv( fname_out, index=False, sep='\t' )
        ok = True
    except Exception as err:
        msg = 'Unable to write output file: %s: %s' % (fname_out + str(err))

    if Verbose:
        print( msg )

    return ok, msg
# ---------------------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------------------
def List_and_Pick_Files( fname, optn, fname_out, subj, ref_time_arg, task_arg, optn2, quiet=False ):
    # ok  = True
    # msg = ''
    Files = pd.DataFrame()
    EPrime_Info = copy.deepcopy( Info_Prot )
    EPrime_Data = pd.DataFrame()

    f_path = fname
    path_to_search = f_path + '/*'
    f_list = glob.glob( path_to_search )

    # Extract "kernel" of subject ID for comparison with file name
    subj = subj.replace('NDAR_INV','').replace('INV','')
    if Verbose:
        print('subj=', subj)

    if not f_list:
        # ok = False
        # msg = 'Error: argument is not a directory. '
        # EPrime_Info = {'fname': fname,
        #                'ok':  ok,
        #                'msg': msg }
        EPrime_Info['fname'] = fname
        EPrime_Info['ok'] = False
        EPrime_Info['msg'] = 'Error: argument is not a directory. '
        # if Verbose or not optn2:
        if Verbose:
            print('f_list:', f_list )
            print( EPrime_Info['msg'] )

        # return ok, msg, Files, EPrime_Info, EPrime_Data
        return Files, EPrime_Info, EPrime_Data


    for j, fname_full in enumerate( f_list ):

        path  = os.path.dirname( fname_full)
        fname = os.path.basename(fname_full)

        if Verbose:
            print('Reading:', fname_full )

        EPrime_Info, EPrime_Data  =  EPrime_Info_and_Data_get( fname_full )

        if EPrime_Info['fname']:

            # If a subject ID was provided: search for a perfect-match substring in the file name then evaluate general string similarity
            pGUIDmatch = 0.0
            if len(subj)>7:
                subloc = fname.find(subj)
                if subloc >= 0:
                    pGUIDmatch = 10.0
            if len(subj)>3:
                pGUIDmatch += round( SequenceMatcher( a=subj, b=fname ).ratio(), 2 )

            record = pd.DataFrame( dict({'file': fname,
                                        'encoding':    EPrime_Info['encoding'],
                                        'sep_tab':     EPrime_Info['sep_tab'],
                                        'quoted_rows': EPrime_Info['quoted_rows'],
                                        'hoffs':       EPrime_Info['hoffs'],
                                        'n_rows':      EPrime_Info['n_rows'],
                                        'n_cols':      EPrime_Info['n_cols'],
                                        'exper':       EPrime_Info['exper'],
                                        'exp_in_file': EPrime_Info['exp_in_file'],
                                        'exp_datime':  EPrime_Info['exp_datime'],
                                        'delay':       EPrime_Info['delay'],
                                        'exp_t0':      EPrime_Info['exp_t0'],
                                        'diagnos':         EPrime_Info['diagnos'],
                                        'ok':              EPrime_Info['ok'],
                                        'naming_ok':       EPrime_Info['naming_ok'],
                                        'fname_exp_match': EPrime_Info['fname_exp_match'],
                                        'msg':             EPrime_Info['msg'],
                                        'pGUIDmatch':      pGUIDmatch,
                                        'path': path }),  index=[j] )

            Files = Files.append( record, ignore_index=True )
        else:
            pass


    EPrime_Info =  copy.deepcopy( Info_Prot )
    EPrime_Data = pd.DataFrame()

    if len(Files) <= 0:
        # ok = False
        # msg = 'Unable to read or interpret files in directory. '
        # EPrime_Info = {'fname': fname,
        #                'ok' : ok,
        #                'msg': msg }
        EPrime_Info['fname'] = fname
        EPrime_Info['ok'] = False
        EPrime_Info['msg'] = 'Unable to read or interpret files in directory. '
        # if Verbose or not optn2:
        if Verbose:
            print( EPrime_Info['msg'] )

        # return ok, msg, Files, EPrime_Info, EPrime_Data
        return Files, EPrime_Info, EPrime_Data


    # Reorder columns
    Files = Files[['file', 'encoding', 'sep_tab', 'quoted_rows', 'hoffs', 'n_rows', 'n_cols', 'exp_in_file',
                    'exper', 'exp_datime', 'delay', 'exp_t0', 'ok', 'pGUIDmatch', 'diagnos', 'naming_ok', 'fname_exp_match', 'msg', 'path']]

    # Sort by ok, date & time, interpreted or not, prefered format
    Files = Files.sort_values( by=['ok','pGUIDmatch','exp_datime', 'exper', 'diagnos', 'naming_ok', 'fname_exp_match'],
                               ascending=[False, False, True, True, True, False, False] )
    Files = Files.reset_index( drop=True )

    if (Verbose or optn == 'ListFiles') and not quiet:
        print('Found files:')
        print( Files )

    if optn == 'PickFile':

        # Filter file list: valid files and, if specified, task.
        # Find prefered-encoding file for which exp_t0 is closest to given reference date & time
        # Calculate time difference in minutes and order files according to the absolute time time differences

        Files = Files[ Files['ok'] ]

        if len(Files) <= 0:
            # ok = False
            # msg = 'No valid files in directory. '
            EPrime_Info['fname'] = fname
            EPrime_Info['ok'] = False
            EPrime_Info['msg'] = 'No valid files in directory. '
            
            # return ok, msg, Files, EPrime_Info, EPrime_Data
            return Files, EPrime_Info, EPrime_Data

        # if subj:
        #     Files = Files[ Files['pGUIDmatch'] >= 10.0 ]

        # if len(Files) <= 0:
        #     ok = False
        #     msg = 'No files in directory satisfy request. '
        #     return ok, msg, Files, EPrime_Info, EPrime_Data

        if task_arg:
            Files = Files[ Files['exper'] == task_arg ]

        if len(Files) <= 0:
            # ok = False
            # msg = 'No files in directory satisfy request. '
            # return ok, msg, Files, EPrime_Info, EPrime_Data
            EPrime_Info['fname'] = fname
            EPrime_Info['ok'] = False
            EPrime_Info['msg'] = 'No files in directory satisfy request. '

            return Files, EPrime_Info, EPrime_Data


        if 13 <= len(ref_time_arg) and len(ref_time_arg) <= 16:
            ref_time = datetime.datetime.strptime( ref_time_arg, '%Y%m%d %H%M%S')

            if Verbose or not optn2:
                print()
                print('Reference time:', ref_time_arg, ':  ', ref_time )

            Files['tdiff']  = [ round( (s-ref_time).total_seconds()/60, 2 )  for s in Files['exp_t0']]
            Files['tdiffa'] = [ abs(s) for s in Files['tdiff'] ]

        #     Files = Files.sort_values( by=['tdiffa', 'pGUIDmatch', 'exp_t0', 'exper', 'diagnos', 'naming_ok'], ascending=[True, False, True, True, True, False] )
        #     Files = Files.drop('tdiffa', axis=1)
        # else:
        #     Files['tdiff']  = float('nan')
        #     Files = Files.sort_values( by=['pGUIDmatch', 'exp_t0', 'exper', 'diagnos', 'naming_ok'], ascending=[False, True, True, True, False] )

            Files = Files.sort_values( by=['tdiffa', 'pGUIDmatch', 'exp_t0', 'exper', 'diagnos', 'naming_ok', 'fname_exp_match'],
                                       ascending=[True, False, True, True, True, False, False] )
            Files = Files.drop('tdiffa', axis=1)
        else:
            Files['tdiff']  = float('nan')
            Files = Files.sort_values( by=['pGUIDmatch', 'exp_t0', 'exper', 'diagnos', 'naming_ok', 'fname_exp_match'], ascending=[False, True, True, True, False, False] )

        if Verbose or not optn2:
            print('Files')
            print( Files )
            print()
            print('Chosen file:')
            print( Files.iloc[0].to_frame().T )
            print()

        # Select top file in list for further processing

        fname = f_path + '/' + Files.iloc[0]['file']
        pGUIDmatch = Files.iloc[0]['pGUIDmatch']
        tdiff      = Files.iloc[0]['tdiff']

        # if optn2 in ['Info', 'ExportFile']:
        #     EPrime_Info, EPrime_Data  =  EPrime_Info_and_Data_get( fname )
        # else:
        #     EPrime_Info, EPrime_Data  =  EPrime_Info_and_Data_get( fname, 'FileNameCheck' )

        EPrime_Info, EPrime_Data  =  EPrime_Info_and_Data_get( fname )

        EPrime_Info['pGUIDmatch'] = pGUIDmatch
        EPrime_Info['tdiff'] = tdiff

        if 'Info' in optn2:
            sep = ',  '
            print( EPrime_Info['diagnos'], sep,
                    1 if EPrime_Info['pGUIDmatch'] > 10.0 else 0,  sep,
                    EPrime_Info['naming_ok'], sep,
                    EPrime_Info['exp_t0'], sep,  EPrime_Info['exper'],   sep,
                    EPrime_Info['tdiff'],  sep,  EPrime_Info['fname'] )

        if optn2 == 'ExportFile':
            ok, msg  =  ExportFile( EPrime_Data, fname_out )
            EPrime_Info['fout_ok']  = ok
            EPrime_Info['fout_msg'] = msg

    # return ok, msg, Files, EPrime_Info, EPrime_Data
    return Files, EPrime_Info, EPrime_Data
# ---------------------------------------------------------------------------------------------------------------
# =============================================================================================================================================



# =============================================================================================================================================
if __name__ == "__main__":
    code = 0   # will be returned by 'sys.exit(diagnos)', whenever the program ends. Check this exit-status code with 'echo $?'

    EPrime_Info = copy.deepcopy( Info_Prot )

    cmnd_syntx_ok, fname, optn, fname_out, ref_time_arg, subj, task_arg, optn2  =  command_line_get_variables()

    if not cmnd_syntx_ok:
        code = 1
        sys.exit( code )


    if optn == 'Summary':
        Verbose = True
        try:
            EPrime_Info, EPrime_Data  =  EPrime_Info_and_Data_get( fname )

            diagnos = EPrime_Info['diagnos']

        except Exception as err:
            print('Unable to interpret file:', fname )
            sys.exit( code )

        print('exp_datime =', EPrime_Info['exp_datime'], ',  t_del =', EPrime_Info['delay'], 's ,  exp_t0 =', EPrime_Info['exp_t0'], '\n')


    elif optn in ['Info', 'InfoCheckName', 'ExportFile']:
        try:
            # if optn == 'Info':
            #     EPrime_Info, EPrime_Data  =  EPrime_Info_and_Data_get( fname )
            # else:
            #     EPrime_Info, EPrime_Data  =  EPrime_Info_and_Data_get( fname, 'FileNameCheck' )

            EPrime_Info, EPrime_Data  =  EPrime_Info_and_Data_get( fname )

            if EPrime_Info['ok']  and  optn == 'InfoCheckName':
                if not EPrime_Info['fname_exp_match']:
                    # Report mismatch by replacing exp_diagnos with its negative
                    if EPrime_Info['diagnos'] > 0:
                        EPrime_Info['diagnos'] = -64 - EPrime_Info['diagnos']
                    else:
                        EPrime_Info['diagnos'] = -64 + EPrime_Info['diagnos']

        except Exception as err:
            msg = 'Unable to read or interpret file. '
            if Verbose:
                print( msg, fname )

        if 'Info' in optn:
            sep = ',  '
            print( EPrime_Info['diagnos'], sep,
                    1 if EPrime_Info['pGUIDmatch'] > 10.0 else 0,  sep,
                    EPrime_Info['naming_ok'], sep,
                    EPrime_Info['exp_t0'], sep,  EPrime_Info['exper'],   sep,  EPrime_Info['fname'] )

        elif optn == 'ExportFile':
            if EPrime_Info['diagnos'] > 0:
                ok, msg  =  ExportFile( EPrime_Data, fname_out )
            else:
                ok = False
                msg += 'Nothing to export. '
                print( msg )

            EPrime_Info['fout_ok']  = ok
            EPrime_Info['fout_msg'] = msg


    elif optn in ['ListFiles', 'PickFile']:
        try:
            Files, EPrime_Info, EPrime_Contents  =  List_and_Pick_Files( fname, optn, fname_out, subj, ref_time_arg, task_arg, optn2 )

        except Exception as err:
            # ok = False
            # msg = 'Unable to get files from: %s' % fname
            # naming_ok =  0
            EPrime_Info['ok'] = False
            EPrime_Info['msg'] = 'Unable to get files from: %s' % fname
            EPrime_Info['pGUIDmatch'] = 0
            EPrime_Info['naming_ok'] = 0
            EPrime_Info['fname'] = ''

        # if ok:
        if EPrime_Info['ok']:
            if optn == 'ListFiles':
                EPrime_Info['diagnos'] = 0
        else:
            print()
        #     sep = ',  '
        #     print( EPrime_Info['diagnos'], sep,
        #             1 if EPrime_Info['pGUIDmatch'] > 10.0 else 0,  sep,
        #             EPrime_Info['naming_ok'], sep,
        #             EPrime_Info['exp_t0'], sep,  EPrime_Info['exper'],   sep,
        #             EPrime_Info['tdiff'],  sep,  EPrime_Info['fname'] )

        if optn2 == 'ExportFile':
            ok, msg  =  ExportFile( EPrime_Data, fname_out )
            EPrime_Info['fout_ok']  = ok
            EPrime_Info['fout_msg'] = msg
    # ---------------------------------------------------------------------------------------------------

    # Unix exit-status code is an unsigned integer 0..255
    if EPrime_Info['diagnos'] > 0:
        code = 100 + EPrime_Info['diagnos']
    sys.exit( code )
# =============================================================================================================================================


                # 'Nback': ['NBACK','WM','REC'],


        # ok, msg, Files, EPrime_Info, EPrime_Contents  =  List_and_Pick_Files( fname, optn, fname_out, subj, ref_time_arg, task_arg, optn2 )
        # print('msg:', msg )
        # print('Files:')
        # print( Files )
        # print('EPrime_Info:')
        # print( EPrime_Info )

        # # print()
        # # print('Files')
        # # print( Files )
        # print( subj, task_arg, ref_time_arg )
        # # if subj:
        # #     print('subj')
        # # if task_arg:
        # #     print( 'task_arg')
        # # if ref_time_arg:
        # #     print('ref_time_arg')
        # print()


    # exp_datime = dateutil.parser.parse( exp_date_str + ' ' + exp_time_str )

        # print('Error (eprime_sprdsht_get.py): unable to extract or interpret experimemt date or time')
        # print('fname:', fname )
        # print('exp_date_str:', exp_date_str )
        # print('exp_time_str:', exp_time_str )
        # print()


        # print( EPrime_Info['diagnos'], sep, naming_ok, sep, EPrime_Info['exp_t0'], sep, EPrime_Info['exper'], sep, tdiff, sep, EPrime_Info['fname'] )

            # fname   = f_path + '/' + Files.iloc[0]['file']
            # naming_ok = Files.iloc[0]['naming_ok'] 
            # EPrime_Info['tdiff'] = Files.iloc[0]['tdiff'] 

            # if Verbose or not optn2:
            #     print('Getting:', fname )

            # # EPrime_Info, EPrime_Data  =  EPrime_Info_and_Data_get( fname )
            # if optn2 == 'Info':
            #     EPrime_Info, EPrime_Data  =  EPrime_Info_and_Data_get( fname )
            # else:
            #     EPrime_Info, EPrime_Data  =  EPrime_Info_and_Data_get( fname, 'FileNameCheck' )

            # # if optn2 == 'Info':
            # if 'Info' in optn2:
            #     sep = ',  '
            #     print( EPrime_Info['diagnos'], sep,  EPrime_Info['naming_ok],  sep,
            #            EPrime_Info['exp_t0'],  sep,  EPrime_Info['exper'],  sep,
            #            EPrime_Info['tdiff'],   sep,  EPrime_Info['fname'] )




            # fname   = f_path + '/' + Files.iloc[0]['file']

            # if Verbose or not optn2:
            #     print('Getting:', fname )

            # if optn2 == 'Info':
            #     EPrime_Info, EPrime_Data  =  EPrime_Info_and_Data_get( fname )
            # else:
            #     EPrime_Info, EPrime_Data  =  EPrime_Info_and_Data_get( fname, 'FileNameCheck' )

            # naming_ok = Files.iloc[0]['naming_ok'] 
            # tdiff   = Files.iloc[0]['tdiff'] 



            # # if optn2 == 'Info':
            # if 'Info' in optn2:
            #     sep = ',  '
            #     print( EPrime_Info['diagnos'], sep, naming_ok, sep, EPrime_Info['exp_t0'], sep, EPrime_Info['exper'], sep, tdiff, sep, EPrime_Info['fname'] )

            # elif optn2 == 'ExportFile':
            #     ok, msg  =  ExportFile( EPrime_Data, fname_out )
            #     EPrime_Info['fout_ok']  = ok
            #     EPrime_Info['fout_msg'] = msg




                # if len(fname):
                #     # Check if experiment matches file name identifier
                #     if isinstance( exp_fnameID_list, list ):
                #         for exp_fn in exp_fnameID_list:
                #             if fname.lower().rfind(exp_fn.lower()) > (len(fname)-10):   # (look for task id near the end of file name)
                #                 fname_exp_match = 1
                #                 break
                #     else:
                #         exp_fn = exp_fnameID_list
                #         if fname.lower().rfind(exp_fn.lower()) > (len(fname)-10):   # (look for task id near the end of file name)
                #             fname_exp_match = 1

                #     else:
                #         ok = False
                #         msg = 'Experiment in spreadsheet does not match file name. '
                #         exp_diagnos = 0
                # else:
                #     msg = ''
                #     exp_diagnos = 0


        # EPrime_Info = copy.deepcopy( Info_Prot )
        # EPrime_Data = pd.DataFrame()
        # record      = pd.DataFrame()

        # path  = os.path.dirname( fname_full)
        # fname = os.path.basename(fname_full)

        # if Verbose:
        #     print('Reading:', fname_full )

        # # try:
        # #     EPrime_Info, EPrime_Data  =  EPrime_Info_and_Data_get( fname_full )
        # # except:
        # #     pass

        # # if len(EPrime_Info) > 0  and  len(EPrime_Data) > 0  and  EPrime_Info['ok']:
        # # if len(EPrime_Info) > 0:
        # # if EPrime_Info['fname']:   # we got info about this file, but we don't know yet if is a valid exp. spreadsheet

    #             Additionaly, if optn = 'FileNameCheck', check that file contents correspond to file's name

    # if optn == 'FileNameCheck':
    #     exper, fname_exp_match, ok, msg, exp_diagnos  =  ExperimentCheck( data, fname )
    # else:
    #     exper, fname_exp_match, ok, msg, exp_diagnos  =  ExperimentCheck( data, '' )
