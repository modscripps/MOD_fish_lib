function matchingFiles = GV_filterFilesByTimestamp(theFiles, inputTimestamp)
% FILTERFILESBYTIMESTAMP Filters files based on a timestamp.
%
%   matchingFiles = FILTERFILESBYTIMESTAMP(inputTimestamp)
%
%   This function looks for files matching the pattern
%   'FCTDYY_MM_DD_hhmmss.modraw' in the current directory and returns a
%   cell array of filenames that have a timestamp *after* the provided
%   inputTimestamp.
%
% Inputs:
%   inputTimestamp - A string in the format 'YYYYMMDD_hhmmss' representing
%                    the cutoff time.
%
% Outputs:
%   matchingFiles - A cell array of strings containing the names of the
%                   files whose timestamps are after the inputTimestamp.

    % --- 1. Define the File Pattern and Find All Matching Files ---
    % The pattern to search for is 'FCTDYY_MM_DD_hhmmss.modraw'
    % filePattern = 'FCTD*_*.modraw';
    % theFiles = dir(filePattern);

    % Pre-allocate a cell array to store the names of the matching files
    matchingFiles = {};
    matchIndex = 1;

    % --- 2. Convert Input Timestamp String to a Datetime Object ---
    % The input is YYYYMMDD_hhmmss (e.g., 20251211_002841)
    try
        % Specify the format to correctly parse the input string
        cutOffTime = datetime(inputTimestamp, 'InputFormat', 'yyyyMMdd_HHmmss');
    catch ME
        error('filterFilesByTimestamp:InvalidInput', ...
              'Input timestamp must be a string in YYYYMMDD_hhmmss format.');
    end


    % --- 3. Iterate Through Files and Compare Timestamps ---
    for k = 1:length(theFiles)
        currentFileName = theFiles{k};

        % --- Extract the Timestamp Substring from the Filename ---
        % The target timestamp is always located from index 5 to index 22.
        % Filename: FCTD25_12_11_002841.modraw
        % Index:    1234567890123456789012
        % Target:     ^-------------^
        timestampString = currentFileName(5:19); % Gets '25_12_11_002841'

        % --- Convert Filename Timestamp to a Datetime Object ---
        % The format is YY_MM_DD_hhmmss (e.g., 25_12_11_002841)
        try
            fileTime = datetime(timestampString, 'InputFormat', 'yy_MM_dd_HHmmss');
        catch ME
            % This step should generally not fail if the file pattern is correct,
            % but it's good practice to handle potential parsing issues.
            warning('filterFilesByTimestamp:ParseError', ...
                    'Could not parse timestamp from file: %s', currentFileName);
            continue;
        end

        % --- Comparison Logic ---
        % We are checking if the fileTime is STRICTLY AFTER the cutOffTime.
        if fileTime > cutOffTime
            matchingFiles{matchIndex} = currentFileName;
            matchIndex = matchIndex + 1;
        end
    end
end