function map_format = map_creator(filename, nCh, precision)
% From file create a format map that memmapfile can use to read the file
%
% INPUTS:
%   filename  - scalar file path, what file to read
%   precision - text scalar, what data dercision.
%               one of:
%               'int8','int16','int32','int64','uint8',
%               'uint16','uint32','uint64','single','double'
%   nCh       - positive int, number of channels in data.
%
% OUTPUT:
%   map_format - cell array, 1 X 3.
%                Formats for memmapfile.
%
% Dependencies:
%       slutsky_ECInVivo\utilities\class2bytes.m
%   See Also memmapfile

arguments
    filename  (1,:) string {mustBeTextScalar,mustBeFile}
    nCh       (1,1) double {mustBePositive,mustBeInteger}
    precision (1,:) string {mustBeTextScalar, mustBeMember(precision,...
        ["int8","int16","int32","int64","uint8",...
        "uint16","uint32","uint64","single","double"])} ...
        = 'int16'
end

% convert precision to bytes
nBytes = class2bytes(precision);

% get number of samples per channel
file_info = dir(filename);
nSamps = file_info.bytes / nBytes / nCh;

% merge all for a format
map_format = {precision, [nCh, nSamps]};
