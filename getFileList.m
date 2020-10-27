function [filelist,names,paths,exts] = getFileList(folder, extension, excludehidden)
%GETFILELIST Returns a list of filenames in the given folder, excluding
%folders and the  "." and ".." listings that dir() returns. Also returns
%filenames, paths, and extensions separately.
%
% Optional arguments:
%   extension: return only list of files matching given extension. If an
%   empty string is passed, getFileList will not match the extension
%   excludehidden: exclude hidden files (true/false). If a list of all
%   files including hidden files is required, use:
%   getFileList(folder,'','true')

matchextension = true;

switch nargin
    case 1
        extension = "";
        excludehidden = true;
    case 2
        excludehidden = true;
end

%if empty string is given as extension, don't match extension
if extension == ""
    matchextension = false;
end

%add . to extension if missing
if ~startsWith(extension,".")
    extension = strcat(".",extension);
end

%get file list
fileliststruct = dir(folder);
fileliststruct = fileliststruct(~ismember({fileliststruct.name},{'.','..','.DS_Store'}));

filelist = [];
fileNum = 1;
for i=1:length(fileliststruct)
    if ~fileliststruct(i).isdir
        filelist = [filelist; string(fileliststruct(i).name)];
        [paths(fileNum),names(fileNum),exts(fileNum)] = fileparts(strcat(folder,'/',filelist(fileNum)));
        extmatch(fileNum) = exts(fileNum) == extension;
        fileNum = fileNum + 1;
    end
end

if isempty(filelist)
    exc = MException('getFileList:noFilesFound', ...
        strcat("No files found in folder ", folder," Check file path for typos."));
    throw(exc)
end

%remove hidden files if required
if excludehidden
    ind = startsWith(names,".");
    
    paths(ind) = [];
    names(ind) = [];
    exts(ind) = [];
end

%extract only those with the correct extension
if matchextension
    paths(~extmatch) = [];
    names(~extmatch) = [];
    exts(~extmatch) = [];
end
filelist = strcat(names,exts);

end

