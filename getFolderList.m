function [filelist] = getFolderList(folder)
%GETFILELIST Returns a list of folders in the given folder, excluding
%files and the  "." and ".." listings that dir() returns
fileliststruct = dir(folder);
fileliststruct = fileliststruct(~ismember({fileliststruct.name},{'.','..','.DS_Store'}));

filelist = [];
for i=1:length(fileliststruct)
    if fileliststruct(i).isdir
        filelist = [filelist; string(fileliststruct(i).name)];
    end
end


end

