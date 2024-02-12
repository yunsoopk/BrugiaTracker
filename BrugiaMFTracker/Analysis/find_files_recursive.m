function [ matches ] = find_files_recursive( folder, pattern )
%FIND_FILES_RECURSIVE will return all the files in the specified folder and
%subfolders that match the given pattern
matches = [];
files = dir(folder);
%eliminate . and ..
files(1:2) = [];
for i = 1:length(files)
    if(files(i).isdir)
        temp = find_files_recursive( fullfile(folder,files(i).name), pattern );
        matches = [matches;temp];
    else
        startIndex = regexp( files(i).name, pattern, 'once' );
        if(~isempty(startIndex))
            matches = [matches; {fullfile(folder,files(i).name)}];
        end
    end
end

end

