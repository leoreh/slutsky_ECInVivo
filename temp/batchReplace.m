% batchReplace

directory = 'D:\Code\slutsky_ECInVivo';

targetWord = 'stateBouts';
replaceWord = 'boutTimes';

files = dir(fullfile(directory, '**', '*.m'));

for k = 1:length(files)
    filePath = fullfile(files(k).folder, files(k).name);

    fid = fopen(filePath, 'r');
    fileContent = fread(fid, '*char')';
    fclose(fid);

    [updatedContent] = regexprep(fileContent, targetWord, replaceWord);

    fid = fopen(filePath, 'w');
    fwrite(fid, updatedContent);
    fclose(fid);
end

