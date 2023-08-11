myFolder = 'C:\Users\Higley402\Desktop\DATA\raw_data\OT calculated';
if ~isdir(myFolder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
  uiwait(warndlg(errorMessage));
  return;
end
filePattern = fullfile(myFolder, '*.ibw');
IgorFiles = dir(filePattern);
for k = 1:length(IgorFiles)
  baseFileName = IgorFiles(k).name;
  fullFileName = fullfile(myFolder, baseFileName);
  [path, name, ext] = fileparts(fullFileName);
  % fprintf(1, 'Now reading %s\n', fullFileName);
  OT.(name)=IBWread(fullFileName);
  % Now do whatever you need to do with the structure you recalled.
end