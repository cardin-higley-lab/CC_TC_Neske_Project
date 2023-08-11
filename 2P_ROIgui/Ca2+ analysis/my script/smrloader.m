function dataout = smrloader( datain )
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here

fid=fopen(datain);
SONImport(fid);
[path, name, ext]=fileparts(datain);
dataout = load(name)


end

