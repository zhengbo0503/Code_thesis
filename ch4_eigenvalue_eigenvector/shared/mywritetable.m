function mywritetable( data, cellName, outputPoisition)
%MYWRITETABLE write a floating point data to table
%	mywritetable( data, cellName, outputPoisition) output
%	the matrix data to a table whose title will be
%	the elements in cellName.
%	outputPosition specify the saved location
%   
%   data : should aligned in columns
%   cellName: A vector of strings contains the names 
%       ["name1","name2",...,"namen"]
%   outputPosition: rel or abs position are all acceptable.
%
%   This is a thin wrapper kept so that existing call sites do not change.
%   It delegates to write_figure_csv, which writes 17 significant digits
%   (round-tripping IEEE double instead of writetable's 15) and additionally
%   publishes the file to the thesis figure directory.

write_figure_csv(outputPoisition, double(data), cellName);

end
