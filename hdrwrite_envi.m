function hdrwrite_envi(filename, path, row, col)
filestring = strcat(path,'\',filename,'.hdr');
fileID = fopen(char(filestring), 'w');
fprintf(fileID, 'ENVI\n');
fprintf(fileID, 'description = {\n');
fprintf(fileID, 'PolSARpro File Imported to ENVI}\n');
samplestring = strcat('samples =  ', num2str(col),'\n');
fprintf(fileID, samplestring);
linesstring = strcat('lines   =  ', num2str(row),'\n');
fprintf(fileID, linesstring);
fprintf(fileID, 'bands   = 1\n');
fprintf(fileID, 'header offset = 0\n');
fprintf(fileID, 'file type = ENVI Standard\n');
fprintf(fileID, 'data type = 4\n');
fprintf(fileID, 'interleave = bsq\n');
fprintf(fileID, 'sensor type = Unknown\n');
fprintf(fileID, 'byte order = 0\n');
fprintf(fileID, 'band names = {\n');
bandstring = strcat(filename,'.bin',' }');
fprintf(fileID, char(bandstring));
fclose(fileID);
end