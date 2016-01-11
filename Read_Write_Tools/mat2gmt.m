function [] = mat2gmt(fld, fname);

tt = '%14.6f ';
[m,n] = size(fld);
frmtwrt = [repmat(tt,1,n), '\n'];

fid = fopen(fname,'w+');
fprintf(fid, 'ncols            %g',n);
fprintf(fid, '\n');
fprintf(fid, 'nrows            %g',m);
fprintf(fid, '\n');
fprintf(fid, 'xllcorner        %g',0);
fprintf(fid, '\n');
fprintf(fid, 'yllcorner        %g',0);
fprintf(fid, '\n');
fprintf(fid, 'cellsize         %2.4f', 1/n);
fprintf(fid, '\n');
fprintf(fid, 'NODATA_value   -9999');
fprintf(fid, '\n');
fprintf(fid, frmtwrt, (flipud(fld))');
fclose(fid);