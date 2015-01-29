bc_table.list
cmd
Av  0.100
Rv  3.100
-------------------------------
aaaff.fff = format
first line gives the filename where the BC tables are specified.  See existing example.

second line is the suffix which will be added to the isochrone filename to create the 
new isochrone file with magnitudes included.  basic format of the new file is similar 
to the original isochrone file except that magnitudes are written at lower resolution.

third and fourth lines list a value of Av and Rv to be used in the interpolation. 
format is (a3,f6.3).