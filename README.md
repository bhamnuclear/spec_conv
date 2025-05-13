# spec_conv.c

**Program `spec_conv`** converts data between common spectrum formats in
used nuclear physics.

Large numbers of spectra can be converted via the list file functionality.
This is can also be used for gain-matching many RadWare-format spectra
by including gain coefficients after each filename in the list file.

Program `spec_conv` accepts the following formats:

- ASCII: 1 (y) and 2 (x,y) column formats;
- RadWare .spe binary format;
- ORTEC (MAESTRO) .Spe (1-column ASCII) and .Chn (binary) formats with headers
and trailers.
- Xtrack .spec format (e.g. AGATA spectra);
- GENIE .IEC spectra.

For ASCII spectra and list files, lines starting with # are
treated as comments and skipped.

The full list of `spec_conv` options is:

1. to convert RadWare (.spe) ==> Ascii (.txt)
2. to convert Ascii (.txt) ==> RadWare (.spe)
3. to convert Ascii (.txt) ==> Xtrack (.spec)
4. to convert Maestro_Chn (.Chn) ==> Ascii (.txt)
5. to convert Maestro_Chn (.Chn) ==> RadWare (.spe)
6. to convert Xtrack (.spec) ==> Ascii (.txt)
7. to convert Xtrack (.spec) ==> RadWare (.spe)
8. to convert GENIE (.IEC) ==> RadWare (.spe)
9. to convert Maestro_Spe (.Spe) ==> RadWare (.spe)
a. to convert Maestro_Spe (.Spe) ==> Ascii (.txt)
g. to gainmatch a RadWare spectrum
0. Quit
