#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <ctype.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <termios.h>

#define CHMAX 	  32768	/*max number of channels in spectra*/
#define MAXCOLS   3     /*max. number data columns in input spectrum*/
#define MXNUMDIG  3     /*max number of digits for the number of multi spectra
                            that can be extracted. i.e. 3 ==> 999 spectra*/
#define NUMOPT    11    /*number of options*/
#define CHLEN     120   /*character length of filename arrays*/
                            
/* Carl Wheldon May 2003 */
/* Lastest up-date May 2025*/

/* To compile:
gcc spec_conv.c -Wall -pedantic -o spec_conv -lm -O2
*/

/*%%%%% A program to convert between different spectra formats %%%%%*/
/* To add a new option:
    add input extension and format to exti[] and fmti[],
    add output extension and format to ext[] and fmt[],
    edit get_mode(), read_spec() and write_spec(),
    and add md (i.e. mode number) to relevant point in main().*/

/*structure of the radware header as written to/read from a spectrum*/
struct radheader {
    unsigned int q1;
    char    	 name[8];
    unsigned int channels;
    unsigned int q2;
    unsigned int q3;
    unsigned int q4;
    unsigned int q5;
    unsigned int size;
} radheader;

/*structure of the radware trailer*/
struct radtrailer {
    unsigned int size;
} radtrailer;

/*structure of the Ortec Maestro header as written to/read from a spectrum*/
struct  maest_header {
    short int q1;           /*must be one*/
    short int q2;           /*MCA/det number*/
    short int q3;           /*segment number*/
    short int q4;           /*ascii second of start time*/
    unsigned int real;      /*real time (increments of 20 ms)*/
    unsigned int lve;       /*live time (increments of 20 ms)*/
    char dt[8];             /*Start date as ASCII DDMMMYY*
                              The * character should be ignored if it is not
                              a "1". If it is a "1", it indicates the data is
                              after the year 2000.*/
    char sttm[4];           /*start time as ASCII HHMM or binary zeros,
                                if not known*/
    short int off;          /*channel offset of data*/
    short int channels;     /*length of data (channels)*/
} maest_header;

/*structure of the Ortec Maestro header as written to/read from a spectrum*/
struct  maest_trailer {
    short int t1;           /*must be 102*/
    short int t2;           /*reserved???*/
    float     g[3];         /*energy calib coeff offset, gain and quadratic term*/  
    char      trailer[496]; /*nothing particularly useful in the rest of the trailer*/
} maest_trailer;

int 	ascii_read(char name[]);
void 	ascii_write(char name[], int numch);
void 	chan_num_ext(char fin[], char fout[], int *numch, char ext[]);
void	check_ext(char fin[], char ext[]);
int  	col_determ(FILE *file);
long 	convert_bytes(char name[]);
int 	cswap4(int decim);
int 	cswap2(int decim);
void	decode_mspec_name(char name[], int *set, int *mxsp, int *numch,
	    int *sz, int bytes);
void 	file_status(char name[], char ext[], int len);
int 	genie_read(char name[]);
void 	get_ans(char ans[], int num);
void 	get_line(char ans[], int len);
void    get_line_file(FILE *file, char ans[], int len);
int 	get_mode(int md);
void 	get_pars(float pars[], int num);
void    get_pars_file(FILE *file, float pars[], int num);
void 	get_val(float *val);
void 	itoa(int n, char s[]);
int 	maestro_read(char name[]);
void 	num_fname(char name[], int num);
int 	rad_read(char name[]);
void 	rad_write(char name[], int numch);
int 	read_lst(char inname[], int lst);
int     read_spec(char name[], int md);
void 	reverse(char s[]);
void 	set_ext(char name[], char ext[]);
void 	skip_hash(FILE *file);
void    skip_lines(FILE *file, int lns);
void    store_colours();
void 	swapb2(char *buf);
void 	swapb4(char *buf);
void    write_spec(char name[], int numch, int md);
void 	xtrack_read(char name[], int *numch, int mxsp, int sz, int nsp,
	    int flg);
void 	xtrack_write(char name[], int numch);
        
float spectrum[CHMAX], gain[3];
int md;
char ext[NUMOPT][11], exti[NUMOPT][11], fmti[NUMOPT][14], fmt[NUMOPT][14];
char clr[10][12];
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* ++++++++++++++++++++++++++++++++++ MAIN ++++++++++++++++++++++++++++++++++ */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
int main(int argc, char *argv[])
{
    extern float spectrum[CHMAX], gain[3];
    extern int	md;
    extern char ext[NUMOPT][11], exti[NUMOPT][11]; 
    extern char fmti[NUMOPT][14], fmt[NUMOPT][14];
    float   spbuf[CHMAX], res = 0.0, cal_chan = 0.0, calib = 2.000, tmpf = -1.0;
    long    bytes = 0;
    int     flg = 1, fn = 0, i = 0, j = 0, lst = 3, mxsp = 0, nsp = -1;
    int     numch = CHMAX, set = 1, sz = 4;
    char    inname[CHLEN] = "", outname[CHLEN] = "", ans[CHLEN] = "";
    struct  stat statbuf;
    FILE    *fl;
    
    /*store colours*/
    store_colours();
    
    /*fill input extension arrays and input format arrays.
        Note strlen+1 is used to terminate strings*/
    strncpy(exti[0], ".spe", 5);         strncpy(fmti[0], "RadWare", 8);
    strncpy(exti[1], ".txt", 5);         strncpy(fmti[1], "Ascii", 6);
    strncpy(exti[2], ".txt", 5);         strncpy(fmti[2], "Ascii", 6);
    strncpy(exti[3], ".Chn", 5);         strncpy(fmti[3], "Maestro_Chn", 12);
    strncpy(exti[4], ".Chn", 5);         strncpy(fmti[4], "Maestro_Chn", 12);
    strncpy(exti[5], ".spec", 6);        strncpy(fmti[5], "Xtrack", 7);
    strncpy(exti[6], ".spec", 6);        strncpy(fmti[6], "Xtrack", 7);
    strncpy(exti[7], ".IEC", 5);         strncpy(fmti[7], "GENIE", 6);
    strncpy(exti[8], ".Spe", 5);         strncpy(fmti[8], "Maestro_Spe", 12);
    strncpy(exti[9], ".Spe", 5);         strncpy(fmti[9], "Maestro_Spe", 12);
    strncpy(exti[NUMOPT-1], ".spe", 5);  strncpy(fmti[NUMOPT-1], "RadWare", 8);
    
    /*fill output extension arrays*/
    strncpy(ext[0], ".txt", 5);                 strncpy(fmt[0], "Ascii", 6);
    strncpy(ext[1], ".spe", 5);                 strncpy(fmt[1], "RadWare", 8);
    strncpy(ext[2], ".spec", 6);                strncpy(fmt[2], "Xtrack", 7);
    strncpy(ext[3], ".txt", 5);                 strncpy(fmt[3], "Ascii", 6);
    strncpy(ext[4], ".spe", 5);                 strncpy(fmt[4], "RadWare", 8);
    strncpy(ext[5], ".txt", 5);                 strncpy(fmt[5], "Ascii", 6);
    strncpy(ext[6], ".spe", 5);                 strncpy(fmt[6], "RadWare", 8);
    strncpy(ext[7], ".spe", 5);                 strncpy(fmt[7], "RadWare", 8);
    strncpy(ext[8], ".spe", 5);                 strncpy(fmt[8], "RadWare", 8);
    strncpy(ext[9], ".txt", 5);                 strncpy(fmt[9], "Ascii", 6);
    strncpy(ext[NUMOPT-1], "_mtchd.spe", 11);   strncpy(fmt[NUMOPT-1], "RadWare", 8);
    
    
    md = 0;
    for (i = 0; i < 3; i++) gain[i] = 0.0;
    i = (int)pow(10,MXNUMDIG) - 1;
    
    printf("\n \t \t     *****Welcome to SPEC_CONV*****\n"
	    "\tThis program converts spectra between RadWare, Ascii,\n"
	    "\tXtrack (GASPWARE) and Ortec (binary Chn & ASCII Spe) formats,\n"
            "\tincluding multiple-spectra (<%d) Xtrack files, e.g. from AGATA.\n"
	    "\tand can gainmatch spectra.\n"
	    "\t(Ascii means (y) or (x y) data starting from channel zero)\n"
	    "\tComment lines starting with # are ignored at the front of\n"
	    "\tascii spectra. The 1 or 2 col. format is auto-detected.\n\n",i); 

/*    printf("    Input ext   Output ext  Format\n");
    for (i = 0; i < NUMOPT; i++)
        printf("    %-11s %-11s %-14s\n",exti[i], ext[i], fmti[i]);
    printf("\n");*/
       
    /*argv[i] is the ith argument, i.e. first is the program name*/
    if ( argc > 2 || (argc == 2 && (stat(argv[1], &statbuf))) )
    {
	printf("\nUnrecognised arguments...usage: spec_conv\n"
		" or: spec_conv SpectrumFileName\n");
	if (argc == 2) printf(" ***File %s does not exist\n",argv[1]);
	return -1;
    }
    else if (argc == 2)
    {
        strcpy(inname,argv[1]);
	/*printf("Filename = %s\n",inname);*/
	lst = 2;
        /*read first non-comment line of file...
            if it's the name of an existing file assume a list file*/
    	if ((fl = fopen(inname, "r" )) == NULL)
    	{
    	    printf("Cannot open file: %s \n", inname);
    	    return -1;			
    	}
        skip_hash(fl);
        if ( (i = fscanf(fl, "%s", ans)) == 1)
        {
           if ( ! stat(ans, &statbuf))
           {
              /*file exists...assuming list file*/
              lst = -1;
              fclose(fl);
              printf("List filename = %s\n",inname);
           }
           else printf("Spectrum filename = %s\n",inname);
        }
    }
    
    if ( (md = get_mode(md)) == 0) return 0;
    
    while (lst == 3)
    {        
	printf("Read spectrum names from list file (y/n) \n");
	get_ans(ans,1);
	if (ans[0] == 'y' || ans[0] == 'Y') lst = 1;
	else if (ans[0] == 'n' || ans[0] == 'N') lst = 0;
    }    

    /* simple spectrum read/write */
    if ( (md >= 1 && md <= 5) || md == 8 || md == 9 || md == 10)
    {
	while (flg == 1)
	{
	    numch = CHMAX;
	    /*if a list file*/
	    if (lst == 1 || lst == -1) fn = read_lst(inname,lst);
	    else if (lst == 0)
	    {
    	    	printf("Type %s filename inc. extension (eg %s):\n",
                        fmti[md-1],exti[md-1]);
	    	get_line(inname, CHLEN);
	    }
	    if (lst != 1 && lst != -1) flg = -1;
            
	    if (fn == -1) return 0;
            
	    /*zero spectrum array*/    	
    	    for (i = 0; i < CHMAX; i++) spectrum[i] = 0.0;
    	    	
    	    /*read spectrum file*/
	    if ( (numch = read_spec(inname, md)) < 0)
            {
		printf("Error, no. channels:%d ...Exiting\n", numch);
                return -1;
            }
	    /* if channels is not multiple of 4096 or <1024, ask for length*/
    	    chan_num_ext(inname, outname, &numch, ext[md-1]);

	    strcpy(outname,inname);
    	    set_ext(outname, ext[md-1]);
    	    /*check file status*/
    	    file_status(outname, ext[md-1], CHLEN);
    	                
    	    printf(" %s", inname);
	    write_spec(outname, numch, md);
	}
    }/*END simple spectrum read/write */        
            
    /*Xtrackn format options*/    
    if ( md == 6 || md == 7 )
    {
	while (flg == 1)
	{
	    /*reset some parameters*/
	    if (mxsp < 1)
	    {
		numch = CHMAX;
		set = 1;
	    }
	    /*if a list file*/
	    if ( (lst == 1 || lst == -1) && mxsp <= 1) fn = read_lst(inname,lst);
	    else if (lst == 0 && nsp == -1)
	    {
    	    	printf("Type %s filename inc. extension (eg %s):\n",
                        fmti[md-1],exti[md-1]);
	    	get_line(inname, CHLEN);
	    }
	    if (fn == -1) return 0;
	    
	    strcpy(outname, inname);
	    if (lst != 1) flg = -1;
	    
	    /*modes 6 and 7 allow for extraction/conversion of
	    	multiple CHMAX channel spectrum in 1 file*/
	    check_ext(inname, exti[md-1]);
    	    /*get and print file size*/
     	    if (mxsp <= 1) bytes = convert_bytes(inname);
	
	    /*check file name for "__" surrounding mult. spec info*/
            if ( mxsp < 1 && strchr(inname,'_') && ! strncmp( strchr(inname,'_'), "__", 2 ) )
	        decode_mspec_name(inname, &set, &mxsp, &numch, &sz, bytes);

	    /*for a standard spectrum check numch is compatible with
	        the filesize*/
    	    while ( mxsp < 1 && 
	        (mxsp = (int)( bytes/(numch*sizeof(unsigned int)) )) < 1
	            && numch > 50) numch /= 2;
    
    	    if ( mxsp > 1 && nsp == -1)
    	    {
	        while (1)
	        {
	            printf("Enter spectrum number you require"
	                " (0-%d, %d for all)\n",mxsp-1,mxsp);
	            get_val(&tmpf);
	            nsp = (int)tmpf;
	            j = 0;
	            if (nsp >= 0 && nsp <= mxsp) break;
	        }
    	    }
    	    else if (mxsp == 1)
    	    {
	        nsp = 0;
	        j = nsp;
    	    }
	
    	    /*set output filename*/
/*    	    if (mxsp > 1 && nsp == mxsp)*/
    	    if (mxsp > 1)
	    {
	        if (nsp != mxsp) j = nsp;
	        flg = 1;
	        strcpy(outname,inname);
                num_fname(outname, j);
	    }
	    
	    /*zero spectrum array*/    	
    	    for (i = 0; i < CHMAX; i++) spectrum[i] = 0.0;
            
    	    if (numch > 50) xtrack_read(inname, &numch, mxsp, sz, j, flg);
            
	    j++;
	    /*reset parameters if last spectrum read*/
	    if ( ( (mxsp > 1) && ((nsp == mxsp && j == mxsp)
	            || (nsp != mxsp))) || mxsp == 1)
	    {
	        mxsp = 0;
	        nsp = -1;
	        /*only break while loop if no more spectra to process*/
	        if (lst != 1) flg = -1;
	    }
	    
	    if (numch <= 0 || numch > CHMAX)
	    {
		printf("Error, no. channels: %d ...Exiting\n", numch);
		return -1;
	    }
	    	    
    	    set_ext(outname, ext[md-1]);
    	    /*check file status*/
    	    file_status(outname, ext[md-1], CHLEN);
    	    
    	    printf(" %s", inname);
	    write_spec(outname, numch, md);
    	}
    } /*END Xtrackn ==> format options*/ 
        
    /* Read in RadWare spectrum, GAINMATCH, and output as RadWare spectrum */
    if (md == NUMOPT)
    {
        printf("md = %d\n",md);
	while (flg == 1)
	{
	    numch = CHMAX;
	    /*if a list file i.e. extension not .spe*/
            /*strcmp returns zero if identical, i.e. if not equal to NULL*/
            if (lst == 2 && strcmp( (strrchr(inname,'.')), exti[NUMOPT-1] ) )
                  lst = -1;
	    
            if (lst == 1 || lst == -1) fn = read_lst(inname,lst);
	    else if (lst == 0 || lst == 2)
	    {
		flg = -1;
		if (lst == 0)
		{
    	    	    printf("Type %s filename inc. extension (eg %s):\n",
                        fmti[md-1],exti[md-1]);
	    	    get_line(inname, CHLEN);
		}
	    	printf("Enter up to 3 gainmatching coeffs. (A0 A1 A2):\n");
		get_pars(gain,3);
	    	printf("A0 = %e, A1 = %e, A2 = %e\n",
			gain[0], gain[1], gain[2]);
	    }
	    if (fn == -1) return 0;
	    
	    if (fn == 1 || lst != 1)
	    {
		if (gain[1] != 0.0) calib = fabs(1.0/gain[1]);
		else calib = 1.0;
	    	printf("You can choose to multiply the coeffs by a"
	    	    " constant factor (e.g. |1/%.2f| = %.2f).\n"
                    "Then your keV/channel for the matched spectra will be"
                    " 1/(factor).\n"
                    "E.g. for a factor of 2.0 the matched spectra will have 0.5 keV/channel\n"
		    "Enter value for factor [<Enter> for 1.0]\n",gain[1],calib);
	    	get_val(&calib);
	    	if (calib <= 0.0)
	    	{
	    	    printf("Mult. factor = 0.0 ==> reset to 1.0\n");
	    	    calib = 1.0;
	    	}
	    }
	    /*zero spectra arrays*/	    	
    	    for (i = 0; i < CHMAX; i++)
	    {
		spectrum[i] = 0.0;
		spbuf[i] = 0.0;
	    }
    	    					
	    gain[0] = calib*gain[0];
	    gain[1] = calib*gain[1];
	    gain[2] = calib*gain[2];
	    if (calib != 1.0) printf("A0 = %e, A1 = %e, A2 = %e\n",
			gain[0], gain[1], gain[2]);
	    
    	    /*read spectrum file*/
	    if ( (numch = read_spec(inname, md)) < 0)
            {
		printf("Error, no. channels:%d ...Exiting\n", numch);
                return -1;
            }
	    
            /*fill spectrum array*/
	    for (j = 0; j < numch; j++)
    	    {
	    	cal_chan = gain[0] + j*gain[1] + j*j*gain[2];
		/*do nothing with counts in channels out of range*/
	    	if ( ( (int)cal_chan ) < 0 || ( (int)cal_chan ) >= numch) ;
	    	else
	    	{
    	    	    res =  cal_chan - (float)( (int)cal_chan );
		    spbuf[(int)cal_chan] += spectrum[j]*(1.0-res);
		    spbuf[((int)cal_chan) + 1] += spectrum[j]*res;	    
	    	}
		spectrum[j] = 0.0;
    	    }
	    
	    /*round counts in spectrum array*/
    	    for (j = 0; j < numch; j++)   	    	
 	    	spectrum[j] = (float)( (int)(spbuf[j] + 0.5) );

	    strcpy(outname,inname);
    	    set_ext(outname, ext[md-1]);
    	    /*check file status*/
    	    file_status(outname, ext[md-1], CHLEN);
    	    
    	    printf(" %s", inname);
	    write_spec(outname, numch, md);
    	}
    } /* END Gainmatch RadWare spectrum and output as RadWare */
                     
    return 0;
} /*END main()*/

/*==========================================================================*/
/* ascii_read: read an ASCII format spectrum	    	    	    	    */
/****************************************************************************/
int ascii_read(char name[])
{
    float   rd = 0, rlt[2];
    int     ascii = 0, chan = 0, res = 0, lchan = CHMAX;
    char    ans[CHLEN] = "";
    FILE    *fsp;
     
    /*opens ascii file*/
    if ((fsp = fopen(name, "r" )) == NULL)
    {
    	printf("Cannot open file: %s \n", name);
    	return -1;
    }
    
    /*for Maestro ASCII format spectrum decode header*/
    if (md == 9 || md == 10)
    {
        /*skip the first 7 header lines*/
        skip_lines(fsp, 7);
        /*get and print date*/
        get_line_file(fsp, ans, CHLEN);
        printf("%sSpectrum date: %s%s\n",clr[2],ans,clr[0]);
        
        /*get and print real/live time*/
        skip_lines(fsp, 1);
        rlt[0] = 0.0; rlt[1] = 0.0;
        get_pars_file(fsp, rlt, 2);
        printf("    %sLive time: %.2f s, Real time %.2f s%s\n",
            clr[2],rlt[0],rlt[1],clr[0]);
        
        /*get channel number*/
        skip_lines(fsp, 1);
        rlt[0] = 0.0; rlt[1] = 0.0;
        get_pars_file(fsp, rlt, 2);
        printf("    %sFirst chan: %d, last chan %d%s\n",
            clr[2],(int)rlt[0],(int)rlt[1],clr[0]);
    }
    
    /*Determine if spectrum is 1 or 2 col. ascii format & set 1 col flag*/
    if ( (ascii = col_determ(fsp)) == 0)
    {
	printf("%s***No suitable data in file %s; Exiting...%s\n\n",
                clr[1],name,clr[0]);
	return -1;
    }  
    printf("Ascii %d column format....", ascii);    
    /*End of deciding if spectrum is 1 or 2 column ascii format*/
    
    for (chan = 0; (((md != 9 && md != 10) && chan < CHMAX) || ((md == 9 || md == 10) && chan < (rlt[1]+1))); chan++)
    { 
/*    	res = fscanf(fsp, "%f %f", &rd, &spectrum[chan]); */
	if (ascii == 2)    /*only for two column data*/
	{
    	    res = fscanf(fsp, "%f", &rd);
	    chan = (int) rd;
	}
    	res = fscanf(fsp, "%f\n", &spectrum[chan]);
	switch (res)
	{
	    case 0:
	    {
	    	printf("\nchan= %d, spec[chan]= %f\n", chan,
		    	spectrum[chan]);
    	    	printf("Read error occurred for file: %s\n", name);
                /*close ascii file*/
	    	fclose(fsp);
	    	return -1;
	    }
	    case EOF:
	    {
	    	if (chan < CHMAX)
		{
		    /*increment chan because chan starts from zero. This is
		    because of line: chan = (int) rd;*/
		    if (ascii == 2) chan += 1;
	    	    printf("Reached EOF after reading chan %d \n", chan);
		    if (chan == 0)
		    {
			printf("\n*******Incorrect file format*******\n");
			printf("....Exiting....\n\n");
			return -1;
		    }
		}
		lchan = chan;
                /*now set chan = CHMAX to exit for() loop*/
	    	chan = CHMAX;
    	    	break;
	    }
	    default:
	    {
/*	    	if (chan == 0) printf("Reading ascii spectrum.......\n");*/
	    	break;
    	    }
    	}
    }
    /*for Maestro EOF not reached as there is a trailer*/
    if (md == 9 || md == 10) lchan = rlt[1] + 1;
    chan = lchan;
    
    /*read trailer for Maestro_Spe formt*/
    if (md == 9 || md == 10)
    {
        printf("Reached EOF after reading chan %d \n", chan);
        /*skip_lines(fsp, 10);*/
        /*trailer is of varying length depending on ROIs etc.
          search for energy calibration indicated by $MCA_CAL
          followed by a line with the integer no. of coefficients
          and a second line with the coefficients themselves*/

        while (strncmp(ans, "$MCA_CAL", 8) )
        {
            get_line_file(fsp, ans, CHLEN);
            /*printf("ans:%s:\n",ans);*/
        }
        /*strcmp returns zero if identical, i.e. if not equal to NULL*/
        if (! strncmp(ans, "$MCA_CAL", 8))
        {
            get_line_file(fsp, ans, CHLEN);
            get_line_file(fsp, ans, CHLEN);
            printf("Maestro calibration coefficients: %s\n\n",ans);
        }
    }
    fclose(fsp); 
    return (chan);
} /*END ascii_read()*/

/*==========================================================================*/
/* ascii_write: write an ASCII format spectrum	    	    	    	    */
/****************************************************************************/
void ascii_write(char name[], int numch)
{
    int j;
    FILE *fasc;
    /*open .txt file*/
    if ( (fasc = fopen(name, "w" )) == NULL)
    {
        printf("Cannot open file: %s \n", name);
        return ; 
    }
   
    /*write .txt file using floats*/
    for (j = 0; j < numch; j++) fprintf(fasc, "%d %.5f\n", j, spectrum[j]);
    
    printf(" ==> %s %d chs.\n", name, numch);

    fclose(fasc);
} /*END ascii_write()*/

/*==========================================================================*/
/*Puts spectrum length in filename before ext 	    	    	    	    */
/****************************************************************************/
void chan_num_ext(char fin[], char fout[], int *numch, char ext[])
{
    float tmpf = 0.0;
    int force = 0, i;
    char buf[20];
    char ans[20];

    strcpy(fout,fin);
    set_ext(fout, ext);
    
    if ( ( (*numch) < 2048) || ( ((*numch))%1024 != 0) )
    {
    	while (force == 0)
    	{
    	    printf("Would you like to force the length of the"
    		" spectrum (y/n)?\n");
	    get_ans(ans,1);
    	    if (ans[0] == 'y' || ans[0] == 'Y')
    	    {
    		printf("Enter length in channels (max.%d) eg 8192\n"
    		"(Enter only multiples of 1024 channels)\n",CHMAX);
    		get_val(&tmpf);
		*numch = (int)tmpf;
    		if ( (*numch%1024 == 0) &&  (*numch < CHMAX) )
    		{
    		    force = 1;
    		    break;
    		}
    		else printf("Length entered is not a multiple of 1024"
    			" or >%d\n",CHMAX);
    	    }
    	    else if (ans[0] == 'n' || ans[0] == 'N')
    	    {
    		force = 1;
    		break;
    	    }
    	}
    }
    /*If not an xtrack spectrum, return current filename*/
    /*enters if statement if equal to NULL..strcmp returns zero if identical*/
    if ( strcmp(ext,".spec") )
	return ;
    
    /*if it is an xtrack spectrum, construct new file extension*/
    i = 0;
    itoa(*numch, buf);
    strncpy(ans+i++, "_", 1);
    strncpy(ans+i++, buf, 1);
    if (*numch >= 10000) strncpy(ans+i++, buf+1, 1 );
    if (*numch < 1000)
    {
    	strncpy(ans+i++, buf+1, 1 );
	strncpy(ans+i++, buf+2, 1 );
    }
    else strncpy(ans+i++, "k", 1 );
    strcpy(ans+i++, ext );
    strcpy( strrchr(fout, '.'), ans );
    return ;
} /*END chan_num_ext()*/

/*==========================================================================*/
/*check_ext: check file extension is ext 	    	    	    	    	    	    */
/****************************************************************************/
void check_ext(char fin[], char ext[])
{
    char ans[3];
    
    /*strcmp returns zero if identical, i.e. if not equal to NULL*/
    if ( ! strcmp( (strrchr(fin,'.')), ext ) ) return ;
    else
    {	
	while (1)
	{
    	    printf("File extension is not '%s'\n"
		    "Are you sure this is the correct filename (y/n) ?\n",ext);    	    	      /*no ext found*/
	    get_ans(ans,1);
	    if (ans[0] == 'y' || ans[0] == 'Y') return ;
	    else if (ans[0] == 'n' || ans[0] == 'N')
	    {
	    	printf("Enter new file name inc. extension (%s):\n",ext);
	    	get_line(fin, CHLEN);
  	    	printf("New filename is: %s \n",fin);
	    	return ;
	    }
	}
    }
} /*END check_ext()*/

/*==========================================================================*/
/* col_determ: find if a file has 1 or 2 column format      	    	    */
/****************************************************************************/
int col_determ(FILE *file)
{
    int     col = 0, hash = 0, blnk = 0;
    fpos_t  pos;
    
    /*skip comment lines starting with #*/
    skip_hash(file);
    
    /*store current file position in pos*/
    fgetpos(file, &pos);
     	   
    while (1)
    {
    	/*read until first digit of first number*/
    	while ( isdigit(hash = fgetc(file)) == 0 ) ;
    	
	col = 1;
        /*found first digit of first number. read rest of digits*/
    	/*ignore decimal points*/
    	while ( (isdigit(hash = fgetc(file)) != 0) || (hash == '.') ) ;
    	
        /*push the last character back on the stream*/
    	ungetc(hash, file);
    	/*now test characters until end of line*/
    	while ( (hash = fgetc(file)) != '\n' )
    	{
    	    /*if blank space,increment blank counter and continue*/
    	    if ( (hash == ' ') || (hash == '\t') || (hash == '\r')) blnk++;
    	    /*check how many columns of numbers are present*/
            else if ( (isdigit(hash)) == 0 )
            {
                printf("%s***Input is not a valid data file."
                        " Illegal characters found.%s\n",clr[1],clr[0]);
                col = 0;
                return 0;
            }
    	    /*if character is a digit*/
    	    /*check how many columns of numbers are present*/
    	    else if ( (isdigit(hash)) != 0 )
    	    {
    		/*if number increment col flag*/
    		if (blnk > 0) col++;
    		
   		/*otherwise if > MAXCOLS columns set 1 col. flag
                    and hope for best*/
    		if ( (blnk != 0) && (col > MAXCOLS) )
    		{
    		    col = 1;
    		    break;
    		}
    		/*reset blank counter to zero*/
    		blnk = 0;
    		/*read rest of digits*/
    		while ( (isdigit(hash = fgetc(file)) != 0) || (hash == '.') ) ;
    		
                /*push the last character back on the stream*/
    		ungetc(hash, file);
    	    }
    	}
    	break;
    }
    /*reset file position to start of data*/
    fsetpos(file, &pos);
    return col;
} /*END col_determ()*/

/*==========================================================================*/
/*convert_bytes: convert bytes into Kb, Mb, etc     	    	    	    */
/****************************************************************************/
long convert_bytes(char name[])
{
    float sz;
    int i;
    char bye[4][2];
    struct stat stbuf;
    
    for (i = 0; i < 3; i++)
    {
	bye[i][0] = (char)0;
	bye[i][1] = (char)0;
    }
    strncpy(bye[0],"b",2);
    strncpy(bye[1],"Kb",2);
    strncpy(bye[2],"Mb",2);
    strncpy(bye[3],"Gb",2);
    
    stat(name, &stbuf);
    sz = (float)stbuf.st_size;
    i = 0;
    while ( sz >= 1024.0 )
    {
	sz /= 1024;
	i++;
    }
    printf("\n File size: %ld bytes (%.1f %c%c)\n",
	    (long)stbuf.st_size,sz,bye[i][0],bye[i][1]);
    return (long)stbuf.st_size;
} /*END convert_bytes()*/

/*==========================================================================*/
/* cswap2: swap bits of a 2 byte number     	    	    	    	    */
/****************************************************************************/
int cswap2(int decim)
{
    int i, j, max = 17;
    char bin[17];
    int swapped = 0;
    
    for (i = 1; i >= 0; --i)
    {
	for (j = 0; j < (max-1)/2; ++j)
	{
    	    if (decim & 0x8000) bin[i*8+j] = '1';
	    else bin[i*8+j] = '0';
	    
	    decim <<= 1;
/*	    printf(" i*8+j = %d \n", i*8+j); */
	}
    }
    bin[max-1] = '\0';
/*    printf("bin = %s \n", bin); */
    for (i = 0; i < max-1; ++i)
    {
	if (bin[max-2-i] == '1')
	{
	    swapped = swapped + (int)pow( (double)2, (double)(i));
	}
    }
    
    return swapped;
} /*END cswap2()*/

/*==========================================================================*/
/* cswap4: swap bits of a 4 byte number	    	    	    	    	    */
/****************************************************************************/
int cswap4(int decim)
{
    int i, j, max = 33;
    char bin[33];
    int swapped = 0;
    
    for (i = 3; i >= 0; --i)
    {
	for (j = 0; j < (max-1)/4; ++j)
	{
    	    if (decim & 0x80000000) bin[i*8+j] = '1';
	    else bin[i*8+j] = '0';
	    
	    decim <<= 1;
/*	    printf(" i*8+j = %d \n", i*8+j); */
	}
    }
    bin[max-1] = '\0';
/*    printf("bin = %s \n", bin); */
    for (i = 0; i < max-1; ++i)
    {
	if (bin[max-2-i] == '1')
	    swapped = swapped + (int)pow( (double)2, (double)(i));
    }
    
    return swapped;
} /*END cswap()*/

/*==========================================================================*/
/* decode_mspec_name: decode multiple spectrum filename     	    	    */
/****************************************************************************/
void decode_mspec_name(char name[], int *set, int *mxsp, int *numch,
	int *sz, int bytes)
{
    float   pars[4];
    int i, j;
    char ans[80] = "", ans0[10];
    
    i = (strrchr(name,'_') - 1) - (strchr(name,'_') + 2);
    memset(ans,'\0',sizeof(ans));
    strncpy(ans,strchr(name,'_')+2, i);

    /*read string BACKWARDS since at start is an optional
     	multiple spectra set indicator*/

    i--; j = 0;
    /*read letter(s) at end*/
    memset(ans0,'\0',sizeof(ans0));
    while ( isalpha(ans[i]) ) ans0[j++] = ans[i--];
    reverse(ans0);
    if (!strncmp(ans0, "UI", j) || !strncmp(ans0, "ui", j))
	*sz = (int)sizeof(unsigned int);
    else if (!strncmp(ans0, "I", j) || !strncmp(ans0, "i", j))
	*sz = (int)sizeof(unsigned int);
    else if (!strncmp(ans0, "S", j) || !strncmp(ans0, "s", j))
	*sz = (int)sizeof(short int);
    else if (!strncmp(ans0, "US", j) || !strncmp(ans0, "us", j))
	*sz = (int)sizeof(short int);
    else if (!strncmp(ans0, "F", j) || !strncmp(ans0, "f", j))
    	*sz = (int)sizeof(float);

    /*read number of channels*/
    while (! isdigit(ans[i]) ) i--;
    j = 0;
    memset(ans0,'\0',sizeof(ans0));
    while ( isdigit(ans[i]) ) ans0[j++] = ans[i--];
    reverse(ans0);
    *numch = atoi(ans0);
    
    /*read number of spectra*/
    while (! isdigit(ans[i]) ) i--;
    j = 0;
    memset(ans0,'\0',sizeof(ans0));
    while ( isdigit(ans[i]) ) ans0[j++] = ans[i--];
    reverse(ans0);
    *mxsp = atoi(ans0);
    
    /*check for optional multiple spectra set indicator*/
    if (i > 0)
    {
	while (! isdigit(ans[i]) ) i--;
    	j = 0;
    	memset(ans0,'\0',sizeof(ans0));
    	while ( isdigit(ans[i]) ) ans0[j++] = ans[i--];
    	reverse(ans0);
    	*set = atoi(ans0);
    }
    
    if ((*set)*(*sz)*(*mxsp)*(*numch) != bytes)
    {
	printf("Enter spectra details: set(s) of spectra, no. of spectra,"
		" channels, bytes/chan.\n   (set, mxsp, numch, sz)\n");
    	get_pars(pars,4);
	*set = (int)pars[0];
	*mxsp = (int)pars[1];
    	*numch = (int)pars[2];
	*sz = (int)pars[3];
	if ((*set)*(*sz)*(*mxsp)*(*numch) != bytes)
	{
	    *numch = 0;
	    return ;
	}
    }
    
    if ((*set)*(*sz)*(*mxsp)*(*numch) == bytes)
	printf("Found %d sets of %d spectra: %d channels"
	    " (%d bytes/channel)\n",*set,*mxsp,*numch,*sz);
    
    *mxsp *= *set;
    return ;
} /*END decode_mspec_name()*/

/*==========================================================================*/
/* file_status: check file status. 0 (file exists), -1 (file doesn't exist) */
/****************************************************************************/
void file_status(char name[], char ext[], int len)
{
    char    ans[3];
    struct  stat statbuf;
    
    while ( (stat(name, &statbuf) == 0) )
    {
    	printf("\n*****Output file %s exists. Overwrite (y/n)?\n", name);
	get_ans(ans,1);
	if (ans[0] == 'y' || ans[0] == 'Y') break;
	else if (ans[0] == 'n' || ans[0] == 'N')
	{
	    printf("Enter new file name inc. extension (eg %s):\n",ext);
	    get_line(name, CHLEN);
    	}
    }
} /*END file_status()*/

/*==========================================================================*/
/* genie_read: read an GENIE IEC format spectrum	    	    	    	    */
/****************************************************************************/
int genie_read(char name[])
{
    int     chan = 0, res = 0, lchan = CHMAX;
    char    jk[10] = "";
    FILE    *fsp;
     
    /*opens ascii file*/
    if ((fsp = fopen(name, "r" )) == NULL)
    {
    	printf("Cannot open file: %s \n", name);
    	return -1;	    	    	    
    }
    
    /*skip the first 58 header lines*/
    skip_lines(fsp, 58);
    
    for (chan = 0; chan < CHMAX; chan++)
    { 
        /*throw away first (ID?) string, e.g. A004 and read channel number*/     
        res = fscanf(fsp, "%s %d",jk,&chan);
        res = fscanf(fsp, "%f %f %f %f %f\n",
            &spectrum[chan],&spectrum[chan+1],&spectrum[chan+2],
            &spectrum[chan+3],&spectrum[chan+4]);
	switch (res)
	{
	    case 0:
	    {
	    	printf("\nchan= %d, spec[chan]= %f\n", chan,
		    	spectrum[chan]);
    	    	printf("Read error occurred for file: %s\n", name);
                /*close genie file*/
	    	fclose(fsp);
	    	return -1;
	    }
	    case EOF:
	    {
	    	if (chan < CHMAX)
		{
                    /*take one from chan that for(chan..) loop added*/
                    chan--;
	    	    printf("Reached EOF after reading chan %d \n", chan );
		    if (chan == 0)
		    {
			printf("\n*******Incorrect file format*******\n");
			printf("....Exiting....\n\n");
			return -1;
		    }
		}
		lchan = chan;
                /*now set chan = CHMAX to exit for() loop*/
	    	chan = CHMAX;
    	    	break;
	    }
	    default:
	    {
/*	    	if (chan == 0) printf("Reading genie spectrum.......\n");
                if (chan < 20)
                {
                  printf("%d %f %f %f %f %f\n",chan-4,
                      spectrum[chan],spectrum[chan+1],spectrum[chan+2],
                      spectrum[chan+3],spectrum[chan+4]);
                }*/
		/*increment chan. Normally starts from zero and
                5 channels are read at a time but loop with also increment chan*/
		chan += 4;
	    	break;
    	    }
    	}
    }
    chan = lchan;
    fclose(fsp); 
    return (chan);
} /*END genie_read()*/

/*==========================================================================*/
/* get_ans: get answer without waiting for carriage return   	    	    */
/****************************************************************************/
void get_ans(char ans[], int num)
{
    int     i;
    struct  termios newt, oldt;
    
    while (1)
    {
    	tcgetattr(0, &oldt);
	newt = oldt;
    	newt.c_lflag &= ~ICANON;
    	newt.c_cc[VMIN] = 1;
    	newt.c_cc[VTIME] = 0;
	/*handle sigs*/
/*    	newt.c_lflag |= ISIG;*/
    	tcsetattr(0, TCSANOW, &newt);
    	i = 0;
    	while( (ans[i++] = (char)getchar()) != '\n' && i < 1) ;
    	
	tcsetattr(0, TCSANOW, &oldt);
    	if (ans[i-1] != '\n') printf("\n");
	else if (ans[0] == '\n') continue;
	
    	ans[i] = '\0';
    	return ;
    }
} /*END get_ans()*/

/*==========================================================================*/
/* get_line: read a line from stdin into s, return a length    	    	    */
/****************************************************************************/
void get_line(char ans[], int len)
{
    int     c, i = 0;

    memset(ans,'\0',sizeof(char)*len);
    /*note that scanf() leaves a carriage return in the keyboard buffer*/
    while ( (c = getchar()) != '\n' && c != EOF && i < len) ans[i++] = c;
    
    ans[i] = '\0';
} /*END get_line()*/

/*==========================================================================*/
/* get_line_file: read a line from a file into s, return a length           */
/****************************************************************************/
void get_line_file(FILE *file, char ans[], int len)
{
    int     c, i = 0;

    memset(ans,'\0',sizeof(char)*len);
    /*note that scanf() leaves a carriage return in the keyboard buffer*/
    while ( (c = fgetc(file)) != '\n' && c != EOF && i < len) ans[i++] = c;
    
    /*remove any trailing '\r' carriage returns that are often prior to
      '\n' line feeds*/
    if (ans[i-1] == '\r') i--;
    
    ans[i] = '\0';
} /*END get_line_file()*/

/*==========================================================================*/
/* get_mode: get mode from user input	    	    	    	    	    */
/****************************************************************************/
int get_mode(int md)
{
    char    ans[10] = "";
    
    while(1)
    {
    	printf(" 1) to convert %s (%s) ==> %s (%s)\n",fmti[0],exti[0],fmt[0],ext[0]);
    	printf(" 2) to convert %s (%s) ==> %s (%s)\n",fmti[1],exti[1],fmt[1],ext[1]);
    	printf(" 3) to convert %s (%s) ==> %s (%s)\n",fmti[2],exti[2],fmt[2],ext[2]);
    	printf(" 4) to convert %s (%s) ==> %s (%s)\n",fmti[3],exti[3],fmt[3],ext[3]);
    	printf(" 5) to convert %s (%s) ==> %s (%s)\n",fmti[4],exti[4],fmt[4],ext[4]);
    	printf(" 6) to convert %s (%s) ==> %s (%s)\n",fmti[5],exti[5],fmt[5],ext[5]);
    	printf(" 7) to convert %s (%s) ==> %s (%s)\n",fmti[6],exti[6],fmt[6],ext[6]);
    	printf(" 8) to convert %s (%s) ==> %s (%s)\n",fmti[7],exti[7],fmt[7],ext[7]);
    	printf(" 9) to convert %s (%s) ==> %s (%s)\n",fmti[8],exti[8],fmt[8],ext[8]);
    	printf(" a) to convert %s (%s) ==> %s (%s)\n",fmti[9],exti[9],fmt[9],ext[9]);
    	printf(" g) to gainmatch a RadWare spectrum\n");
    	printf(" 0) Quit\n");
	get_ans(ans,1);
 	if (ans[0] - '0' >= 0 && ans[0] - '0' <= NUMOPT)
	{
	    md = ans[0] - '0';
	    break;
	}
	else if (ans[0] == 'a' || ans[0] == 'A')
	{
	    md = 10;
	    break;
	}
	else if (ans[0] == 'g' || ans[0] == 'G')
	{
	    md = NUMOPT;
	    break;
	}
    }
    return md;
} /*END get_mode()*/

/*==========================================================================*/
/* get_pars: extract comma or space separated numbers from string ans0	    */
/****************************************************************************/
void get_pars(float pars[], int num)
{
    int     i, j = 0, k, minus;
    char    ans0[CHLEN] = "", ans1[40] = "";

    get_line (ans0, CHLEN);
    
    for (i = 0; i < num; i++)
    {
	k = 0;
	minus = 0;
    	memset(ans1,'\0',sizeof(ans1));
	while ( ! isdigit(ans0[j]) )
	{
	    if (ans0[j] == '-') minus = 1;
	    j++;
	}
	if (minus == 1)
	{
	    j -= 1;
	    ans0[j] = '-';
	}	
	while( ans0[j] != '\n' && ans0[j] != ' ' && ans0[j] != ','
		&& j < strlen(ans0) ) ans1[k++] = ans0[j++];
	
	j++;
	pars[i] = atof(ans1);
    }
    printf("\n");
} /*END get_pars()*/

/*==========================================================================*/
/* get_pars_file: extract separated numbers from string ans0 read from file */
/****************************************************************************/
void get_pars_file(FILE *file, float pars[], int num)
{
    int     i, j = 0, k, minus;
    char    ans0[CHLEN] = "", ans1[40] = "";

    get_line_file(file, ans0, CHLEN);
    
    for (i = 0; i < num; i++)
    {
	k = 0;
	minus = 0;
    	memset(ans1,'\0',sizeof(ans1));
	while ( ! isdigit(ans0[j]) )
	{
	    if (ans0[j] == '-') minus = 1;
	    j++;
	}
	if (minus == 1)
	{
	    j -= 1;
	    ans0[j] = '-';
	}	
	while( ans0[j] != '\n' && ans0[j] != ' ' && ans0[j] != ','
		&& j < strlen(ans0) ) ans1[k++] = ans0[j++];
	
	j++;
	pars[i] = atof(ans1);
    }
} /*END get_pars_file()*/

/*==========================================================================*/
/* get_val: extract a number from string ans0	    	    	    	    */
/****************************************************************************/
void get_val(float *val)
{
    int     j = 0, k, minus;
    char    ans0[CHLEN] = "", ans1[40] = "";

    get_line(ans0, CHLEN);
    
    k = 0;
    minus = 0;
    memset(ans1,'\0',sizeof(ans1));
    while ( ! isdigit(ans0[j]) )
    {
    	if (ans0[j] == '-') minus = 1;
    	j++;
    }
    if (minus == 1)
    {
    	j -= 1;
    	ans0[j] = '-';
    }	    
    while( ans0[j] != '\n' && ans0[j] != ' ' && ans0[j] != ','
    	    && j < strlen(ans0) ) ans1[k++] = ans0[j++];

    j++;
    *val = atof(ans1);
    printf("\n");
} /*END get_val()*/

/*===========================================================================*/
/* itoa: convert integer n to string s*/
/*****************************************************************************/
void itoa(int n, char s[])
{
    int i = 0, sign;
    
    /*record sign and make n positive*/
    if ((sign = n) < 0) n = -n;

    /* generates digits in reverse order */
    do {
        /*get next digit*/
    	s[i++] = n % 10 + '0';
      /* delete it*/
    } while ((n /= 10) > 0);
    if (sign < 0) s[i++] = '-';
    
    s[i] = '\0';
    reverse(s);
} /*END itoa()*/

/*===========================================================================*/
/* maestro_read: read the maestro format spectrum                            */
/*****************************************************************************/
int maestro_read(char name[])
{
    int *counts, i = 0;
    char dt[10] = "";
    FILE *fsp;
        
    /*opens read only Maestro file*/
    if ( (fsp = fopen(name, "r")) == NULL)
    {
    	printf("Cannot open file: %s\n", name);
	return -1;
    }    
    /*clear the header*/
    memset(&maest_header, 0, sizeof(maest_header));
    /*construct Maestro header*/
    fread(&maest_header, sizeof(maest_header), 1, fsp);
        
/*    printf(" maest_header.q1 = %d\n maest_header.q2 = %d \n"
	   " maest_header.q3 = %d\n maest_header.q4 = %d \n"
	   " maest_header.real = %d\n maest_header.lve = %d \n"
	   " maest_header.dt = %s\n maest_header.sttm = %s \n"
	   " maest_header.off = %d\n maest_header.channels = %d\n",
	    maest_header.q1, maest_header.q2, maest_header.q3,
	    maest_header.q4, maest_header.real, maest_header.lve,
            maest_header.dt, maest_header.sttm, maest_header.off,
            maest_header.channels);*/
    
    /*unix byte swapping option*/
    if (maest_header.channels > CHMAX && cswap4(maest_header.channels) > CHMAX )
    {
    	printf("Unrecognised format. Exiting.....\n");
	fclose(fsp);
    	return -1;
    }
    if (maest_header.channels > CHMAX || 
            maest_header.channels < (i = cswap2(CHMAX)) )
    {
	/*swap the bytes in the headers*/
    	maest_header.channels = cswap2(maest_header.channels);
	maest_header.real = cswap4(maest_header.real);
	maest_header.lve = cswap4(maest_header.lve);
        
	printf(".......SWAPPING BYTES read from file.......\n");
    	/*allocate sufficient memory for spectrum*/
    	counts = (int *) malloc( maest_header.channels*sizeof(int) );
    	/*read the data*/
	fread(counts, maest_header.channels*sizeof(int), 1, fsp);
	/*read the trailer*/
	fread(&maest_trailer, sizeof(maest_trailer), 1, fsp);
    	maest_trailer.g[0] = cswap4(maest_trailer.g[0]);
    	maest_trailer.g[1] = cswap4(maest_trailer.g[1]);
    	maest_trailer.g[2] = cswap4(maest_trailer.g[2]);

      	/*fill spectrum array*/
    	for (i = 0; i < (int) maest_header.channels; i++)
	    swapb4( (char *) (counts + i) );	
    } /*end of byte swapping loop for unix*/    
    else
    {
    	/*allocate sufficient memory for spectrum*/
    	counts = (int *) malloc( maest_header.channels*sizeof(int) );
	/*read the data*/
	fread(counts, maest_header.channels*sizeof(int), 1, fsp);
	/*read the trailer*/
	fread(&maest_trailer, sizeof(maest_trailer), 1, fsp);
    }
    /*fill spectrum array*/
    for (i = 0; i < maest_header.channels; i++)
        spectrum[i] = (float)*(counts + i);
    
    
    /*print real and live times to screen and convert from 20ms units to sec*/
    /*Copy day and month*/
    strncpy(dt, maest_header.dt, 5);
    /*check year, if maest_header.dt[7] == 1 year is 2000+YY else 1900+YY*/
    if (maest_header.dt[7] == '1') strncpy(dt+5, "20", 2);
    else strncpy(dt+5, "19", 2);
    strncpy(dt+7, maest_header.dt+5, 2);
    dt[9] = '\0';
    printf("Spectrum info: %s at %c%c:%s\n"
           "               Real time: %d s\n"
           "               Live time: %d s\n",dt,maest_header.sttm[0],
            maest_header.sttm[1],maest_header.sttm+2,
            (int)(maest_header.real*0.02),(int)(maest_header.lve*0.02));
    
    /*print potentially useful trailer information*/
    printf("Maestro energy calibration coeffs: %f %f %f\n",
          maest_trailer.g[0],maest_trailer.g[1],maest_trailer.g[2]);
    
    free(counts);
    fclose(fsp);
    return maest_header.channels;
} /*END maestro_read()*/

/*==========================================================================*/
/* num_fname: create numbered filenames	    	    	    	    	    */
/****************************************************************************/
void num_fname(char name[], int num)
{
    int i = 0, j, len;
    char ans[MXNUMDIG+1] = "", buf[MXNUMDIG+3] = "";
    
    len = strlen(name);
    strncpy(buf+i++, "_", 1);
    itoa(num, ans);
    
    /*add zeroes to pad out file name if num is less than the maximum
        number of digits*/
    for (j = 0; j < (MXNUMDIG - (int)strlen(ans)); j++)
        strncpy(buf+i++, "0", 1);
    
    strncpy(buf+i, ans, (int)strlen(ans));
    i += (int)strlen(ans);
    strncpy(buf+i++, ".", 1);
    buf[i] = '\0';
    if ( (strrchr(name,'.')) )
	strncpy( (strrchr(name,'.')), buf, i);
    else strncat(name, buf, i);
    
    /*length without ext.*/
    i = strrchr(name,'.') - &name[0];
    for (j = 1; j < (len - i); j++) name[i+j] = (char)0;
    
    return;
} /*END num_fname()*/

/*===========================================================================*/
/* rad_read: read the radware format spectrum */
/*****************************************************************************/
int rad_read(char name[])
{
    float *counts;
    int i = 0;
    FILE *fsp;
        
    /*opens read only RadWare file*/
    if ( (fsp = fopen(name, "r")) == NULL)
    {
    	printf("Cannot open file: %s\n", name);
	return -1;
    }
    
    /*clear the header*/
    memset(&radheader, 0, sizeof(radheader));
    /*construct radware header*/
    fread(&radheader, sizeof(radheader), 1, fsp);
        
/*    printf(" radheader.channels = %d\n radheader.q1 = %d \n"
	   " radheader.q3 = %d\n radheader.q4 = %d \n"
	   " radheader.q5 = %d\n radheader.size = %d \n",
	    radheader.channels, radheader.q1, radheader.q3,
	    radheader.q4, radheader.q5, radheader.size);*/
       
    /*unix byte swapping option*/
    if (radheader.channels > CHMAX && cswap4(radheader.channels) > CHMAX )
    {
    	printf("Unrecognised format. Exiting.....\n");
	fclose(fsp);
    	return -1;
    }
    if (radheader.channels > CHMAX)
    {
	/*swap the bytes in the headers*/
    	radheader.channels = cswap4(radheader.channels);
	radheader.q1 = cswap4(radheader.q1);
	radheader.q2 = cswap4(radheader.q2);
    	radheader.q3 = cswap4(radheader.q3);
    	radheader.q4 = cswap4(radheader.q4);
    	radheader.q5 = cswap4(radheader.q5);
	radheader.size = cswap4(radheader.size);
		
	printf(".......SWAPPING BYTES read from file.......\n");
    	/*allocate sufficient memory for spectrum*/
    	counts = (float *) malloc( radheader.channels*sizeof(float) );
	
    	/*read the data*/
	fread(counts, radheader.size, 1, fsp);
	/*read the trailer*/
	fread(&radtrailer.size, 4, 1, fsp);
	radtrailer.size = cswap4(radtrailer.size);
    	
    	/*fill spectrum array*/
    	for (i = 0; i < (int) radheader.channels; i++)
	    swapb4( (char *) (counts + i) );	
    }
    /*end of byte swapping loop for unix*/    
    else
    {
    	/*allocate sufficient memory for spectrum*/
    	counts = (float *) malloc( radheader.channels*sizeof(float) );
	/*read the data*/
	fread(counts, radheader.size, 1, fsp);
	/*read the trailer*/
	fread(&radtrailer.size, 4, 1, fsp);    
    }
    /*fill spectrum array*/
    for (i = 0; i < radheader.channels; i++)
    	spectrum[i] = (float)*(counts + i);
    
    free(counts);
    fclose(fsp);
    return radheader.channels;
} /*END rad_read()*/

/*==========================================================================*/
/* rad_write: prepare and write the radware format spectrum      	    */
/****************************************************************************/
void rad_write(char name[], int numch)
{
    int     j;
    FILE    *fsp;
        
    /*open .spe write only file*/
    if ( (fsp = fopen(name, "w" )) == NULL)
    {
    	printf("Cannot open file: %s \n", name);
	return ; 
    } 	    	    
    	  
    /*clear the header*/
    memset(&radheader, 0, sizeof(radheader));
    /*clear the trailer*/
    memset(&radtrailer, 0, sizeof(radtrailer));
    
    /*construct radware header*/
    radheader.q1 = 24;
    
    /*length without ext.*/
    j = strrchr(name,'.') - &name[0];
    /*copy max. 8 bytes to radheader.name*/
    strncpy(radheader.name, name, 8);
    /*set any extra characters so spaces*/
    if (j < 8) memset(&radheader.name[j], ' ', 8-j);
    
    radheader.channels = numch;
    radheader.q2 = 1;
    radheader.q3 = 1;
    radheader.q4 = 1;
    radheader.q5 = 24;
    radheader.size = numch * sizeof(float);
    
    radtrailer.size = numch * sizeof(float);
  
    fwrite(&radheader, sizeof(radheader), 1, fsp);
    fwrite(spectrum, radheader.size, 1, fsp);
    fwrite(&radtrailer, sizeof(float), 1, fsp);

    printf(" ==> %s %d chs.\n", name, numch);
    
    fclose(fsp);
} /*END rad_write()*/

/*==========================================================================*/
/* read_lst: read next spectrum name from list file   	    	    	    */
/****************************************************************************/
int read_lst(char inname[], int lst)
{
    static int  fn = 0;
    int         res;
    char        listname[CHLEN] = "";
    static FILE *flst;
    
    /*open list file to read spec names*/
    if (fn == 0 && lst == 1)
    {
        if (md == NUMOPT)
    	    printf("Type filename containing list of spectrum file names"
		    " and coefficients: \n");
	else printf("Type filename containing list of spectrum file names:\n");
	
	get_line(listname, CHLEN);	
    }
    /*open file on first time in function*/
    if (fn == 0)
    {
        if (lst == -1) strcpy(listname,inname);
    	if ((flst = fopen(listname, "r" )) == NULL)
    	{
    	    printf("Cannot open file: %s \n", listname);
    	    return -1;			
    	}
    }
    
    /*read and store ascii file names*/
    /*skip comments lines starting with #*/
    skip_hash(flst);
    if (md == NUMOPT) res = fscanf(flst, "%s %f %f %f", inname,
    		&gain[0], &gain[1], &gain[2]);
    else res = fscanf(flst, "%s", inname);
    
    switch (res)
    {
    	case 0:
    	{
    	    printf("Read error occurred for file: %s\n", listname);
    	    fclose(flst);
    	    return -1;
    	}
    	case EOF:
    	{
    	    printf("\n\tRead %d spectrum names\n\n", fn);
    	    fclose(flst);	    
    	    return -1;
    	}
    	default:
    	{
    	    fn++;
    	    printf("Read filename %d from list: %s\n",fn, inname);
	    return fn;
    	}
    }
} /*END read_lst()*/

/*==========================================================================*/
/* read_spec: call appropriate spectrum_read function based on mode    	    */
/****************************************************************************/
int read_spec(char name[], int md)
{
    int i = 0;
    
    if (md == 1) i = rad_read(name);
    else if (md == 2) i = ascii_read(name);
    else if (md == 3) i = ascii_read(name);
    else if (md == 4) i = maestro_read(name);
    else if (md == 5) i = maestro_read(name);
    else if (md == 8) i = genie_read(name);
    else if (md == 9) i = ascii_read(name);
    else if (md == 10) i = ascii_read(name);
    else if (md == NUMOPT) i = rad_read(name);
    else return -1;
    
    return i;
} /*END read_spec()*/

/*==========================================================================*/
/* reverse: reverse string s in place	    	    	    	    	    */
/****************************************************************************/
void reverse(char s[])
{
    int c, i, j;

    for (i = 0, j = strlen(s)-1; i < j; i++, j--)
    {
        c = s[i];
        s[i] = s[j];
        s[j] = c;
    }
} /*END reverse()*/

/*==========================================================================*/
/* set_ext: set file extension of string name[] to ext[]   	    	    */
/****************************************************************************/
void set_ext(char name[], char ext[])
{    
    /*if not equal to NULL, find '.' then copy ext with +1 for term char*/
    if ( (strrchr(name,'.')) )
        strncpy( (strrchr(name,'.')), ext, (int)(strlen(ext)+1) );
    /*if equal to NULL just add ext to the end*/
    else strcat( name, ext );
} /*END set_ext()*/

/*==========================================================================*/
/* skip_hash: skip comment lines starting with hash (#) at start of file    */
/****************************************************************************/
void skip_hash(FILE *file)
{
    int     i = 0, hash = 0;
    fpos_t  pos;
    
    for (i = 0; i < 1000; i++)
    {
    	/*store position corresponding to start of line*/
    	fgetpos(file, &pos);
        /*check for carriage return in the keyboard buffer
            due to scanf function*/
        if ( (hash = fgetc(file)) == '\n' ) fgetpos(file, &pos);
    	else fsetpos(file, &pos);
    	if ( (hash = fgetc(file)) != '#' )
    	{
    	    /*reset position to start of line*/
    	    fsetpos(file, &pos);
    	    break;
    	}
    	/*read rest of line*/
    	while ( ( (hash = fgetc(file)) != '\n' ) && (hash != EOF) )
            ;
    }
} /*END skip_hash()*/

/*==========================================================================*/
/*skip_lines: skip lns lines in file                                        */
/****************************************************************************/
void skip_lines(FILE *file, int lns)
{
    int cnt = 0, hash = 0;
    
    while (cnt < lns)
    {
        while ( ( (hash = fgetc(file)) != '\n' ) && (hash != EOF) )
           ;
        
        if (hash == EOF)
        {
            printf("Found EOF...returning\n");
            return ;
        }
        /*increment counter on each carriage return*/
        cnt++;
    }
} /*END skip_lines()*/

/*==========================================================================*/
/* store_colours: store colours in clr[][] array                            */
/****************************************************************************/
void store_colours()
{
    /*none*/
    strcpy(clr[0], "\033[0;0m");
    /*white on red background*/
    strcpy(clr[1], "\033[0;41m");
    /*green no background*/
    strcpy(clr[2], "\033[0;32m");
    /*white on blue background*/
    strcpy(clr[3], "\033[0;44m");
    /*bold*/
    strcpy(clr[4], "\033[1m");
    /*underline*/
    strcpy(clr[5], "\033[0;4m");
} /*END store_colours()*/  

/*==========================================================================*/
/* swapb2: swap 2 array elements (a 2 byte number) in place 	    	    */
/****************************************************************************/
void swapb2(char *buf)
{
    char c;
    c = buf[1]; buf[1] = buf[0]; buf[0] = c;
} /*END swapb2()*/

/*==========================================================================*/
/* swapb4: swap 4 array elements (a 4 byte number) in place 	    	    */
/****************************************************************************/
void swapb4(char *buf)
{
    char c;
    c = buf[3]; buf[3] = buf[0]; buf[0] = c;
    c = buf[2]; buf[2] = buf[1]; buf[1] = c;    
} /*END swapb4()*/

/*==========================================================================*/
/* write_spec: call appropriate spectrum_write function based on mode       */
/****************************************************************************/
void write_spec(char name[], int numch, int md)
{
    if (md == 1) ascii_write(name, numch);
    else if (md == 2) rad_write(name, numch);
    else if (md == 3) xtrack_write(name, numch);
    else if (md == 4) ascii_write(name, numch);
    else if (md == 5) rad_write(name, numch);
    else if (md == 6) ascii_write(name, numch);
    else if (md == 7) rad_write(name, numch);
    else if (md == 8) rad_write(name, numch);
    else if (md == 9) rad_write(name, numch);
    else if (md == 10) ascii_write(name, numch);
    else if (md == NUMOPT) rad_write(name, numch);

    return ;
} /*END write_spec()*/

/*==========================================================================*/
/* xtrack_read: read an xtrack (GASPWARE) format spectrum   	    	    */
/****************************************************************************/
void xtrack_read(char name[], int *numch, int mxsp, int sz, int nsp, int flg)
{
    int i = 0, last_nonzero_channel = 0, mxcnts = 0;
    unsigned int *xtrack_spec;
    FILE *fp;
         
    /*open xtrack file*/
    if ((fp = fopen(name, "r" )) == NULL)
    {
    	printf("Cannot open file: %s \n", name);
    	*numch = -1;
	return ;	    	    	    
    }
    
    if (mxsp > 1 && nsp > 0) fseek(fp, (int)(nsp*(*numch)*sz), 0);
    
    xtrack_spec = (unsigned int *) malloc(*numch*sz);
    
    while ( fread(xtrack_spec, *numch*sz, 1, fp) != 1 )
    {
    	*numch /= 2;	    
	printf("Trying spectrum length: %d channels\n", *numch);
	if (*numch <= 1024)
	{
    	    printf("Error reading file: %s \n", name);
	    *numch = -1;
	    return ;
	}
    }
    printf(" Length = %d channels was successful\n", *numch);
	        
    for (i = 0; i < *numch; i++)
    {
	/*Note that in the conversion int is necessary first to get the
	    	sign correct*/
	spectrum[i] = (float) (int) *(xtrack_spec + i);
/*	printf("i = %d spectrum[i] = %f\n", i, spectrum[i]); */
	
	if (spectrum[i] != 0) last_nonzero_channel = i;
	
	if (spectrum[i] > mxcnts) mxcnts = (int) spectrum[i];
	
    }
    
    if (mxcnts > 10000000)    /*Probably not an Xtrack format spectrum*/
    {
	printf("***WRONG FORMAT. NOT AN XTRACK SPECTRUM***\n");
	*numch = -1;
	return ;
    }
    
    /*only remove zeroes if not multiple spectrum file or list*/
    if (flg != 1)
    {
    	while (last_nonzero_channel < (*numch/2) && *numch >= 1024) *numch /= 2;
    	printf("Real length of spectrum (after removing zeros) = %d channels\n",
	    *numch);
    }
/*    else printf("Length of spectrum = %d channels\n",*numch);*/
       
    fclose(fp);
    return ; 
    
} /*END xtrack_read()*/

/*==========================================================================*/
/* xtrack_write: write a (GASPWARE) format spectrum    	    	    	    */
/****************************************************************************/
void xtrack_write(char name[], int numch)
{
    int i = 0;
    unsigned int tmp_spec[CHMAX];
    unsigned int numbytes;
    FILE *fsp;
         
    /*open .spec write only file*/
    if ( (fsp = fopen(name, "w" )) == NULL)
    {
    	printf("Cannot open file: %s \n", name);
	return ; 
    } 	    	    

    numbytes = numch*sizeof(unsigned int);
        
    for (i = 0; i < numch; i++) tmp_spec[i] = (unsigned int)spectrum[i];
    
    fwrite(&tmp_spec, numbytes, 1, fsp);

    printf(" ==> %s %d chs.\n", name, numch);
    
    fclose(fsp);
} /*END xtrack_write()*/
