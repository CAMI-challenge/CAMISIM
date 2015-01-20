/* Draw secondary structure diagrams given a postscript SS template.
 * Initial development of this program was for SSU rRNA structures
 * with templates derived from Gutell's CRW. 
 *
 * EPN, Mon Jun 23 14:46:05 2008
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <time.h>

#include "easel.h"
#include "esl_distance.h"
#include "esl_dmatrix.h"
#include "esl_fileparser.h"
#include "esl_getopts.h"
#include "esl_keyhash.h"
#include "esl_msa.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_stack.h"
#include "esl_tree.h"
#include "esl_vectorops.h"
#include "esl_wuss.h"

#define SSDRAWINFINITY 987654321
#define ERRBUFSIZE 1024

#define ALIMODE 0
#define INDIMODE 1
#define SIMPLEMASKMODE 2
#define DRAWFILEMODE 3

#define RAINBOWRHSCHEME 0
#define RAINBOWRLSCHEME 1
#define NRAINBOWRHSCHEME 11
#define NRAINBOWRLSCHEME 11

#define RBSIXRHSCHEME 2
#define RBSIXRLSCHEME 3
#define NRBSIXRHSCHEME 6
#define NRBSIXRLSCHEME 6

#define RBFIVERHSCHEME 4
#define RBFIVERLSCHEME 5
#define NRBFIVERHSCHEME 5
#define NRBFIVERLSCHEME 5


#define NOC 10
#define CYANOC 0
#define MAGENTAOC 1
#define YELLOWOC 2
#define BLACKOC 3
#define LIGHTGREYOC 4
#define DARKGREYOC 5
#define REDOC 6
#define PURPLEOC 7
#define ORANGEOC 8
#define WHITEOC 9

#define LEGTEXTNCHARS 60
#define NCMYK 4
#define ICYAN 0
#define IMAGENTA 1
#define IYELLOW 2
#define IBLACK 3

/* define the color for blank cells where no value is appropriate */
#define BLANKCYAN    0.0
#define BLANKMAGENTA 0.0
#define BLANKYELLOW  0.0
#define BLANKBLACK   0.5
/*if I extend to allow RED, BLUE, GREEN with combos of CMYK, but it's a pain right now due to implementation 
  #define IRED 4
  #define IBLUE 5
  #define IGREEN 6
  #define NCOLORS 7
*/
#define LEG_NBOXES  11
#define LEG_BOXSIZE 24.
#define LEG_MINFONTSIZE 10
#define LEGX_OFFSET 24.
#define LEGY_OFFSET -24.
#define SPECIAL_FONT "Courier-BoldOblique"
#define LEG_FONT "Courier-Bold"
#define LEG_EXTRA_COLUMNS 12 /* how many extra columns we need for printing stats in the legend */

#define DEFAULT_FONT "Courier-Bold"
#define RESIDUES_FONT "Helvetica-Bold"
#define HUNDREDS_FONT "Helvetica"

#define SS_BOXSIZE 8.

#define RESIDUES_FONTSIZE 8.
#define HUNDREDS_FONTSIZE 8.
#define LEG_FONTSIZE_UNSCALED 9.6
/*#define HEADER_FONTSIZE_UNSCALED 14.4*/
#define HEADER_FONTSIZE_UNSCALED 12
#define HEADER_MODELNAME_MAXCHARS 20
#define TICKS_LINEWIDTH 2.
#define BP_LINEWIDTH 1.

#define POSTSCRIPT_PAGEWIDTH 612.
#define POSTSCRIPT_PAGEHEIGHT 792.
#define PAGE_TOPBUF 12.
#define PAGE_SIDEBUF 12.
#define PAGE_BOTBUF 12.
#define COURIER_HEIGHT_WIDTH_RATIO 1.65

/* Structure: scheme_color_legend
 * Incept:    EPN, Thu Jun 25 20:20:38 2009
 *
 * Parameters describing a one-dimensional legend of colors
 * from a preset scheme for use in a SSPostscript_t data structure.
 */
typedef struct scheme_color_legend_s {
  int    scheme;            /* preset color scheme index */
  int    nbins;             /* number of colors (bins) in this scheme */
  char   *text1;            /* first line of text for legend, a single string */
  char   *text2;            /* second line of text for legend, a single string */
  float *limits;            /* [nbins+1] limits for each bin, limits[0] is min value we would expect to see, limits[nbins] is max */
  int   *counts;            /* [nbins] number of cells we've painted each color */
  int   *counts_masked;     /* [nbins] number of cells within mask ('1's) that we've painted each color */
  int    ints_only_flag;    /* TRUE if possible values are only integers, legend values will be drawn differently in this case */
} SchemeColorLegend_t;

/* Structure: onecell_color_legend
 * Incept:    EPN, Tue Sep 30 13:06:15 2008
 *
 * Parameters describing a single colored cell legend for a
 * SSPostscript_t data structure.
 */
typedef struct onecell_color_legend_s {
  float  col[NCMYK];        /* [CMYK] color value for the cell */
  char   *text;             /* text for legend */
  int    nres;              /* number of residues colored by the color in col[NCMYK] */
  int    nres_masked;       /* number of residues within a mask colored by the color in col[NCMYK] */
} OneCellColorLegend_t;

/* Structure: ss_postscript
 * Incept:    EPN, Mon Jun 23 15:50:30 2008
 *
 * A clumsy data structure for storing the information that will
 * become a postscript secondary structure diagram based on a 
 * template created by Robin Gutell and colleagues.
 *
 */
typedef struct ss_postscript_s {
  int     npage;        /* number of pages in eventual postscript */
  char   *modelname;    /* name of model, read from template file */
  int    *modeA;        /* [0..npage-1] page mode, ALIMODE, INDIMODE, or SIMPLEMASKMODE */
  char  **descA;        /* [0..npage-1] description for each page */
  int     desc_max_chars; /* max num characters for a page description */
  float   headerx;      /* x coordinate (bottom left corner) of header area */
  float   headery;      /* y coordinate (bottom left corner) of header area */
  float   headerx_charsize;/* size of a character in x-dimension in the header */
  float   headery_charsize;/* size of a character in y-dimension in the header */
  float   headerx_desc; /* x coordinate (bottom left corner) of header area */
  float   legx;         /* x coordinate (bottom left corner) of legend area */
  float   legy;         /* y coordinate (bottom left corner) of legend area */
  float   cur_legy;     /* y coordinate of current line in legend */
  float   legx_charsize;/* size of a character in x-dimension in the legend */
  float   legy_charsize;/* size of a character in y-dimension in the legend */
  int     legx_max_chars; /* max num residues in x direction we can print in legend before running off page */
  int     legy_max_chars; /* max num residues in y direction we can print in legend before running off page */
  int     legx_stats;   /* x position for printing stats in the legend */
  float   pagex_max;    /* max x position on page */
  float   pagey_max;    /* max y position on page */
  float   scale;        /* scale parameter, read from template file */
  char  **regurgA;      /* [0..nregurg-1][] lines from the template file to regurgitate, these are unchanged. */
  int     nregurg;      /* number of lines (char *'s) in the regurg_textAA 2D array */
  float  *hundredsxA;   /* [0..nhundreds-1] x value for hundreds (el 0 is for '100', 1 is for '200', etc.) */
  float  *hundredsyA;   /* [0..nhundreds-1] y value for hundreds (el 0 is for '100', 1 is for '200', etc.) */
  int     nhundreds;    /* number of elements in hundredsx and hundredsy */
  float  *ticksx1A;     /* [0..nticks-1] x begin value for ticks */
  float  *ticksx2A;     /* [0..nticks-1] x end   value for ticks */
  float  *ticksy1A;     /* [0..nticks-1] y begin value for ticks */
  float  *ticksy2A;     /* [0..nticks-1] x end   value for ticks */
  int     nticks;       /* number of ticks */
  float  *bpx1A;        /* [0..nbp-1] x begin value for bp connect line */
  float  *bpx2A;        /* [0..nbp-1] x end   value for bp connect line */
  float  *bpy1A;        /* [0..nbp-1] y begin value for bp connect line */
  float  *bpy2A;        /* [0..nbp-1] x end   value for bp connect line */
  int     nbp;          /* number of bp */
  float  *rxA;          /* [0..rflen-1] x coordinate for each residue in the eventual postscript */
  float  *ryA;          /* [0..rflen-1] y coordinate for each residue in the eventual postscript */
  int     rflen;         /* the number of residues in the template file */
  char  **rrAA;         /* [0..npage-1][0..rflen-1] residue character in the eventual postscript */
  float ***rcolAAA;     /* [0..npage-1][0..rflen-1][0..3] color for block on page p, position c, CMYK in the eventual postscript */
  OneCellColorLegend_t ***occlAAA;/* [0..npage-1][0..l..nocclA[p]  ptr to one cell color legend l for page p */
  int     *nocclA;      /* [0..npage-1] number of one cell color legends for each page */
  SchemeColorLegend_t  **sclAA;/* [0..npage-1]  ptr to scheme color legend l for page p, NULL if none */
  char    *mask;        /* mask for this postscript, columns which are '0' get drawn differently */
  int      nalloc;      /* number of elements to add to arrays when reallocating */
  int      msa_nseq;    /* number of sequences in the msa, impt b/c msa->nseq will be 0 if --small */
  int     *msa_ct;      /* [1..ps->rflen] CT array for msa this postscript corresponds to, 
			 * msa_ct[i] is the position that consensus residue i base pairs to, or 0 if i is unpaired. */
  int      msa_nbp;     /* number of bps read from current MSA (in msa_ct), should equal nbp, but only if bps read from template file */
  int     *msa_rf2a_map;/* [0..ps->rflen-1]     = apos, apos is the alignment position (0..msa->alen-1) that is non-gap RF position rfpos (for rfpos in 0..rflen-1) */
  int     *msa_a2rf_map;/* [0..ps->msa->alen-1] = rfpos, rfpos is the non-gap RF position (0..ps->rflen) that is alignment position apos (for apos in 0..ps->msa_alen-1) */
  int     *uaseqlenA;   /* [0..ps->msa->nseq-1] unaligned sequence length for all sequences in the MSA, only computed if --indi */
  int     *seqidxA;     /* [0..ps->npage-1] the sequence index in the MSA each page corresponds to, only valid if --indi */
  ESL_MSA *msa;         /* pointer to MSA this object corresponds to */
} SSPostscript_t;

static SSPostscript_t *create_sspostscript();
static int  setup_sspostscript(SSPostscript_t *ps, char *errbuf);
static OneCellColorLegend_t *create_onecell_colorlegend(float *cmykA, int nres, int nres_masked);
static SchemeColorLegend_t  *create_scheme_colorlegend(int scheme, int ncols, float *limits, int ints_only_flag);
static int  add_text_to_scheme_colorlegend(SchemeColorLegend_t *scl, char *text, int legx_max_chars, char *errbuf);
static int  add_text_to_onecell_colorlegend(SSPostscript_t *ps, OneCellColorLegend_t *occl, char *text, int legx_max_chars, char *errbuf);
static int  add_page_desc_to_sspostscript(SSPostscript_t *ps, int page, char *text, char *errbuf);
static int  add_diffmask_page_desc_to_sspostscript(SSPostscript_t *ps, int page, char *mask1, char *mask2, char *errbuf);
static int  draw_sspostscript(FILE *fp, const ESL_GETOPTS *go, char *errbuf, char *command, char *date, float ***hc_scheme, SSPostscript_t *ps, int nused);
static int  draw_legend_column_headers(FILE *fp, SSPostscript_t *ps, char *errbuf);
static int  draw_onecell_colorlegend(FILE *fp, OneCellColorLegend_t *occl, SSPostscript_t *ps, int occl_idx);
static int  draw_scheme_colorlegend(const ESL_GETOPTS *go, FILE *fp, SchemeColorLegend_t *scl, float **hc_scheme, SSPostscript_t *ps, int page);
static void free_sspostscript(SSPostscript_t *ps);
static int  add_pages_sspostscript(SSPostscript_t *ps, int ntoadd, int page_mode);
static int  parse_template_file(char *filename, const ESL_GETOPTS *go, char *errbuf, int msa_rflen, SSPostscript_t **ret_ps);
static int  parse_template_page(ESL_FILEPARSER *efp, const ESL_GETOPTS *go, char *errbuf, SSPostscript_t **ret_ps);
static int  parse_modelname_section(ESL_FILEPARSER *efp, char *errbuf, SSPostscript_t *ps);
static int  parse_scale_section(ESL_FILEPARSER *efp, char *errbuf, SSPostscript_t *ps);
static int  parse_ignore_section(ESL_FILEPARSER *efp, char *errbuf, int *ret_read_showpage);
static int  parse_regurgitate_section(ESL_FILEPARSER *efp, char *errbuf, SSPostscript_t *ps);
static int  parse_text_section(ESL_FILEPARSER *efp, char *errbuf, SSPostscript_t *ps);
static int  parse_lines_section(ESL_FILEPARSER *efp, char *errbuf, SSPostscript_t *ps);
static int  validate_justread_sspostscript(SSPostscript_t *ps, char *errbuf);
static int  validate_and_update_sspostscript_given_msa(const ESL_GETOPTS *go, SSPostscript_t *ps, ESL_MSA *msa, int msa_nseq, char *errbuf);
static int  read_mask_file(char *filename, char *errbuf, char **ret_mask, int *ret_masklen, int *ret_mask_has_internal_zeroes);
static void PairCount(const ESL_ALPHABET *abc, double *counters, ESL_DSQ syml, ESL_DSQ symr, double wt);
static int  get_command(const ESL_GETOPTS *go, char *errbuf, char **ret_command);
static int  get_date(char *errbuf, char **ret_date);
static int  set_scheme_values(char *errbuf, float *vec, int ncolvals, float **scheme, float val, SchemeColorLegend_t *scl, int within_mask, int *ret_bi);
static int  set_onecell_values(char *errbuf, float *vec, int ncolvals, float *onecolor);
static int  add_mask_to_ss_postscript(SSPostscript_t *ps, char *mask);
static int  draw_masked_block(FILE *fp, float x, float y, float *colvec, int do_circle_mask, int do_square_mask, int do_x_mask, int do_border, float boxsize);
static int  draw_header_and_footer(FILE *fp, const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, int page, int pageidx2print);
static int  read_seq_list_file_bigmem  (char *filename, ESL_MSA *msa, int **ret_useme, int *ret_nused);
static int  read_seq_list_file_smallmem(char *filename, ESL_KEYHASH **ret_useme_keyhash, int *ret_nused);
static void get_insert_info_from_msa(ESL_MSA *msa, int rflen, int **ret_nseq_with_ins_ct, int ***ret_ins_ct);
static void get_insert_info_from_abc_ct(double **abc_ct, ESL_ALPHABET *abc, char *msa_rf, int64_t msa_alen, int rflen, int **ret_nseq_with_ins_ct);
static void get_insert_info_from_ifile(char *ifile, int rflen, int msa_nseq, ESL_KEYHASH *useme_keyhash, int **ret_nseq_with_ins_ct, int ***ret_ins_ct);
static int  count_msa(ESL_MSA *msa, char *errbuf, int *a2rf_map, int rflen, double ***ret_abc_ct, double ****ret_bp_ct, int ***ret_pp_ct, int **ret_srfpos_ct, int **ret_erfpos_ct);
static int  infocontent_sspostscript(const ESL_GETOPTS *go, ESL_ALPHABET *abc, char *errbuf, SSPostscript_t *ps, double **abc_ct, int msa_nseq, float ***hc_scheme, int hc_scheme_idx, int hc_nbins, float **hc_onecell, int hc_onecell_idx, FILE *tabfp);
static int  mutual_information_sspostscript(const ESL_GETOPTS *go, ESL_ALPHABET *abc, char *errbuf, SSPostscript_t *ps, double ***bp_ct, int msa_nseq, float ***hc_scheme, int hc_scheme_idx, int hc_nbins, float **hc_onecell, int ss_idx, int zerores_idx, FILE *tabfp);
static int  delete_sspostscript(const ESL_GETOPTS *go, ESL_ALPHABET *abc, char *errbuf, SSPostscript_t *ps, double **abc_ct, int *srfpos_ct, int *erfpos_ct, int msa_nseq, int do_all, float ***hc_scheme, int hc_scheme_idx, int hc_nbins, float **hc_onecell, int hc_onecell_idx, FILE *tabfp);
static int  avg_posteriors_sspostscript(const ESL_GETOPTS *go, ESL_ALPHABET *abc, char *errbuf, SSPostscript_t *ps, int **pp_ct, int msa_nseq, float ***hc_scheme, int hc_scheme_idx, int hc_nbins, float **hc_onecell, int hc_onecell_idx, FILE *tabfp);
static int  insert_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, int *nseq_with_ins_ct, int msa_nseq, float ***hc_scheme, int hc_scheme_idx, int hc_nbins, float **hc_onecell, int hc_onecell_idx, FILE *tabfp);
static int  span_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, int *srfpos_ct, int *erfpos_ct, int msa_nseq, float ***hc_scheme, int hc_scheme_idx, int hc_nbins, float **hc_onecell, int zercov_idx, int maxcov_idx, FILE *tabfp);
static int  individual_seqs_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, int **ins_ct, int *useme, int nused, float ***hc_scheme, int hc_scheme_idx, int hc_nbins, float **hc_onecell, int zeroins_idx, int extdel_idx);
static int  rf_seq_sspostscript (const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa);
static int  individual_posteriors_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, int *useme, int nused, float ***hc_scheme, int hc_scheme_idx, int hc_nbins, float **hc_onecell, int hc_onecell_idx);
static int  colormask_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, char *mask, float **hc_onecell, int incmask_idx, int excmask_idx);
static int  diffmask_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, char *mask1, char *mask2, float **hc_onecell, int incboth_idx, int inc1_idx, int inc2_idx, int excboth_idx);
static int  drawfile2sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps);
static int  get_pp_idx(ESL_ALPHABET *abc, char ppchar);
static int  get_span_ct(int rflen, int nseq, int *srfpos_ct, int *erfpos_ct, int **ret_span_ct);

static char banner[] = "draw postscript secondary structure diagrams.";
static char usage[]  = "[options] <msafile> <SS postscript template> <output postscript file name>\n\
The <msafile> must be in Stockholm format.";

#define MASKTYPEOPTS "-d,-c,-x" /* exclusive choice for mask types */
#define INCOMPATWITHSINGLEOPTS "--prob,--ins,--dall,--dint,--mutinfo,--indi,--list,--all,--dfile" /* incompatibile with --mask and --mask-col */
#define INCOMPATWITHDFILEOPTS "-q,--prob,--ins,--dall,--dint,--mutinfo,--indi,--list,--all,--mask-col,--mask-diff" /* incompatible with --dfile */

static ESL_OPTIONS options[] = {
  /* name       type        default env   range togs  reqs  incomp      help                                                   docgroup */
  { "-h",       eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "help; show brief info on version and usage",              1 },
  { "-q",       eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "do not draw info content diagram (or RF sequence if --indi)", 1 },
  { "--mask",   eslARG_INFILE, NULL, NULL, NULL, NULL,NULL, NULL,            "for all diagrams, mark masked ('0') columns from mask in <f>", 1 },
  { "--prob",   eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "draw posterior probability diagram(s)", 1 },
  { "--small",  eslARG_NONE,   NULL, NULL, NULL, NULL,NULL, "--all",         "operate in small memory mode (aln must be 1 line/seq Pfam format)", 1 },

  { "--ins",    eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL,"--indi",         "draw insert diagram", 2 },
  { "--dall",   eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL,"--indi",         "draw delete diagram w/all deletions (incl. terminal deletes)", 2 },
  { "--dint",   eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL,"--indi",         "draw delete diagram w/only internal (non-terminal) deletions", 2 },
  { "--mutinfo",eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL,"--indi",         "draw base pair mutual information diagram", 2 },
  { "--span",   eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL,"--indi",         "draw diagram showing fraction of seqs that span each posn", 2 },

  { "--indi",   eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "draw diagrams for individual sequences instead of the aln", 3 },
  { "--list",   eslARG_INFILE,FALSE, NULL, NULL, NULL,"--indi", "--all",     "with --indi, draw individual diagrams for seqs listed in <f>", 3 },
  { "--all",    eslARG_NONE,  FALSE, NULL, NULL, NULL,"--indi", "--list,--small","with --indi, draw individual diagrams of all sequences", 3 },
  { "--keep-list",eslARG_OUTFILE,FALSE,NULL,NULL,NULL,"--indi,--list,--small","--all",   "w/--list,--indi & --small, save aln of seqs in list to <f>", 3 },

  { "--mask-u", eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "with --mask, mark masked columns as squares", 4 },
  { "--mask-x", eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "with --mask, mark masked columns as x's", 4 },
  { "--mask-a", eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,            "with --mask-u or --mask-x, draw alternative mask style", 4 },

  { "--mask-col",eslARG_NONE, NULL, NULL, NULL, NULL,"--mask",  INCOMPATWITHSINGLEOPTS, "w/--mask draw black/cyan diagram denoting masked columns", 5 },
  { "--mask-diff",eslARG_INFILE,NULL, NULL, NULL, NULL,"--mask",INCOMPATWITHSINGLEOPTS, "with --mask-col <f1>, compare mask in <f1> to mask in <f>", 5 },

  { "--dfile",   eslARG_INFILE, NULL, NULL, NULL, NULL,NULL, INCOMPATWITHDFILEOPTS, "read 'draw' file specifying >=1 diagrams", 6 },
  { "--ifile",   eslARG_INFILE, NULL, NULL, NULL, NULL,NULL, NULL,            "read insert information from cmalign insert file <f>", 6 },

  { "--tabfile", eslARG_OUTFILE, NULL, NULL, NULL, NULL, NULL, "--indi",  "output per position data in tabular format to file <f>", 7 },

  { "--no-leg", eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,          "do not draw legend", 8 },
  { "--no-head",eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,          "do not draw header", 8 },
  { "--no-foot",eslARG_NONE,  FALSE, NULL, NULL, NULL,NULL, NULL,          "do not draw footer", 8 },
  { 0,0,0,0,0,0,0,0,0,0 },
};

int
main(int argc, char **argv)
{
  ESL_GETOPTS    *go      = NULL;	/* application configuration       */
  ESL_ALPHABET   *abc     = NULL;	/* biological alphabet             */
  char           *alifile = NULL;	/* alignment file name             */
  char           *outfile = NULL;	/* output ps file name             */
  char           *templatefile = NULL;  /* template file, specifying >= 1 SS diagrams 
		  		         * (each must have a unique consensus length) */
  int             fmt;		        /* format code for alifile         */
  ESL_MSAFILE    *afp     = NULL;	/* open alignment file             */
  ESL_MSA        *msa     = NULL;	/* one multiple sequence alignment */
  int             status;		/* easel return code               */
  int             rflen;                 /* non-gap RF (consensus) length of each alignment */
  int             apos;                 /* counter over alignment positions */
  char            errbuf[ERRBUFSIZE];   /* for printing error messages */
  FILE           *ofp       = NULL;     /* output file for postscript */
  FILE           *tabfp    = NULL;     /* output file for text data, only set to non-null if --tabfile enabled */
  SSPostscript_t *ps        = NULL;     /* the postscript data structure we create */
  int            *hc_nbins  = NULL;     /* hard-coded number of bins (colors) for each possible scheme */
  float        ***hc_scheme = NULL;     /* colors for each pre-defined scheme */
  float         **hc_onecell = NULL;    /* single colors used for binary color schemes */
  int             z;                    /* counter */
  char           *command = NULL;       /* string for printing command to stdout */
  char           *date = NULL;          /* date command was executed */
  int             master_mode;          /* which mode we're operating in */ 
  char           *mask = NULL;          /* first mask we'll use, string of '0's and '1's */
  int             masklen;              /* length of mask */
  char           *mask2 = NULL;         /* second mask we'll use, string of '0's and '1's */
  int             mask2len;             /* length of second mask */
  int             mask_has_internal_zeroes = FALSE;  /* does mask have any '0's with at least '1' on both sides? */
  int             mask2_has_internal_zeroes = FALSE; /* does mask have any '0's with at least '1' on both sides? */
  int            *useme = NULL;         /* only relevant if --list, [0..i..msa->nseq] TRUE to include indi diagram for seq i, FALSE not to */
  int             nused = 0;            /* only relevant if --list, number of TRUEs in useme */
 /* counts of relevant values from the msa */
  int            *nseq_with_ins_ct = NULL; /* [0..ps->rflen] number of sequences with >=1 insert after each consensus position, only used if --ins */
  int           **ins_ct = NULL;        /* [0..msa->nseq-1][0..ps->rflen] for each sequence, the number of inserts after each consensus position */
  double        **abc_ct = NULL;        /* [0..msa->alen-1][0..abc->K], count of each residue at each position, over all sequences, missing and nonresidues are *not counted* */
  double       ***bp_ct = NULL;         /* [0..msa->alen-1][0..abc->Kp][0..abc->Kp], count of each possible base pair at each position, over all sequences, missing and nonresidues are *not counted* 
                                           base pairs are indexed by 'i' for a base pair between positions i and j, where i < j. */
  int           **pp_ct = NULL;         /* [0..msa->alen-1][0..11], count of reach posterior probability (PP) code, over all sequences, gap is 11 */
  int            *srfpos_ct = NULL;     /* [0..msa->alen-1] per position count of first non-gap position, over all seqs */
  int            *erfpos_ct = NULL;     /* [0..msa->alen-1] per position count of final non-gap position, over all seqs */
  int64_t         msa_alen;             /* msa->alen */
  int             msa_nseq;             /* msa->nseq */
  /* variables related to small memory mode */
  int             do_small = TRUE;      /* TRUE to operate in small memory mode, FALSE not to */
  /* small memory mode variables used only if --list and --indi */
  char            tmp_indi_alifile[32] = "esltmpXXXXXX"; /* the name of the indi alignment file */
  ESL_MSAFILE    *indi_afp = NULL;      /* MSA file pointer for newly created indi alignment in */
  FILE           *indi_fp = NULL;       /* file pointer for outputting indi alignment */
  ESL_MSA        *indi_msa = NULL;      /* new indi msa */
  int           **indi_ins_ct = NULL;   /* [0..indi_msa->nseq-1][0..ps->rflen] for each sequence, the number of inserts after each consensus position */
  ESL_KEYHASH    *useme_keyhash;        /* keyhash of sequence names listed in list file, only used if --list and --indi enabled */
  int             i;                    /* counter of sequences */

  /***********************************************
   * Parse command line
   ***********************************************/

  go = esl_getopts_Create(options);
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK ||
      esl_opt_VerifyConfig(go)               != eslOK)
    {
      printf("Failed to parse command line: %s\n", go->errbuf);
      esl_usage(stdout, argv[0], usage);
      printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
      exit(1);
    }

  if (esl_opt_GetBoolean(go, "-h") )
    {
      esl_banner(stdout, argv[0], banner);
      esl_usage (stdout, argv[0], usage);
      puts("\n where basic options are:");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
      puts("\noptions for alignment summary diagrams (incompatible with --indi):");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80); 
      puts("\noptions for individual mode (require --indi):");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80); 
      puts("\noptions controlling style of masked positions:");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80); 
      puts("\noptions for drawing simple block-only diagrams of masks:");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 80); 
      puts("\noptions related to optional input files:");
      esl_opt_DisplayHelp(stdout, go, 6, 2, 80); 
      puts("\noptions related to optional output files:");
      esl_opt_DisplayHelp(stdout, go, 7, 2, 80); 
      puts("\noptions for omitting parts of the diagram:");
      esl_opt_DisplayHelp(stdout, go, 8, 2, 80); 
      exit(0);
    }

  if (esl_opt_ArgNumber(go) != 3) 
    {
      printf("Incorrect number of command line arguments.\n");
      esl_usage(stdout, argv[0], usage);
      printf("\nTo see more help on available options, do %s -h\n\n", argv[0]);
      exit(1);
    }

  alifile      = esl_opt_GetArg(go, 1);
  templatefile = esl_opt_GetArg(go, 2);
  outfile      = esl_opt_GetArg(go, 3);

  if((status = get_command(go, errbuf, &command)) != eslOK) esl_fatal(errbuf);
  if((status = get_date(errbuf, &date))           != eslOK) esl_fatal(errbuf);

  /****************/
  /* Premlinaries */
  /****************/
  /* allocate and fill predefined one-cell colors, these are hardcoded */
  ESL_ALLOC(hc_onecell, sizeof(float *) * NOC);
  for(z = 0; z < NOC; z++) hc_onecell[z] = NULL;
  for(z = 0; z < NOC; z++) { ESL_ALLOC(hc_onecell[z], sizeof(float) * NCMYK); }
  
  hc_onecell[CYANOC][0] = 1.0;
  hc_onecell[CYANOC][1] = 0.0;
  hc_onecell[CYANOC][2] = 0.0;
  hc_onecell[CYANOC][3] = 0.0;

  hc_onecell[MAGENTAOC][0] = 0.0;
  hc_onecell[MAGENTAOC][1] = 1.0;
  hc_onecell[MAGENTAOC][2] = 0.0;
  hc_onecell[MAGENTAOC][3] = 0.0;

  hc_onecell[YELLOWOC][0] = 0.0;
  hc_onecell[YELLOWOC][1] = 0.0;
  hc_onecell[YELLOWOC][2] = 1.0;
  hc_onecell[YELLOWOC][3] = 0.0;

  hc_onecell[BLACKOC][0] = 0.0;
  hc_onecell[BLACKOC][1] = 0.0;
  hc_onecell[BLACKOC][2] = 0.0;
  hc_onecell[BLACKOC][3] = 1.0;

  hc_onecell[LIGHTGREYOC][0] = 0.0;
  hc_onecell[LIGHTGREYOC][1] = 0.0;
  hc_onecell[LIGHTGREYOC][2] = 0.0;
  hc_onecell[LIGHTGREYOC][3] = 0.2;

  hc_onecell[DARKGREYOC][0] = 0.0;
  hc_onecell[DARKGREYOC][1] = 0.0;
  hc_onecell[DARKGREYOC][2] = 0.0;
  hc_onecell[DARKGREYOC][3] = 0.5;

  hc_onecell[REDOC][0] = 0.0;
  hc_onecell[REDOC][1] = 1.0;
  hc_onecell[REDOC][2] = 1.0;
  hc_onecell[REDOC][3] = 0.0;

  hc_onecell[PURPLEOC][0] = 1.0;
  hc_onecell[PURPLEOC][1] = 1.0;
  hc_onecell[PURPLEOC][2] = 0.0;
  hc_onecell[PURPLEOC][3] = 0.0;

  hc_onecell[ORANGEOC][0] = 0.0;
  hc_onecell[ORANGEOC][1] = 0.5;
  hc_onecell[ORANGEOC][2] = 1.0;
  hc_onecell[ORANGEOC][3] = 0.0;

  hc_onecell[WHITEOC][0] = 0.0;
  hc_onecell[WHITEOC][1] = 0.0;
  hc_onecell[WHITEOC][2] = 0.0;
  hc_onecell[WHITEOC][3] = 0.0;

  /***********************************/
  /* allocate and fill predefined color schemes, these are hardcoded */
  ESL_ALLOC(hc_scheme, sizeof(float **) * 6);
  for (z = 0; z < 6; z++) hc_scheme[z] = NULL;
  ESL_ALLOC(hc_scheme[0], sizeof(float *) * 11); 
  for(z = 0; z < 11; z++) hc_scheme[0][z] = NULL;
  for(z = 0; z < 11; z++) { ESL_ALLOC(hc_scheme[0][z], sizeof(float) * NCMYK); }
  ESL_ALLOC(hc_scheme[1], sizeof(float *) * 11); 
  for(z = 0; z < 11; z++) hc_scheme[1][z] = NULL;
  for(z = 0; z < 11; z++) { ESL_ALLOC(hc_scheme[1][z], sizeof(float) * NCMYK); }
  ESL_ALLOC(hc_scheme[2], sizeof(float *) * 6); 
  for(z = 0; z < 6; z++) hc_scheme[2][z] = NULL;
  for(z = 0; z < 6; z++) { ESL_ALLOC(hc_scheme[2][z], sizeof(float) * NCMYK); }
  ESL_ALLOC(hc_scheme[3], sizeof(float *) * 6); 
  for(z = 0; z < 6; z++) hc_scheme[3][z] = NULL;
  for(z = 0; z < 6; z++) { ESL_ALLOC(hc_scheme[3][z], sizeof(float) * NCMYK); }
  ESL_ALLOC(hc_scheme[4], sizeof(float *) * 5); 
  for(z = 0; z < 5; z++) hc_scheme[4][z] = NULL;
  for(z = 0; z < 5; z++) { ESL_ALLOC(hc_scheme[4][z], sizeof(float) * NCMYK); }
  ESL_ALLOC(hc_scheme[5], sizeof(float *) * 5); 
  for(z = 0; z < 5; z++) hc_scheme[5][z] = NULL;
  for(z = 0; z < 5; z++) { ESL_ALLOC(hc_scheme[5][z], sizeof(float) * NCMYK); }

  ESL_ALLOC(hc_nbins, sizeof(int) * 6);
  hc_nbins[0] = NRAINBOWRHSCHEME;
  hc_nbins[1] = NRAINBOWRLSCHEME;
  hc_nbins[2] = NRBSIXRHSCHEME;
  hc_nbins[3] = NRBSIXRLSCHEME;
  hc_nbins[4] = NRBFIVERHSCHEME;
  hc_nbins[5] = NRBFIVERLSCHEME;

  /***********************************/
  /*Scheme 0 and 1: Rainbow(red high) 11 is 0, Rainbow (red low) 11 is 1 */
  hc_scheme[0][ 0][0] = 0.92; hc_scheme[0][ 0][1] = 0.84; hc_scheme[0][ 0][2] = 0.00; hc_scheme[0][ 0][3] = 0.08; /*blue*/
  hc_scheme[1][10][0] = 0.92; hc_scheme[1][10][1] = 0.84; hc_scheme[1][10][2] = 0.00; hc_scheme[1][10][3] = 0.08; /*blue*/

  hc_scheme[0][ 1][0] = 0.78; hc_scheme[0][ 1][1] = 0.56; hc_scheme[0][ 1][2] = 0.00; hc_scheme[0][ 1][3] = 0.22;
  hc_scheme[1][ 9][0] = 0.78; hc_scheme[1][ 9][1] = 0.56; hc_scheme[1][ 9][2] = 0.00; hc_scheme[1][ 9][3] = 0.22;

  hc_scheme[0][ 2][0] = 0.50; hc_scheme[0][ 2][1] = 0.00; hc_scheme[0][ 2][2] = 0.00; hc_scheme[0][ 2][3] = 0.50;
  hc_scheme[1][ 8][0] = 0.50; hc_scheme[1][ 8][1] = 0.00; hc_scheme[1][ 8][2] = 0.00; hc_scheme[1][ 8][3] = 0.50;

  hc_scheme[0][ 3][0] = 0.61; hc_scheme[0][ 3][1] = 0.00; hc_scheme[0][ 3][2] = 0.56; hc_scheme[0][ 3][3] = 0.22;
  hc_scheme[1][ 7][0] = 0.61; hc_scheme[1][ 7][1] = 0.00; hc_scheme[1][ 7][2] = 0.56; hc_scheme[1][ 7][3] = 0.22;

  hc_scheme[0][ 4][0] = 0.42; hc_scheme[0][ 4][1] = 0.00; hc_scheme[0][ 4][2] = 1.00; hc_scheme[0][ 4][3] = 0.00;
  hc_scheme[1][ 6][0] = 0.42; hc_scheme[1][ 6][1] = 0.00; hc_scheme[1][ 6][2] = 1.00; hc_scheme[1][ 6][3] = 0.00;

  hc_scheme[0][ 5][0] = 0.00; hc_scheme[0][ 5][1] = 0.00; hc_scheme[0][ 5][2] = 1.00; hc_scheme[0][ 5][3] = 0.00;
  hc_scheme[1][ 5][0] = 0.00; hc_scheme[1][ 5][1] = 0.00; hc_scheme[1][ 5][2] = 1.00; hc_scheme[1][ 5][3] = 0.00;

  hc_scheme[0][ 6][0] = 0.00; hc_scheme[0][ 6][1] = 0.21; hc_scheme[0][ 6][2] = 1.00; hc_scheme[0][ 6][3] = 0.00;
  hc_scheme[1][ 4][0] = 0.00; hc_scheme[1][ 4][1] = 0.21; hc_scheme[1][ 4][2] = 1.00; hc_scheme[1][ 4][3] = 0.00;

  hc_scheme[0][ 7][0] = 0.00; hc_scheme[0][ 7][1] = 0.42; hc_scheme[0][ 7][2] = 1.00; hc_scheme[0][ 7][3] = 0.00;
  hc_scheme[1][ 3][0] = 0.00; hc_scheme[1][ 3][1] = 0.42; hc_scheme[1][ 3][2] = 1.00; hc_scheme[1][ 3][3] = 0.00;

  hc_scheme[0][ 8][0] = 0.00; hc_scheme[0][ 8][1] = 0.63; hc_scheme[0][ 8][2] = 1.00; hc_scheme[0][ 8][3] = 0.00;
  hc_scheme[1][ 2][0] = 0.00; hc_scheme[1][ 2][1] = 0.63; hc_scheme[1][ 2][2] = 1.00; hc_scheme[1][ 2][3] = 0.00;

  hc_scheme[0][ 9][0] = 0.00; hc_scheme[0][ 9][1] = 0.84; hc_scheme[0][ 9][2] = 1.00; hc_scheme[0][ 9][3] = 0.00;
  hc_scheme[1][ 1][0] = 0.00; hc_scheme[1][ 1][1] = 0.84; hc_scheme[1][ 1][2] = 1.00; hc_scheme[1][ 1][3] = 0.00;

  hc_scheme[0][10][0] = 0.00; hc_scheme[0][10][1] = 0.94; hc_scheme[0][10][2] = 1.00; hc_scheme[0][10][3] = 0.00; /*red*/
  hc_scheme[1][ 0][0] = 0.00; hc_scheme[1][ 0][1] = 0.94; hc_scheme[1][ 0][2] = 1.00; hc_scheme[1][ 0][3] = 0.00; /*red*/
  /***********************************/
  /*Scheme 2 and 3: Rainbow(red high) 6 is 2, Rainbow (red low) 6 is 3 */
  hc_scheme[2][0][0] = 0.92; hc_scheme[2][0][1] = 0.84; hc_scheme[2][0][2] = 0.00; hc_scheme[2][0][3] = 0.08; /*blue*/
  hc_scheme[3][5][0] = 0.92; hc_scheme[3][5][1] = 0.84; hc_scheme[3][5][2] = 0.00; hc_scheme[3][5][3] = 0.08; /*blue*/

  hc_scheme[2][1][0] = 0.50; hc_scheme[2][1][1] = 0.00; hc_scheme[2][1][2] = 0.00; hc_scheme[2][1][3] = 0.50;
  hc_scheme[3][4][0] = 0.50; hc_scheme[3][4][1] = 0.00; hc_scheme[3][4][2] = 0.00; hc_scheme[3][4][3] = 0.50;

  hc_scheme[2][2][0] = 0.42; hc_scheme[2][2][1] = 0.00; hc_scheme[2][2][2] = 1.00; hc_scheme[2][2][3] = 0.00;
  hc_scheme[3][3][0] = 0.42; hc_scheme[3][3][1] = 0.00; hc_scheme[3][3][2] = 1.00; hc_scheme[3][3][3] = 0.00;

  hc_scheme[2][3][0] = 0.00; hc_scheme[2][3][1] = 0.21; hc_scheme[2][3][2] = 1.00; hc_scheme[2][3][3] = 0.00;
  hc_scheme[3][2][0] = 0.00; hc_scheme[3][2][1] = 0.21; hc_scheme[3][2][2] = 1.00; hc_scheme[3][2][3] = 0.00;

  hc_scheme[2][4][0] = 0.00; hc_scheme[2][4][1] = 0.63; hc_scheme[2][4][2] = 1.00; hc_scheme[2][4][3] = 0.00;
  hc_scheme[3][1][0] = 0.00; hc_scheme[3][1][1] = 0.63; hc_scheme[3][1][2] = 1.00; hc_scheme[3][1][3] = 0.00;

  hc_scheme[2][5][0] = 0.00; hc_scheme[2][5][1] = 0.94; hc_scheme[2][5][2] = 1.00; hc_scheme[2][5][3] = 0.00; /*red*/
  hc_scheme[3][0][0] = 0.00; hc_scheme[3][0][1] = 0.94; hc_scheme[3][0][2] = 1.00; hc_scheme[3][0][3] = 0.00; /*red*/
  /***************************************************************/
  /***********************************/
  /*Scheme 4 and 5: Rainbow(red high) 5 is 4, Rainbow (red low) 5 is 5 */
  /*Same as schemes 3 and 4 (rainbow 6s, except no final blue, which makes black text difficult overlaid on it difficult to read */
  hc_scheme[4][0][0] = 0.50; hc_scheme[4][0][1] = 0.00; hc_scheme[4][0][2] = 0.00; hc_scheme[4][0][3] = 0.50; /*teal*/
  hc_scheme[5][4][0] = 0.50; hc_scheme[5][4][1] = 0.00; hc_scheme[5][4][2] = 0.00; hc_scheme[5][4][3] = 0.50; /*teal*/

  hc_scheme[4][1][0] = 0.42; hc_scheme[4][1][1] = 0.00; hc_scheme[4][1][2] = 1.00; hc_scheme[4][1][3] = 0.00;
  hc_scheme[5][3][0] = 0.42; hc_scheme[5][3][1] = 0.00; hc_scheme[5][3][2] = 1.00; hc_scheme[5][3][3] = 0.00;

  hc_scheme[4][2][0] = 0.00; hc_scheme[4][2][1] = 0.21; hc_scheme[4][2][2] = 1.00; hc_scheme[4][2][3] = 0.00;
  hc_scheme[5][2][0] = 0.00; hc_scheme[5][2][1] = 0.21; hc_scheme[5][2][2] = 1.00; hc_scheme[5][2][3] = 0.00;

  hc_scheme[4][3][0] = 0.00; hc_scheme[4][3][1] = 0.63; hc_scheme[4][3][2] = 1.00; hc_scheme[4][3][3] = 0.00;
  hc_scheme[5][1][0] = 0.00; hc_scheme[5][1][1] = 0.63; hc_scheme[5][1][2] = 1.00; hc_scheme[5][1][3] = 0.00;

  hc_scheme[4][4][0] = 0.00; hc_scheme[4][4][1] = 0.94; hc_scheme[4][4][2] = 1.00; hc_scheme[4][4][3] = 0.00; /*red*/
  hc_scheme[5][0][0] = 0.00; hc_scheme[5][0][1] = 0.94; hc_scheme[5][0][2] = 1.00; hc_scheme[5][0][3] = 0.00; /*red*/
  /***************************************************************/
  master_mode = ALIMODE;
  if(esl_opt_GetBoolean (go, "--indi"))      master_mode = INDIMODE;
  if(esl_opt_GetBoolean (go, "--mask-col"))  master_mode = SIMPLEMASKMODE;
  if(! esl_opt_IsDefault(go, "--mask-diff")) master_mode = SIMPLEMASKMODE;
  if(! esl_opt_IsDefault(go, "--dfile"))     master_mode = DRAWFILEMODE;

  /*****************************************
   * Open the MSA file; determine alphabet; 
   *****************************************/
  
  do_small = esl_opt_GetBoolean(go, "--small") ? TRUE : FALSE;
  fmt = do_small ? eslMSAFILE_PFAM : eslMSAFILE_STOCKHOLM;
  status = esl_msafile_Open(alifile, fmt, NULL, &afp);
  if      (status == eslENOTFOUND) esl_fatal("Alignment file %s doesn't exist or is not readable\n", alifile);
  else if (status == eslEFORMAT)   esl_fatal("Couldn't determine format of alignment %s\n", alifile);
  else if (status != eslOK)        esl_fatal("Alignment file open failed with error %d\n", status);

  /* Assert RNA, it's the ribosome */
  abc = esl_alphabet_Create(eslRNA);

  /* Read the mask files, if nec */
  if(! esl_opt_IsDefault(go, "--mask")) { 
    if((status = read_mask_file(esl_opt_GetString(go, "--mask"), errbuf, &mask, &masklen, &mask_has_internal_zeroes)) != eslOK) esl_fatal(errbuf);
  }
  if(! esl_opt_IsDefault(go, "--mask-diff")) { 
    if((status = read_mask_file(esl_opt_GetString(go, "--mask-diff"), errbuf, &mask2, &mask2len, &mask2_has_internal_zeroes)) != eslOK) esl_fatal(errbuf);
    if(masklen != mask2len) esl_fatal("Mask in %f length (%d) differs from mask in %f (%d)!", esl_opt_GetString(go, "--mask"), masklen, esl_opt_GetString(go, "--mask-diff"), mask2len);
  }

  /* Open output files, if necessary */
  if (esl_opt_IsUsed(go, "--tabfile")) { 
    if((tabfp = fopen(esl_opt_GetString(go, "--tabfile"), "w")) == NULL) esl_fatal("Failed to open output file %s\n", esl_opt_GetString(go, "--tabfile"));
  }

  /**********************
   * Read the alignment *
   **********************/
  status = (do_small) ? 
    esl_msa_ReadNonSeqInfoPfam(afp, abc, -1, NULL, NULL, &msa, &msa_nseq, &msa_alen, NULL, NULL, NULL, NULL, NULL, &abc_ct, &pp_ct, NULL, NULL, NULL) :
    esl_msa_Read              (afp, &msa); /* if ! do_small, we read full aln into memory */
  if      (status == eslEFORMAT) esl_fatal("Alignment file parse error:\n%s\n", afp->errbuf);
  else if (status == eslEINVAL)  esl_fatal("Alignment file parse error:\n%s\n", afp->errbuf);
  else if (status == eslEOF)     esl_fatal("No alignments found in file %s\n", alifile);
  else if (status != eslOK)      esl_fatal("Alignment file read failed with error code %d\n%s", status, afp);

  if(do_small) { 
    msa->alen = msa_alen; }
  else         
    { 
    msa_nseq = msa->nseq; 
    msa_alen = msa->alen; 
  }

  if(do_small && (esl_opt_GetBoolean(go, "--dint") || esl_opt_GetBoolean(go, "--span") || esl_opt_GetBoolean(go, "--mutinfo"))) { 
    /* If we're in small mem mode, now that we know msa->rf and msa->ss_cons, 
     * we do a second read of the msa to get bpct, srfpos_ct and erfpos_ct, 
     * but only if we need them (if --span, --dint, or --mutinfo). 
     * To do this, we close the alifile and reopen it, so we can read the 1st
     * alignment again.
     */
    esl_msafile_Close(afp);
    status = esl_msafile_Open(alifile, fmt, NULL, &afp);
    if      (status == eslENOTFOUND) esl_fatal("2nd pass, alignment file %s doesn't exist or is not readable\n", alifile);
    else if (status == eslEFORMAT)   esl_fatal("2nd pass, couldn't determine format of alignment %s\n", alifile);
    else if (status != eslOK)        esl_fatal("2nd pass, alignment file open failed with error %d\n", status);
    
    status = esl_msa_ReadNonSeqInfoPfam(afp, abc, msa_alen, msa->rf, msa->ss_cons, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, &bp_ct, &srfpos_ct, &erfpos_ct);
    if      (status == eslEFORMAT) esl_fatal("2nd pass, Alignment file parse error:\n%s\n", afp->errbuf);
    else if (status == eslEINVAL)  esl_fatal("2nd pass, Alignment file parse error:\n%s\n", afp->errbuf);
    else if (status == eslEOF)     esl_fatal("2nd pass, No alignments found in file %s\n", alifile);
    else if (status != eslOK)      esl_fatal("2nd pass, Alignment file read failed with error code %d\n%s", status, afp);
  }

  msa->abc = abc;
  if(msa->rf == NULL) esl_fatal("First MSA in %s does not have RF annotation.", alifile);
  /* determine non-gap RF length (consensus length) */
  rflen = 0;
  for(apos = 0; apos < msa->alen; apos++) { 
    if((! esl_abc_CIsGap(msa->abc, msa->rf[apos])) && 
       (! esl_abc_CIsMissing(msa->abc, msa->rf[apos])) && 
       (! esl_abc_CIsNonresidue(msa->abc, msa->rf[apos]))) { 
      rflen++;
    }
  }
  /* We've read the alignment, now read the template postscript file (we do this second b/c the RF len of the alignment tells us which postscript template to use) */
  if((status = parse_template_file(templatefile, go, errbuf, rflen, &ps) != eslOK)) esl_fatal(errbuf);
  
  /* if we get here, the postscript file has been successfully read; now we open the output file */
  if(ofp == NULL) { 
    /* open postscript output file for writing */
    if ((ofp = fopen(outfile, "w")) == NULL)
      ESL_FAIL(eslFAIL, errbuf, "Failed to open output postscript file %s\n", esl_opt_GetArg(go, 2));
  }
  
  /* determine position for header and legend */
  if((status = setup_sspostscript(ps, errbuf) != eslOK)) esl_fatal(errbuf);
  if(ps->rflen == 0)     esl_fatal("MSA has consensus (non-gap RF) length of %d which != template file consensus length of %d.", rflen, ps->rflen);
  if(rflen != ps->rflen) esl_fatal("MSA has consensus (non-gap RF) length of %d which != template file consensus length of %d.", rflen, ps->rflen);
  
  /* add the mask if there is one */
  if(mask != NULL && (master_mode != SIMPLEMASKMODE)) add_mask_to_ss_postscript(ps, mask);
  if(mask != NULL && ps->rflen != masklen) esl_fatal("MSA has consensus (non-gap RF) length of %d which != lane mask length of %d from mask file %s.", rflen, masklen, esl_opt_GetString(go, "--mask"));
  
  if((status = validate_and_update_sspostscript_given_msa(go, ps, msa, msa_nseq, errbuf)) != eslOK) esl_fatal(errbuf);

  /* get the information we need from the alignment and create the SS postscript files as we go */
  if(master_mode == ALIMODE) { 
    if(! do_small) { /* derive counts from the msa for postscript diagrams, we do this our functions for drawing diagrams work in small mem or big mem mode */
      if((status = count_msa(msa, errbuf, ps->msa_a2rf_map, ps->rflen, &(abc_ct), &(bp_ct), (esl_opt_GetBoolean(go, "--prob") ? &(pp_ct) : NULL), &(srfpos_ct), &(erfpos_ct))) != eslOK) esl_fatal(errbuf);
    }    
    if(! esl_opt_GetBoolean(go, "-q")) { 
      if((status = infocontent_sspostscript(go, abc, errbuf, ps, abc_ct, msa_nseq, hc_scheme, RBSIXRLSCHEME, hc_nbins[RBSIXRLSCHEME], hc_onecell, LIGHTGREYOC, tabfp)) != eslOK) esl_fatal(errbuf);
    }
    if(esl_opt_GetBoolean(go, "--mutinfo")) { 
      if((status = mutual_information_sspostscript(go, abc, errbuf, ps, bp_ct, msa_nseq, hc_scheme, RBSIXRHSCHEME, hc_nbins[RBSIXRHSCHEME], hc_onecell, DARKGREYOC, LIGHTGREYOC, tabfp)) != eslOK) esl_fatal(errbuf);
    }
    if(esl_opt_GetBoolean(go, "--ins")) { /* make a new postscript page marking insertions */
      /* first, determine number of sequences with inserts after each position, 3 different ways */
      if(! esl_opt_IsDefault(go, "--ifile")) { /* read the insert file from cmalign */
	get_insert_info_from_ifile((esl_opt_GetString(go, "--ifile")), ps->rflen, msa_nseq, NULL, &(nseq_with_ins_ct), NULL); /* dies with esl_fatal() upon an error */
      }
      else if (do_small) { /* use abc_ct to derive nseq_with_ins_ct */
	get_insert_info_from_abc_ct(abc_ct, abc, msa->rf, msa_alen, ps->rflen, &(nseq_with_ins_ct)); /* dies with esl_fatal() upon an error */
      }
      else { 
	get_insert_info_from_msa(msa, ps->rflen, &(nseq_with_ins_ct), NULL); /* dies with esl_fatal() upon an error */
      }
      /* now draw the insert diagram */
      if((status = insert_sspostscript(go, errbuf, ps, nseq_with_ins_ct, msa_nseq, hc_scheme, RBSIXRHSCHEME, hc_nbins[RBSIXRHSCHEME], hc_onecell, LIGHTGREYOC, tabfp)) != eslOK) esl_fatal(errbuf);
    }
    if(esl_opt_GetBoolean(go, "--dall")) { /* make a new postscript page marking all deletes */
      if((status = delete_sspostscript(go, abc, errbuf, ps, abc_ct, srfpos_ct, erfpos_ct, msa_nseq, TRUE, hc_scheme, RBSIXRHSCHEME, hc_nbins[RBSIXRHSCHEME], hc_onecell, LIGHTGREYOC, tabfp)) != eslOK) esl_fatal(errbuf);
    }
    if(esl_opt_GetBoolean(go, "--dint")) { /* make a new postscript page marking internal deletes */
      if((status = delete_sspostscript(go, abc, errbuf, ps, abc_ct, srfpos_ct, erfpos_ct, msa_nseq, FALSE, hc_scheme, RBSIXRHSCHEME, hc_nbins[RBSIXRHSCHEME], hc_onecell, LIGHTGREYOC, tabfp)) != eslOK) esl_fatal(errbuf);
    }
    if(esl_opt_GetBoolean(go, "--prob")) { 
      if((status = avg_posteriors_sspostscript(go, abc, errbuf, ps, pp_ct, msa_nseq, hc_scheme, RBSIXRLSCHEME, hc_nbins[RBSIXRLSCHEME], hc_onecell, LIGHTGREYOC, tabfp)) != eslOK) esl_fatal(errbuf);
    }
    if(esl_opt_GetBoolean(go, "--span")) { 
      if((status = span_sspostscript(go, errbuf, ps, srfpos_ct, erfpos_ct, msa_nseq, hc_scheme, RBSIXRLSCHEME, hc_nbins[RBSIXRLSCHEME], hc_onecell, LIGHTGREYOC, BLACKOC, tabfp)) != eslOK) esl_fatal(errbuf);
    }
  }
  else if(master_mode == INDIMODE) { 
    if(! esl_opt_GetBoolean(go, "-q")) { 
      /* this will work in either small memory or normal memory mode */
      if((status = rf_seq_sspostscript(go, errbuf, ps, msa)) != eslOK) esl_fatal(errbuf);
    }
    /* if nec, draw individual sequence diagrams */
    if((esl_opt_GetBoolean(go, "--all")) || (! esl_opt_IsDefault(go, "--list"))) { 
      if(esl_opt_IsOn(go, "--all")) { 
	
	/* get insert info we'll use in individual_seqs_sspostscript() */
	get_insert_info_from_msa(msa, ps->rflen, NULL, &ins_ct); /* dies with esl_fatal() upon an error */

	if(do_small) esl_fatal("Operating in small memory mode, --all is not allowed.");
	ESL_ALLOC(useme, sizeof(int) * msa_nseq); 
	esl_vec_ISet(useme, msa_nseq, TRUE);
	nused = msa_nseq;
      }
      else if(esl_opt_IsOn(go, "--list")) {
	/* first read the list file */
	if(! do_small) { 
	  read_seq_list_file_bigmem(esl_opt_GetString(go, "--list"), msa, &useme, &nused);    /* this will die with esl_fatal() upon an error */
	}
	else { /* do_small == TRUE */
	  read_seq_list_file_smallmem(esl_opt_GetString(go, "--list"), &useme_keyhash, &nused); /* this will die with esl_fatal() upon an error */
	}

	/* At this point, we know which sequences we're going to draw indi diagrams for.
	 * Next step, get insert info. We need this so we can allow --ifile to work (which reads insert info from a file)
	 * instead of relying on always reading insert info from the msa itself.
	 * The way we get insert info varies, depending on if --small and --ifile. 
	 * If ! --small, we read all insert info, either from the msa (if ! --ifile) or from the ifile (if --ifile). 
	 * If   --small, we read only the insert info for the seqs we're going to draw diagrams for.
	 * When --small and ! --ifile, first we have to create the smaller alignment containing only those
	 * seqs we'll draw diagrams for, then we get the insert info from it.
	 */
	if(! do_small) { 
	  if(! esl_opt_IsDefault(go, "--ifile")) { /* read the insert file from cmalign, with info from all seqs  */
	    get_insert_info_from_ifile((esl_opt_GetString(go, "--ifile")), ps->rflen, msa_nseq, NULL, NULL, &(ins_ct)); /* dies with esl_fatal() upon an error */
	  }
	  else { /* get insert info from the msa */
	    get_insert_info_from_msa(msa, ps->rflen, NULL, &ins_ct); /* dies with esl_fatal() upon an error */
	  }
	}
	else { /* do_small */ 
	  if(! esl_opt_IsDefault(go, "--ifile")) { /* read the insert file and get insert info on only those seqs we'll draw diagrams for */
	    get_insert_info_from_ifile((esl_opt_GetString(go, "--ifile")), ps->rflen, msa_nseq, useme_keyhash, NULL, &(indi_ins_ct)); /* dies with esl_fatal() upon an error */
	  }
	  /* the 'else' half of this statement must wait until we've created indi_msa (see notes in comment block immediately
	   * below). We put the if() part above so that if there's an error in the ifile we find out before we do the 
	   * full msa regurgitation below.
	   */

	  /* We're in small memory mode, which means we never read the full alignment into memory. Instead we will just read
	   * the sequences that we're going to draw indi diagrams for into memory by creating a new msa: indi_msa.
	   * The way we do this is convoluted: open the msa file, and regurgitate it to a temp file, but taking care only to 
	   * regurgitate the seqs we want to draw diagrams for. Finally we read that temp file into memory as indi_msa.
	   */
	  esl_msafile_Close(afp);
	  status = esl_msafile_Open(alifile, fmt, NULL, &afp);

	  /* small memory mode, write a temporary alignment with only the sequences we want individual diagrams for */
	  if(esl_opt_IsOn(go, "--keep-list")) { 
	    if ((indi_fp  = fopen(esl_opt_GetString(go, "--keep-list"), "w")) == NULL) esl_fatal("Failed to open temporary output file %s for --indi and --list", esl_opt_GetString(go, "--keep-list"));
	  }
	  else { 
	    if ((esl_tmpfile_named(tmp_indi_alifile, &indi_fp)) != eslOK) esl_fatal("Failed to open temporary output file %s for --indi and --list");
	  }
	  if      (status == eslENOTFOUND) esl_fatal("Final pass, alignment file %s doesn't exist or is not readable\n", alifile);
	  else if (status == eslEFORMAT)   esl_fatal("Final pass, couldn't determine format of alignment %s\n", alifile);
	  else if (status != eslOK)        esl_fatal("Final pass, alignment file open failed with error %d\n", status);
	  status = esl_msa_RegurgitatePfam(afp, indi_fp, 
					   -1, -1, -1, -1, /* don't care about max width of fields */
					   TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, /* regurgitate all non-seq info */
					   useme_keyhash, /* only regurgitate seqs in useme_keyhash */
					   NULL,          /* no list of seqs to skip */
					   NULL, NULL, -1, '.'); 
	  if(status == eslEOF) esl_fatal("Final pass, no alignments in file");
	  if(status != eslOK)  esl_fatal("Final pass, error reading alignment");
	  fclose(indi_fp); 
	  indi_fp = NULL;
	  
	  /* now read in that small alignment */
	  if(esl_opt_IsOn(go, "--keep-list")) { 
	    status = esl_msafile_Open(esl_opt_GetString(go, "--keep-list"), fmt, NULL, &indi_afp);
	  }
	  else { 
	    status = esl_msafile_Open(tmp_indi_alifile, fmt, NULL, &indi_afp);
	  }
	  esl_msa_Read(indi_afp, &indi_msa); /* read the full thing into memory */
	  esl_msafile_Close(indi_afp);
	  indi_msa->abc = abc;

	  if(esl_opt_IsOn(go, "--keep-list")) { 
	    printf("# Alignment with the %d sequences from %s saved to file %s.\n", nused, esl_opt_GetString(go, "--list"), esl_opt_GetString(go, "--keep-list"));
	  }
	  else { 
	    remove(tmp_indi_alifile);
	  }

	  /* check to make sure all the sequences from the list file were in the alignment */
	  if(indi_msa->nseq != nused) { 
	    for(i = 0; i < nused; i++) { 
	      if((status = esl_key_Lookup(indi_msa->index, (esl_keyhash_Get(useme_keyhash, i)), NULL)) == eslENOTFOUND) { 
		esl_fatal("Error with list file %s, sequence %s does not exist in the alignment.", esl_opt_GetString(go, "--list"), esl_keyhash_Get(useme_keyhash, i));
	      }
	    }
	    /* we should never get here, but just in case */
	    esl_fatal("Error, couldn't find all the sequences from the list file %s in the alignment (%d expected, %d found).", esl_opt_GetString(go, "--list"), nused, indi_msa->nseq);
	  }

	  /* Now, the complement to the if statement above that begins: if( ! esl_opt_IsDefault(go, "--ifile")) *
	   * we need to get insert info from indi_msa in the event --ifile was not invoked */
	  if(esl_opt_IsDefault(go, "--ifile")) { 
	    get_insert_info_from_msa(indi_msa, ps->rflen, NULL, &indi_ins_ct); /* dies with esl_fatal() upon an error */
	  }

	  /* now indi_msa includes exactly the seqs listed in the list file, rewrite useme and nused */
	  ESL_ALLOC(useme, sizeof(int) * indi_msa->nseq); 
	  esl_vec_ISet(useme, indi_msa->nseq, TRUE);
	  nused = indi_msa->nseq;
	}
      }
      if((status = individual_seqs_sspostscript(go, errbuf, ps, (do_small ? indi_msa : msa), (do_small ? indi_ins_ct : ins_ct), useme, nused, hc_scheme, RBFIVERHSCHEME, hc_nbins[RBFIVERHSCHEME], hc_onecell, WHITEOC, LIGHTGREYOC)) != eslOK) esl_fatal(errbuf);

      if(esl_opt_GetBoolean(go, "--prob")) { 
	if((status = individual_posteriors_sspostscript(go, errbuf, ps, (do_small ? indi_msa : msa), useme, nused, hc_scheme, RBSIXRLSCHEME, hc_nbins[RBSIXRLSCHEME], hc_onecell, LIGHTGREYOC)) != eslOK) esl_fatal(errbuf);
      }
    }
  }
  else if(master_mode == SIMPLEMASKMODE) { 
    if(esl_opt_GetBoolean(go, "--mask-col")) { 
      if(ps->rflen != masklen) esl_fatal("MSA has consensus (non-gap RF) length of %d which != lane mask length of %d.", rflen, masklen);
      /* Paint positions excluded by the mask magenta, unless the mask has zero internal exclusions.
       * Such a mask is a 'truncating' mask, that excludes only a 5' contiguous set of columns, and a 3' contiguous set of columns, in this case paint excluded positions light grey. */
      if((status = colormask_sspostscript(go, errbuf, ps, msa, mask, hc_onecell, BLACKOC, (mask_has_internal_zeroes ? MAGENTAOC : LIGHTGREYOC))) != eslOK) esl_fatal(errbuf);
    }
    if(! esl_opt_IsDefault(go, "--mask-diff")) { 
      if((status = diffmask_sspostscript(go, errbuf, ps, msa, mask, mask2, hc_onecell, BLACKOC, CYANOC, MAGENTAOC, LIGHTGREYOC)) != eslOK) esl_fatal(errbuf);
    }
  }
  else if(master_mode == DRAWFILEMODE) { 
    if((status = drawfile2sspostscript(go, errbuf, ps)) != eslOK) esl_fatal(errbuf);
  }
  
  if((status = draw_sspostscript(ofp, go, errbuf, command, date, hc_scheme, ps, nused)) != eslOK) esl_fatal(errbuf);
  free(command);
  fclose(ofp);
  printf("# %d page postscript saved to file %s.\n", ps->npage, outfile);

  /* Cleanup, normal return
   */
  if(tabfp != NULL) { 
    fclose(tabfp); 
    printf("# Per position data saved to tab-delimited text file %s.\n", esl_opt_GetString(go, "--tabfile"));
  }

  if(abc_ct != NULL)  esl_Free2D((void **) abc_ct, msa->alen);
  if(pp_ct != NULL)   esl_Free2D((void **) pp_ct, msa->alen);
  if(bp_ct  != NULL)  esl_Free3D((void ***) bp_ct, msa->alen, abc->Kp);
  if(srfpos_ct != NULL) free(srfpos_ct);
  if(erfpos_ct != NULL) free(erfpos_ct);
  if(nseq_with_ins_ct != NULL) free(nseq_with_ins_ct);
  if(ins_ct != NULL) esl_Free2D((void **) ins_ct, msa_nseq);
  if(indi_ins_ct != NULL) esl_Free2D((void **) indi_ins_ct, indi_msa->nseq);
  if(mask != NULL) free(mask);
  if(date != NULL) free(date);
  if(useme != NULL) free(useme);
  free_sspostscript(ps);
  esl_alphabet_Destroy(abc);
  esl_msafile_Close(afp);
  esl_getopts_Destroy(go);
  free(hc_nbins);
  for(z = 0; z < NOC; z++) free(hc_onecell[z]);
  free(hc_onecell);
  for(z = 0; z < 11; z++) free(hc_scheme[0][z]);
  for(z = 0; z < 11; z++) free(hc_scheme[1][z]);
  for(z = 0; z < 6; z++) free(hc_scheme[2][z]);
  for(z = 0; z < 6; z++) free(hc_scheme[3][z]);
  for(z = 0; z < 5; z++) free(hc_scheme[4][z]);
  for(z = 0; z < 5; z++) free(hc_scheme[5][z]);
  free(hc_scheme[0]);
  free(hc_scheme[1]);
  free(hc_scheme[2]);
  free(hc_scheme[3]);
  free(hc_scheme[4]);
  free(hc_scheme[5]);
  free(hc_scheme);
  if(msa != NULL) esl_msa_Destroy(msa);
  if(indi_msa != NULL) esl_msa_Destroy(indi_msa);
  return 0;

  ERROR: 
  esl_fatal("Memory allocation error in main().");
  }

/* Function: create_sspostscript()
 * 
 * Purpose:  Create and initialize a SS postscript data structure.
 * Return:   ps
 */
static SSPostscript_t *
create_sspostscript()
{
  int status;
  SSPostscript_t *ps;

  ESL_ALLOC(ps, sizeof(SSPostscript_t));

  ps->npage    = 0;
  ps->modelname = NULL;
  ps->modeA = NULL;
  ps->descA = NULL;
  ps->headerx = 0.;
  ps->headery = 0.;
  ps->headerx_desc = 0.;
  ps->headerx_charsize = 0.;
  ps->headery_charsize = 0.;
  ps->desc_max_chars = 0;
  ps->legx = 0.;
  ps->legy = 0.;
  ps->cur_legy = 0.;
  ps->legx_charsize = 0.;
  ps->legy_charsize = 0.;
  ps->legx_max_chars = 0;
  ps->legx_stats = 0.;
  ps->pagex_max = 0.;
  ps->pagey_max = 0.;
  ps->scale = 0.;
  ps->regurgA  = NULL;
  ps->nregurg  = 0;
  ps->hundredsxA = ps->hundredsyA = NULL;
  ps->nhundreds = 0;
  ps->ticksx1A = ps->ticksx2A = ps->ticksy1A = ps->ticksy2A = NULL;
  ps->nticks = 0;
  ps->bpx1A = ps->bpx2A = ps->bpy1A = ps->bpy2A = NULL;
  ps->nbp = 0;
  ps->rxA = ps->ryA = NULL;
  ps->rflen = 0;
  ps->rrAA        = NULL;
  ps->rcolAAA     = NULL;
  ps->occlAAA     = NULL;
  ps->nocclA      = NULL;
  ps->sclAA       = NULL;
  ps->mask        = NULL;
  ps->nalloc      = 50;
  ps->msa_ct      = NULL;
  ps->msa_nbp     = 0;
  ps->msa_rf2a_map = NULL;
  ps->msa_a2rf_map = NULL;
  ps->uaseqlenA   = NULL;
  ps->seqidxA     = NULL;
  ps->msa         = NULL;
  return ps;

 ERROR: esl_fatal("create_sspostscript(): memory allocation error.");
  return NULL; /* NEVERREACHED */
}

/* Function: setup_sspostscript()
 * 
 * Purpose:  Determine positions for header and legend in a SSPostscript_t()
 * Return:   eslOK
 */
static int
setup_sspostscript(SSPostscript_t *ps, char *errbuf)
{
  if(ps->rflen == 0) ESL_FAIL(eslEINVAL, errbuf, "Failed to ready any residues in template file.");

  /* set up legx, legy, this is a hack (takes advantage of position of 3' residue in all SSU models) */
  ps->legx = ps->rxA[ps->rflen-1] + LEGX_OFFSET;
  ps->legy = ps->ryA[ps->rflen-1] + LEGY_OFFSET;
  ps->cur_legy = ps->legy;

  ps->pagex_max = POSTSCRIPT_PAGEWIDTH / ps->scale;
  ps->pagey_max = POSTSCRIPT_PAGEHEIGHT / ps->scale;

  ps->headerx = 0. + PAGE_SIDEBUF;
  ps->headery = ps->pagey_max - PAGE_TOPBUF - ((HEADER_FONTSIZE_UNSCALED) / ps->scale);

  /* determine max number of residues we can print before we run off the page in the legend section */
  float xroom, yroom;
  xroom  = ps->pagex_max - ps->legx;
  yroom  = (ps->legy - ps->pagey_max) * -1;
  ps->legx_charsize = (LEG_FONTSIZE_UNSCALED / COURIER_HEIGHT_WIDTH_RATIO) / ps->scale; 
  ps->legy_charsize = (LEG_FONTSIZE_UNSCALED) / ps->scale; 
  ps->legx_max_chars = (int) (xroom / ps->legx_charsize);
  ps->legy_max_chars = (int) (yroom / ps->legy_charsize);
  ps->legx_stats     = ps->pagex_max - PAGE_SIDEBUF - (LEG_EXTRA_COLUMNS * ps->legx_charsize);

  /* determine max size of description that will fit in header */
  float header_fontwidth, header_max_chars;
  header_fontwidth      = (HEADER_FONTSIZE_UNSCALED / COURIER_HEIGHT_WIDTH_RATIO) / ps->scale; 
  ps->headerx_charsize  = (HEADER_FONTSIZE_UNSCALED / COURIER_HEIGHT_WIDTH_RATIO) / ps->scale; 
  header_max_chars      = (int) (ps->pagex_max / ps->headerx_charsize) - 2;
  ps->headery_charsize  = (HEADER_FONTSIZE_UNSCALED) / ps->scale; 
  ps->desc_max_chars    = header_max_chars - (HEADER_MODELNAME_MAXCHARS + 6 + 6 + 8 +2); /*6,6,8 for #res,#bps,#seq plus 2 spaces each, plus 2 for after name */
  ps->headerx_desc      = ps->pagex_max - PAGE_SIDEBUF - (ps->desc_max_chars * ps->headerx_charsize);

  return eslOK;
}

/* Function: free_sspostscript()
 * 
 * Purpose:  Free a SS postscript data structure.
 * Return:   (void)
 */
static void
free_sspostscript(SSPostscript_t *ps)
{
  int i, p, c, l;

  if(ps->modelname != NULL) free(ps->modelname);

  if(ps->modeA != NULL)  free(ps->modeA);

  if(ps->descA != NULL) {
    for(i = 0; i < ps->npage; i++) { 
      if(ps->descA[i] != NULL) {
	free(ps->descA[i]);
      }
    }
    free(ps->descA);
  }

  if(ps->regurgA != NULL) {
    for(i = 0; i < ps->nregurg; i++) { 
      if(ps->regurgA[i] != NULL) {
	free(ps->regurgA[i]);
      }
    }
    free(ps->regurgA);
  }

  if(ps->hundredsxA != NULL) free(ps->hundredsxA);
  if(ps->hundredsyA != NULL) free(ps->hundredsyA);
  if(ps->ticksx1A != NULL) free(ps->ticksx1A);
  if(ps->ticksy1A != NULL) free(ps->ticksy1A);
  if(ps->ticksx2A != NULL) free(ps->ticksx2A);
  if(ps->ticksy2A != NULL) free(ps->ticksy2A);
  if(ps->bpx1A != NULL) free(ps->bpx1A);
  if(ps->bpy1A != NULL) free(ps->bpy1A);
  if(ps->bpx2A != NULL) free(ps->bpx2A);
  if(ps->bpy2A != NULL) free(ps->bpy2A);
  if(ps->rxA != NULL) free(ps->rxA);
  if(ps->ryA != NULL) free(ps->ryA);

  if(ps->rrAA != NULL) { 
    for(p = 0; p < ps->npage; p++) 
      if(ps->rrAA[p] != NULL) free(ps->rrAA[p]);
    free(ps->rrAA);
  }

  if(ps->rcolAAA != NULL) { 
    for(p = 0; p < ps->npage; p++) { 
      if(ps->rcolAAA[p] != NULL) { 
	for(c = 0; c < ps->rflen; c++) free(ps->rcolAAA[p][c]);
	free(ps->rcolAAA[p]); 
      }
    }
    free(ps->rcolAAA); 
  }

  if(ps->occlAAA != NULL) { 
    for(p = 0; p < ps->npage; p++) { 
      if(ps->occlAAA[p] != NULL) { 
	for(l = 0; l < ps->nocclA[p]; l++) { 
	  if(ps->occlAAA[p][l]->text != NULL) free(ps->occlAAA[p][l]->text);
	  free(ps->occlAAA[p][l]); /* rest is statically allocated memory */
	}
      free(ps->occlAAA[p]);
      }
    }
    free(ps->occlAAA);
  }

  if(ps->sclAA != NULL) { 
    for(p = 0; p < ps->npage; p++) { 
      if(ps->sclAA[p] != NULL) { 
	if(ps->sclAA[p]->limits != NULL) free(ps->sclAA[p]->limits);
	if(ps->sclAA[p]->counts != NULL) free(ps->sclAA[p]->counts);
	if(ps->sclAA[p]->counts_masked != NULL) free(ps->sclAA[p]->counts_masked);
	if(ps->sclAA[p]->text1 != NULL) free(ps->sclAA[p]->text1);
	if(ps->sclAA[p]->text2 != NULL) free(ps->sclAA[p]->text2);
	free(ps->sclAA[p]); /* statically allocated memory */
      }
    }
    free(ps->sclAA);
  }

  if(ps->nocclA != NULL) free(ps->nocclA);
  if(ps->msa_ct != NULL) free(ps->msa_ct);
  if(ps->msa_rf2a_map != NULL) free(ps->msa_rf2a_map);
  if(ps->msa_a2rf_map != NULL) free(ps->msa_a2rf_map);
  if(ps->mask != NULL) free(ps->mask);
  if(ps->seqidxA != NULL) free(ps->seqidxA);

  free(ps);
  return;
}

/* Function: create_onecell_colorlegend()
 * 
 * Purpose:  Create and initialize a one cell color legend data structure.
 * Return:   occl
 */
static OneCellColorLegend_t *
create_onecell_colorlegend(float *col, int nres, int nres_masked)
{
  int status;
  OneCellColorLegend_t *occl;

  ESL_ALLOC(occl, sizeof(OneCellColorLegend_t));

  /* initialize */
  esl_vec_FSet(occl->col, NCMYK, 0.);
  occl->text = NULL;

  /* set caller specified values */
  esl_vec_FCopy(col, NCMYK, occl->col);

  occl->nres = nres;
  occl->nres_masked = nres_masked;
  return occl;

 ERROR: esl_fatal("create_onecell_colorlegend(): memory allocation error.");
  return NULL; /* NEVERREACHED */
}


/* Function: create_scheme_colorlegend()
 * 
 * Purpose:  Create and initialize a scheme color legend data structure.
 * Return:   scl
 */
static SchemeColorLegend_t *
create_scheme_colorlegend(int scheme, int nbins, float *limits, int ints_only_flag)
{
  int status;
  SchemeColorLegend_t *scl;
  int i;
  int *counts;
  int *counts_masked;
  ESL_ALLOC(scl, sizeof(SchemeColorLegend_t));

  /* initialize */
  scl->text1 = NULL;
  scl->text2 = NULL;
  
  /* set caller specified values */
  scl->scheme = scheme;
  scl->nbins = nbins;
  ESL_ALLOC(scl->limits, sizeof(float) * (nbins+1));
  for(i = 0; i <= nbins; i++) { scl->limits[i] = limits[i]; }

  ESL_ALLOC(counts, sizeof(int) * nbins);
  ESL_ALLOC(counts_masked, sizeof(int) * nbins);
  esl_vec_ISet(counts, nbins, 0);
  esl_vec_ISet(counts_masked, nbins, 0);
  scl->counts = counts;
  scl->counts_masked = counts_masked;
  scl->ints_only_flag = ints_only_flag;
  return scl;

 ERROR: esl_fatal("create_scheme_colorlegend(): memory allocation error.");
  return NULL; /* NEVERREACHED */
}

/* Function: add_text_to_scheme_colorlegend()
 * 
 * Purpose:  Add text to an existing scheme color legend data structure.
 * Throws:   Exception if the text is too long.
 */
int
add_text_to_scheme_colorlegend(SchemeColorLegend_t *scl, char *text, int legx_max_chars, char *errbuf)
{
  int status;
  int i, idx;
  int max_chars_per_line;

  if(scl->text1 != NULL) esl_fatal("add_text_to_scheme_colorlegend(), text already exists!\n"); 
  if(scl->text2 != NULL) esl_fatal("add_text_to_scheme_colorlegend(), text already exists!\n"); 
  if(text == NULL) esl_fatal("add_text_to_scheme_colorlegend(), passed in text is NULL!\n"); 

  max_chars_per_line = legx_max_chars - LEG_EXTRA_COLUMNS -2;
  if(((int) strlen(text)) <= max_chars_per_line) {
    /* case 1, entire text can fit in one line */
    if((status = esl_strdup(text, -1, &(scl->text1))) != eslOK) esl_fatal("add_text_to_scheme_colorlegend(), error copying text");
    return eslOK;
  }
  else if(((int) strlen(text)) > ((2 * max_chars_per_line) - 6)) { 
    /* case 2, entire text can't even fit in two lines, 
     * (this is inexact doesn't account for size of words which is an issue b/c we break at newline,
     * this is why I have the extra '- 6');
     */
    ESL_FAIL(eslEINVAL, errbuf, "add_text_to_scheme_colorlegend(), text is %d chars, max allowed is %d (%s)\n", (int) strlen(text), ((2 * max_chars_per_line) - 6), text);
  }
  else { /* split it up into two lines */
    idx = max_chars_per_line - 1;
    while(text[idx] != ' ') { 
      idx--;
      if(idx < 0) ESL_FAIL(eslEINVAL, errbuf, "add_text_to_scheme_colorlegend(), couldn't find a breakpoint for splitting the string (%s)\n", text);
    }
    /* copy first string into text1 */
    ESL_ALLOC(scl->text1, sizeof(char) * (idx+1));
    for(i = 0; i < idx; i++) scl->text1[i] = text[i];
    scl->text1[idx] = '\0';

    /* copy remainder into text2 */
    int len = (int) strlen(text);
    idx++;
    ESL_ALLOC(scl->text2, sizeof(char) * (len - idx+1));
    for(i = idx; i < len; i++) scl->text2[i-idx] = text[i];
    scl->text2[len-idx] = '\0';
  }
  return eslOK;

 ERROR: ESL_FAIL(status, errbuf, "Error adding text to scheme_colorlegend probably out of memory");
}


/* Function: add_text_to_onecell_colorlegend()
 * 
 * Purpose:  Add text to an existing one cell color legend data structure.
 * Throws:   Exception if the text is too long.
 */
int
add_text_to_onecell_colorlegend(SSPostscript_t *ps, OneCellColorLegend_t *occl, char *text, int legx_max_chars, char *errbuf)
{
  int status;
  int max_chars_per_line;
  if(occl->text != NULL) esl_fatal("add_text_to_onecell_colorlegend(), text already exists!\n"); 
  if(text == NULL) esl_fatal("add_text_to_onecell_colorlegend(), passed in text is NULL!\n"); 

  max_chars_per_line = legx_max_chars - LEG_EXTRA_COLUMNS - 2 - ((int) ((LEG_BOXSIZE * 1.5) / ps->legx_charsize));
  /*printf("max: %d cur: %d text: %s\n", max_chars_per_line, (int) strlen(text), text);*/
  if(((int) strlen(text)) > (max_chars_per_line)) { 
    ESL_FAIL(eslEINVAL, errbuf, "add_text_to_onecell_colorlegend(), text is %d chars, max allowed is %d (%s)\n", (int) strlen(text), max_chars_per_line, text);
  }
  if((status = esl_strdup(text, -1, &(occl->text))) != eslOK) esl_fatal("add_text_to_onecell_colorlegend(), error copying text");
  return eslOK;
}

/* Function: add_page_desc_to_sspostscript()
 * 
 * Purpose:  Add text describing a particular page of a postscript object.
 * Throws:   Exception if the text is too long.
 */
int
add_page_desc_to_sspostscript(SSPostscript_t *ps, int page, char *text, char *errbuf)
{
  int status;
  int i, j;
  int max_both_lines;

  if(ps->descA[page] != NULL) ESL_FAIL(eslEINVAL, errbuf, "add_page_desc_to_sspostscript(), description for page %d already exists!\n", page); 
  if(text == NULL)            ESL_FAIL(eslEINVAL, errbuf, "add_page_desc_to_sspostscript(), passed in text is NULL!\n"); 

  max_both_lines =(2. * ps->desc_max_chars);
  if(ps->modeA[page] == INDIMODE || ps->modeA[page] == SIMPLEMASKMODE) 
    max_both_lines--; /* b/c we have to add a '-' to split up the larger strings onto 2 lines */

  /* check to see if we can fit the text onto two lines of max width ps->desc_max_chars */
  int textlen = (int) strlen(text);
  if(textlen <= ps->desc_max_chars) { /* fine, this will fit on one line */
    if((status = esl_strdup(text, -1, &(ps->descA[page]))) != eslOK) ESL_FAIL(eslEINVAL, errbuf, "add_page_desc_to_sspostscript(), error copying text");
  }
  else if (textlen <= max_both_lines) { 
    if(ps->modeA[page] == ALIMODE) { 
      /* maybe fine, make sure there's a break point (space ' ') that will break this string into two strings of length <= ps->desc_max_chars */
      i = ps->desc_max_chars;
      while(text[i] != ' ') { 
	i--; 
	if(i < 0) ESL_FAIL(eslEINVAL, errbuf, "add_page_desc_to_sspostscript(), first word of text (%s) is more than max allowed of %d chars", text, ps->desc_max_chars);
      }
      /* found last space before max width, make sure it breaks the line up into two chunks of valid size */
      if((textlen - (i+1)) <= ps->desc_max_chars) { /* we're good */
	if((status = esl_strdup(text, -1, &(ps->descA[page]))) != eslOK) ESL_FAIL(eslEINVAL, errbuf, "add_page_desc_to_sspostscript(), error copying text");
	ps->descA[page][i] = '\n'; /* so we can remember where the break is */
      }
      else ESL_FAIL(eslEINVAL, errbuf, "add_page_desc_to_sspostscript(), couldn't find break point (' ') for partitioning text into two valid size chunks (%s)", text);
    }
    else { /* INDIMODE or SIMPLEMASKMODE, sequence/mask name bigger than 1 line, but not 2, we put a '-' in it at the end of line 1 and add a '\n' so we remember where it was */
      ESL_ALLOC(ps->descA[page], sizeof(char) * (textlen + 3)); /* +3 so we have space for the extra '-' and '\n' */
      for(i = 0; i < ps->desc_max_chars; i++) ps->descA[page][i] = text[i];
      i = ps->desc_max_chars;
      ps->descA[page][i] = '-';
      i++; 
      ps->descA[page][i] = '\n';
      for(i = ps->desc_max_chars; i < textlen; i++) ps->descA[page][i+2] = text[i];
      ps->descA[page][textlen+2] = '\0';
    }
  }
  else /* the text won't fit on two lines */
    if(ps->modeA[page] != INDIMODE) { /* not fine, this won't fit on two lines */
      ESL_FAIL(eslEINVAL, errbuf, "add_page_desc_to_sspostscript(), text is %d chars, max allowed is %d (%s)\n", textlen, max_both_lines, text);
    }
    else { /* INDIMODE or SIMPLEMASKMODE, sequence/mask name exceeds max, we put a '-' in it at the end of line 1 and truncate it */
      ESL_ALLOC(ps->descA[page], sizeof(char) * (max_both_lines + 2)); /* +2 so we have space for the extra '\n' (which we won't print), and the '\0' */
      /* first look to see if there's a space before desc_max_chars, this will never happen for a seq name, but may for a mask description */
      j = ps->desc_max_chars;
      while(text[j] != ' ' && j > 0) { 
	j--; 
      }
      if(j == 0) { j = ps->desc_max_chars; } /* no space */
      for(i = 0; i < j; i++) ps->descA[page][i] = text[i];
      i = j;
      if(j == ps->desc_max_chars && text[j] != ' ') { 
	ps->descA[page][i] = '-';
      }
      i++; 
      ps->descA[page][i] = '\n';
      for(i = j; i < j + ps->desc_max_chars; i++) ps->descA[page][i+2] = text[i];
      ps->descA[page][i+2] = '\0';
  }
  return eslOK;

 ERROR: ESL_FAIL(status, errbuf, "add_page_desc_to_sspostscript() error, probably out of memory.");
}


/* Function: add_diffmask_page_desc_to_sspostscript()
 * 
 * Purpose:  Add text describing a diff mask page of a postscript object.
 * Throws:   Exception if the text is too long.
 */
int
add_diffmask_page_desc_to_sspostscript(SSPostscript_t *ps, int page, char *mask1, char *mask2, char *errbuf)
{
  int status;
  int i;
  char *mask1desc = NULL;
  char *mask2desc = NULL;
  int len2copy;
  int mask1len;
  int mask2len;
  void *tmp;

  if(ps->descA[page] != NULL) ESL_FAIL(eslEINVAL, errbuf, "add_diffmask_page_desc_to_sspostscript(), description for page %d already exists!\n", page); 
  if(mask1 == NULL)            ESL_FAIL(eslEINVAL, errbuf, "add_diffmask_page_desc_to_sspostscript(), passed in mask1 is NULL!\n"); 
  if(mask2 == NULL)            ESL_FAIL(eslEINVAL, errbuf, "add_diffmask_page_desc_to_sspostscript(), passed in mask2 is NULL!\n"); 

  /* check to see if we can fit the text onto two lines of max width ps->desc_max_chars */
  mask1len = (int) strlen(mask1);
  mask2len = (int) strlen(mask2);
  if((status = esl_strcat(&(mask1desc), -1, "mask 1: ", -1)) != eslOK) ESL_FAIL(eslEINVAL, errbuf, "add_diffmask_page_desc_to_sspostscript(), error copying text");
  if((mask1len + 8) <= ps->desc_max_chars) { /* this will fit on one line */
    if((status = esl_strcat(&(mask1desc), -1, mask1, -1)) != eslOK) ESL_FAIL(eslEINVAL, errbuf, "add_diffmask_page_desc_to_sspostscript(), error copying text");
  }
  else { /* won't fit on one line, include as much as we can */
    len2copy = ps->desc_max_chars - 8 - 3; /* 8 is for "mask 1: " at beginning, 3 is for "..." at end */
    ESL_RALLOC(mask1desc, tmp, sizeof(char) * ps->desc_max_chars+1);
    for(i = 0; i < len2copy; i++) { 
      mask1desc[8+i] = mask1[i];
    }
    mask1desc[8+len2copy]   = '.';
    mask1desc[8+len2copy+1] = '.';
    mask1desc[8+len2copy+2] = '.';
    mask1desc[8+len2copy+3] = '\0'; 
  }

  /* repeat for mask 2 */
  if((status = esl_strcat(&(mask2desc), -1, "mask 2: ", -1)) != eslOK) ESL_FAIL(eslEINVAL, errbuf, "add_diffmask_page_desc_to_sspostscript(), error copying text");
  if((mask2len + 8) <= ps->desc_max_chars) { /* this will fit on one line */
    if((status = esl_strcat(&(mask2desc), -1, mask2, -1)) != eslOK) ESL_FAIL(eslEINVAL, errbuf, "add_diffmask_page_desc_to_sspostscript(), error copying text");
  }
  else { /* won't fit on one line, include as much as we can */
    len2copy = ps->desc_max_chars - 8 - 3; /* 8 is for "mask 1: " at beginning, 3 is for "..." at end */
    ESL_RALLOC(mask2desc, tmp, sizeof(char) * ps->desc_max_chars+1);
    for(i = 0; i < len2copy; i++) { 
      mask2desc[8+i] = mask2[i];
    }
    mask2desc[8+len2copy]   = '.';
    mask2desc[8+len2copy+1] = '.';
    mask2desc[8+len2copy+2] = '.';
    mask2desc[8+len2copy+3] = '\0'; 
  }

  /* concatenate them in ps->descA[page] */
  if((status = esl_strcat(&(ps->descA[page]), -1, mask1desc, -1)) != eslOK) ESL_FAIL(eslEINVAL, errbuf, "add_diffmask_page_desc_to_sspostscript(), error copying text");
  if((status = esl_strcat(&(ps->descA[page]), -1, "\n", -1)) != eslOK) ESL_FAIL(eslEINVAL, errbuf, "add_diffmask_page_desc_to_sspostscript(), error copying text");
  if((status = esl_strcat(&(ps->descA[page]), -1, mask2desc, -1)) != eslOK) ESL_FAIL(eslEINVAL, errbuf, "add_diffmask_page_desc_to_sspostscript(), error copying text");

  free(mask1desc);
  free(mask2desc);
  return eslOK;

 ERROR: ESL_FAIL(status, errbuf, "add_page_desc_to_sspostscript() error, probably out of memory.");
}

/* Function: add_mask_to_sspostscript
 * 
 * Purpose:  Add a mask to a sspostscript object.
 */
int
add_mask_to_ss_postscript(SSPostscript_t *ps, char *mask)
{
  int status;
  if(ps->mask != NULL) { esl_fatal("add_mask_to_ss_postscript(), mask is non-null!\n"); }
  if(mask == NULL) esl_fatal("add_mask_to_ss_postscript(), passed in mask is NULL!\n"); 
  if((status = esl_strdup(mask, -1, &(ps->mask))) != eslOK) esl_fatal("add_mask_to_ss_postscript(), error copying mask");
  return eslOK;
}


/* Function: draw_legend_column_headers()
 * 
 * Purpose:  Draw the legend column headers.
 * Return:   eslOK
 */
static int 
draw_legend_column_headers(FILE *fp, SSPostscript_t *ps, char *errbuf)
{
  int status;
  int i;
  float x, y, legend_fontsize;
  char *cur_string;
  int cur_width = 0;

  legend_fontsize = LEG_FONTSIZE_UNSCALED / ps->scale;

  x = ps->legx;
  y = ps->cur_legy;
  if(ps->mask != NULL) { 
    y -= 0.625 * LEG_BOXSIZE;
  }
  /*fprintf(fp, "/%s findfont %f scalefont setfont\n", SPECIAL_FONT, legend_fontsize);*/
  fprintf(fp, "(%s) %.4f %.4f moveto show\n", "LEGEND", x, (y + (LEG_BOXSIZE * .25)));
  /*fprintf(fp, "/%s findfont %f scalefont setfont\n", LEG_FONT, legend_fontsize);*/

  x = ps->legx_stats;
  y = ps->cur_legy;
  cur_width = ps->legx_max_chars - LEG_EXTRA_COLUMNS -2;

  ESL_ALLOC(cur_string, sizeof(char) * (cur_width+1));
  for(i = 0; i < cur_width; i++) cur_string[i] = '-'; 
  cur_string[cur_width] = '\0';

  if(ps->mask != NULL) { 
    fprintf(fp, "(%4s  %4s) %.4f %.4f moveto show\n", "", " in ", x, (y + (LEG_BOXSIZE * .25)));
    y -= 0.625 * LEG_BOXSIZE;
    fprintf(fp, "(%4s  %4s) %.4f %.4f moveto show\n", "all", "mask", x, (y + (LEG_BOXSIZE * .25)));
    y -= 0.625 * LEG_BOXSIZE;
    fprintf(fp, "(%s) %.4f %.4f moveto show\n", cur_string, ps->legx, (y + (LEG_BOXSIZE * .25)));
    fprintf(fp, "(----  ----) %.4f %.4f moveto show\n", x, (y + (LEG_BOXSIZE * .25)));
  }
  else { 
    fprintf(fp, "(%5s) %.4f %.4f moveto show\n", "count", x, (y + (LEG_BOXSIZE * .25)));
    y -= 0.625 * LEG_BOXSIZE;
    fprintf(fp, "(%s) %.4f %.4f moveto show\n", cur_string, ps->legx, (y + (LEG_BOXSIZE * .25)));
    fprintf(fp, "(-----) %.4f %.4f moveto show\n", x, (y + (LEG_BOXSIZE * .25)));
  }
  ps->cur_legy = y - (1.0 * LEG_BOXSIZE);
  
  free(cur_string);

  return eslOK;

 ERROR: ESL_FAIL(status, errbuf, "ERROR drawing legend column headers, probably out of memory.");
}

/* Function: draw_onecell_colorlegend()
 * 
 * Purpose:  Print a one cell color legend to an open file.
 * Return:   eslOK
 */
static int 
draw_onecell_colorlegend(FILE *fp, OneCellColorLegend_t *occl, SSPostscript_t *ps, int occl_idx)
{
  float x, y;
  int cp;
  float fontsize;

  /* object is valid, print it */
  x = ps->legx;
  y = ps->cur_legy;

  fontsize = LEG_FONTSIZE_UNSCALED / ps->scale;

  /* print cell */
  fprintf(fp, "newpath\n");
  fprintf(fp, "  %.2f %.2f moveto", x, y);
  fprintf(fp, "  0 %.3f rlineto %.3f 0 rlineto 0 %.3f rlineto closepath\n", LEG_BOXSIZE, LEG_BOXSIZE, (-1 * LEG_BOXSIZE));
  fprintf(fp, "  ");
  for(cp = 0; cp < NCMYK; cp++) fprintf(fp, "%.4f ", occl->col[cp]);
  fprintf(fp, "setcmykcolor\n");
  fprintf(fp, "  fill\n");
  
  x += LEG_BOXSIZE * 1.5;

  /* print text for this legend */
  if(occl->text != NULL) { 
    /* back to black */
    fprintf(fp, "  0.00 0.00 0.00 1.00 setcmykcolor\n");
    fprintf(fp, "/%s findfont %f scalefont setfont\n", LEG_FONT, fontsize);
    fprintf(fp, "(%s) %.4f %.4f moveto show\n", occl->text, x, y + (LEG_BOXSIZE * 0.25));

    /* print stats */
    x = ps->legx_stats;
    if(ps->mask != NULL) { 
      fprintf(fp, "(%4d  %4d) %.4f %.4f moveto show\n", occl->nres, occl->nres_masked, x, y + (LEG_BOXSIZE * 0.25));
    }
    else { 
      fprintf(fp, "(%5d) %.4f %.4f moveto show\n", occl->nres, x, y + (LEG_BOXSIZE * 0.25));
    }
  }

  /* reset color to black */ 
  fprintf(fp, "  %.4f %.4f %.4f %.4f setcmykcolor\n", 0., 0., 0., 1.);
  y -= LEG_BOXSIZE * 1.5;
  ps->cur_legy = y;
  
  return eslOK;
}


/* Function: draw_scheme_colorlegend()
 * 
 * Purpose:  Print a scheme color legend to an open file.
 * Return:   eslOK
 */
static int 
draw_scheme_colorlegend(const ESL_GETOPTS *go, FILE *fp, SchemeColorLegend_t *scl, float **hc_scheme, SSPostscript_t *ps, int page)
{
  float x, y;
  int cp;
  int c,i,n1s;
  float fontsize;
  int do_circle_mask, do_square_mask, do_x_mask, do_border;
  int do_mask;
  float old_x;

  do_mask = (ps->mask == NULL) ? FALSE : TRUE;
  do_border = (!esl_opt_GetBoolean(go, "--mask-a"));
  do_circle_mask = do_square_mask = do_x_mask = FALSE;
  if(esl_opt_GetBoolean(go, "--mask-u")) { do_square_mask = TRUE; }
  else if(esl_opt_GetBoolean(go, "--mask-x")) { do_x_mask = TRUE; }
  else do_circle_mask = TRUE;

  x = ps->legx;
  y = ps->cur_legy;
  //y = ps->legy - (ps->nocclA[page] * (LEG_BOXSIZE * 1.5));
  fontsize = LEG_FONTSIZE_UNSCALED / ps->scale;
  fprintf(fp, "/%s findfont %f scalefont setfont\n", LEG_FONT, fontsize);
  fprintf(fp, "  0.00 0.00 0.00 1.00 setcmykcolor\n");

  float colvec[NCMYK];
  colvec[0] = colvec[1] = colvec[2] = 0.;
  colvec[3] = 1.0;
  if(do_mask) { /* print cells showing difference between masked and unmasked */
    /*x -= LEG_BOXSIZE;*/
    /*y -= LEG_BOXSIZE;*/
    fprintf(fp, "%.1f setlinewidth\n", LEG_BOXSIZE/4.);
    fprintf(fp, "newpath\n");
    fprintf(fp, "  %.2f %.2f moveto", x, y);
    fprintf(fp, "  0 %.3f rlineto %.3f 0 rlineto 0 %.3f rlineto closepath\n", LEG_BOXSIZE, LEG_BOXSIZE, (-1 * LEG_BOXSIZE));
    fprintf(fp, "  ");
    for(cp = 0; cp < NCMYK; cp++) { 
      fprintf(fp, "%.4f ", colvec[cp]);
    }
    fprintf(fp, "setcmykcolor\n");
    fprintf(fp, "  fill\n");

    /* print label */
    x += LEG_BOXSIZE * 1.5;
    y += LEG_BOXSIZE * 0.625;
    fprintf(fp, "(included by mask) %.4f %.4f moveto show\n", x, y);
    y -= LEG_BOXSIZE * 0.625;
    fprintf(fp, "((all colors)) %.4f %.4f moveto show\n", x, y);
    x -= LEG_BOXSIZE * 1.5;

    /* print stats for included by mask */
    old_x = x;
    n1s = 0;
    for(i = 0; i < ps->rflen; i++) if(ps->mask[i] == '1') n1s++; 
    x = ps->legx_stats;
    y += LEG_BOXSIZE * 0.3125;
    fprintf(fp, "(%4s  %4d) %.4f %.4f moveto show\n", "-", n1s, x, y);
    y -= LEG_BOXSIZE * 0.3125;

    x = old_x;
    y -= LEG_BOXSIZE * 1.5;
    draw_masked_block(fp, x, y, colvec, do_circle_mask, do_square_mask, do_x_mask, do_border, LEG_BOXSIZE);

    x += LEG_BOXSIZE * 1.5;
    y += LEG_BOXSIZE * 0.625;
    fprintf(fp, "(excluded by mask) %.4f %.4f moveto show\n", x, y);
    y -= LEG_BOXSIZE * 0.625;
    fprintf(fp, "((all colors)) %.4f %.4f moveto show\n", x, y);

    /* print stats for excluded by mask */
    old_x = x;
    x = ps->legx_stats;
    y += LEG_BOXSIZE * 0.3125;
    fprintf(fp, "(%4s  %4d) %.4f %.4f moveto show\n", "-", ps->rflen-n1s, x, y);

    y -= LEG_BOXSIZE * 1.8125;
    x = ps->legx;
  }

  /* print text for this legend */
  if(scl->text1 != NULL) { 
    if(scl->text2 == NULL) { 
      fprintf(fp, "(%s:) %.4f %.4f moveto show\n", scl->text1, x, (y + (LEG_BOXSIZE * .25)));
    }
    else { 
      fprintf(fp, "(%s) %.4f %.4f moveto show\n", scl->text1, x, (y + (LEG_BOXSIZE * .25)));
      y -= LEG_BOXSIZE * 0.625;
      fprintf(fp, "(%s:) %.4f %.4f moveto show\n", scl->text2, x, (y + (LEG_BOXSIZE * .25)));
    }
  }
  y -= LEG_BOXSIZE;
  
  /* print masked scheme color cells */
  /*if(do_mask) { 
    fprintf(fp, "%.1f setlinewidth\n", LEG_BOXSIZE/4.);
    for(c = 0; c < scl->nbins; c++) { 
    draw_masked_block(fp, x, y, hc_scheme[c], do_circle_mask, do_square_mask, do_x_mask, do_border, LEG_BOXSIZE);
    y -= LEG_BOXSIZE;
    }
    y += (LEG_BOXSIZE * scl->nbins);
    x += 1.5 * LEG_BOXSIZE;
    fprintf(fp, "1.0 setlinewidth\n");
  }
  */

  /* print scheme color cells and labels next to them */
  for(c = 0; c < scl->nbins; c++) { 
    fprintf(fp, "newpath\n");
    fprintf(fp, "  %.2f %.2f moveto", x, y);
    fprintf(fp, "  0 %.3f rlineto %.3f 0 rlineto 0 %.3f rlineto closepath\n", LEG_BOXSIZE, LEG_BOXSIZE, (-1 * LEG_BOXSIZE));
    fprintf(fp, "  ");
    for(cp = 0; cp < NCMYK; cp++) { 
      fprintf(fp, "%.4f ", hc_scheme[c][cp]);
    }
    fprintf(fp, "setcmykcolor\n");
    fprintf(fp, "  fill\n");

    /* print label */
    x += LEG_BOXSIZE * 1.5;
    y += LEG_BOXSIZE * 0.25;
    fprintf(fp, "  0.00 0.00 0.00 1.00 setcmykcolor\n");
    if(esl_FCompare(scl->limits[c+1], SSDRAWINFINITY, eslSMALLX1) == eslOK) { /* max value is infinity, special case */
      if(c != scl->nbins-1) esl_fatal("ERROR when drawing color legend, limits[%d] is INFINITY, but this is reserved only for the max limit", c+1);
      if(scl->ints_only_flag) fprintf(fp, "(>=%d) %.4f %.4f moveto show\n",   (int) scl->limits[c], x, y);
      else                    fprintf(fp, "(>=%3.f) %.4f %.4f moveto show\n", scl->limits[c], x, y);
    }
    else if(scl->ints_only_flag) { 
      if(c == scl->nbins-1) { 
	fprintf(fp, "(\\[%d-%d\\]) %.4f %.4f moveto show\n", (int) scl->limits[c], (int) scl->limits[c+1], x, y);
      }
      else if(esl_FCompare(scl->limits[c], scl->limits[c+1]-1, eslSMALLX1) == eslOK) { /* next limit is exactly 1 plus cur limit, don't do range, define single int */
	fprintf(fp, "(%d) %.4f %.4f moveto show\n", (int) scl->limits[c], x, y);
      }
      else { 
	fprintf(fp, "(\\[%d-%d\\]) %.4f %.4f moveto show\n", (int) scl->limits[c], (int) scl->limits[c+1]-1, x, y);
      }
    }
    else { 
      if(c == scl->nbins-1) fprintf(fp, "(\\[%.3f-%.3f\\]) %.4f %.4f moveto show\n", scl->limits[c], scl->limits[c+1], x, y);
      else                  fprintf(fp, "(\\[%.3f-%.3f\\)) %.4f %.4f moveto show\n", scl->limits[c], scl->limits[c+1], x, y);
    }
    /* print stats */
    old_x = x;
    x = ps->legx_stats;
    if(ps->mask != NULL) { 
      fprintf(fp, "(%4d  %4d) %.4f %.4f moveto show\n", scl->counts[c], scl->counts_masked[c], x, y);
    }
    else { 
      fprintf(fp, "(%5d) %.4f %.4f moveto show\n", scl->counts[c], x, y);
    }

    x = old_x - LEG_BOXSIZE * 1.5;
    y -= LEG_BOXSIZE * 0.25;
    y -= LEG_BOXSIZE;
  }

  /* reset color to black */ 
  fprintf(fp, "  %.4f %.4f %.4f %.4f setcmykcolor\n", 0., 0., 0., 1.);
  
  ps->cur_legy = y;
  return eslOK;
}

/* Function: draw_sspostscript()
 * 
 * Purpose:  Print a SS postscript data structure.
 * Return:   eslOK on success;
 *           eslEINCOMPAT if ps->npage == 0
 */
static int
draw_sspostscript(FILE *fp, const ESL_GETOPTS *go, char *errbuf, char *command, char *date, float ***hc_scheme, SSPostscript_t *ps, int nused)
{
  int status;
  int p, pi, i, c, l, si;
  int do_circle_mask, do_square_mask, do_x_mask, do_border;
  int *page_orderA;
  int rfoffset;

  do_border = (!esl_opt_GetBoolean(go, "--mask-a"));
  do_circle_mask = do_square_mask = do_x_mask = FALSE;
  if(esl_opt_GetBoolean(go, "--mask-u")) { do_square_mask = TRUE; }
  else if(esl_opt_GetBoolean(go, "--mask-x")) { do_x_mask = TRUE; }
  else do_circle_mask = TRUE;

  if(ps->npage == 0) ESL_FAIL(eslEINCOMPAT, errbuf, "draw_sspostscript, ps->npage == 0\n");

  /* determine print order or pages, this is just 0..npage-1 UNLESS:
   *  --indi and --prob and either (--all or --list) enabled, in which case, we put the 
   *    posteriors immediately after the sequences
   */
  ESL_ALLOC(page_orderA, sizeof(int) * ps->npage);
  if(esl_opt_GetBoolean(go, "--indi") && 
     esl_opt_GetBoolean(go, "--prob") &&
     (esl_opt_GetBoolean(go, "--all") || (! esl_opt_IsDefault(go, "--list")))) {
    pi = 0;
    rfoffset = 0;
    if(! esl_opt_GetBoolean(go, "-q")) { 
      page_orderA[0] = 0; /* print consensus sequence first */ 
      pi++;
      rfoffset = 1;
    } /* otherwise, we didn't print the consensus sequence */
    for(si = 0; si < nused; si++) { 
      page_orderA[pi] = si + rfoffset; /* the indi sequence page */
      pi++;
      page_orderA[pi] = (si + nused) + rfoffset; /* the posterior page */
      pi++;
    }
  }
  else { 
    for(pi = 0; pi < ps->npage; pi++) page_orderA[pi] = pi;
  }

  /* draw the pages */
  for(pi = 0; pi < ps->npage; pi++) { 
    p = page_orderA[pi];
    ps->cur_legy = ps->legy;

    /* scale section */
    fprintf(fp, "%% begin scale\n");
    fprintf(fp, "%.2f %.2f scale\n", ps->scale, ps->scale);
    fprintf(fp, "%% end scale\n\n");
      
    /* header section */
    if((status = draw_header_and_footer(fp, go, errbuf, ps, p, pi+1)) != eslOK) return status;

    /* regurgitated section */
    if(ps->regurgA != NULL) {
      fprintf(fp, "%% begin regurgitate\n");
      for(i = 0; i < ps->nregurg; i++)
	fprintf(fp, "%s", ps->regurgA[i]);
      fprintf(fp, "%% end regurgitate\n\n");
    }

    /* 'text hundreds' section */
    for(i = 0; i < ps->nhundreds; i++) { 
      if(i == 0) { 
	fprintf(fp, "%% begin text hundreds\n");
	fprintf(fp, "/%s findfont %.2f scalefont setfont\n", HUNDREDS_FONT, HUNDREDS_FONTSIZE);
	fprintf(fp, "0.00 0.00 0.00 1.00 setcmykcolor\n"); /* black */
      }
      fprintf(fp, "(%d) %.2f %.2f moveto show\n", (i+1) * 100, ps->hundredsxA[i], ps->hundredsyA[i]); 
      if(i == (ps->nhundreds-1)) { 
	fprintf(fp, "%% end text hundreds\n\n");
      }
    }
    
    /* 'lines ticks' section */
    for(i = 0; i < ps->nticks; i++) { 
      if(i == 0) { 
	fprintf(fp, "%% begin lines ticks\n");
	fprintf(fp, "%.2f setlinewidth\n", TICKS_LINEWIDTH);
	fprintf(fp, "0.00 0.00 0.00 1.00 setcmykcolor\n"); /* black */
      }
      fprintf(fp, "%.2f %.2f %.2f %.2f newpath moveto lineto stroke\n", ps->ticksx1A[i], ps->ticksy1A[i], ps->ticksx2A[i], ps->ticksy2A[i]);
      if(i == (ps->nticks-1)) { 
	fprintf(fp, "%% end lines ticks\n\n");
      }
    }

    /* 'lines bpconnects' section */
    for(i = 0; i < ps->nbp; i++) { 
      if(i == 0) { 
	fprintf(fp, "%% begin lines bpconnects\n");
	fprintf(fp, "%.2f setlinewidth\n", BP_LINEWIDTH);
	fprintf(fp, "0.00 0.00 0.00 1.00 setcmykcolor\n"); /* black */
      }
      fprintf(fp, "%.2f %.2f %.2f %.2f newpath moveto lineto stroke\n", ps->bpx1A[i], ps->bpy1A[i], ps->bpx2A[i], ps->bpy2A[i]);
      if(i == (ps->nbp-1)) { 
	fprintf(fp, "%% end lines bpconnects\n\n");
      }
    }

    /* 'text residues' section */
    /* NOTE: I only print this out so that this file could possibly be used as a template */
    fprintf(fp, "%% begin text residues\n");
    fprintf(fp, "/%s findfont %.2f scalefont setfont\n", RESIDUES_FONT, RESIDUES_FONTSIZE);
    fprintf(fp, "0.00 0.00 0.00 1.00 setcmykcolor\n"); /* black */
    for(i = 0; i < ps->rflen; i++) { 
      fprintf(fp, "() %.2f %.2f moveto show\n", ps->rxA[i], ps->ryA[i]);
    }
    fprintf(fp, "%% end text residues\n");

    /* the rest of the text will be ignored by ssu-draw if the output
     * file we're creating is read in as a template file later on
     */
    fprintf(fp, "%% begin ignore\n");
    fprintf(fp, "0.00 0.00 0.00 1.00 setcmykcolor\n"); /* set to black */
    fprintf(fp, "/%s findfont %f scalefont setfont\n", LEG_FONT, LEG_FONTSIZE_UNSCALED / ps->scale);

    /* draw legend headers, if we have a legend */
    if((ps->nocclA[p] > 0) || (ps->sclAA != NULL && ps->sclAA[p] != NULL)) { 
      if(! (esl_opt_GetBoolean(go, "--no-leg"))) { 
	if((status = draw_legend_column_headers(fp, ps, errbuf)) != eslOK) return status;
      }
    }

    /* print one cell color legends, if any */
    if(ps->occlAAA != NULL && ps->occlAAA[p] != NULL) { 
      for(l = 0; l < ps->nocclA[p]; l++) 
      if(! (esl_opt_GetBoolean(go, "--no-leg"))) { 
	draw_onecell_colorlegend(fp, ps->occlAAA[p][l], ps, l);
      }
    }
    /* print scheme color legends, if any */
    if(ps->sclAA != NULL && ps->sclAA[p] != NULL) { 
      if(! (esl_opt_GetBoolean(go, "--no-leg"))) { 
	draw_scheme_colorlegend(go, fp, ps->sclAA[p], hc_scheme[ps->sclAA[p]->scheme], ps, p);
      }
    }

    if(ps->rcolAAA != NULL && ps->rcolAAA[p] != NULL) { 
      if(ps->mask != NULL) { 
	fprintf(fp, "2.0 setlinewidth\n");
	if(do_border && do_x_mask)      { fprintf(fp, "1.0 setlinewidth\n"); }
	if(do_border && do_square_mask) { fprintf(fp, "2.0 setlinewidth\n"); }
	if(do_border && do_circle_mask) { fprintf(fp, "2.5 setlinewidth\n"); }
	for(c = 0; c < ps->rflen; c++) { 
	  fprintf(fp, "%sresidue %d\n", "%", c+1);
	  if(ps->mask[c] == '0') { 
	    draw_masked_block(fp, ps->rxA[c]-1., ps->ryA[c]-1., ps->rcolAAA[p][c], do_circle_mask, do_square_mask, do_x_mask, do_border, SS_BOXSIZE);
	  }
	  else { /* cell is within mask, ps->mask[c] == '1' */
	    fprintf(fp, "newpath\n");
	    fprintf(fp, "  %.2f %.2f moveto", ps->rxA[c] - 1., ps->ryA[c] -1.);
	    fprintf(fp, "  0 8 rlineto 8 0 rlineto 0 -8 rlineto closepath\n");
	    fprintf(fp, "  %.4f %.4f %.4f %.4f setcmykcolor\n", ps->rcolAAA[p][c][0], ps->rcolAAA[p][c][1], ps->rcolAAA[p][c][2], ps->rcolAAA[p][c][3]);
	    fprintf(fp, "  fill\n"); 
	  }
	}
	fprintf(fp, "1.00 setlinewidth\n");
      }
      else { /* no mask, all cells are printed the same */
	for(c = 0; c < ps->rflen; c++) { 
	  fprintf(fp, "%sresidue %d\n", "%", c+1);
	  fprintf(fp, "newpath\n");
	  fprintf(fp, "  %.2f %.2f moveto", ps->rxA[c] - 1., ps->ryA[c] -1.);
	  fprintf(fp, "  0 8 rlineto 8 0 rlineto 0 -8 rlineto closepath\n");
	  fprintf(fp, "  %.4f %.4f %.4f %.4f setcmykcolor\n", ps->rcolAAA[p][c][0], ps->rcolAAA[p][c][1], ps->rcolAAA[p][c][2], ps->rcolAAA[p][c][3]);
	  fprintf(fp, "  fill\n");
	}
      }
      /* back to black */
      fprintf(fp, "  0.00 0.00 0.00 1.00 setcmykcolor\n");
    }

    if(ps->rrAA[p] != NULL) { 
      fprintf(fp, "/%s findfont %f scalefont setfont\n", RESIDUES_FONT, RESIDUES_FONTSIZE);
      for(c = 0; c < ps->rflen; c++) { 
	fprintf(fp, "(%c) %.2f %.2f moveto show\n", ps->rrAA[p][c], ps->rxA[c], ps->ryA[c]);
      }
    }
    fprintf(fp, "grestore\nshowpage\n");
    fprintf(fp, "%% end ignore\n\n");
  }
  free(page_orderA);
  return eslOK;

 ERROR: ESL_FAIL(eslEINVAL, errbuf, "draw_sspostscript() error, probably out of memory.");
}


/* parse_template_file
 *
 * Read secondary structure templates from a postscript 
 * template file derived from the Gutell CRW website until
 * the one that corresponds to our alignment is found, 
 * (defined as the first template structure that has the
 * same consensus length as the passed in <msa_rflen>)
 * or we run out structure templates. This function
 * repeatedly calls parse_template() which actually parses
 * the templates.
 * 
 * If we run out of structure templates of anything is 
 * invalid in any of the templates in the file, we 
 * return a non-eslOK status code to caller and fill
 * errbuf with the error message.
 * 
 * Returns:  eslOK on success.
 */
int
parse_template_file(char *filename, const ESL_GETOPTS *go, char *errbuf, int msa_rflen, SSPostscript_t **ret_ps)
{
  int             status;
  ESL_FILEPARSER *efp;
  SSPostscript_t *ps;
  int             found_match = FALSE;

  if (esl_fileparser_Open(filename, NULL, &efp) != eslOK) esl_fatal("ERROR, failed to open template file %s in parse_template_file\n", filename);
  esl_fileparser_SetCommentChar(efp, '#');

  status = eslOK;
  while((found_match == FALSE) && (status == eslOK)) {
    status = parse_template_page(efp, go, errbuf, &ps);
    if((status != eslOK) && (status != eslEOF)) { 
      esl_fileparser_Close(efp);
      return status;
    }
    if(ps->rflen == msa_rflen) { found_match = TRUE; }
    else                     { free_sspostscript(ps); }
  }
  if(found_match == FALSE) { 
    esl_fileparser_Close(efp);
    esl_fatal("ERROR, did not find template structure to match alignment consensus length of %d in:\n%s\n", msa_rflen, filename);
  }

  /* if we get here, we've found a match */

  /* validate the template we just read */
  if((status = validate_justread_sspostscript(ps, errbuf)) != eslOK) return status;

  esl_fileparser_Close(efp);
  *ret_ps = ps;

  return eslOK;
}
  

/* parse_template_page
 *
 * Read a single secondary structure template page for a
 * single CM from a postscript template file derived from the 
 * Gutell CRW website. This function is called repeatedly
 * on the same file to read multiple structure templates.
 * The logic is to keep reading until we find one with the
 * same consensus length as our input alignment. 
 * 
 * The structure template is read in sections.
 * Each section begins with a line like this: 
 * % begin <type1> <type2> 
 * 
 * list of valid tokens for <type1>:
 * modelname
 * scale
 * regurgitate
 * ignore 
 * lines
 * text
 * 
 * if <type1> is lines or text, then <type2> is read, 
 * valid tokens for <type2> if <type1> is 'text'
 * hundreds
 * residues
 * 
 * valid tokens for <type2> if <type1> is 'lines'
 * ticks
 * bpconnects
 * 
 * The 'regurgitate' lines are stored, but never changed.
 * The 'ignore' lines are not even stored.
 * All other <type1> lines are stored in data structures that
 * can be manipulated (though not all are at this point).
 * 
 * If anything is invalid, we return a non-eslOK status code
 * to caller. 
 * 
 * Returns:  eslOK on success.
 */
int
parse_template_page(ESL_FILEPARSER *efp, const ESL_GETOPTS *go, char *errbuf, SSPostscript_t **ret_ps)
{
  int             status;
  char           *tok;
  int             toklen;
  SSPostscript_t *ps = NULL;
  int            read_showpage = FALSE;
  int            reached_eof = FALSE;

  /* Create the postscript object */
  ps = create_sspostscript();

  while ((read_showpage == FALSE) && ((status = esl_fileparser_GetToken(efp, &tok, &toklen))  == eslOK)) {
    if(strcmp(tok, "%") == 0) { 
      if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) == eslOK) { 
	if(strcmp(tok, "begin") == 0) { 
	  if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) == eslOK) { 
	    if(strcmp(tok, "modelname") == 0) { 
	      if((status = parse_modelname_section(efp, errbuf, ps)) != eslOK) return status;
	    }
	    else if(strcmp(tok, "scale") == 0) { 
	      /*printf("parsing scale\n");*/
	      if((status = parse_scale_section(efp, errbuf, ps)) != eslOK) return status;
	    }
	    else if(strcmp(tok, "ignore") == 0) { 
	      /*printf("parsing ignore\n");*/
	      if((status = parse_ignore_section(efp, errbuf, &read_showpage)) != eslOK) return status;
	    }	
	    else if(strcmp(tok, "regurgitate") == 0) { 
	      /*printf("parsing regurgitate\n");*/
	      if((status = parse_regurgitate_section(efp, errbuf, ps)) != eslOK) return status;
	    }	
	    else if(strcmp(tok, "text") == 0) { 
	      /*printf("parsing text\n");*/
	      if((status = parse_text_section(efp, errbuf, ps)) != eslOK) return status;		   
	    }
	    else if(strcmp(tok, "lines") == 0) { 
	      /*printf("parsing lines\n");*/
	      if((status = parse_lines_section(efp, errbuf, ps)) != eslOK) return status;		   
	    }	
	    else { 
	      ESL_FAIL(eslEINVAL, errbuf, "parse_template_page(), error, unknown section type %s.", tok);
	    }
	  }
	  else { 
	    ESL_FAIL(eslEINVAL, errbuf, "parse_template_page(), error last read line number %d.", efp->linenumber);
	  }
	}
	else { 
	  ESL_FAIL(eslEINVAL, errbuf, "parse_template_page(), expected line beginning with %% begin, but read tok: %s instead of begin, last read line number %d.", tok, efp->linenumber);
	}
      }
      else { 
	ESL_FAIL(eslEINVAL, errbuf, "parse_template_page(), ran out of tokens early, error last read line number %d.", efp->linenumber);
      }
    }
    else { 
      ESL_FAIL(eslEINVAL, errbuf, "parse_template_page(), expected line beginning with %%, read tok: %s, last read line number %d.", tok, efp->linenumber);
    }
  }
  if(read_showpage == FALSE && status != eslEOF) { 
    ESL_FAIL(status, errbuf, "parse_template_page(), error, ran out of tokens, but not at end of file?, last read line number %d.", efp->linenumber);
  }
  if(status == eslEOF) reached_eof = TRUE;

  *ret_ps = ps;

  if(reached_eof) return eslEOF;
  else            return eslOK;
}

/* parse_modelname_section
 *
 * Parse the modelname section of a template postscript file.
 * If anything is invalid, we return a non-eslOK status code
 * to caller. 
 * 
 * Returns:  eslOK on success.
 */
int
parse_modelname_section(ESL_FILEPARSER *efp, char *errbuf, SSPostscript_t *ps)
{
  int status;
  char *tok;
  int   toklen;
  char *curstr = NULL;
  int   curlen = 0;
  int  ntok = 0;
  /* this section should be exactly 3 lines, one of which we've already read,
   * here's an example, next token should be the first 0.65
   * % begin modelname
   * % archaea
   * % end scale
   */
  if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing modelname section, reading token 1 of 3"); 
  if (strcmp(tok, "%") != 0)  ESL_FAIL(eslEINVAL, errbuf, "Error, parsing modelname section, middle line token 1 should be a percent sign but it's %s", tok); 
  /* read remainder of line, this is the model name */
  while ((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen))  == eslOK) { 
    if(ntok > 0) {
      if((status = esl_strcat(&curstr, curlen, " ",  1))      != eslOK) ESL_FAIL(status, errbuf, "parse_modelname_section(), error parsing model name.");
      curlen += 1;
    }
    if((status = esl_strcat(&curstr, curlen, tok,  toklen)) != eslOK) ESL_FAIL(status, errbuf, "parse_modelname_section(), error parsing model name.");
    curlen += toklen;
    ntok++;
  }
  ESL_ALLOC(ps->modelname, sizeof(char) * (curlen+1));
  strcpy(ps->modelname, curstr);

  /* next line should be '% end modelname' */
  if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing modelname section, reading end line token 1 of 3"); 
  if (strcmp(tok, "%") != 0)  ESL_FAIL(eslEINVAL, errbuf, "Error, parsing modelname section, end line token 1 of 3 should be a percent sign but it's %s", tok); 
  if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing modelname section, reading end line token 2 of 3"); 
  if (strcmp(tok, "end") != 0)  ESL_FAIL(eslEINVAL, errbuf, "Error, parsing modelname section, end line token 2 of 3 should be 'end' but it's %s", tok); 
  if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing modelname section, reading end line token 3 of 3"); 
  if (strcmp(tok, "modelname") != 0)  ESL_FAIL(eslEINVAL, errbuf, "Error, parsing modelname section, end line token 3 of 3 should be 'modelname' but it's %s", tok); 

  if(curstr != NULL) free(curstr);
  return eslOK;

 ERROR: ESL_FAIL(status, errbuf, "Error, parsing modelname section, memory error?");
}


/* parse_scale_section
 *
 * Parse the scale section of a template postscript file.
 * If anything is invalid, we return a non-eslOK status code
 * to caller. 
 * 
 * Returns:  eslOK on success.
 */
int
parse_scale_section(ESL_FILEPARSER *efp, char *errbuf, SSPostscript_t *ps)
{
  int status;
  char *tok;
  int   toklen;

  /* this section should be exactly 3 lines, one of which we've already read,
   * here's an example, next token should be the first 0.65
   * % begin scale
   * 0.65 0.65 scale
   * % end scale
   */
  if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing scale section, reading token 1 of 3"); 
  ps->scale = atof(tok);
  if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing scale section, reading token 2 of 3"); 
  if(esl_FCompare(ps->scale, atof(tok), eslSMALLX1) != eslOK) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing scale section, x and y scales are not equal %.2f != %.2f", ps->scale, atof(tok)); 
  if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK)  ESL_FAIL(status, errbuf, "Error, parsing scale section, reading token 3 of 3"); 
  if (strcmp(tok, "scale") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing scale section, token 3 of 3 should be 'scale' but it's %s", tok); 

  if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing scale section, reading end line token 1 of 3"); 
  if (strcmp(tok, "%") != 0)  ESL_FAIL(eslEINVAL, errbuf, "Error, parsing scale section, end line token 1 of 3 should be a percent sign but it's %s", tok); 
  if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing scale section, reading end line token 2 of 3"); 
  if (strcmp(tok, "end") != 0)  ESL_FAIL(eslEINVAL, errbuf, "Error, parsing scale section, end line token 2 of 3 should be 'end' but it's %s", tok); 
  if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing scale section, reading end line token 3 of 3"); 
  if (strcmp(tok, "scale") != 0)  ESL_FAIL(eslEINVAL, errbuf, "Error, parsing scale section, end line token 3 of 3 should be 'scale' but it's %s", tok); 

  return eslOK;
}


/* parse_ignore_section
 *
 * Parse an ignore section of a template postscript file.
 * We ignore this data. This function's purpose is to read
 * tokens until we see the "% end ignore" line signalling
 * the end of the ignore section.
 * 
 * Returns:  eslOK on success.
 */
int
parse_ignore_section(ESL_FILEPARSER *efp, char *errbuf, int *ret_read_showpage)
{
  int status;
  char *tok;
  int   toklen;
  int   keep_reading = TRUE;
  int   read_showpage = FALSE;
  while((keep_reading) && (status = esl_fileparser_GetToken(efp, &tok, &toklen)) == eslOK) { 
    if (strcmp(tok, "%") == 0) { 
      if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing ignore section, read %% prefixed line without ' end ignore' after it"); 
      if (strcmp(tok, "end") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing ignore section, read %% prefixed line without ' end ignore' after it"); 
      if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing ignore section, read %% prefixed line without ' end ignore' after it"); 
      if (strcmp(tok, "ignore") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing ignore section, read %% prefixed line without ' end ignore' after it"); 
      keep_reading = FALSE;
      status = eslOK;
    }
    else if(strcmp(tok, "showpage") == 0) { 
      read_showpage = TRUE;
    }
  }
  if(status == eslEOF) ESL_FAIL(status, errbuf, "Error, parsing ignore section, finished file looking for '%% end ignore' line");
  if(status != eslOK)  ESL_FAIL(status, errbuf, "Error, parsing ignore section, last line number read %d", efp->linenumber);

  *ret_read_showpage = read_showpage;
  return eslOK;
}


/* parse_regurgitate_section
 *
 * Parse a regurgitate section of a template postscript file.
 * If anything is invalid, we return a non-eslOK status code
 * to caller. 
 * 
 * Returns:  eslOK on success.
 */
int
parse_regurgitate_section(ESL_FILEPARSER *efp, char *errbuf, SSPostscript_t *ps)
{
  int status;
  char *tok;
  int   toklen;
  int   seen_end = FALSE;
  char *curstr = NULL;
  char *newstr;
  int   nalloc = ps->nregurg;
  int   curlen;
  void *tmp;
  int   ntok = 0;

  while (((status = esl_fileparser_NextLine(efp)) == eslOK) && (!seen_end))
  {
    curlen = ntok = 0;
    while ((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen))  == eslOK) { 
      if (strcmp(tok, "%") == 0) { /* should be the end, make sure it's properly formatted */
	if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing regurgitate section, read %% prefixed line without ' end regurgitate' after it"); 
	if (strcmp(tok, "end") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing regurgitate section, read %% prefixed line without ' end regurgitate' after it"); 
	if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing regurgitate section, read %% prefixed line without ' end regurgitate' after it"); 
	if (strcmp(tok, "regurgitate") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing regurgitate section, read %% prefixed line without ' end regurgitate' after it"); 
	seen_end = TRUE;
	break;
      }
      else { 
	if(ntok > 0) { 
	  if((status = esl_strcat(&curstr, curlen, " ",  1))  != eslOK) ESL_FAIL(status, errbuf, "parse_regurgitate_section(), error (2).");
	  curlen += 1;
	}
	if((status = esl_strcat(&curstr, curlen, tok,  toklen)) != eslOK) ESL_FAIL(status, errbuf, "parse_regurgitate_section(), error (1).");
	curlen += toklen;
	ntok++;
      }
    }
    if(seen_end) break;
    if((status = esl_strcat(&curstr, curlen, "\n", 1)) != eslOK) ESL_FAIL(status, errbuf, "parse_regurgitate_section(), error (3).");
    curlen += 1;
    ESL_ALLOC(newstr, sizeof(char) * (curlen+1));
    strcpy(newstr, curstr);
    if(ps->nregurg == nalloc) { 
      nalloc += ps->nalloc; ESL_RALLOC(ps->regurgA, tmp, sizeof(char *) * nalloc); 
    }
    ps->regurgA[ps->nregurg++] = newstr;
    free(curstr);
    curstr = NULL;
  }
  if(status == eslEOF) ESL_FAIL(status, errbuf, "Error, parsing regurgitate section, finished file looking for '%% end regurgitate' line");
  if(status != eslOK)  ESL_FAIL(status, errbuf, "Error, parsing regurgitate section, last line number read %d", efp->linenumber);

  return eslOK;

 ERROR: ESL_FAIL(status, errbuf, "Memory error parsing regurgitate section");
}


/* parse_text_section
 *
 * Parse a text section of a template postscript file.
 * If anything is invalid, we return a non-eslOK status code
 * to caller. 
 * 
 * Returns:  eslOK on success.
 */
int
parse_text_section(ESL_FILEPARSER *efp, char *errbuf, SSPostscript_t *ps)
{
  int status;
  char *tok;
  int   toklen;
  int   seen_end = FALSE;
  int   nalloc;
  void *tmp;
  int do_hundreds = FALSE;
  int do_residues = FALSE;

  /* find out which section we're in, 'hundreds' or 'residues' */
  if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing text section, last line %d\n", efp->linenumber);
  if      (strcmp(tok, "hundreds") == 0) { do_hundreds = TRUE; nalloc = ps->nhundreds; }
  else if (strcmp(tok, "residues") == 0) { do_residues = TRUE; nalloc = ps->rflen; }

  /* read the first two special lines, should be a 5-token line ending with setfont and a 5-token line ending with setcmykcolor,
   * we don't store these, but we require that they're there. */
  if((status = esl_fileparser_NextLine(efp) != eslOK)) ESL_FAIL(status, errbuf, "Error, parsing text section, last line %d\n", efp->linenumber);
  if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing text section first line should be 5-tokens ending with 'setfont'");
  if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing text section first line should be 5-tokens ending with 'setfont'");
  if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing text section first line should be 5-tokens ending with 'setfont'");
  if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing text section first line should be 5-tokens ending with 'setfont'");
  if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing text section first line should be 5-tokens ending with 'setfont'");
  if(strcmp(tok, "setfont") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing text section first line should be 5-tokens ending with 'setfont'");

  if((status = esl_fileparser_NextLine(efp) != eslOK)) ESL_FAIL(status, errbuf, "Error, parsing text section, last line %d\n", efp->linenumber);
  if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing text section second line should be 5-tokens ending with 'setcmykcolor'");
  if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing text section second line should be 5-tokens ending with 'setcmykcolor'");
  if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing text section second line should be 5-tokens ending with 'setcmykcolor'");
  if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing text section second line should be 5-tokens ending with 'setcmykcolor'");
  if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing text section second line should be 5-tokens ending with 'setcmykcolor'");
  if(strcmp(tok, "setcmykcolor") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing text section second line should be 5-tokens ending with 'setcmykcolor'");

  while (((status = esl_fileparser_NextLine(efp)) == eslOK) && (!seen_end))
  {
    if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing text main section should include 5-tokens ending with 'show'");
    if (strcmp(tok, "%") == 0) { /* should be the end, make sure it's properly formatted */
      if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing text section, read %% prefixed line without ' end text' after it"); 
      if (strcmp(tok, "end") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing text section, read %% prefixed line without ' end text' after it"); 
      if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing text section, read %% prefixed line without ' end text' after it"); 
      if (strcmp(tok, "text") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing text section, read %% prefixed line without ' end text' after it"); 
      if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing text section, read %% prefixed line without ' end text' after it"); 
      if(do_hundreds) { 
	if (strcmp(tok, "hundreds") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing text section, read %% prefixed line without ' end text hundreds' after it"); 
      }
      if(do_residues) {
	if (strcmp(tok, "residues") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing text section, read %% prefixed line without ' end text residues' after it"); 
      }
      seen_end = TRUE;
      break;
    }
    /* if we get here, we haven't seen the end, we're reading a normal line, tok is the string, we discard this */
    if(do_hundreds && ps->nhundreds == nalloc) { 
      nalloc += ps->nalloc; 
      ESL_RALLOC(ps->hundredsxA, tmp, sizeof(float) * nalloc); 
      ESL_RALLOC(ps->hundredsyA, tmp, sizeof(float) * nalloc); 
    }
    if(do_residues && ps->rflen == nalloc) { 
      nalloc += ps->nalloc; 
      ESL_RALLOC(ps->rxA, tmp, sizeof(float) * nalloc); 
      ESL_RALLOC(ps->ryA, tmp, sizeof(float) * nalloc); 
    }
    /* get x */
    if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing text main section should include 5-tokens ending with 'show'");
    if(do_hundreds) ps->hundredsxA[ps->nhundreds] = atof(tok);
    if(do_residues) ps->rxA[ps->rflen] = atof(tok);
    /* get y */
    if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing text main section should include 5-tokens ending with 'show'");
    if(do_hundreds) ps->hundredsyA[ps->nhundreds] = atof(tok);
    if(do_residues) ps->ryA[ps->rflen] = atof(tok);
    
    /* verify moveto */
    if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing text main section should include 5-tokens ending with 'show'");
    if (strcmp(tok, "moveto") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing text main section, fourth token should be 'moveto', line %d", efp->linenumber);
    /* verify show */
    if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing text main section should include 5-tokens ending with 'show'");
    if (strcmp(tok, "show") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing text main section, fifth token should be 'show', line %d", efp->linenumber);
    
    if(do_hundreds) ps->nhundreds++;
    if(do_residues) ps->rflen++;
  }
  if(!seen_end) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing text section, didn't see end! line: %d\n", efp->linenumber);
  if(status == eslEOF && do_hundreds) 
    ESL_FAIL(status, errbuf, "Error, parsing text section, finished file looking for '%% end text hundreds' line");
  if(status == eslEOF && do_residues) 
    ESL_FAIL(status, errbuf, "Error, parsing text section, finished file looking for '%% end text residues' line");
  if(status != eslOK)  ESL_FAIL(status, errbuf, "Error, parsing text section, last line number read %d", efp->linenumber);

  return eslOK;

 ERROR: ESL_FAIL(status, errbuf, "Memory error parsing text section");
}

/* parse_lines_section
 *
 * Parse a lines section of a template postscript file.
 * If anything is invalid, we return a non-eslOK status code
 * to caller. 
 * 
 * Returns:  eslOK on success.
 */
int
parse_lines_section(ESL_FILEPARSER *efp, char *errbuf, SSPostscript_t *ps)
{
  int status;
  char *tok;
  int   toklen;
  int   seen_end = FALSE;
  int   nalloc;
  void *tmp;
  int do_ticks = FALSE;
  int do_bpconnects = FALSE;

  /* find out which section we're in, 'ticks' or 'bpconnects' */
  if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing lines section, last line %d\n", efp->linenumber);
  if      (strcmp(tok, "ticks") == 0)      { do_ticks = TRUE; nalloc = ps->nticks; }
  else if (strcmp(tok, "bpconnects") == 0) { do_bpconnects = TRUE; nalloc = ps->nbp; }

  /* read the first two special lines, should be a 2-token line ending with setlinewidth and a 5-token line ending with setcmykcolor,
   * we don't store these, but we require that they're there. */
  if((status = esl_fileparser_NextLine(efp) != eslOK)) ESL_FAIL(status, errbuf, "Error, parsing lines section, last line %d\n", efp->linenumber);
  if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing lines section first line should be 2-tokens ending with 'setlinewidth'");
  if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing lines section first line should be 2-tokens ending with 'setlinewidth'");
  if(strcmp(tok, "setlinewidth") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing lines section first line should be 2-tokens ending with 'setlinewidth'");

  if((status = esl_fileparser_NextLine(efp) != eslOK)) ESL_FAIL(status, errbuf, "Error, parsing lines section, last line %d\n", efp->linenumber);
  if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing lines section second line should be 5-tokens ending with 'setcmykcolor'");
  if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing lines section second line should be 5-tokens ending with 'setcmykcolor'");
  if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing lines section second line should be 5-tokens ending with 'setcmykcolor'");
  if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing lines section second line should be 5-tokens ending with 'setcmykcolor'");
  if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing lines section second line should be 5-tokens ending with 'setcmykcolor'");
  if(strcmp(tok, "setcmykcolor") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing lines section second line should be 5-tokens ending with 'setcmykcolor'");

  while (((status = esl_fileparser_NextLine(efp)) == eslOK) && (!seen_end))
  {
    if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing lines main section should include 5-tokens ending with 'show'");
    if (strcmp(tok, "%") == 0) { /* should be the end, make sure it's properly formatted */
      if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing lines section, read %% prefixed line without ' end lines' after it"); 
      if (strcmp(tok, "end") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing lines section, read %% prefixed line without ' end lines' after it"); 
      if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing lines section, read %% prefixed line without ' end lines' after it"); 
      if (strcmp(tok, "lines") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing lines section, read %% prefixed line without ' end lines' after it"); 
      if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing lines section, read %% prefixed line without ' end lines' after it"); 
      if(do_ticks) { 
	if (strcmp(tok, "ticks") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing lines section, read %% prefixed line without ' end lines ticks' after it"); 
      }
      if(do_bpconnects) {
	if (strcmp(tok, "bpconnects") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing lines section, read %% prefixed line without ' end lines bpconnects' after it"); 
      }
      seen_end = TRUE;
      break;
    }
    /* if we get here, we haven't seen the end, we're reading a normal line, tok is the first x coord, we record this */
    /* first we expand our arrays if nec */
    if(do_ticks && ps->nticks == nalloc) { 
      nalloc += ps->nalloc; 
      ESL_RALLOC(ps->ticksx1A, tmp, sizeof(float) * nalloc); 
      ESL_RALLOC(ps->ticksy1A, tmp, sizeof(float) * nalloc); 
      ESL_RALLOC(ps->ticksx2A, tmp, sizeof(float) * nalloc); 
      ESL_RALLOC(ps->ticksy2A, tmp, sizeof(float) * nalloc); 
    }
    if(do_bpconnects && ps->nbp == nalloc) { 
      nalloc += ps->nalloc; 
      ESL_RALLOC(ps->bpx1A, tmp, sizeof(float) * nalloc); 
      ESL_RALLOC(ps->bpy1A, tmp, sizeof(float) * nalloc); 
      ESL_RALLOC(ps->bpx2A, tmp, sizeof(float) * nalloc); 
      ESL_RALLOC(ps->bpy2A, tmp, sizeof(float) * nalloc); 
    }
    /* store x1 */
    if(do_ticks) ps->ticksx1A[ps->nticks] = atof(tok);
    if(do_bpconnects) ps->bpx1A[ps->nbp] = atof(tok);
    /* get y1 */
    if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing lines main section should include 8-tokens ending with 'stroke'");
    if(do_ticks) ps->ticksy1A[ps->nticks] = atof(tok);
    if(do_bpconnects) ps->bpy1A[ps->nbp] = atof(tok);
    /* get x2 */
    if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing lines main section should include 8-tokens ending with 'stroke'");
    if(do_ticks) ps->ticksx2A[ps->nticks] = atof(tok);
    if(do_bpconnects) ps->bpx2A[ps->nbp] = atof(tok);
    /* get y2 */
    if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing lines main section should include 8-tokens ending with 'stroke'");
    if(do_ticks) ps->ticksy2A[ps->nticks] = atof(tok);
    if(do_bpconnects) ps->bpy2A[ps->nbp] = atof(tok);
    
    /* verify newpath */
    if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing lines main section should include 8-tokens ending with 'stroke'");
    if (strcmp(tok, "newpath") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing lines main section, fifth token should be 'newpath', line %d", efp->linenumber);
    /* verify moveto */
    if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing lines main section should include 8-tokens ending with 'stroke'");
    if (strcmp(tok, "moveto") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing lines main section, sixth token should be 'moveto', line %d", efp->linenumber);
    /* verify lineto */
    if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing lines main section should include 8-tokens ending with 'stroke'");
    if (strcmp(tok, "lineto") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing lines main section, seventh token should be 'lineto', line %d", efp->linenumber);
    /* verify stroke */
    if((status = esl_fileparser_GetTokenOnLine(efp, &tok, &toklen)) != eslOK) ESL_FAIL(status, errbuf, "Error, parsing lines main section should include 8-tokens ending with 'stroke'");
    if (strcmp(tok, "stroke") != 0) ESL_FAIL(eslEINVAL, errbuf, "Error, parsing lines main section, eigth token should be 'stroke', line %d", efp->linenumber);
    
    if(do_ticks) ps->nticks++;
    if(do_bpconnects) ps->nbp++;
  }
  if(!seen_end) ESL_FAIL(status, errbuf, "Error, parsing lines section, didn't see end! line: %d\n", efp->linenumber);
  if(status == eslEOF && do_ticks) 
    ESL_FAIL(status, errbuf, "Error, parsing lines section, finished file looking for '%% end lines ticks' line");
  if(status == eslEOF && do_bpconnects) 
    ESL_FAIL(status, errbuf, "Error, parsing lines section, finished file looking for '%% end lines bpconnects' line");

  if(status != eslOK)  ESL_FAIL(status, errbuf, "Error, parsing lines section, last line number read %d", efp->linenumber);

  return eslOK;
 ERROR: ESL_FAIL(status, errbuf, "Memory error parsing lines section");
}

/* Function: individual_seqs_sspostscript()
 * 
 * Purpose:  Fill a postscript data structure with info for individual seqs in the MSA 
 * Return:   eslOK on success.
 * 
 * 
 * ins_ct - [0..i..msa->nseq-1][0..rflen] number of inserts 
 *          insert after each position per sequence. 
 */
static int
individual_seqs_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, int **ins_ct, int *useme, int nused, float ***hc_scheme, int hc_scheme_idx, int hc_nbins, float **hc_onecell, int zeroins_idx, int extdel_idx)
{
  int status;
  int p, i, pp, ai;
  int rfpos, apos;
  int orig_npage = ps->npage;
  int nzeroins = 0;  /* number of positions with 0 inserts after them */
  int nextdel = 0;   /* number of external deletes (3' and 5' flush deletes) */
  float *limits;     /* [nbins+1] limits for each bin, limits[0] is min value we would expect to see, limits[nbins] is max */
  int nins;          /* number of inserts after the current rfpos for current sequence */
  int spos, epos;    /* first/final nongap position for cur sequence */

  if((status = add_pages_sspostscript(ps, nused, INDIMODE)) != eslOK) ESL_FAIL(status, errbuf, "memory error adding pages to the postscript object.");

  for(p = orig_npage; p < ps->npage; p++) { 
    ps->rrAA[p]    = NULL;
    ps->rcolAAA[p] = NULL;
    ps->sclAA[p]   = NULL;
    ps->occlAAA[p] = NULL;
  }

  for(p = orig_npage; p < ps->npage; p++) { 
    ESL_ALLOC(ps->rrAA[p], sizeof(char) *  (ps->rflen+1));
    ESL_ALLOC(ps->rcolAAA[p], sizeof(float *) * ps->rflen);
    ESL_ALLOC(ps->sclAA[p],    sizeof(SchemeColorLegend_t *) * 1);
    ESL_ALLOC(ps->occlAAA[p], sizeof(OneCellColorLegend_t **) * 2);
    for(rfpos = 0; rfpos < ps->rflen; rfpos++) 
      ps->rcolAAA[p][rfpos] = NULL;
    for(rfpos = 0; rfpos < ps->rflen; rfpos++) { 
      ESL_ALLOC(ps->rcolAAA[p][rfpos], sizeof(float) * NCMYK); /* CMYK colors */
    }
  }

  ESL_ALLOC(limits, sizeof(float) * (hc_nbins+1)); 
  limits[0] = 1;
  limits[1] = 2;
  limits[2] = 4;
  limits[3] = 6;
  limits[4] = 10;
  limits[5] = SSDRAWINFINITY;

  /* allocate the uaseqlenA data structure to store unaligned seq lengths */
  if(ps->uaseqlenA != NULL) { free(ps->uaseqlenA); ps->uaseqlenA = NULL; }
  ESL_ALLOC(ps->uaseqlenA, sizeof(int) * msa->nseq);
  esl_vec_ISet(ps->uaseqlenA, msa->nseq, 0);

  /* fill ps->rrAA with residues and gaps */
  ai = 0;
  for(i = 0; i < msa->nseq; i++) {
    if(useme[i]) { 
      pp = orig_npage + ai;

      spos = epos = -1;
      /* determine first and final non-gap position */
      for(apos = 0; apos < msa->alen; apos++) { /* find first non-gap RF position */
	if((! esl_abc_CIsGap(msa->abc, msa->aseq[i][apos]))) { 
	  spos = apos;
	  break;
	}
      }
      for(apos = msa->alen-1; apos >= 0; apos--) { 
	if((! esl_abc_CIsGap(msa->abc, msa->aseq[i][apos]))) { 
	  epos = apos;
	  break;
	}
      }
      ps->sclAA[pp] = create_scheme_colorlegend(hc_scheme_idx, hc_nbins, limits, TRUE);
      nextdel = 0;
      nzeroins = 0;
      ai++;
      for(rfpos = 0; rfpos < ps->rflen; rfpos++) { 
	apos = ps->msa_rf2a_map[rfpos];
	nins = ins_ct[i][rfpos];

	if(! esl_abc_CIsGap(msa->abc, msa->aseq[i][apos])) ps->uaseqlenA[i]++;
	ps->uaseqlenA[i] += nins;
	ps->rrAA[pp][rfpos] = msa->aseq[i][apos];
	/* printf("ps->rrAA[%3d][%4d]: %c\n", pp, rfpos, ps->rrAA[pp][rfpos]);  */
	if((spos != -1 && epos != -1) && /* this should always be true, unless seq has length 0! */
	   (apos < spos || apos > epos)) { 
	  if((status = set_onecell_values(errbuf, ps->rcolAAA[pp][rfpos], NCMYK, hc_onecell[extdel_idx])) != eslOK) return status;
	  nextdel++;
	}
	else if(nins == 0) { 
	    if((status = set_onecell_values(errbuf, ps->rcolAAA[pp][rfpos], NCMYK, hc_onecell[zeroins_idx])) != eslOK) return status;
	    nzeroins++;
	}
	else { 
	  if((status = set_scheme_values(errbuf, ps->rcolAAA[pp][rfpos], NCMYK, hc_scheme[hc_scheme_idx], nins, ps->sclAA[pp], TRUE, NULL)) != eslOK) return status;
	}	  
      }

      ps->rrAA[pp][ps->rflen] = '\0';
      ps->seqidxA[pp] = i;

      /* add one-cell color legend for zero inserts */
      ps->occlAAA[pp][0] = create_onecell_colorlegend(hc_onecell[zeroins_idx], nzeroins, -1);
      if((status = add_text_to_onecell_colorlegend(ps, ps->occlAAA[pp][0], "(blank) zero inserts", ps->legx_max_chars, errbuf)) != eslOK) return status;

      /* add one-cell color legend for external gaps (deletes) */
      ps->occlAAA[pp][1] = create_onecell_colorlegend(hc_onecell[extdel_idx], nextdel, -1);
      if((status = add_text_to_onecell_colorlegend(ps, ps->occlAAA[pp][1], "5'/3'-flush gaps", ps->legx_max_chars, errbuf)) != eslOK) return status;
      ps->nocclA[pp] = 2;

      /* add description to ps */
      if((status = add_text_to_scheme_colorlegend(ps->sclAA[pp], "# inserted residues after each consensus position", ps->legx_max_chars, errbuf)) != eslOK) return status;
      if((status = add_page_desc_to_sspostscript(ps, pp, msa->sqname[i], errbuf)) != eslOK) return status;
    }
  }
  free(limits);
  return eslOK;

 ERROR: ESL_FAIL(status, errbuf, "individual_seqs_sspostscript(): memory allocation error.");
  return status; /* NEVERREACHED */
}

/* Function: rf_seq_sspostscript()
 * 
 * Purpose:  Fill a postscript data structure with 1 new page, the RF sequence.
 * Return:   eslOK on success.
 */
static int
rf_seq_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa)
{
  int status;
  int p, pp;
  int cpos, apos;
  int orig_npage = ps->npage;

  if((status = add_pages_sspostscript(ps, 1, INDIMODE)) != eslOK) ESL_FAIL(status, errbuf, "memory error adding pages to the postscript object.");

  for(p = orig_npage; p < ps->npage; p++) { 
    ESL_ALLOC(ps->rrAA[p], sizeof(char) *  ps->rflen);
  }

  /* fill ps->rrAA with residues and gaps for RF sequence */
  pp = orig_npage;
  cpos = 0;
  for(apos = 0; apos < msa->alen; apos++) {
    if((! esl_abc_CIsGap(msa->abc, msa->rf[apos])) && 
       (! esl_abc_CIsMissing(msa->abc, msa->rf[apos])) && 
       (! esl_abc_CIsNonresidue(msa->abc, msa->rf[apos]))) { 
      ps->rrAA[pp][cpos] = msa->rf[apos];
      /* printf("ps->rrAA[%3d][%4d]: %c\n", pp, cpos, ps->rrAA[pp][cpos]); */
      cpos++;
    }
  }

  /* add description to ps */
  if((status = add_page_desc_to_sspostscript(ps, pp, "*CONSENSUS*", errbuf)) != eslOK) return status;

  return eslOK;

 ERROR: ESL_FAIL(status, errbuf, "rf_seq_sspostscript(): memory allocation error.");
  return status; /* NEVERREACHED */
}

/* count_msa()
 *                   
 * Given an msa, count residues, post probs, basepairs, and consensus start and end 
 * positions and store them in <ret_abc_ct>, <ret_pp_ct>, <ret_bp_ct>, 
 * <ret_srfpos_ct>, and <ret_erfpos_ct>.
 * 
 * <ret_abc_ct> [0..apos..alen-1][0..abc->K]:
 * - per position count of each symbol in alphabet over all seqs.
 * 
 * <ret_bp_ct>  [0..apos..alen-1][0..abc->Kp-1][0..abc->Kp-1] 
 * - per (non-pknotted) consensus basepair count of each possible basepair 
 *   over all seqs basepairs are indexed by 'i' the minimum of 'i:j' for a 
 *   pair between i and j, where i < j. Note that non-canonicals and 
 *   gaps and the like are all stored independently.
 *
 * <ret_pp_ct> [0..apos..alen-1][0..10]
 * - per position count of each posterior probability code over all seqs.
 * 
 * <ret_srfpos_ct> [0..rfpos..rflen-1]
 * - per position count of first non-gap RF position over all seqs.
 *                                
 * <ret_erfpos_ct> [0..rfpos..rflen-1]
 * - per position count of final non-gap RF position over all seqs.
 *
 * A 'gap' has a looser defintion than in esl_abc here, esl_abc's gap, 
 * missing residues and nonresidues are all considered 'gaps' here.
 * 
 * If we encounter an error, we return non-eslOK status and fill
 * errbuf with error message.
 * 
 * Returns eslOK upon success.
 */
int count_msa(ESL_MSA *msa, char *errbuf, int *a2rf_map, int rflen, double ***ret_abc_ct, double ****ret_bp_ct, int ***ret_pp_ct, int **ret_srfpos_ct, int **ret_erfpos_ct)
{
  int status;
  double  **abc_ct = NULL;
  double ***bp_ct = NULL;
  int      *srfpos_ct = NULL;
  int      *erfpos_ct = NULL;
  int       apos, i, j, x, epos, rfpos;
  ESL_DSQ  *tmp_dsq = NULL;
  int       seen_start = FALSE;
  /* variables related to getting bp counts */
  int      *ct = NULL;            /* 0..alen-1 base pair partners array for current sequence */
  char     *ss_nopseudo = NULL;   /* no-pseudoknot version of structure */
  int       nppvals = 12;         /* '0'-'9' = 0-9, '*' = 10, gap = '11' */
  int     **pp_ct = NULL;         /* [0..alen-1][0..nppvals-1] per position count of each possible PP char over all seqs */
  int       ppidx; 

  /* contract check, msa should be in text mode and have ss_cons */
  if(msa->flags & eslMSA_DIGITAL) ESL_FAIL(eslEINVAL, errbuf, "count_msa() contract violation, MSA is digitized");
  if(msa->ss_cons == NULL) ESL_FAIL(eslEINVAL, errbuf, "the alignment lacks SS_cons annotation");
  if(ret_pp_ct != NULL && msa->pp == NULL) ESL_FAIL(eslEINVAL, errbuf, "--prob requires all sequences in the alignment have PP, but none do.");

  /* allocate pp_ct array, if nec */
  if(ret_pp_ct != NULL) { 
    ESL_ALLOC(pp_ct, sizeof(int *) * msa->alen);
    for(apos = 0; apos < msa->alen; apos++) { 
      ESL_ALLOC(pp_ct[apos], sizeof(int) * nppvals);
      esl_vec_ISet(pp_ct[apos], nppvals, 0);
    }
  }

  /* get ct array which defines the consensus base pairs */
  ESL_ALLOC(ct,  sizeof(int)  * (msa->alen+1));
  ESL_ALLOC(ss_nopseudo, sizeof(char) * (msa->alen+1));
  esl_wuss_nopseudo(msa->ss_cons, ss_nopseudo);
  if ((status = esl_wuss2ct(ss_nopseudo, msa->alen, ct)) != eslOK) ESL_FAIL(status, errbuf, "Consensus structure string is inconsistent.");

  ESL_ALLOC(tmp_dsq, (msa->alen+2) * sizeof(ESL_DSQ));
  ESL_ALLOC(abc_ct, sizeof(double *) * msa->alen); 
  ESL_ALLOC(bp_ct,  sizeof(double **) * msa->alen); 
  for(apos = 0; apos < msa->alen; apos++) { 
    ESL_ALLOC(abc_ct[apos], sizeof(double) * (msa->abc->K+1));
    esl_vec_DSet(abc_ct[apos], (msa->abc->K+1), 0.);
    /* careful ct is indexed 1..alen, not 0..alen-1 */
    if(ct[(apos+1)] > (apos+1)) { /* apos+1 is an 'i' in an i:j pair, where i < j */
      ESL_ALLOC(bp_ct[apos], sizeof(double *) * (msa->abc->Kp));
      for(x = 0; x < msa->abc->Kp; x++) { 
	ESL_ALLOC(bp_ct[apos][x], sizeof(double) * (msa->abc->Kp));
	esl_vec_DSet(bp_ct[apos][x], msa->abc->Kp, 0.);
      }
    }
    else { /* apos+1 is not an 'i' in an i:j pair, where i < j, set to NULL */
      bp_ct[apos] = NULL;
    }
  }
  ESL_ALLOC(srfpos_ct, sizeof(int) * rflen);
  ESL_ALLOC(erfpos_ct, sizeof(int) * rflen);
  esl_vec_ISet(srfpos_ct, rflen, 0);
  esl_vec_ISet(erfpos_ct, rflen, 0);

  for(i = 0; i < msa->nseq; i++) { 
    seen_start = FALSE;
    if((status = esl_abc_Digitize(msa->abc, msa->aseq[i], tmp_dsq)) != eslOK) ESL_FAIL(status, errbuf, "problem digitizing sequence %d", i);
    rfpos = 0;
    for(apos = 0; apos < msa->alen; apos++) { /* update appropriate abc count, careful, tmp_dsq ranges from 1..msa->alen (not 0..msa->alen-1) */
      if((status = esl_abc_DCount(msa->abc, abc_ct[apos], tmp_dsq[apos+1], 1.0)) != eslOK) ESL_FAIL(status, errbuf, "problem counting residue %d of seq %d", apos, i);
      if(a2rf_map[apos] != -1) { 
	if(! esl_abc_XIsGap(msa->abc, tmp_dsq[apos+1])) { 
	  if(! seen_start) { 
	    srfpos_ct[rfpos]++; 
	    seen_start = TRUE;
	  }
	  epos = rfpos;
	}
	rfpos++;
      }
      /* get bp count, if nec */
      if(bp_ct[apos] != NULL) { /* our flag for whether position (apos+1) is an 'i' in an i:j pair where i < j */
	j = ct[apos+1] - 1; /* ct is indexed 1..alen */
	bp_ct[apos][tmp_dsq[(apos+1)]][tmp_dsq[(j+1)]]++;
      }
    }
    erfpos_ct[epos]++;
    /* get PP counts, if nec  */
    if(ret_pp_ct != NULL) { 
      if(msa->pp[i] == NULL) ESL_FAIL(eslEINVAL, errbuf, "--prob requires all sequences in the alignment have PP, seq %d does not.", i+1);
      for(apos = 0; apos < msa->alen; apos++) { /* update appropriate pp count, careful, tmp_dsq ranges from 1..msa->alen (not 0..msa->alen-1) */
	if((ppidx = get_pp_idx(msa->abc, msa->pp[i][apos])) == -1) ESL_FAIL(eslEFORMAT, errbuf, "bad #=GR PP char: %c", msa->pp[i][apos]);
	pp_ct[apos][ppidx]++;
      }
    }
  }

  *ret_abc_ct  = abc_ct;
  *ret_bp_ct   = bp_ct;
  *ret_srfpos_ct = srfpos_ct;
  *ret_erfpos_ct = erfpos_ct;
  if(ret_pp_ct != NULL) *ret_pp_ct = pp_ct; /* we only allocated pp_ct if ret_pp_ct != NULL */

  if(tmp_dsq != NULL) free(tmp_dsq);
  if(ss_nopseudo != NULL) free(ss_nopseudo);
  if(ct != NULL) free(ct);
  return eslOK;

 ERROR:
  if(abc_ct != NULL)  esl_Free2D((void **) abc_ct, msa->alen);
  if(pp_ct != NULL)   esl_Free2D((void **) pp_ct, msa->alen);
  if(bp_ct  != NULL)  esl_Free3D((void ***) bp_ct, msa->alen, msa->abc->Kp);
  if(srfpos_ct != NULL) free(srfpos_ct);
  if(erfpos_ct != NULL) free(erfpos_ct);
  if(tmp_dsq != NULL) free(tmp_dsq);
  ESL_FAIL(status, errbuf, "Error, out of memory while counting important values in the msa.");
  return status; /* NEVERREACHED */
}

/* Function: infocontent_sspostscript()
 * 
 * Purpose:  Fill a postscript data structure with 1 new page, colored squares indicating
 *           the information content of each consensus column.
 *           
 * Return:   eslOK on success.
 */
static int
infocontent_sspostscript(const ESL_GETOPTS *go, ESL_ALPHABET *abc, char *errbuf, SSPostscript_t *ps, double **abc_ct, int msa_nseq, float ***hc_scheme, int hc_scheme_idx, int hc_nbins, float **hc_onecell, int hc_onecell_idx, FILE *tabfp)
{
  int status;
  int p, pp, c;
  int rfpos, apos;
  int orig_npage = ps->npage;
  double  *tmp_obs = NULL;
  double  *ent   = NULL;
  double  *bg    = NULL;
  float   *limits = NULL;
  int zero_obs;
  int nonecell = 0;
  int nonecell_masked = 0;
  int within_mask;
  int bi;  /* bin index for current position */
  int l;

  if((status = add_pages_sspostscript(ps, 1, ALIMODE)) != eslOK) ESL_FAIL(status, errbuf, "memory error adding pages to the postscript object.");

  for(p = orig_npage; p < ps->npage; p++) { 
    ps->rrAA[p]    = NULL;
    ps->rcolAAA[p] = NULL;
    ps->sclAA[p]   = NULL;
    ps->occlAAA[p] = NULL;
  }
  for(p = orig_npage; p < ps->npage; p++) { 
    ESL_ALLOC(ps->rrAA[p], sizeof(char) *  ps->rflen);
    ESL_ALLOC(ps->rcolAAA[p], sizeof(float *) * ps->rflen);
    ESL_ALLOC(ps->sclAA[p],    sizeof(SchemeColorLegend_t *) * 1);
    ESL_ALLOC(ps->occlAAA[p], sizeof(OneCellColorLegend_t **) * 1);
    for(c = 0; c < ps->rflen; c++) 
      ps->rcolAAA[p][c] = NULL;
    for(c = 0; c < ps->rflen; c++) { 
      ESL_ALLOC(ps->rcolAAA[p][c], sizeof(float) * NCMYK); /* CMYK colors */
    }
  }

  ESL_ALLOC(ent, sizeof(double) * ps->rflen);
  ESL_ALLOC(bg, sizeof(double) * abc->K);
  esl_vec_DSet(bg, abc->K, 1./(abc->K));
  ESL_ALLOC(tmp_obs, sizeof(double) * abc->K);
  esl_vec_DSet(tmp_obs, abc->K, 0.);

  pp = orig_npage;

  /* add color legend */
  ESL_ALLOC(limits, sizeof(float) * (hc_nbins+1)); 
  limits[0] = 0.0;
  limits[1] = 0.4;
  limits[2] = 0.8;
  limits[3] = 1.2;
  limits[4] = 1.6;
  limits[5] = 1.99;
  limits[6] = 2.00;
  ps->sclAA[pp] = create_scheme_colorlegend(hc_scheme_idx, hc_nbins, limits, FALSE);

  if(tabfp != NULL) { 
    fprintf(tabfp, "# ------------------------\n");
    fprintf(tabfp, "# Information content data\n");
    fprintf(tabfp, "# ------------------------\n");
    fprintf(tabfp, "# This section includes %d non #-prefixed lines, one for each consensus position\n", ps->rflen);
    fprintf(tabfp, "# in the alignment and corresponding template.\n");
    fprintf(tabfp, "# Each line includes %d tokens, separated by whitespace:\n", ps->mask == NULL ? 5 : 6);
    fprintf(tabfp, "# \ttoken 1: 'infocontent' (tag defining line type to ease parsing)\n");
    fprintf(tabfp, "# \ttoken 2: consensus position (starting at 1)\n");
    fprintf(tabfp, "# \ttoken 3: information content for position (bits)\n");
    fprintf(tabfp, "# \ttoken 4: number of non-gap residues in position (max possible is %d (num seqs in aln))\n", msa_nseq); 
    fprintf(tabfp, "# \ttoken 5: bin index this positions falls in (see bin values below)\n");
    if(ps->mask != NULL) { 
      fprintf(tabfp, "# \ttoken 6: '1' if position is included by mask, '0' if not\n");
    }
    fprintf(tabfp, "#\n");
    fprintf(tabfp, "# Information content is calculated as 2.0 - H, where\n");
    fprintf(tabfp, "# H = - \\sum_x p_x \\log_2 p_x, for x in {A, C, G, U}\n");
    fprintf(tabfp, "# p_x is the frequency of x for *non-gap* residues at the position.\n");
    fprintf(tabfp, "# For example, p_A in a column that includes 4 As, 3 Cs, 2 Gs, 1 U and 5 gaps\n");
    fprintf(tabfp, "# would be 4/10 = 0.4.\n");
    fprintf(tabfp, "# Maximum possible value for token 3 is %d, the number of sequences in the file.\n", msa_nseq); 
    fprintf(tabfp, "#\n");
    fprintf(tabfp, "# Value ranges for bins:\n");
    fprintf(tabfp, "# \tbin  0: special case, 0 non-gap residues in this position\n");
    for(l = 0; l < hc_nbins; l++) { 
      fprintf(tabfp, "# \tbin %2d: [%.3f-%.3f%s information per position (bits)\n", l+1, limits[l], limits[l+1], (l == hc_nbins-1) ? "]" : ")");
    }
    fprintf(tabfp, "#\n");
    fprintf(tabfp, "# %11s  %6s  %8s  %10s  %3s", "type", "cpos", "info", "nongap", "bin");
    if(ps->mask != NULL) fprintf(tabfp, "  %4s", "mask");
    fprintf(tabfp, "\n");
    fprintf(tabfp, "# %11s  %6s  %8s  %10s  %3s", "-----------", "------", "--------", "----------", "---");
    if(ps->mask != NULL) fprintf(tabfp, "  %4s", "----");
    fprintf(tabfp, "\n");
  }
  
  if(ps->mask == NULL) nonecell_masked = -1; /* special flag */
  for(rfpos = 0; rfpos < ps->rflen; rfpos++) { 
    apos = ps->msa_rf2a_map[rfpos];
    esl_vec_DCopy(abc_ct[apos], abc->K, tmp_obs); /* only copy first abc->K values, don't copy gaps */
    zero_obs = (esl_DCompare(esl_vec_DSum(tmp_obs, abc->K), 0., eslSMALLX1) == eslOK) ? TRUE : FALSE;
    esl_vec_DNorm(tmp_obs, abc->K);
    ent[rfpos] = esl_vec_DEntropy(bg, abc->K) - esl_vec_DEntropy(tmp_obs, abc->K);

    if(zero_obs) { 
      if((status = set_onecell_values(errbuf, ps->rcolAAA[pp][rfpos], NCMYK, hc_onecell[hc_onecell_idx])) != eslOK) return status;
      nonecell++;
      if(ps->mask != NULL && ps->mask[rfpos] == '1') nonecell_masked++; 
      bi = -1;
    }
    else { 
      within_mask = (ps->mask != NULL && ps->mask[rfpos] == '1') ? TRUE : FALSE;
      if((status = set_scheme_values(errbuf, ps->rcolAAA[pp][rfpos], NCMYK, hc_scheme[hc_scheme_idx], ent[rfpos], ps->sclAA[pp], within_mask, &bi)) != eslOK) return status;
    }

    ps->rrAA[pp][rfpos] = ' ';
    /*ps->rrAA[pp][rfpos] = (esl_FCompare(ent[rfpos], 0., eslSMALLX1) == eslOK) ? '-' : ' ';*/

    if(tabfp != NULL) { 
      fprintf(tabfp, "  infocontent  %6d  %8.5f  %10d  %3d", rfpos+1, ent[rfpos], (int) esl_vec_DSum(abc_ct[apos], abc->K), bi+1);
      if(ps->mask != NULL) fprintf(tabfp, "  %4d\n", ps->mask[rfpos] == '1' ? 1 : 0);
      fprintf(tabfp, "\n");
    }
  }

  /* add one-cell color legend */
  ps->occlAAA[pp][0] = create_onecell_colorlegend(hc_onecell[hc_onecell_idx], nonecell, nonecell_masked);
  if((status = add_text_to_onecell_colorlegend(ps, ps->occlAAA[pp][0], "100% gaps", ps->legx_max_chars, errbuf)) != eslOK) return status;
  ps->nocclA[pp] = 1;

  /* add text to legend */
  /*sprintf(text, "information content (bits) (total: %.2f bits)", esl_vec_DSum(ent, ps->rflen));*/
  if((status = add_text_to_scheme_colorlegend(ps->sclAA[pp], "information content (bits)", ps->legx_max_chars, errbuf)) != eslOK) return status;

  /* add description to ps */
  if((status = add_page_desc_to_sspostscript(ps, pp, "information content per position", errbuf)) != eslOK) return status;

  free(ent);
  free(tmp_obs);
  free(bg);
  free(limits);

  if(tabfp != NULL) fprintf(tabfp, "//\n");
  return eslOK;
  
 ERROR: ESL_FAIL(status, errbuf, "infocontent_sspostscript(): memory allocation error.");
  return status; /* NEVERREACHED */
}

/* Function: delete_sspostscript()
 * 
 * Purpose:  Fill a postscript data structure with a new page w/colored squares indicating
 *           the number of sequences with gaps (deletions) at each consensus column. 
 *           If do_all is TRUE the page shows all deletions. If false only 'internal' deletions,
 *           those that come after the first occupied consensus column of each sequence and
 *           before the final occupied consensus column for each sequence, are shown.
 *           
 * Return:   eslOK on success.
 */
static int
delete_sspostscript(const ESL_GETOPTS *go, ESL_ALPHABET *abc, char *errbuf, SSPostscript_t *ps, double **abc_ct, int *srfpos_ct, int *erfpos_ct, int msa_nseq, int do_all, float ***hc_scheme, int hc_scheme_idx, int hc_nbins, float **hc_onecell, int hc_onecell_idx, FILE *tabfp)
{
  int status;
  int p, pp, c, l, bi;
  int rfpos, apos;
  int orig_npage = ps->npage;
  double dfreq;
  int nonecell = 0;
  int nonecell_masked = 0;
  int within_mask;
  float *limits = NULL;
  int *span_ct = NULL;   /* [0..rfpos..ps->msa->rflen-1] number of sequences that 'span' each position rfpos
			  * (have >= 1 residue at a consensus position <= rfpos, and >= 1 residue at a consensus position >= rfpos),
			  * msa_nseq - span_ct[rfpos] gives the number of 'external' deletions at position apos,
			  * these are deletes that come before the first consensus residue of a seq, or after the final consensus residue 
			  */
  float n_ext_del; /* msa_nseq - span_ct[rfpos], number of external deletions */

  if((status = add_pages_sspostscript(ps, 1, ALIMODE)) != eslOK) ESL_FAIL(status, errbuf, "memory error adding pages to the postscript object.");

  for(p = orig_npage; p < ps->npage; p++) { 
    ESL_ALLOC(ps->rrAA[p], sizeof(char) *  ps->rflen);
    ESL_ALLOC(ps->rcolAAA[p], sizeof(int *) * ps->rflen);
    ESL_ALLOC(ps->sclAA[p],    sizeof(SchemeColorLegend_t *) * 1);
    ESL_ALLOC(ps->occlAAA[p], sizeof(OneCellColorLegend_t **) * 1);
    for(c = 0; c < ps->rflen; c++) { 
      ESL_ALLOC(ps->rcolAAA[p][c], sizeof(int) * NCMYK); /* CMYK colors */
    }
  }

  pp = orig_npage;

  /* add color legend */
  ESL_ALLOC(limits, sizeof(float) * (hc_nbins+1)); 
  limits[0] = 0.0;
  limits[1] = 0.167;
  limits[2] = 0.333;
  limits[3] = 0.500;
  limits[4] = 0.667;
  limits[5] = 0.833;
  limits[6] = 1.000;
  ps->sclAA[pp] = create_scheme_colorlegend(hc_scheme_idx, hc_nbins, limits, FALSE);

  if(tabfp != NULL) { 
    if(do_all) { 
      fprintf(tabfp, "# -----------\n");
      fprintf(tabfp, "# Delete data\n");
      fprintf(tabfp, "# -----------\n");
      fprintf(tabfp, "# This section includes %d non #-prefixed lines, one for each consensus position\n", ps->rflen);
      fprintf(tabfp, "# in the alignment and corresponding template.\n");
      fprintf(tabfp, "# Each line includes %d tokens, separated by whitespace:\n", ps->mask == NULL ? 4 : 5);
      fprintf(tabfp, "# \ttoken 1: 'deleteall' (tag defining line type to ease parsing)\n");
      fprintf(tabfp, "# \ttoken 2: consensus position (starting at 1)\n");
      fprintf(tabfp, "# \ttoken 3: frequency of deletions (gaps) for position\n");
      fprintf(tabfp, "# \ttoken 4: bin index this positions falls in (see bin values below)\n");
      if(ps->mask != NULL) { 
	fprintf(tabfp, "# \ttoken 5: '1' if position is included by mask, '0' if not\n");
      }
      fprintf(tabfp, "#\n");
      fprintf(tabfp, "# A sequence s has a 'delete' at consensus position x if position\n");
      fprintf(tabfp, "# x is a gap for aligned sequence s.\n");
      fprintf(tabfp, "# Total number of sequences in the file is %d\n", msa_nseq);
      fprintf(tabfp, "#\n");
      fprintf(tabfp, "# Value ranges for bins:\n");
      fprintf(tabfp, "# \tbin  0: special case, 0 sequences have a delete at position\n");
      for(l = 0; l < hc_nbins; l++) { 
	fprintf(tabfp, "# \tbin %2d: [%.3f-%.3f%s frequency of deletes per position\n", l+1, limits[l], limits[l+1], (l == hc_nbins-1) ? "]" : ")");
      }
      fprintf(tabfp, "#\n");
      fprintf(tabfp, "# %9s  %6s  %8s  %3s", "type", "cpos", "dfreq", "bin");
      if(ps->mask != NULL) fprintf(tabfp, "  %4s", "mask");
      fprintf(tabfp, "\n");
      fprintf(tabfp, "# %9s  %6s  %8s  %3s", "---------", "------", "--------", "---");
      if(ps->mask != NULL) fprintf(tabfp, "  %4s", "----");
      fprintf(tabfp, "\n");
    }
    else { /* ! do_all (internal deletes only) */
      fprintf(tabfp, "# --------------------\n");
      fprintf(tabfp, "# Internal delete data\n");
      fprintf(tabfp, "# --------------------\n");
      fprintf(tabfp, "# This section includes %d non #-prefixed lines, one for each consensus position\n", ps->rflen);
      fprintf(tabfp, "# in the alignment and corresponding template.\n");
      fprintf(tabfp, "# Each line includes %d tokens, separated by whitespace:\n", ps->mask == NULL ? 5 : 6);
      fprintf(tabfp, "# \ttoken 1: 'deleteint' (tag defining line type to ease parsing)\n");
      fprintf(tabfp, "# \ttoken 2: consensus position (starting at 1)\n");
      fprintf(tabfp, "# \ttoken 3: frequency of internal deletions (gaps) for position\n");
      fprintf(tabfp, "# \ttoken 4: number of sequences that span (begin at or prior to and end at or after) position (max is %d)\n", msa_nseq);
      fprintf(tabfp, "# \ttoken 5: bin index this positions falls in (see bin values below)\n");
      if(ps->mask != NULL) fprintf(tabfp, "# \ttoken 6: '1' if position is included by mask, '0' if not\n");
      fprintf(tabfp, "#\n");
      fprintf(tabfp, "# A sequence s has an 'internal delete' at consensus position x if position\n");
      fprintf(tabfp, "# x is a gap for aligned sequence s, and s has at least one non-gap residue aligned\n");
      fprintf(tabfp, "# to a consensus position a <= x, and at least one non-gap residue aligned to a \n");
      fprintf(tabfp, "# consensus position b >= x.\n");
      fprintf(tabfp, "#\n");
      fprintf(tabfp, "# Value ranges for bins:\n");
      fprintf(tabfp, "# \tbin  0: special case, 0 sequences have an internal delete at position\n");
      for(l = 0; l < hc_nbins; l++) { 
	fprintf(tabfp, "# \tbin %2d: [%.3f-%.3f%s frequency of internal deletes per position\n", l+1, limits[l], limits[l+1], (l == hc_nbins-1) ? "]" : ")");
      }
      fprintf(tabfp, "#\n");
      fprintf(tabfp, "# %9s  %6s  %8s  %10s  %3s", "type", "cpos", "dfreq", "nspan", "bin");
      if(ps->mask != NULL) fprintf(tabfp, "  %4s", "mask");
      fprintf(tabfp, "\n");
      fprintf(tabfp, "# %9s  %6s  %8s  %10s  %3s", "---------", "------", "--------", "----------", "---");
      if(ps->mask != NULL) fprintf(tabfp, "  %4s", "----");
      fprintf(tabfp, "\n");
    }
  }

  /* determine number of sequences that span each position */
  if((status = get_span_ct(ps->rflen, msa_nseq, srfpos_ct, erfpos_ct, &span_ct)) != eslOK) ESL_FAIL(eslEMEM, errbuf, "Out of memory, getting span_ct array.");

  if(ps->mask == NULL) nonecell_masked = -1; /* special flag */
  /* draw delete page */
  for(rfpos = 0; rfpos < ps->rflen; rfpos++) { 
    ps->rrAA[pp][rfpos] = ' ';
    apos = ps->msa_rf2a_map[rfpos];
    n_ext_del = (float) (msa_nseq - span_ct[rfpos]); /* num external deletes is num seqs minus number of seqs that 'span' apos, see function header comments for explanation */
    if((( do_all) && ( abc_ct[apos][abc->K] < eslSMALLX1)) ||               /* abc_ct[apos][abc->K] == 0 */
       ((!do_all) && ((abc_ct[apos][abc->K] - n_ext_del) < eslSMALLX1))) {  /* abc_ct[apos][abc->K] == n_ext_del (all deletes are external) */
      if((status = set_onecell_values(errbuf, ps->rcolAAA[pp][rfpos], NCMYK, hc_onecell[hc_onecell_idx])) != eslOK) return status; 
      nonecell++;
      if(ps->mask != NULL && ps->mask[rfpos] == '1') nonecell_masked++; 
      bi = -1;
      dfreq = 0.;
    }
    else {
      within_mask = (ps->mask != NULL && ps->mask[rfpos] == '1') ? TRUE : FALSE;
      dfreq = do_all ? 
	( abc_ct[apos][abc->K]              / (float) msa_nseq) : 
	((abc_ct[apos][abc->K] - n_ext_del) / (float) msa_nseq);
      /* printf("do_all: %d rfpos: %d dfreq: %.4f\n", do_all, rfpos, dfreq); */
      if((status = set_scheme_values(errbuf, ps->rcolAAA[pp][rfpos], NCMYK, hc_scheme[hc_scheme_idx], dfreq, ps->sclAA[pp], within_mask, &bi)) != eslOK) return status;
    }
    if(tabfp != NULL) { 
      if(do_all) fprintf(tabfp, "  deleteall  %6d  %8.5f  %3d",       rfpos+1, dfreq, bi+1);
      else       fprintf(tabfp, "  deleteint  %6d  %8.5f  %10d  %3d", rfpos+1, dfreq, span_ct[rfpos], bi+1);
      if(ps->mask != NULL) fprintf(tabfp, "  %4d", ps->mask[rfpos] == '1' ? 1 : 0); 
      fprintf(tabfp, "\n");
    }
  }
    
  /* add one-cell color legend */
  ps->occlAAA[pp][0] = create_onecell_colorlegend(hc_onecell[hc_onecell_idx], nonecell, nonecell_masked);
  ps->nocclA[pp] = 1;

  if(do_all) { 
    if((status = add_text_to_onecell_colorlegend(ps, ps->occlAAA[pp][0], "zero deletions", ps->legx_max_chars, errbuf)) != eslOK) return status;
  }
  else {
    if((status = add_text_to_onecell_colorlegend(ps, ps->occlAAA[pp][0], "zero internal deletions", ps->legx_max_chars, errbuf)) != eslOK) return status; 
  }

  /* add color legend and description */
  if(do_all) { 
    /*sprintf(text, "fraction seqs w/deletes ('-'=0 deletes; avg/seq: %.2f)", (float) esl_vec_ISum(dct, ps->rflen) / (float) msa_nseq);*/
    if((status = add_text_to_scheme_colorlegend(ps->sclAA[pp], "fraction of seqs with deletes", ps->legx_max_chars, errbuf)) != eslOK) return status;
    if((status = add_page_desc_to_sspostscript(ps, ps->npage-1, "frequency of deletions at each position", errbuf)) != eslOK) return status;
  }
  else { /* !do_all, only internal deletes counted */
    /*sprintf(text, "fraction seqs w/internal deletes ('-'=0; avg/seq: %.2f)", (float) esl_vec_ISum(dct_internal, ps->rflen) / (float) msa_nseq);*/
    if((status = add_text_to_scheme_colorlegend(ps->sclAA[pp], "fraction of seqs w/internal deletions", ps->legx_max_chars, errbuf)) != eslOK) return status;
    if((status = add_page_desc_to_sspostscript(ps, ps->npage-1, "frequency of internal \\(non-terminal\\) deletions in each position", errbuf)) != eslOK) return status;
  }
  
  if(span_ct != NULL)   free(span_ct);
  free(limits);
  
  if(tabfp != NULL) fprintf(tabfp, "//\n");
  return eslOK;
  
 ERROR: ESL_FAIL(status, errbuf, "delete_sspostscript(): memory allocation error.");
  return status; /* NEVERREACHED */
}

/* Function: insert_sspostscript()
 * 
 * Purpose:  Fill a postscript data structure with 1 new page, with colors 
 *           indicating the fraction of seqs with inserts after each
 *           position. This function differs from the others in that the number
 *           of inserts per column is passed in. This is because insert information
 *           can either be extracted from the msa (in get_insert_info_from_msa()) or
 *           from an insert information file <f> (in get_insert_info_from_ifile(), 
 *           where <f> was created by cmalign --matchonly --iinfo <f>).
 *           
 * Return:   eslOK on success.
 */
static int
insert_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, int *nseq_with_ins_ct, int msa_nseq, float ***hc_scheme, int hc_scheme_idx, int hc_nbins, float **hc_onecell, int hc_onecell_idx, FILE *tabfp)
{
  int status;
  int p, pp, c, l;
  int rfpos;
  int orig_npage = ps->npage;
  int apos;
  int nonecell = 0;
  int nonecell_masked = 0;
  float *limits;
  int within_mask;
  float ifreq;
  int bi;

  if(ps->mask == NULL) nonecell_masked = -1; /* special flag */

  if((status = add_pages_sspostscript(ps, 1, ALIMODE)) != eslOK) ESL_FAIL(status, errbuf, "memory error adding pages to the postscript object.");

  for(p = orig_npage; p < ps->npage; p++) { 
    ESL_ALLOC(ps->rrAA[p], sizeof(char) *  (ps->rflen+1));
    ESL_ALLOC(ps->rcolAAA[p], sizeof(float *) * ps->rflen);
    ESL_ALLOC(ps->sclAA[p],   sizeof(SchemeColorLegend_t *) * 1);
    ESL_ALLOC(ps->occlAAA[p], sizeof(OneCellColorLegend_t **) * 1);
    for(c = 0; c < ps->rflen; c++) { 
      ESL_ALLOC(ps->rcolAAA[p][c], sizeof(float) * NCMYK); /* CMYK colors */
    }
  }

  pp = orig_npage;

  /* add color legend */
  ESL_ALLOC(limits, sizeof(float) * (hc_nbins+1)); 
  limits[0] = 0.0;
  limits[1] = 0.167;
  limits[2] = 0.333;
  limits[3] = 0.500;
  limits[4] = 0.667;
  limits[5] = 0.833;
  limits[6] = 1.00;
  ps->sclAA[pp] = create_scheme_colorlegend(hc_scheme_idx, hc_nbins, limits, FALSE);

  if(tabfp != NULL) { 
    fprintf(tabfp, "# -----------\n");
    fprintf(tabfp, "# Insert data\n");
    fprintf(tabfp, "# -----------\n");
    fprintf(tabfp, "# This section includes %d non #-prefixed lines, one for each possible insert position\n", ps->rflen+1);
    fprintf(tabfp, "# after each of the %d consensus positions and one more for inserts prior to the first consensus position.\n", ps->rflen);
    fprintf(tabfp, "# Each line includes %d tokens, separated by whitespace:\n", ps->mask == NULL ? 4 : 5);
    fprintf(tabfp, "# \ttoken 1: 'insert' (tag defining line type to ease parsing)\n");
    fprintf(tabfp, "# \ttoken 2: consensus position after which inserts occur ('0' == before posn 1)\n");
    fprintf(tabfp, "# \ttoken 3: fraction of sequences with an >= 1 inserted residues after position\n");
    fprintf(tabfp, "# \ttoken 4: bin index this positions falls in (see bin values below)\n");
    if(ps->mask != NULL) { 
      fprintf(tabfp, "# \ttoken 5: '1' if position is included by mask, '0' if not\n");
    }
    fprintf(tabfp, "#\n");
    fprintf(tabfp, "# Value ranges for bins:\n");
    fprintf(tabfp, "# \tbin -1: special case, reserved for inserts before position 1,\n");
    fprintf(tabfp, "# \t        these are NOT SHOWN in the postscript diagram (!)\n");
    fprintf(tabfp, "# \tbin  0: special case, 0 sequences have inserts after this position\n");
    for(l = 0; l < hc_nbins; l++) { 
      fprintf(tabfp, "# \tbin %2d: [%.3f-%.3f%s fraction of sequences with >= 1 inserts after each position\n", l+1, limits[l], limits[l+1], (l == hc_nbins-1) ? "]" : ")");
    }
    fprintf(tabfp, "#\n");
    fprintf(tabfp, "# %6s  %6s  %8s  %3s", "type", "cpos", "ifreq", "bin");
    if(ps->mask != NULL) fprintf(tabfp, "  %4s", "mask");
    fprintf(tabfp, "\n");
    fprintf(tabfp, "# %6s  %6s  %8s  %3s", "------", "------", "--------", "---");
    if(ps->mask != NULL) fprintf(tabfp, "  %4s", "----");
    fprintf(tabfp, "\n");
  }

  /* print info on inserts before rfpos 1 to tabfile, if nec */
  if(tabfp != NULL) { 
    apos = ps->msa_rf2a_map[0];
    ifreq = (float) nseq_with_ins_ct[0] / (float) msa_nseq;
    fprintf(tabfp, "  insert  %6d  %8.5f  %3d", 0, ifreq, -1);
    if(ps->mask != NULL) fprintf(tabfp, "  %4d", 0);
    fprintf(tabfp, "\n");
  }

  for(rfpos = 0; rfpos < ps->rflen; rfpos++) { 
    ps->rrAA[pp][rfpos] = ' ';
    apos = ps->msa_rf2a_map[rfpos]; 
    if(nseq_with_ins_ct[(rfpos+1)] == 0) {  /* careful, nseq_with_ins_ct goes from 1..rflen, its off-by-one with other arrays */
      if((status = set_onecell_values(errbuf, ps->rcolAAA[pp][rfpos], NCMYK, hc_onecell[hc_onecell_idx])) != eslOK) return status;
      nonecell++;
      if(ps->mask != NULL && ps->mask[rfpos] == '1') nonecell_masked++; 
      ifreq = 0.;
      bi = -1;
    }
    else {
      within_mask = (ps->mask != NULL && ps->mask[rfpos] == '1') ? TRUE : FALSE;
      ifreq = (float) nseq_with_ins_ct[rfpos+1] / (float) msa_nseq;
      if((status = set_scheme_values(errbuf, ps->rcolAAA[pp][rfpos], NCMYK, hc_scheme[hc_scheme_idx], ifreq, ps->sclAA[pp], within_mask, &bi)) != eslOK) return status;
    }
    /* printf("rfpos: %5d ifreq: %.3f\n", rfpos, ifreq); */
    if(tabfp != NULL) { 
      fprintf(tabfp, "  insert  %6d  %8.5f  %3d", rfpos, ifreq, bi+1);
      if(ps->mask != NULL) fprintf(tabfp, "  %4d", ps->mask[rfpos] == '1' ? 1 : 0);
      fprintf(tabfp, "\n");
    }
  }

  /* add one-cell color legend */
  ps->occlAAA[pp][0] = create_onecell_colorlegend(hc_onecell[hc_onecell_idx], nonecell, nonecell_masked);
  ps->nocclA[pp] = 1;
  if((status = add_text_to_onecell_colorlegend(ps, ps->occlAAA[pp][0], "zero insertions", ps->legx_max_chars, errbuf)) != eslOK) return status;

  /* add color legend */
  if((status = add_text_to_scheme_colorlegend(ps->sclAA[pp], "fraction of seqs w/insertions", ps->legx_max_chars, errbuf)) != eslOK) return status;
  if((status = add_page_desc_to_sspostscript(ps, ps->npage-1, "frequency of insertions after each position", errbuf)) != eslOK) return status;

  free(limits);

  if(tabfp != NULL) fprintf(tabfp, "//\n");
  return eslOK;
  
 ERROR: ESL_FAIL(status, errbuf, "insert_sspostscript(): memory allocation error.");
  return status; /* NEVERREACHED */
}


/* Function: span_sspostscript()
 * 
 * Purpose:  Fill a postscript data structure with 1 new page, with colors 
 *           indicating the fraction of seqs that 'span' each consensus 
 *           position. A consensus position cpos is spanned by a sequence x if
 *           if at least 1 residue in x is aligned to a consensus position <= cpos
 *           *and*  at least 1 residue in x is aligned to a consensus position >= cpos.
 *           
 * Return:   eslOK on success.
 */
static int
span_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, int *srfpos_ct, int *erfpos_ct, int msa_nseq, float ***hc_scheme, int hc_scheme_idx, int hc_nbins, float **hc_onecell, int zerocov_idx, int maxcov_idx, FILE *tabfp)
{
  int status;
  int p, pp, c, l;
  int rfpos;
  int orig_npage = ps->npage;
  int nzerocov = 0;
  int nzerocov_masked = 0;
  int nmaxcov = 0;
  int nmaxcov_masked = 0;
  float *limits;
  int within_mask;
  float cfract;
  int bi;
  int *span_ct; /* [0..rfpos..rflen-1], number of sequences that 'span' position rfpos */

  if(ps->mask == NULL) nzerocov_masked = nmaxcov_masked = -1; /* special flag */

  if((status = add_pages_sspostscript(ps, 1, ALIMODE)) != eslOK) ESL_FAIL(status, errbuf, "memory error adding pages to the postscript object.");

  for(p = orig_npage; p < ps->npage; p++) { 
    ESL_ALLOC(ps->rrAA[p], sizeof(char) *  (ps->rflen+1));
    ESL_ALLOC(ps->rcolAAA[p], sizeof(float *) * ps->rflen);
    ESL_ALLOC(ps->sclAA[p],   sizeof(SchemeColorLegend_t *) * 1);
    ESL_ALLOC(ps->occlAAA[p], sizeof(OneCellColorLegend_t **) * 2);
    for(c = 0; c < ps->rflen; c++) { 
      ESL_ALLOC(ps->rcolAAA[p][c], sizeof(float) * NCMYK); /* CMYK colors */
    }
  }

  /* determine number of sequences that span each position */
  if((status = get_span_ct(ps->rflen, msa_nseq, srfpos_ct, erfpos_ct, &span_ct)) != eslOK) ESL_FAIL(eslEMEM, errbuf, "Out of memory, getting span_ct array.");

  pp = orig_npage;

  /* add color legend */
  ESL_ALLOC(limits, sizeof(float) * (hc_nbins+1)); 
  limits[0] = 0.0;
  limits[1] = 0.167;
  limits[2] = 0.333;
  limits[3] = 0.500;
  limits[4] = 0.667;
  limits[5] = 0.833;
  limits[6] = 1.00;
  ps->sclAA[pp] = create_scheme_colorlegend(hc_scheme_idx, hc_nbins, limits, FALSE);

  if(tabfp != NULL) { 
    fprintf(tabfp, "# ---------\n");
    fprintf(tabfp, "# Span data\n");
    fprintf(tabfp, "# ---------\n");
    fprintf(tabfp, "# This section includes %d non #-prefixed lines, one for each consensus position\n", ps->rflen);
    fprintf(tabfp, "# in the alignment and corresponding template.\n");
    fprintf(tabfp, "# Each line includes %d tokens, separated by whitespace:\n", ps->mask == NULL ? 4 : 5);
    fprintf(tabfp, "# \ttoken 1: 'span' (tag defining line type to ease parsing)\n");
    fprintf(tabfp, "# \ttoken 2: consensus position (starting at 1)\n");
    fprintf(tabfp, "# \ttoken 3: fraction of sequences that 'span' position\n");
    fprintf(tabfp, "# \ttoken 4: bin index this positions falls in (see bin values below)\n");
    if(ps->mask != NULL) { 
      fprintf(tabfp, "# \ttoken 5: '1' if position is included by mask, '0' if not\n");
    }
    fprintf(tabfp, "#\n");
    fprintf(tabfp, "# A sequence s spans consensus position x if s has at least one\n");
    fprintf(tabfp, "# non-gap residue aligned to a consensus position a <= x and at least one\n");
    fprintf(tabfp, "# non-gap residue aligned to a consensus position b >= x..\n");
    fprintf(tabfp, "#\n");
    fprintf(tabfp, "# Value ranges for bins:\n");
    fprintf(tabfp, "# \tbin  0: special case, 0 sequences span this position\n");
    for(l = 0; l < hc_nbins; l++) { 
      fprintf(tabfp, "# \tbin %2d: [%.3f-%.3f%s fraction of sequences that span each position\n", l+1, limits[l], limits[l+1], (l == hc_nbins-1) ? "]" : ")");
    }
    fprintf(tabfp, "#\n");
    fprintf(tabfp, "# %8s  %6s  %8s  %3s", "type", "cpos", "span", "bin");
    if(ps->mask != NULL) fprintf(tabfp, "  %4s", "mask");
    fprintf(tabfp, "\n");
    fprintf(tabfp, "# %8s  %6s  %8s  %3s", "------", "------", "--------", "---");
    if(ps->mask != NULL) fprintf(tabfp, "  %4s", "----");
    fprintf(tabfp, "\n");
  }

  for(rfpos = 0; rfpos < ps->rflen; rfpos++) { 
    ps->rrAA[pp][rfpos] = ' ';
    if(span_ct[rfpos] == 0) {  
      if((status = set_onecell_values(errbuf, ps->rcolAAA[pp][rfpos], NCMYK, hc_onecell[zerocov_idx])) != eslOK) return status;
      nzerocov++;
      if(ps->mask != NULL && ps->mask[rfpos] == '1') nzerocov_masked++;
      cfract = 0.;
      bi = -1;
    }
    else if(span_ct[rfpos] == msa_nseq) {  
      if((status = set_onecell_values(errbuf, ps->rcolAAA[pp][rfpos], NCMYK, hc_onecell[maxcov_idx])) != eslOK) return status;
      nmaxcov++;
      if(ps->mask != NULL && ps->mask[rfpos] == '1') nmaxcov_masked++; 
      cfract = 0.;
      bi = -1;
    }
    else {
      within_mask = (ps->mask != NULL && ps->mask[rfpos] == '1') ? TRUE : FALSE;
      cfract = (float) span_ct[rfpos] / (float) msa_nseq;
      if((status = set_scheme_values(errbuf, ps->rcolAAA[pp][rfpos], NCMYK, hc_scheme[hc_scheme_idx], cfract, ps->sclAA[pp], within_mask, &bi)) != eslOK) return status;
    }
    if(tabfp != NULL) { 
      fprintf(tabfp, "  span  %6d  %8.5f  %3d", rfpos, cfract, bi+1); 
      if(ps->mask != NULL) fprintf(tabfp, "  %4d", ps->mask[rfpos] == '1' ? 1 : 0);
      fprintf(tabfp, "\n");
    }
  }

  /* add one-cell color legend for zero span positions */
  ps->occlAAA[pp][0] = create_onecell_colorlegend(hc_onecell[zerocov_idx], nzerocov, nzerocov_masked);
  if((status = add_text_to_onecell_colorlegend(ps, ps->occlAAA[pp][0], "no sequences span", ps->legx_max_chars, errbuf)) != eslOK) return status;

  /* add one-cell color legend for maximum span positions */
  ps->occlAAA[pp][1] = create_onecell_colorlegend(hc_onecell[maxcov_idx], nmaxcov, nmaxcov_masked);
  if((status = add_text_to_onecell_colorlegend(ps, ps->occlAAA[pp][1], "100% of seqs span", ps->legx_max_chars, errbuf)) != eslOK) return status;
  ps->nocclA[pp] = 2;

  /* add color legend */
  if((status = add_text_to_scheme_colorlegend(ps->sclAA[pp], "fraction of seqs that span each position", ps->legx_max_chars, errbuf)) != eslOK) return status;
  if((status = add_page_desc_to_sspostscript(ps, ps->npage-1, "fraction of sequences that span each position", errbuf)) != eslOK) return status;

  free(span_ct);
  free(limits);

  if(tabfp != NULL) fprintf(tabfp, "//\n");
  return eslOK;
  
 ERROR: ESL_FAIL(status, errbuf, "insert_sspostscript(): memory allocation error.");
  return status; /* NEVERREACHED */
}


/* Function: avg_posteriors_sspostscript()
 * 
 * Purpose:  Fill a postscript data structure with info average posterior probabilities in the MSA.
 * Return:   eslOK on success.
 */
static int
avg_posteriors_sspostscript(const ESL_GETOPTS *go, ESL_ALPHABET *abc, char *errbuf, SSPostscript_t *ps, int **pp_ct, int msa_nseq, float ***hc_scheme, int hc_scheme_idx, int hc_nbins, float **hc_onecell, int hc_onecell_idx, FILE *tabfp)
{
  int status;
  int p;
  int rfpos;
  int pp;
  float *limits; /* bin limits for the color scheme */
  int nonecell_allgap = 0;
  int nonecell_allgap_masked = 0;
  int within_mask;
  int orig_npage = ps->npage;
  float ppavgA[11];
  int l;
  int apos;
  int nnongap;
  float ppsum;
  int ppidx;
  float ppavg;
  int bi;

  ppavgA[0]  = 0.025;
  ppavgA[1]  = 0.10;
  ppavgA[2]  = 0.20;
  ppavgA[3]  = 0.30;
  ppavgA[4]  = 0.40;
  ppavgA[5]  = 0.50;
  ppavgA[6]  = 0.60;
  ppavgA[7]  = 0.70;
  ppavgA[8]  = 0.80;
  ppavgA[9]  = 0.90;
  ppavgA[10] = 0.975;

  if(ps->mask == NULL) { nonecell_allgap_masked = -1; } /* special flag */

  if((status = add_pages_sspostscript(ps, 1, ALIMODE)) != eslOK) ESL_FAIL(status, errbuf, "memory error adding pages to the postscript object.");

  for(p = orig_npage; p < ps->npage; p++) { 
    ESL_ALLOC(ps->rrAA[p], sizeof(char) *  (ps->rflen+1));
    ESL_ALLOC(ps->rcolAAA[p], sizeof(float *) * ps->rflen);
    ESL_ALLOC(ps->sclAA[p],   sizeof(SchemeColorLegend_t *) * 1);
    ESL_ALLOC(ps->occlAAA[p], sizeof(OneCellColorLegend_t **) * 1);
    for(rfpos = 0; rfpos < ps->rflen; rfpos++) { 
      ESL_ALLOC(ps->rcolAAA[p][rfpos], sizeof(float) * NCMYK); /* CMYK colors */
    }
  }

  ESL_ALLOC(limits, sizeof(float) * (hc_nbins+1)); 
  limits[0] = 0.0;
  limits[1] = 0.70;
  limits[2] = 0.80;
  limits[3] = 0.85;
  limits[4] = 0.90;
  limits[5] = 0.95;
  limits[6] = 1.00;

  if(tabfp != NULL) { 
    fprintf(tabfp, "# ----------------------------------\n");
    fprintf(tabfp, "# Average posterior probability data\n");
    fprintf(tabfp, "# ----------------------------------\n");
    fprintf(tabfp, "# This section includes %d non #-prefixed lines, one for each consensus position\n", ps->rflen);
    fprintf(tabfp, "# in the alignment and corresponding template.\n");
    fprintf(tabfp, "# Each line includes %d tokens, separated by whitespace:\n", ps->mask == NULL ? 5 : 6);
    fprintf(tabfp, "# \ttoken 1: 'avgpostprob' (tag defining line type to ease parsing)\n");
    fprintf(tabfp, "# \ttoken 2: consensus position (starting at 1)\n");
    fprintf(tabfp, "# \ttoken 3: average posterior probability of non-gap residues for position\n");
    fprintf(tabfp, "# \ttoken 4: number of non-gap residues in position (max possible is %d (num seqs in aln))\n", msa_nseq); 
    fprintf(tabfp, "# \ttoken 5: bin index this positions falls in (see bin values below)\n");
    if(ps->mask != NULL) fprintf(tabfp, "# \ttoken 6: '1' if position is included by mask, '0' if not\n");
    fprintf(tabfp, "#\n");
    fprintf(tabfp, "# Posterior probability (PP) values in the alignment file can have 12 possible values,\n");
    fprintf(tabfp, "# the average per position is calculated by defining each as the average of its range given\n");
    fprintf(tabfp, "# below. (For example, a '8' which indicates PP between 0.75 and 0.85 is treated as 0.8).\n");
    fprintf(tabfp, "# \t'.': gap, corresponds to a gap in the sequence (not counted)\n");
    fprintf(tabfp, "# \t'0': posterior probability of between 0.00 and 0.05\n");
    fprintf(tabfp, "# \t'1': posterior probability of between 0.05 and 0.15\n");
    fprintf(tabfp, "# \t'2': posterior probability of between 0.15 and 0.25\n");
    fprintf(tabfp, "# \t'3': posterior probability of between 0.25 and 0.35\n");
    fprintf(tabfp, "# \t'4': posterior probability of between 0.35 and 0.45\n");
    fprintf(tabfp, "# \t'5': posterior probability of between 0.45 and 0.55\n");
    fprintf(tabfp, "# \t'6': posterior probability of between 0.55 and 0.65\n");
    fprintf(tabfp, "# \t'7': posterior probability of between 0.65 and 0.75\n");
    fprintf(tabfp, "# \t'8': posterior probability of between 0.75 and 0.85\n");
    fprintf(tabfp, "# \t'9': posterior probability of between 0.85 and 0.95\n");
    fprintf(tabfp, "# \t'*': posterior probability of between 0.95 and 1.00\n");
    fprintf(tabfp, "#\n");
    fprintf(tabfp, "# Value ranges for bins:\n");
    fprintf(tabfp, "# \tbin  0: special case, 0 sequences have a non-gap residue at position\n");
    for(l = 0; l < hc_nbins; l++) { 
      fprintf(tabfp, "# \tbin %2d: [%.3f-%.3f%s average posterior probability per position\n", l+1, limits[l], limits[l+1], (l == hc_nbins-1) ? "]" : ")");
    }
    fprintf(tabfp, "#\n");
    fprintf(tabfp, "# %11s  %6s  %8s  %10s  %3s", "type", "cpos", "avgpp", "nongap", "bin");
    if(ps->mask != NULL) fprintf(tabfp, "  %4s", "mask");
    fprintf(tabfp, "\n");
    fprintf(tabfp, "# %11s  %6s  %8s  %10s  %3s", "-----------", "------", "--------", "----------", "---");
    if(ps->mask != NULL) fprintf(tabfp, "  %4s", "----");
    fprintf(tabfp, "\n");
  }

  /* step through each sequence and each column, collecting stats */
  pp = orig_npage;
  ps->sclAA[pp] = create_scheme_colorlegend(hc_scheme_idx, hc_nbins, limits, FALSE);
  for(rfpos = 0; rfpos < ps->rflen; rfpos++) { 
    apos = ps->msa_rf2a_map[rfpos]; 
    nnongap = esl_vec_ISum(pp_ct[apos], 11);
    if(nnongap == 0) { /* all gaps */
      if((status = set_onecell_values(errbuf, ps->rcolAAA[pp][rfpos], NCMYK, hc_onecell[hc_onecell_idx])) != eslOK) return status;
      nonecell_allgap++;
      if(ps->mask != NULL && ps->mask[rfpos] == '1') nonecell_allgap_masked++; 
      bi = -1;
      ppavg = 0.;
    }
    else { /* at least 1 non-gap residue in this position */
      if(pp_ct[apos][10] == nnongap) { /* all nongap residues have highest possible posterior probability */
      }
      else { 
      }
      ppsum = 0.;
      for(ppidx = 0; ppidx < 11; ppidx++) {
	ppsum += pp_ct[apos][ppidx] * ppavgA[ppidx]; /* Note: PP value is considered average of range, not minimum ('9' == 0.90 (0.95-0.85/2) */
      }
      ppavg = (float) ppsum / (float) nnongap;
      within_mask = (ps->mask != NULL && ps->mask[rfpos] == '1') ? TRUE : FALSE;
      if((status = set_scheme_values(errbuf, ps->rcolAAA[pp][rfpos], NCMYK, hc_scheme[hc_scheme_idx], ppavg, ps->sclAA[pp], within_mask, &bi)) != eslOK) return status;
    }
    if(tabfp != NULL) { 
      fprintf(tabfp, "  avgpostprob  %6d  %8.5f  %10d  %3d", rfpos+1, ppavg, nnongap, bi+1);
      if(ps->mask != NULL) fprintf(tabfp, "  %4d", ps->mask[rfpos] == '1' ? 1 : 0);
      fprintf(tabfp, "\n");
    }

    ps->rrAA[pp][rfpos] = ' ';
  }

  /* add one-cell color legend for all gap positions */
  ps->occlAAA[pp][0] = create_onecell_colorlegend(hc_onecell[hc_onecell_idx], nonecell_allgap, nonecell_allgap_masked);
  if((status = add_text_to_onecell_colorlegend(ps, ps->occlAAA[pp][0], "100% gaps", ps->legx_max_chars, errbuf)) != eslOK) return status;

  ps->nocclA[pp] = 1;
  
  /* add color legend */
  if((status = add_text_to_scheme_colorlegend(ps->sclAA[pp], "average posterior probability \\(alnment confidence\\)", ps->legx_max_chars,errbuf)) != eslOK) return status;

  /* add description to ps */
  if((status = add_page_desc_to_sspostscript(ps, pp, "average posterior probability \\(confidence\\) per position", errbuf)) != eslOK) return status;

  free(limits);

  if(tabfp != NULL) fprintf(tabfp, "//\n");
  return eslOK;

 ERROR:
  return status;
}


/* Function: individual_posteriors_sspostscript()
 * 
 * Purpose:  Fill a postscript data structure with posterior probabilities from single sequences in an MSA.
 * Return:   eslOK on success.
 */
static int
individual_posteriors_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, int *useme, int nused, float ***hc_scheme, int hc_scheme_idx, int hc_nbins, float **hc_onecell, int hc_onecell_idx)
{
  int    status;
  int    i,rfpos, apos;       /* counters over sequences, alignment posns, RF positions of MSA */
  int    p;
  int orig_npage = ps->npage;
  int new_npage = 0;
  int pp;
  float *limits; /* bin limits for the color scheme */
  int nonecell_seq = 0;
  int nonecell_seq_masked = 0;
  int within_mask;
  float ppavgA[11];
  int ppidx;

  ppavgA[0]  = 0.025;
  ppavgA[1]  = 0.10;
  ppavgA[2]  = 0.20;
  ppavgA[3]  = 0.30;
  ppavgA[4]  = 0.40;
  ppavgA[5]  = 0.50;
  ppavgA[6]  = 0.60;
  ppavgA[7]  = 0.70;
  ppavgA[8]  = 0.80;
  ppavgA[9]  = 0.90;
  ppavgA[10] = 0.975;

  if(ps->mask == NULL) { nonecell_seq_masked = -1; } /* special flag */
  if(msa->rf == NULL) esl_fatal("No RF annotation in alignment");

  new_npage = nused; 

  if((status = add_pages_sspostscript(ps, new_npage, INDIMODE)) != eslOK) ESL_FAIL(status, errbuf, "memory error adding pages to the postscript object.");

  for(p = orig_npage; p < ps->npage; p++) { 
    ESL_ALLOC(ps->rrAA[p], sizeof(char) *  (ps->rflen+1));
    ESL_ALLOC(ps->rcolAAA[p], sizeof(float *) * ps->rflen);
    ESL_ALLOC(ps->sclAA[p],   sizeof(SchemeColorLegend_t *) * 1);
    ESL_ALLOC(ps->occlAAA[p], sizeof(OneCellColorLegend_t **) * 1);
    for(rfpos = 0; rfpos < ps->rflen; rfpos++) { 
      ESL_ALLOC(ps->rcolAAA[p][rfpos], sizeof(float) * NCMYK); /* CMYK colors */
    }
  }

  ESL_ALLOC(limits, sizeof(float) * (hc_nbins+1)); 
  limits[0] = 0.0;
  limits[1] = 0.35;
  limits[2] = 0.55;
  limits[3] = 0.75;
  limits[4] = 0.85;
  limits[5] = 0.95;
  limits[6] = 1.00;

  /* step through each sequence and each column, collecting stats */
  pp = orig_npage;
  for(i = 0; i < msa->nseq; i++) { 
    if(useme[i]) { /* add color legend for this sequence */
      ps->sclAA[pp] = create_scheme_colorlegend(hc_scheme_idx, hc_nbins, limits, FALSE);
      nonecell_seq = 0;
      nonecell_seq_masked = (ps->mask == NULL) ? -1 : 0;

      for(rfpos = 0; rfpos < ps->rflen; rfpos++) { 
	apos = ps->msa_rf2a_map[rfpos];
	if(! esl_abc_CIsGap(msa->abc, msa->aseq[i][apos])) {
	  if((ppidx = get_pp_idx(msa->abc, msa->pp[i][apos])) == -1) ESL_FAIL(eslEFORMAT, errbuf, "bad #=GR PP char: %c", msa->pp[i][apos]);
	  if(ppidx == 11) ESL_FAIL(eslEFORMAT, errbuf, "nongap residue: %c, annotated with gap #=GR PP char: %c", msa->aseq[i][apos], msa->pp[i][apos]);
	  within_mask = (ps->mask != NULL && ps->mask[rfpos] == '1') ? TRUE : FALSE;
	  if((status = set_scheme_values(errbuf, ps->rcolAAA[pp][rfpos], NCMYK, hc_scheme[hc_scheme_idx], ppavgA[ppidx], ps->sclAA[pp], within_mask, NULL)) != eslOK) return status;
	  ps->rrAA[pp][rfpos] = ' ';
	}
	else {
	  if((status = set_onecell_values(errbuf, ps->rcolAAA[pp][rfpos], NCMYK, hc_onecell[hc_onecell_idx])) != eslOK) return status;
	  nonecell_seq++;
	  if(ps->mask != NULL && ps->mask[rfpos] == '1') nonecell_seq_masked++; 
	  ps->rrAA[pp][rfpos] = ' ';
	}
      }

      /* add one-cell color legend */
      ps->occlAAA[pp][0] = create_onecell_colorlegend(hc_onecell[hc_onecell_idx], nonecell_seq, nonecell_seq_masked);
      if((status = add_text_to_onecell_colorlegend(ps, ps->occlAAA[pp][0], "gap", ps->legx_max_chars, errbuf)) != eslOK) return status;
      ps->nocclA[pp] = 1;

      if((status = add_text_to_scheme_colorlegend(ps->sclAA[pp], "posterior probability \\(alignment confidence\\)", ps->legx_max_chars, errbuf)) != eslOK) return status;
      ps->seqidxA[pp] = i;
      if((status = add_page_desc_to_sspostscript(ps, pp, msa->sqname[i], errbuf)) != eslOK) return status;
      pp++;
    }
  } /* done with all sequences */


  free(limits);

  return eslOK;

 ERROR:
  return status;
}

/* Function: colormask_sspostscript()
 * 
 * Purpose:  Fill a postscript data structure with 1 new page based on a lanemask, each column
 *           is either black (if included, a '1' in the mask) or pink (not included, a '0' in the
 *           mask.
 *
 * Return:   eslOK on success.
 */
static int
colormask_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, char *mask, float **hc_onecell, int incmask_idx, int excmask_idx)
{
  int status;
  int p, pp, c;
  int cpos;
  int orig_npage = ps->npage;
  int ncols_inside_mask = 0;
  int ncols_outside_mask = 0;
  char *mask_desc = NULL;
  char *mask_file = NULL;

  if((status = add_pages_sspostscript(ps, 1, SIMPLEMASKMODE)) != eslOK) ESL_FAIL(status, errbuf, "memory error adding pages to the postscript object.");

  for(p = orig_npage; p < ps->npage; p++) { 
    ESL_ALLOC(ps->rrAA[p], sizeof(char) *  ps->rflen);
    ESL_ALLOC(ps->rcolAAA[p], sizeof(float *) * ps->rflen);
    ESL_ALLOC(ps->occlAAA[p], sizeof(OneCellColorLegend_t **) * 2);
    for(c = 0; c < ps->rflen; c++) { 
      ESL_ALLOC(ps->rcolAAA[p][c], sizeof(float) * NCMYK); /* CMYK colors */
    }
  }
  pp = orig_npage;

  for(cpos = 0; cpos < ps->rflen; cpos++) { 
    if(mask[cpos] == '1') { /* included */
      if((status = set_onecell_values(errbuf, ps->rcolAAA[pp][cpos], NCMYK, hc_onecell[incmask_idx])) != eslOK) return status; 
      ncols_inside_mask++;
    }
    else if(mask[cpos] == '0') {
      if((status = set_onecell_values(errbuf, ps->rcolAAA[pp][cpos], NCMYK, hc_onecell[excmask_idx])) != eslOK) return status; 
      ncols_outside_mask++;
    }
    else ESL_FAIL(eslEINVAL, errbuf, "--mask mask char number %d is not a 1 nor a 0, but a %c\n", cpos, mask[cpos]);
    ps->rrAA[pp][cpos] = ' ';
  }

  /* add color legend */
  ps->occlAAA[pp][0] = create_onecell_colorlegend(hc_onecell[incmask_idx], ncols_inside_mask, -1);
  if((status = add_text_to_onecell_colorlegend(ps, ps->occlAAA[pp][0], "columns included by mask", ps->legx_max_chars, errbuf)) != eslOK) return status;

  ps->occlAAA[pp][1] = create_onecell_colorlegend(hc_onecell[excmask_idx], ncols_outside_mask, -1);
  if((status = add_text_to_onecell_colorlegend(ps, ps->occlAAA[pp][1], "columns excluded by mask", ps->legx_max_chars, errbuf)) != eslOK) return status;

  if((status = esl_strcat(&(mask_desc), -1, "mask file: ", -1)) != eslOK) ESL_FAIL(status, errbuf, "error copying mask file name string");;
  if((status = esl_FileTail(esl_opt_GetString(go, "--mask"), FALSE, &mask_file)) != eslOK) ESL_FAIL(status, errbuf, "error copying mask file name string (probably out of memory)."); 
  if((strlen(mask_file) + strlen(mask_desc)) > (ps->desc_max_chars*2 - 2)) { /* desc would be too long, shorten mask_file so desc is legal */
    /* the -5 below is so we can add '...' to end */
    if((status = esl_strcat(&(mask_desc), -1, mask_file, ((ps->desc_max_chars*2) - strlen(mask_desc) - 5))) != eslOK) ESL_FAIL(status, errbuf, "error copying mask file name string");
    if((status = esl_strcat(&(mask_desc), -1, "...", 3)) != eslOK) ESL_FAIL(status, errbuf, "error copying mask file name string");
  }
  else { /* desc will not be too long */
    if((status = esl_strcat(&(mask_desc), -1, mask_file, -1)) != eslOK) ESL_FAIL(status, errbuf, "error copying mask file name string");
  }
  free(mask_file);
  if((status = add_page_desc_to_sspostscript(ps, pp, mask_desc, errbuf)) != eslOK) return status;

  ps->nocclA[pp] = 2;

  return eslOK;
  
 ERROR: ESL_FAIL(status, errbuf, "colormask_sspostscript(): memory allocation error.");
  return status; /* NEVERREACHED */
}



/* Function: diffmask_sspostscript()
 * 
 * Purpose:  Fill a postscript data structure with 1 new page based on a comparison between
 *           two masks, each column is either black (if included (a '1') in both masks), red
 *           (if a '1' in mask 1 and a '0' in mask 2), cyan (a '0' in mask 1 and a '1' in mask 2)
 *           or grey (a '0' in both masks).
 *
 * Return:   eslOK on success.
 */
static int
diffmask_sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, ESL_MSA *msa, char *mask1, char *mask2, float **hc_onecell, int incboth_idx, int inc1_idx, int inc2_idx, int excboth_idx)
{
  int status;
  int p, pp, c;
  int cpos;
  int orig_npage = ps->npage;
  int ncols_in_both = 0;
  int ncols_out_both = 0;
  int ncols_in_1_out_2 = 0;
  int ncols_out_1_in_2 = 0;

  if((status = add_pages_sspostscript(ps, 1, SIMPLEMASKMODE)) != eslOK) ESL_FAIL(status, errbuf, "memory error adding pages to the postscript object.");

  for(p = orig_npage; p < ps->npage; p++) { 
    ESL_ALLOC(ps->rrAA[p], sizeof(char) *  ps->rflen);
    ESL_ALLOC(ps->rcolAAA[p], sizeof(float *) * ps->rflen);
    ESL_ALLOC(ps->occlAAA[p],  sizeof(OneCellColorLegend_t **) * 4);
    for(c = 0; c < ps->rflen; c++) { 
      ESL_ALLOC(ps->rcolAAA[p][c], sizeof(float) * NCMYK); /* CMYK colors */
    }
  }
  pp = orig_npage;

  for(cpos = 0; cpos < ps->rflen; cpos++) { 
    if(mask1[cpos] == '1' && mask2[cpos] == '1') { 
      if((status = set_onecell_values(errbuf, ps->rcolAAA[pp][cpos], NCMYK, hc_onecell[incboth_idx])) != eslOK) return status; 
      ncols_in_both++;
    }
    else if(mask1[cpos] == '1' && mask2[cpos] == '0') {
      if((status = set_onecell_values(errbuf, ps->rcolAAA[pp][cpos], NCMYK, hc_onecell[inc1_idx])) != eslOK) return status; 
      ncols_in_1_out_2++;
    }
    else if(mask1[cpos] == '0' && mask2[cpos] == '1') {
      if((status = set_onecell_values(errbuf, ps->rcolAAA[pp][cpos], NCMYK, hc_onecell[inc2_idx])) != eslOK) return status; 
      ncols_out_1_in_2++;
    }
    else if(mask1[cpos] == '0' && mask2[cpos] == '0') {
      if((status = set_onecell_values(errbuf, ps->rcolAAA[pp][cpos], NCMYK, hc_onecell[excboth_idx])) != eslOK) return status; 
      ncols_out_both++;
    }
    else if(mask1[cpos] != '0' && mask1[cpos] != '1') ESL_FAIL(eslEINVAL, errbuf, "--mask-col char number %d is not a 1 nor a 0, but a %c\n", cpos, mask1[cpos]);
    else if(mask2[cpos] != '0' && mask2[cpos] != '1') ESL_FAIL(eslEINVAL, errbuf, "--mask-diff char number %d is not a 1 nor a 0, but a %c\n", cpos, mask2[cpos]);
    ps->rrAA[pp][cpos] = ' ';
  }

  /* add color legend */
  ps->occlAAA[pp][0] = create_onecell_colorlegend(hc_onecell[incboth_idx], ncols_in_both, -1);
  if((status = add_text_to_onecell_colorlegend(ps, ps->occlAAA[pp][0], "included by both masks", ps->legx_max_chars, errbuf)) != eslOK) return status;

  ps->occlAAA[pp][1] = create_onecell_colorlegend(hc_onecell[inc1_idx], ncols_in_1_out_2, -1);
  if((status = add_text_to_onecell_colorlegend(ps, ps->occlAAA[pp][1], "incl. mask 1, excl. mask 2", ps->legx_max_chars, errbuf)) != eslOK) return status;

  ps->occlAAA[pp][2] = create_onecell_colorlegend(hc_onecell[inc2_idx], ncols_out_1_in_2, -1);
  if((status = add_text_to_onecell_colorlegend(ps, ps->occlAAA[pp][2], "excl. mask 1, incl. mask 1", ps->legx_max_chars, errbuf)) != eslOK) return status;

  ps->occlAAA[pp][3] = create_onecell_colorlegend(hc_onecell[excboth_idx], ncols_out_both, -1);
  if((status = add_text_to_onecell_colorlegend(ps, ps->occlAAA[pp][3], "excluded by both masks", ps->legx_max_chars, errbuf)) != eslOK) return status;
  ps->nocclA[pp] = 4;

  if((status = add_diffmask_page_desc_to_sspostscript(ps, pp, esl_opt_GetString(go, "--mask"), esl_opt_GetString(go, "--mask-diff"), errbuf)) != eslOK) return status;

  return eslOK;
  
 ERROR: ESL_FAIL(status, errbuf, "diffmask_sspostscript(): memory allocation error.");
  return status; /* NEVERREACHED */
}
/* Function: add_pages_sspostscript()
 * 
 * Purpose:  Add and initialize blank pages to a postscript object.
 */ 
int 
add_pages_sspostscript(SSPostscript_t *ps, int ntoadd, int page_mode)
{
  int status;
  void *tmp;
  int p;

  if(ps->npage == 0) { 
    assert(ps->rrAA    == NULL);
    assert(ps->rcolAAA == NULL);
    ESL_ALLOC(ps->rrAA,    sizeof(char *)   * (ps->npage + ntoadd));
    ESL_ALLOC(ps->rcolAAA, sizeof(float **) * (ps->npage + ntoadd));
    ESL_ALLOC(ps->occlAAA, sizeof(OneCellColorLegend_t ***) * (ps->npage + ntoadd));
    ESL_ALLOC(ps->nocclA,  sizeof(int) * (ps->npage + ntoadd));
    ESL_ALLOC(ps->sclAA,   sizeof(SchemeColorLegend_t **) * (ps->npage + ntoadd));
    ESL_ALLOC(ps->descA,   sizeof(char *) * (ps->npage + ntoadd));
    ESL_ALLOC(ps->modeA,   sizeof(int) * (ps->npage + ntoadd));
    ESL_ALLOC(ps->seqidxA, sizeof(int) * (ps->npage + ntoadd));
  }
  else { 
    assert(ps->rrAA    != NULL);
    assert(ps->rcolAAA != NULL);
    ESL_RALLOC(ps->rrAA,   tmp, sizeof(char *)   * (ps->npage + ntoadd));
    ESL_RALLOC(ps->rcolAAA,tmp, sizeof(float **) * (ps->npage + ntoadd));
    ESL_RALLOC(ps->occlAAA,tmp, sizeof(OneCellColorLegend_t ***) * (ps->npage + ntoadd));
    ESL_RALLOC(ps->nocclA, tmp, sizeof(int) * (ps->npage + ntoadd));
    ESL_RALLOC(ps->sclAA,  tmp, sizeof(SchemeColorLegend_t **) * (ps->npage + ntoadd));
    ESL_RALLOC(ps->descA,  tmp, sizeof(char *) * (ps->npage + ntoadd));
    ESL_RALLOC(ps->modeA,  tmp, sizeof(int) * (ps->npage + ntoadd));
    ESL_RALLOC(ps->seqidxA,tmp, sizeof(int) * (ps->npage + ntoadd));
  }
  for(p = ps->npage; p < (ps->npage + ntoadd); p++) { 
    ps->rrAA[p]    = NULL;
    ps->rcolAAA[p] = NULL;
    ps->occlAAA[p] = NULL;
    ps->nocclA[p]  = 0;
    ps->sclAA[p]   = NULL;
    ps->descA[p]   = NULL;
    ps->modeA[p]   = page_mode;
    ps->seqidxA[p] = -1;
  }
  ps->npage += ntoadd;

  return eslOK;

 ERROR: 
  return status;
}

/* read_mask_file
 *
 * Given an open file pointer, read the first token of the
 * file and return it as *ret_mask. Also return length of
 * mask in *ret_masklen, and a flag indicating whether or
 * not the mask has any internal zeroes in *ret_mask_has_internal_zeroes.
 * An internal '0' is one that occurs between at least one
 * 5' and one 3' '1'
 *
 * Returns:  eslOK on success.
 */
int
read_mask_file(char *filename, char *errbuf, char **ret_mask, int *ret_masklen, int *ret_mask_has_internal_zeroes)
{
  int             status;
  ESL_FILEPARSER *efp;
  char           *tok;
  char           *mask;
  int             toklen;
  int             n;
  /* for determining if we have internal zeroes */
  int             seen_1 = FALSE;                /* becomes TRUE when we see first (5'-most) '1' */
  int             seen_1_then_0 = FALSE;         /* becomes TRUE when we see first (5'-most) '0' that is 3' of first '1' */
  int             seen_1_then_0_then_1 = FALSE;  /* becomes TRUE when we see first (5'-most) '1' that is 3' of first '0' that is 3' of first '1' */

  if (esl_fileparser_Open(filename, NULL, &efp) != eslOK) ESL_FAIL(eslFAIL, errbuf, "failed to open %s in read_mask_file\n", filename);
  esl_fileparser_SetCommentChar(efp, '#');
  
  if((status = esl_fileparser_GetToken(efp, &tok, &toklen)) != eslOK) ESL_FAIL(eslFAIL, errbuf, "failed to read a single token from %s\n", filename);

  ESL_ALLOC(mask, sizeof(char) * (toklen+1));

  for(n = 0; n < toklen; n++) { 
    mask[n] = tok[n];
    if(mask[n] == '0') { 
      if((seen_1) && (!seen_1_then_0)) { seen_1_then_0 = TRUE; }
    }
    else if (mask[n] == '1') { 
      if(!seen_1) { seen_1 = TRUE; }
      if((seen_1) && (seen_1_then_0) && (!seen_1_then_0_then_1)) { seen_1_then_0_then_1 = TRUE; }
    }
    else { ESL_FAIL(eslEINVAL, errbuf, "character %d of mask file is invalid: %c (must be a '1' or a '0')\n", n, mask[n]); }

    mask[n] = tok[n];
  }
  mask[n] = '\0';

  *ret_mask = mask;
  *ret_masklen= toklen;
  *ret_mask_has_internal_zeroes = seen_1_then_0_then_1;

  esl_fileparser_Close(efp);
  return eslOK;
  
 ERROR:
  return eslEMEM;
}

/* Function: drawfile2sspostscript()
 * 
 * Purpose:  Fill a postscript data structure with >= 1 new page(s), with colors described
 *           in an input 'draw' file, with >= 1 sets of <x> lines of data, each set 
 *           is separated by a line with only "//". <x> must be equal to the consensus
 *           ps->rflen. Each line has at least 4 floats explaining 
 *           the CMYK values for the color to use at each position of the SS diagram,
 *           and optionally contains an extra single character which is the residue
 *           to put at that position.
 *           
 * Return:   eslOK on success.
 */
static int
drawfile2sspostscript(const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps)
{
  int status;
  int p, pp;
  int cpos, c;
  int orig_npage = ps->npage;
  ESL_FILEPARSER *efp;
  char           *s;
  char *dfile = esl_opt_GetString(go, "--dfile");

  if (esl_fileparser_Open(dfile, NULL, &efp) != eslOK) ESL_FAIL(eslFAIL, errbuf, "failed to open %s in draw_file2sspostscript\n", dfile);
  esl_fileparser_SetCommentChar(efp, '#');

  pp = orig_npage - 1;
  cpos = 0;

  while (esl_fileparser_NextLine(efp) == eslOK)
    {
      /* example line without residue markup:
       * 0.000 0.000 0.000 0.500
       *
       * example line with residue markup:
       * 0.000 0.000 0.000 0.500 A
       */

      cpos++;
      if(cpos == 1) { /* add a new page */
	if((status = add_pages_sspostscript(ps, 1, SIMPLEMASKMODE)) != eslOK) ESL_FAIL(status, errbuf, "memory error adding pages to the postscript object.");
	
	for(p = (ps->npage-1); p < ps->npage; p++) { 
	  ESL_ALLOC(ps->rrAA[p], sizeof(char) *  (ps->rflen+1));
	  ESL_ALLOC(ps->rcolAAA[p], sizeof(float *) * ps->rflen);
	  for(c = 0; c < ps->rflen; c++) { 
	    ESL_ALLOC(ps->rcolAAA[p][c], sizeof(float) * NCMYK); /* CMYK colors */
	  }
	}
	pp++; /* if first page, pp == orig_npage now */
      }
      if(cpos == (ps->rflen+1)) { /* should be a single token, a "\\" on this line */ 
	if (esl_fileparser_GetTokenOnLine(efp, &s, NULL) != eslOK)
	  esl_fatal("Failed to read a final token at the end of description of draw page %d on line %d of drawfile %s\n", (pp - orig_npage + 1), efp->linenumber, dfile);
	if (strcmp(s, "//") != 0) 
	  esl_fatal("Failed to read a final \"//\" token (read %s) at the end of description of draw page %d on line %d of drawfile %s\n", s, (pp - orig_npage + 1), efp->linenumber, dfile);
	cpos = 0;
      }
      else { 
	/* get C value */
	if (esl_fileparser_GetTokenOnLine(efp, &s, NULL) != eslOK)
	  esl_fatal("Failed to read C of CMYK value on line %d of drawfile %s\n", efp->linenumber, dfile);
	ps->rcolAAA[pp][(cpos-1)][0] = atof(s);

	/* get M value */
	if (esl_fileparser_GetTokenOnLine(efp, &s, NULL) != eslOK)
	  esl_fatal("Failed to read M of CMYK value on line %d of drawfile %s\n", efp->linenumber, dfile);
	ps->rcolAAA[pp][(cpos-1)][1] = atof(s);

	/* get Y value */
	if (esl_fileparser_GetTokenOnLine(efp, &s, NULL) != eslOK)
	esl_fatal("Failed to read Y of CMYK value on line %d of drawfile %s\n", efp->linenumber, dfile);
	ps->rcolAAA[pp][(cpos-1)][2] = atof(s);

	/* get K value */
	if (esl_fileparser_GetTokenOnLine(efp, &s, NULL) != eslOK)
	  esl_fatal("Failed to read K of CMYK value on line %d of drawfile %s\n", efp->linenumber, dfile);
	ps->rcolAAA[pp][(cpos-1)][3] = atof(s);

	/* optionally read a residue value */
	if (esl_fileparser_GetTokenOnLine(efp, &s, NULL) == eslOK) {
	  if(((int) strlen(s)) != 1) esl_fatal("Read multi-character string (%s) for consensus residue %d on line %d of drawfile %s\n", s, cpos, efp->linenumber, dfile);
	  ps->rrAA[pp][(cpos-1)] = s[0];
	}
	else ps->rrAA[pp][(cpos-1)] = ' ';
      }
    }
  if(pp == (orig_npage - 1)) { /* no new pages were read, this is an error */
    esl_fatal("Failed to read a single page from drawfile %s\n", dfile);
  }

  esl_fileparser_Close(efp);
  return eslOK;

 ERROR: ESL_FAIL(status, errbuf, "drawfile2sspostscript(): memory allocation error.");
  return status; /* NEVERREACHED */
}


/* Function: mutual_information_sspostscript()
 * 
 * Purpose:  Fill a postscript data structure with 1 new page, colored squares indicating
 *           the mutual information of each base paired consensus column.
 *           mutual information is the extra information gained from modelling
 *           the pair together (info of vector of bps, size 16) versus separately (sum
 *           of info of the two independent vector of singlets, size 4).
 *           
 * Return:   eslOK on success.
 */
static int
mutual_information_sspostscript(const ESL_GETOPTS *go, ESL_ALPHABET *abc, char *errbuf, SSPostscript_t *ps, double ***bp_ct, int msa_nseq, float ***hc_scheme, int hc_scheme_idx, int hc_nbins, float **hc_onecell, int ss_idx, int zerores_idx, FILE *tabfp)
{
  int status;
  int p, pp, c, i;
  int rfpos, apos; 
  int orig_npage = ps->npage;
  double *bg, *bg_pair;
  double bg_ent, bg_pair_ent;
  double *obs_left, *obs_right, *obs_pair;
  double ent_left, ent_right, ent_pair;
  int j;
  ESL_DSQ lres;
  ESL_DSQ rres;
  double nres;
  int nss = 0;
  int nzerores = 0;
  int nss_masked = 0;
  int nzerores_masked = 0;
  int i_within_mask, j_within_mask;
  int i_bi, j_bi;
  float *limits;
  int l;
  int idx =1;

  if(ps->mask == NULL) nss_masked = nzerores_masked = -1; /* special flag */
  if((status = add_pages_sspostscript(ps, 1, ALIMODE)) != eslOK) ESL_FAIL(status, errbuf, "memory error adding pages to the postscript object.");

  for(p = orig_npage; p < ps->npage; p++) { 
    ESL_ALLOC(ps->rrAA[p], sizeof(char) *  ps->rflen);
    ESL_ALLOC(ps->rcolAAA[p], sizeof(float *) * ps->rflen);
    ESL_ALLOC(ps->sclAA[p],   sizeof(SchemeColorLegend_t *) * 1);
    ESL_ALLOC(ps->occlAAA[p], sizeof(OneCellColorLegend_t **) * 2);
    for(c = 0; c < ps->rflen; c++) { 
      ESL_ALLOC(ps->rcolAAA[p][c], sizeof(float) * NCMYK); /* CMYK colors */
    }
  }
  pp = orig_npage;

  /* add color legend */
  ESL_ALLOC(limits, sizeof(float) * (hc_nbins+1)); 
  limits[0] = 0.0;
  limits[1] = 0.167;
  limits[2] = 0.333;
  limits[3] = 0.500;
  limits[4] = 0.667;
  limits[5] = 0.833;
  limits[6] = 1.000;
  ps->sclAA[pp] = create_scheme_colorlegend(hc_scheme_idx, hc_nbins, limits, FALSE);

  if(tabfp != NULL) { 
    fprintf(tabfp, "# -----------------------\n");
    fprintf(tabfp, "# Mutual information data\n");
    fprintf(tabfp, "# -----------------------\n");
    fprintf(tabfp, "# This section includes %d non #-prefixed lines, one for each consensus position\n", ps->rflen);
    fprintf(tabfp, "# in the alignment and corresponding template.\n");
    fprintf(tabfp, "# Each line includes %d tokens, separated by whitespace:\n", ps->mask == NULL ? 9 : 11);
    fprintf(tabfp, "# \ttoken  1: 'mutualinfo' (tag defining line type to ease parsing)\n");
    fprintf(tabfp, "# \ttoken  2: base pair index\n");
    fprintf(tabfp, "# \ttoken  3: 5' consensus position of base pair (starting at 1)\n");
    fprintf(tabfp, "# \ttoken  4: 3' consensus position of base pair (starting at 1)\n");
    fprintf(tabfp, "# \ttoken  5: sequence information content at 5' position (bits)\n");
    fprintf(tabfp, "# \ttoken  6: sequence information content at 3' position (bits)\n");
    fprintf(tabfp, "# \ttoken  7: mutual information of the base pair (bits)\n");
    fprintf(tabfp, "# \ttoken  8: number of sequences with non-gap at 5' and 3' posn (max possible is %d)\n", msa_nseq); 
    fprintf(tabfp, "# \ttoken  8: number of sequences with non-gap at 5' and 3' position\n");
    fprintf(tabfp, "# \ttoken  9: bin index this positions falls in (see bin values below).\n");
    if(ps->mask != NULL) { 
      fprintf(tabfp, "# \ttoken 10: '1' if 5' position is included by mask, '0' if not\n");
      fprintf(tabfp, "# \ttoken 11: '1' if 3' position is included by mask, '0' if not\n");
    }
    fprintf(tabfp, "#\n");
    fprintf(tabfp, "# Information content is calculated as 2.0 - H, where\n");
    fprintf(tabfp, "# H = - \\sum_x p_x \\log_2 p_x, for x in {A, C, G, U}\n");
    fprintf(tabfp, "# p_x is the frequency of x for *non-gap* residues at the position.\n");
    fprintf(tabfp, "# Only residues for sequences which have a non-gap residue at both\n"); 
    fprintf(tabfp, "# the 5' and 3' positions of the pair are counted.\n"); 
    fprintf(tabfp, "# Mutual information is calculated as\n");
    fprintf(tabfp, "# \\sum_{x,y} p_{x,y} \\log_2 ((p_x * p_y) / p_{x,y}\n");
    fprintf(tabfp, "# Value ranges for bins:\n");
    fprintf(tabfp, "# \tbin  0: special case, 0 sequences have non-gaps at both 5' and 3' position of pair\n");
	    for(l = 0; l < hc_nbins; l++) { 
      fprintf(tabfp, "# \tbin %2d: [%.3f-%.3f%s mutual information per position (bits)\n", l+1, limits[l], limits[l+1], (l == hc_nbins-1) ? "]" : ")");
    }
    fprintf(tabfp, "#\n");
    fprintf(tabfp, "# %10s  %4s  %5s  %5s  %8s  %8s  %9s  %10s  %3s", "type", "idx", "5'pos", "3'pos", "5'info", "3'info", "mutinfo/2", "nongap", "bin");
    if(ps->mask != NULL) fprintf(tabfp, "  %6s  %6s", "5'mask", "3'mask");
    fprintf(tabfp, "\n");
    fprintf(tabfp, "# %10s  %4s  %5s  %5s  %8s  %8s  %9s  %10s  %3s", "----------", "----", "-----", "-----", "--------", "--------", "---------", "----------", "---");
    if(ps->mask != NULL) fprintf(tabfp, "  %6s  %6s", "------", "------");
    fprintf(tabfp, "\n");
  }

  /* determine background entropy */
  ESL_ALLOC(bg,        sizeof(double) * abc->K);
  ESL_ALLOC(bg_pair,   sizeof(double) * abc->K*abc->K);
  esl_vec_DSet(bg, abc->K, 1./(abc->K));
  esl_vec_DSet(bg_pair, abc->K*abc->K, 1./(abc->K*abc->K));
  bg_pair_ent = esl_vec_DEntropy(bg_pair, abc->K*abc->K);
  bg_ent      = esl_vec_DEntropy(bg, abc->K);
  free(bg);
  free(bg_pair);

  ESL_ALLOC(obs_left,  sizeof(double) * (abc->K));
  ESL_ALLOC(obs_right, sizeof(double) * (abc->K));
  ESL_ALLOC(obs_pair,  sizeof(double) * (abc->K*abc->K));

  /* get observed residues and base pairs at each rfpos */
  for(rfpos = 0; rfpos < ps->rflen; rfpos++) {
    esl_vec_DSet(obs_left,  abc->K, 0.);
    esl_vec_DSet(obs_right, abc->K, 0.);
    esl_vec_DSet(obs_pair,  abc->K*abc->K, 0.);
    i = rfpos;
    nres = 0.;
    apos = ps->msa_rf2a_map[rfpos];
    /* check if we're base paired */
    if(ps->msa_ct[rfpos+1] != 0) { 
      /* printf("msa_ct rfpos+1: %d %d\n", rfpos+1, ps->msa_ct[rfpos+1]); */
      if(ps->msa_ct[rfpos+1] > (rfpos+1)) { /* rfpos is left half of base pair */
	/* add up contributions of all possible base pairs,
	 * we have to be careful to skip residues that are gaps, missing or nonresidues in 
	 * EITHER left or right half of the bp.
	 * We use the raw counts from msa_count() as the wt, (bp_ct[apos][lres][rres]) */
	j = ps->msa_ct[i+1] - 1; 
	for(lres = 0; lres < abc->K; lres++) { 
	  for(rres = 0; rres < abc->K; rres++) { 
	    /* printf("apos: %d rfpos: %d lres: %d rres: %d bp_ct[][][]: %.5f\n", apos, rfpos, lres, rres, bp_ct[apos][lres][rres]);*/
	    esl_abc_DCount(abc, obs_left,  lres, bp_ct[apos][lres][rres]);
	    esl_abc_DCount(abc, obs_right, rres, bp_ct[apos][lres][rres]);
	    PairCount(abc, obs_pair, lres, rres, bp_ct[apos][lres][rres]);
	    nres += bp_ct[apos][lres][rres];
	  }
	}
	for(lres = abc->K+1; lres < abc->Kp-2; lres++) { /* do all degenerates and 'any' (N) */
	  for(rres = abc->K+1; rres < abc->Kp-2; rres++) { /* do all degenerates and 'any' (N) */
	    esl_abc_DCount(abc, obs_left,  lres, bp_ct[apos][lres][rres]);
	    esl_abc_DCount(abc, obs_right, rres, bp_ct[apos][lres][rres]);
	    PairCount(abc, obs_pair, lres, rres, bp_ct[apos][lres][rres]);
	    nres += bp_ct[apos][lres][rres];
	  }
	}
	/* esl_vec_DDump(stdout, obs_left, abc->K, NULL);
	   esl_vec_DDump(stdout, obs_right, abc->K, NULL);
	   esl_vec_DDump(stdout, obs_pair, abc->K*abc->K, NULL);
	*/
	esl_vec_DNorm(obs_left,  abc->K);
	esl_vec_DNorm(obs_right, abc->K);
	esl_vec_DNorm(obs_pair, abc->K*abc->K);
	ent_left  = bg_ent - esl_vec_DEntropy(obs_left, abc->K);      
	ent_right = bg_ent - esl_vec_DEntropy(obs_right, abc->K);      
	ent_pair =  bg_pair_ent - esl_vec_DEntropy(obs_pair, abc->K*abc->K);      
	/* printf("lpos: %5d  rpos: %5d  entP: %8.3f  entL: %8.3f  entR: %8.3f  nres: %.4f  ",  i+1, j+1, ent_pair, ent_left, ent_right, nres); */
	ent_pair -= ent_left + ent_right;
	ent_pair /= 2.;
	/* printf("Final: %8.3f\n", ent_pair);  */
	
	/* To verify that the ent_pair is mutual information, calculated a different way, uncomment the following block */
	/* double mi = 0;
	   for(lres = 0; lres < abc->K; lres++) { 
	   for(rres = 0; rres < abc->K; rres++) { 
	   if(obs_pair[lres*abc->K+rres] > eslSMALLX1) { 
	   mi += obs_pair[lres*abc->K+rres] * (1.44269504 * log((obs_pair[lres*abc->K+rres])/(obs_left[lres] * obs_right[rres])));
	   printf("mi: %.4f obs_left %.4f  obs_right: %.4f obs_pair %.4f \n", mi, obs_left[lres], obs_right[rres], obs_pair[lres*abc->K+rres]);  
	   }
	   }
	   }
	   printf("MI/2: %.4f  EP: %.4f  %.4f\n", mi/2., ent_pair, (mi/2.) - ent_pair); 
	*/

	if(ent_pair < (-1. * eslSMALLX1)) { 
	  ESL_FAIL(eslEINCONCEIVABLE, errbuf, "pair information < 0.: %f (lpos: %d rpos: %d)\n", ent_pair, i, j);
	}
	if(esl_DCompare(nres, 0., eslSMALLX1) == eslOK) { /* nres is 0 */
	  if((status = set_onecell_values(errbuf, ps->rcolAAA[pp][i], NCMYK, hc_onecell[zerores_idx])) != eslOK) return status;
	  if((status = set_onecell_values(errbuf, ps->rcolAAA[pp][j], NCMYK, hc_onecell[zerores_idx])) != eslOK) return status;
	  nzerores += 2;
	  if(ps->mask != NULL && ps->mask[i] == '1') nzerores_masked++; 
	  if(ps->mask != NULL && ps->mask[j] == '1') nzerores_masked++; 
	  if(tabfp != NULL) { 
	    fprintf(tabfp, "  mutualinfo  %4d  %5d  %5d  %8.5f  %8.5f  %9.5f  %10d  %3d", idx++, i+1, j+1, 0., 0., 0., 0, 0);
	    if(ps->mask != NULL) fprintf(tabfp, "  %6d  %6d", ((ps->mask == NULL || ps->mask[i] == '1') ? 1 : 0), ((ps->mask == NULL || ps->mask[j] == '1') ? 1 : 0));
	    fprintf(tabfp, "\n");
	  }
	}
	else { 
	  i_within_mask = (ps->mask != NULL && ps->mask[i] == '1') ? TRUE : FALSE;
	  j_within_mask = (ps->mask != NULL && ps->mask[j] == '1') ? TRUE : FALSE;
	  if((status = set_scheme_values(errbuf, ps->rcolAAA[pp][i], NCMYK, hc_scheme[hc_scheme_idx], ent_pair, ps->sclAA[pp], i_within_mask, &i_bi)) != eslOK) return status;
	  if((status = set_scheme_values(errbuf, ps->rcolAAA[pp][j], NCMYK, hc_scheme[hc_scheme_idx], ent_pair, ps->sclAA[pp], j_within_mask, &j_bi)) != eslOK) return status;
	  if(tabfp != NULL) { 
	    fprintf(tabfp, "  mutualinfo  %4d  %5d  %5d  %8.5f  %8.5f  %9.5f  %10d  %3d", 
		    idx++, i+1, j+1, 
		    ent_left, 
		    ent_right, 
		    ent_pair,
		    (int) nres, 
		    i_bi+1);
	    if(ps->mask != NULL) fprintf(tabfp, "  %6d  %6d", ((ps->mask == NULL || ps->mask[i] == '1') ? 1 : 0), ((ps->mask == NULL || ps->mask[j] == '1') ? 1 : 0));
	    fprintf(tabfp, "\n");
	  }
	}
      } /* end of if(ps->msa_ct[rfpos+1] > (rfpos+1)) { */
    }
    else { /* single stranded, paint grey  */
      nss++;
      if(ps->mask != NULL && ps->mask[rfpos] == '1') nss_masked++; 
      if((status = set_onecell_values(errbuf, ps->rcolAAA[pp][rfpos], NCMYK, hc_onecell[ss_idx])) != eslOK) return status;
    }
    ps->rrAA[pp][rfpos] = ' '; /* leave residue blank, no text, just color (rcolAAA) */
  }

  /* add text to the one cell legend */
  ps->occlAAA[pp][0] = create_onecell_colorlegend(hc_onecell[ss_idx], nss, nss_masked);
  if((status = add_text_to_onecell_colorlegend(ps, ps->occlAAA[pp][0], "single-stranded", ps->legx_max_chars, errbuf)) != eslOK) return status;
  
  /* add text to the second one cell legend */
  ps->occlAAA[pp][1] = create_onecell_colorlegend(hc_onecell[zerores_idx], nzerores, nzerores_masked);
  if((status = add_text_to_onecell_colorlegend(ps, ps->occlAAA[pp][1], "0 complete basepairs", ps->legx_max_chars, errbuf)) != eslOK) return status;
  ps->nocclA[pp] = 2;
  
  /* add text to the scheme legend */
  if((status = add_text_to_scheme_colorlegend(ps->sclAA[pp], "mutual information per position (bits)", ps->legx_max_chars, errbuf)) != eslOK) return status;
  
  /* add description to ps */
  if((status = add_page_desc_to_sspostscript(ps, pp, "mutual information per basepaired position", errbuf)) != eslOK) return status;
  
  free(limits);
  free(obs_left);
  free(obs_right);
  free(obs_pair);
  
  if(tabfp != NULL) fprintf(tabfp, "//\n");
  return eslOK;
  
 ERROR: ESL_FAIL(status, errbuf, "mutual_information_sspostscript(): memory allocation error.");
  return status; /* NEVERREACHED */
}

/* Function: PairCount()
 * Date:     SRE, Tue Aug  1 10:34:20 2000 [St. Louis]
 *
 * Purpose:  Given a possibly degenerate symbol code for left
 *           and right symbols in a pair, increment a symbol
 *           counter array appropriately.
 *           
 * Args:     abc      - pointer to the internal alphabet
 *           counters - vector to count into [0..abc->K^2-1]
 *           syml     - index of left symbol  [0..abc->sym_iupac-1]
 *           symr     - index of right symbol [0..abc->sym_iupac-1]
 *           wt       - weight to use for the count (often 1.0).          
 *
 * Returns:  void
 */
static void
PairCount(const ESL_ALPHABET *abc, double *counters, ESL_DSQ syml, ESL_DSQ symr, double wt)
{
  int status;
  if (syml < abc->K && symr < abc->K) {
    counters[(int) (syml * abc->K + symr)] += wt;
    return;
  }
  else {
    int   l,r;
    double *left = NULL;
    double *right = NULL;
    ESL_ALLOC(left,  sizeof(double) * abc->K);
    ESL_ALLOC(right, sizeof(double) * abc->K);
    
    esl_vec_DSet(left,  abc->K, 0.);
    esl_vec_DSet(right, abc->K, 0.);
    esl_abc_DCount(abc, left,  syml, 1.0);
    esl_abc_DCount(abc, right, symr, 1.0);

    for (l = 0; l < abc->K; l++)
      for (r = 0; r < abc->K; r++)
	counters[l*abc->K +r] += left[l] * right[r] * wt;
    free(left);
    free(right);
  }
  return;

 ERROR:
  esl_fatal("Memory error");
}

/* Function: get_command
 * Date:     EPN, Fri Jan 25 13:56:10 2008
 *
 * Purpose:  Return the command used to call ssu-draw
 *           in <ret_command>.
 *
 * Returns:  eslOK on success; eslEMEM on allocation failure.
 */
static int 
get_command(const ESL_GETOPTS *go, char *errbuf, char **ret_command)
{
  int status;
  int i;
  char *command = NULL;

  for (i = 0; i < go->argc; i++) { /* copy all command line options and args */
    if((status = esl_strcat(&(command),  -1, go->argv[i], -1)) != eslOK) goto ERROR;
    if(i < (go->argc-1)) if((status = esl_strcat(&(command), -1, " ", 1)) != eslOK) goto ERROR;
  }
  *ret_command = command;

  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "get_command(): memory allocation error.");
  return status;
}

/* Function: get_date
 * Date:     EPN, Fri Jan 25 13:59:22 2008
 *
 * Purpose:  Return a string that gives the current date.
 *
 * Returns:  eslOK on success; eslEMEM on allocation failure.
 */
static int 
get_date(char *errbuf, char **ret_date)
{
  int    status;
  time_t date = time(NULL);
  char  *sdate = NULL;

  if((status = esl_strdup(ctime(&date), -1, &sdate)) != eslOK) goto ERROR;
  esl_strchop(sdate, -1); /* doesn't return anything but eslOK */

  *ret_date = sdate;
  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "get_date() error status: %d, probably out of memory.", status);
  return status; 
}

/* Function: set_scheme_values()
 * 
 * Purpose:  Set color values from a predefined scheme given min, max, 
 *           value and number of colors.
 */ 
static int 
set_scheme_values(char *errbuf, float *vec, int ncolvals, float **scheme, float val, SchemeColorLegend_t *scl, int within_mask, int *ret_bi)
{
  float min, max;
  int ci, bi;
  min = scl->limits[0];
  max = scl->limits[scl->nbins];
  if((min-val) > eslSMALLX1) { ESL_FAIL(eslEINVAL, errbuf, "set_scheme_values(), val: %.4f < min: %.4f\n", val, min); }
  if((val-max) > eslSMALLX1) { ESL_FAIL(eslEINVAL, errbuf, "set_scheme_values(), val: %.4f > max: %.4f\n", val, max); }

  bi = 0;
  while((bi < (scl->nbins-1)) &&
	((val > scl->limits[bi+1]) ||                                    /* val exceeds limits[bi+1] OR */
	 (esl_FCompare(val, scl->limits[bi+1], eslSMALLX1) == eslOK))) { /* val equals limits[bi+1] */
    bi++; 
  }    
  /* printf("%.3f %d (%.3f)\n", val, bi, scl->limits[bi+1]);  */
  scl->counts[bi]++;
  if(within_mask) scl->counts_masked[bi]++;
  for(ci = 0; ci < ncolvals; ci++) { 
    vec[ci] = scheme[bi][ci]; 
  }

  if(ret_bi != NULL) *ret_bi = bi;
  return eslOK;
}

/* Function: set_onecell_values()
 * 
 * Purpose:  Set color values as a  predefined single
 *           color.
 */ 
static int 
set_onecell_values(char *errbuf, float *vec, int ncolvals, float *onecolor)
{
  int ci;
  for(ci = 0; ci < ncolvals; ci++) { vec[ci] = onecolor[ci]; }
  return eslOK;
}
  


/* Function: draw_masked_block()
 * 
 * Purpose:  Given coords, color, and mask style options draw a masked block.
 */ 
static int 
draw_masked_block(FILE *fp, float x, float y, float *colvec, int do_circle_mask, int do_square_mask, int do_x_mask, int do_border, float boxsize)
{
  if (do_circle_mask) { 
    if(do_border) { 
      fprintf(fp, "newpath\n");
      fprintf(fp, " %.2f %.2f %.1f 0 360 arc closepath\n", x + (boxsize/2.), y + (boxsize/2.), boxsize * (3./8.));
      fprintf(fp, "  %.4f %.4f %.4f %.4f setcmykcolor\n", colvec[0], colvec[1], colvec[2], colvec[3]);
      fprintf(fp, "  stroke\n"); 
    }
    else { 
      fprintf(fp, "newpath\n");
      fprintf(fp, " %.2f %.2f %.1f 0 360 arc closepath\n", x + (boxsize/2.), y + (boxsize/2.), boxsize * (3./8.));
      fprintf(fp, "  %.4f %.4f %.4f %.4f setcmykcolor\n", colvec[0], colvec[1], colvec[2], colvec[3]);
      fprintf(fp, "  fill\n"); 
    }
  }
  else if(do_square_mask) { 
    if(do_border) { 
      fprintf(fp, "newpath\n");
      fprintf(fp, "  %.2f %.2f moveto", x +1., y +1.);
      fprintf(fp, "  0 %.1f rlineto %.1f 0 rlineto 0 -%.1f rlineto closepath\n", boxsize*0.75, boxsize*0.75, boxsize*0.75);
      fprintf(fp, "  %.4f %.4f %.4f %.4f setcmykcolor\n", colvec[0], colvec[1], colvec[2], colvec[3]);
      fprintf(fp, "  stroke\n");
    }
    else { 
      fprintf(fp, "newpath\n");
      fprintf(fp, "  %.2f %.2f moveto", x + 1.5, y + 1.5);
      fprintf(fp, "  0 %.1f rlineto %.1f 0 rlineto 0 -%.1f rlineto closepath\n", boxsize*(5./8.), boxsize*(5./8.), boxsize*(5./8.));
      fprintf(fp, "  %.4f %.4f %.4f %.4f setcmykcolor\n", colvec[0], colvec[1], colvec[2], colvec[3]);
      fprintf(fp, "  fill\n"); 
    }
  }
  else if (do_x_mask) { 
    if(do_border) { 
      fprintf(fp, "newpath\n");
      fprintf(fp, "  %.2f %.2f moveto", x, y);
      fprintf(fp, "  0 %.1f rlineto %.1f 0 rlineto 0 -%.1f rlineto closepath\n", boxsize, boxsize, boxsize);
      fprintf(fp, "  %.4f %.4f %.4f %.4f setcmykcolor\n", colvec[0], colvec[1], colvec[2], colvec[3]);
      fprintf(fp, "  fill\n"); 
      
      fprintf(fp, "  %.4f %.4f %.4f %.4f setcmykcolor\n", 0., 0., 0., 0.);
      fprintf(fp, "newpath\n");
      fprintf(fp, "  %.2f %.2f moveto", x, y);
      fprintf(fp, "  %.1f %.1f rlineto closepath\n", boxsize, boxsize);
      fprintf(fp, "  stroke\n"); 
      fprintf(fp, "  %.2f %.2f moveto", x + boxsize, y);
      fprintf(fp, "  -%.1f %.1f rlineto closepath\n", boxsize, boxsize);
      fprintf(fp, "  stroke\n"); 
    }
    else { 
      fprintf(fp, "newpath\n");
      fprintf(fp, "  %.4f %.4f %.4f %.4f setcmykcolor\n", colvec[0], colvec[1], colvec[2], colvec[3]);
      fprintf(fp, "  %.2f %.2f moveto", x, y);
      fprintf(fp, "  %.1f %.1f rlineto closepath\n", boxsize, boxsize);
      fprintf(fp, "  stroke\n"); 
      fprintf(fp, "newpath\n");
      fprintf(fp, "  %.2f %.2f moveto", x + boxsize, y);
      fprintf(fp, "  -%.1f %.1f rlineto closepath\n", boxsize, boxsize);
      fprintf(fp, "  stroke\n"); 
    }
  }
  return eslOK;
}

/* Function: validate_justread_sspostscript()
 * 
 * Purpose:  Validate a sspostscript just created by parsing
 *           a template file. Nothing fancy here, just make
 *           sure all we've written everything we expect.
 */ 
static int 
validate_justread_sspostscript(SSPostscript_t *ps, char *errbuf)
{
  if(ps->modelname == NULL) ESL_FAIL(eslEINVAL, errbuf, "validate_justread_sspostscript(), failed to read modelname from template file.");
  if(ps->nbp == 0) ESL_FAIL(eslEINVAL, errbuf, "validate_justread_sspostscript(), failed to read 'lines bpconnects' section from template file.");
  if(ps->rflen == 0) ESL_FAIL(eslEINVAL, errbuf, "validate_justread_sspostscript(), failed to read 'text residues' section from template file.");

  /* Stuff we don't currently require, but we may want to eventually */
  /*if(ps->nhundreds == 0) ESL_FAIL(eslEINVAL, errbuf, "validate_justread_sspostscript(), failed to read 'text hundreds' section from template file.");*/
  /*if(ps->nticks == 0) ESL_FAIL(eslEINVAL, errbuf, "validate_justread_sspostscript(), failed to read 'lines ticks' section from template file.");*/
  /*if(ps->nbp == 0) ESL_FAIL(eslEINVAL, errbuf, "validate_justread_sspostscript(), failed to read 'lines bpconnects' section from template file.");*/

  return eslOK;
}


/* Function: validate_and_update_sspostscript_given_msa()
 * 
 * Purpose:  Validate that a sspostscript works with a MSA.
 */ 
static int 
validate_and_update_sspostscript_given_msa(const ESL_GETOPTS *go, SSPostscript_t *ps, ESL_MSA *msa, int msa_nseq, char *errbuf)
{
  int status;
  int *msa_ct;
  int msa_nbp = 0;
  int *tmp_ct;
  int apos, rfpos_i, rfpos_j, rfpos;
  int *rf2a_map = NULL;
  int *a2rf_map = NULL;
  int rflen = 0;

  ps->msa = msa;

  /* contract check */
  if(msa->rf == NULL)      ESL_FAIL(eslEINVAL, errbuf, "Error, msa does not have RF annotation.");
  if(msa->ss_cons == NULL) ESL_FAIL(eslEINVAL, errbuf, "Error, msa does not have SS_cons annotation.");

  /* count non-gap RF columns */
  for(apos = 0; apos < msa->alen; apos++) { 
    if((! esl_abc_CIsGap(msa->abc, msa->rf[apos])) && 
       (! esl_abc_CIsMissing(msa->abc, msa->rf[apos])) && 
       (! esl_abc_CIsNonresidue(msa->abc, msa->rf[apos])))
      { 
	rflen++;
	/* I don't use esl_abc_CIsResidue() b/c that would return FALSE for 'x' with RNA and DNA */
      }
  }
  if(ps->rflen != rflen) ESL_FAIL(eslEINVAL, errbuf, "validate_and_update_sspostscript_given_msa(), expected consensus length of %d in MSA, but read %d\n", ps->rflen, rflen);

  /* build map of non-gap RF positions to alignment positions and vice versa */
  ESL_ALLOC(rf2a_map, sizeof(int) * rflen);
  ESL_ALLOC(a2rf_map, sizeof(int) * msa->alen);
  esl_vec_ISet(a2rf_map, msa->alen, -1); 
  rfpos = 0;
  for(apos = 0; apos < msa->alen; apos++) {
    if((! esl_abc_CIsGap(msa->abc, msa->rf[apos])) && 
       (! esl_abc_CIsMissing(msa->abc, msa->rf[apos])) && 
       (! esl_abc_CIsNonresidue(msa->abc, msa->rf[apos]))) { 
      rf2a_map[rfpos] = apos;
      a2rf_map[apos]  = rfpos;
      rfpos++;
    }
  }

  /* get the CT array for this msa */
  ESL_ALLOC(tmp_ct, sizeof(int) * (msa->alen+1));
  if (esl_wuss2ct(msa->ss_cons, msa->alen, tmp_ct) != eslOK) ESL_FAIL(status, errbuf, "Problem getting ct from SS_cons, does first alignment of MSA file have SS_cons annotation?");

  ESL_ALLOC(msa_ct, sizeof(int) * (rflen+1));
  esl_vec_ISet(msa_ct, rflen+1, 0);
  
  /* convert tmp_ct which is in alignment coords [1..alen] to non-gap RF coords [1..rflen]*/
  for(apos = 0; apos < msa->alen; apos++) {
    if(tmp_ct[apos+1] > (apos+1)) { 
      rfpos_i = a2rf_map[apos];
      rfpos_j = a2rf_map[tmp_ct[apos+1]-1];
      if(rfpos_i != -1 && rfpos_j != -1) { /* a consensus basepair */
	msa_ct[rfpos_i+1] = rfpos_j+1;
	msa_ct[rfpos_j+1] = rfpos_i+1;
	msa_nbp++;
      }
    }
  }
  free(tmp_ct);
  if(ps->nbp != 0 && ps->nbp != msa_nbp) ESL_FAIL(eslEINVAL, errbuf, "validate_and_update_sspostscript_given_msa(), expected %d basepairs in MSA's SS_cons, but read %d\n", ps->nbp, msa_nbp);

  /* for(rfpos = 0; rfpos < rflen; rfpos++) printf("ct %5d %5d\n", rfpos, msa_ct[rfpos+1]-1); */

  if(ps->msa_ct != NULL) { free(ps->msa_ct); ps->msa_ct = NULL;}
  if(ps->msa_rf2a_map != NULL) { free(ps->msa_rf2a_map); ps->msa_rf2a_map = NULL; }
  if(ps->msa_a2rf_map != NULL) { free(ps->msa_a2rf_map); ps->msa_a2rf_map = NULL; }
  ps->msa_ct = msa_ct;
  ps->msa_nbp = msa_nbp;
  ps->msa_ct = msa_ct;
  ps->msa_rf2a_map = rf2a_map;
  ps->msa_a2rf_map = a2rf_map;
  ps->msa_nseq = msa_nseq;

  return eslOK;

 ERROR:
  ESL_FAIL(status, errbuf, "validate_and_update_sspostscript_given_msa(), error status %d, probably out of memory.\n", status);
  return status; 
}


/* Function: draw_header_and_footer()
 * 
 * Purpose:  Draw header for a page
 */ 
static int 
draw_header_and_footer(FILE *fp, const ESL_GETOPTS *go, char *errbuf, SSPostscript_t *ps, int page, int pageidx2print)
{
  int status;
  int i, split_idx;
  float x, y;
  float header_fontsize;
  int model_width, desc_width, desc_column_width;
  char *model_dashes, *desc_dashes, *desc2print;
  char *desc_string = NULL;
  float xmodel;
  char *model2print = NULL;

  header_fontsize = HEADER_FONTSIZE_UNSCALED / ps->scale; 

  fprintf(fp, "%% begin ignore\n");
  fprintf(fp, "/%s findfont %.2f scalefont setfont\n", DEFAULT_FONT, header_fontsize);
  fprintf(fp, "0.00 0.00 0.00 1.00 setcmykcolor\n"); /* black */

  if(! (esl_opt_GetBoolean(go, "--no-head"))) { 
    model_width = ESL_MAX(strlen("model"), (int) strlen(ps->modelname));
    if(model_width > HEADER_MODELNAME_MAXCHARS) { 
      ESL_ALLOC(model2print, sizeof(char) * (HEADER_MODELNAME_MAXCHARS+1));
      for(i = 0; i < (HEADER_MODELNAME_MAXCHARS-3); i++) model2print[i] = ps->modelname[i];
      model2print[HEADER_MODELNAME_MAXCHARS-3] = '.';
      model2print[HEADER_MODELNAME_MAXCHARS-2] = '.';
      model2print[HEADER_MODELNAME_MAXCHARS-1] = '.';
      model2print[HEADER_MODELNAME_MAXCHARS] = '\0';
    }
    else { 
      if((status = esl_strdup(ps->modelname, -1, &(model2print))) != eslOK) ESL_FAIL(eslEINVAL, errbuf, "draw_header_and_footer(), error copying modelname");
    }
    
    model_width = ESL_MIN(model_width, HEADER_MODELNAME_MAXCHARS);
    ESL_ALLOC(model_dashes, sizeof(char) * (model_width+1));
    for(i = 0; i < model_width; i++) model_dashes[i] = '-'; 
    model_dashes[model_width] = '\0';
    
    if(ps->modeA[page] == ALIMODE || ps->modeA[page] == SIMPLEMASKMODE) { 
      if((status = esl_strdup("description", -1, &(desc_string))) != eslOK) ESL_FAIL(eslEINVAL, errbuf, "draw_header_and_footer(), error copying description");
    }
    else if(ps->modeA[page] == INDIMODE) { 
      if((status = esl_strdup("sequence name", -1, &(desc_string))) != eslOK) ESL_FAIL(eslEINVAL, errbuf, "draw_header_and_footer(), error copying description");
    }
    
    xmodel = ps->headerx_desc - (ps->headerx_charsize * (model_width  + 6 + 6 + 8 + 2)); /*6,6,8 for #res,#bps,#seq|seqlen) plus 2 spaces each, first 2 for 2 spaces before desc*/
    x = xmodel;
    y = ps->headery;
    
    fprintf(fp, "(%-*s  %4s  %4s) %.2f %.2f moveto show\n", model_width, "model", "#res", "#bps", x, y);
    y -= header_fontsize * 0.75;
    fprintf(fp, "(%-*s  %4s  %4s) %.2f %.2f moveto show\n", model_width, model_dashes, "----", "----", x, y);
    y -= header_fontsize * 0.75;
    fprintf(fp, "(%-*s  %4d  %4d) %.2f %.2f moveto show", model_width, model2print, ps->rflen, ps->msa_nbp, x, y);
    free(model_dashes);
    x += (model_width + 6 + 6 +2) * ps->headerx_charsize;
    
    if(ps->modeA[page] == ALIMODE) { 
      y+= header_fontsize * 1.5;
      fprintf(fp, "(%6s) %.2f %.2f moveto show\n", "#seqs", x, y);
      y -= header_fontsize * 0.75;
      fprintf(fp, "(%6s) %.2f %.2f moveto show\n", "------", x, y);
      y -= header_fontsize * 0.75;
      fprintf(fp, "(%6d) %.2f %.2f moveto show", ps->msa_nseq, x, y);
    }
    else if(ps->modeA[page] == INDIMODE && (ps->seqidxA[page] != -1)) { /* ps->seqidxA[page] == -1 if we're printing the consensus sequence */
      y+= header_fontsize * 1.5;
      fprintf(fp, "(%6s) %.2f %.2f moveto show\n", "seqlen", x, y);
      y -= header_fontsize * 0.75;
      fprintf(fp, "(%6s) %.2f %.2f moveto show\n", "------", x, y);
      y -= header_fontsize * 0.75;
      fprintf(fp, "(%6d) %.2f %.2f moveto show", ps->uaseqlenA[ps->seqidxA[page]], x, y);
    }
    
    if(ps->descA[page] != NULL) { 
      x =  ps->headerx_desc;
      y += 2. * header_fontsize * 0.75;
      desc_width = ESL_MAX((int) strlen(desc_string), (int) strlen(ps->descA[page]));
      if(desc_width > ps->desc_max_chars) { 
	/* split into two lines, the add_page_desc_to_sspostscript() function added a '\n' where the split should be */
	i = 0; 
	while(ps->descA[page][i] != '\n') { 
	  i++; 
	  if(i >= desc_width) ESL_FAIL(eslEINVAL, errbuf, "drawing header, failed to find split point from add_page_desc_to_() in two-line description (%s)", ps->descA[page]);
	}
	split_idx = i;
	desc_column_width = split_idx;
      }
      else desc_column_width = desc_width;
      
      ESL_ALLOC(desc_dashes, sizeof(char) * (desc_column_width+1));
      for(i = 0; i < desc_column_width; i++) desc_dashes[i] = '-'; 
      desc_dashes[desc_column_width] = '\0';
      
      fprintf(fp, "(%-*s) %.2f %.2f moveto show\n", desc_column_width, desc_string, x, y);
      y -= header_fontsize * 0.75;
      free(desc_string);
      fprintf(fp, "(%-*s) %.2f %.2f moveto show\n", desc_column_width, desc_dashes, x, y);
      y -= header_fontsize * 0.75;
      free(desc_dashes);
      
      if(desc_width > ps->desc_max_chars) {
	ESL_ALLOC(desc2print, sizeof(char) * (split_idx+1));
	for(i = 0; i < split_idx; i++) desc2print[i] = ps->descA[page][i];
	desc2print[split_idx] = '\0';
	fprintf(fp, "(%-*s) %.2f %.2f moveto show\n", desc_column_width, desc2print, x, y);
	free(desc2print);
	
	x = ps->headerx_desc;
	y-= ps->headery_charsize * 1;
	ESL_ALLOC(desc2print, sizeof(char) * ((desc_width - split_idx) -1+1));
	for(i = split_idx+1; i < desc_width; i++) desc2print[(i-(split_idx+1))] = ps->descA[page][i];
	desc2print[(desc_width-(split_idx+1))] = '\0';
	fprintf(fp, "(%-*s) %.2f %.2f moveto show\n", desc_column_width, desc2print, x, y);
	free(desc2print);
      }
      else { 
	fprintf(fp, "(%-*s) %.2f %.2f moveto show\n", desc_width, ps->descA[page], x, y);
      }
    }
    /* masked row of header goes here if desired */
  }

  /* draw footer */
  float footer_fontsize, footerx_charsize;
  footer_fontsize = LEG_FONTSIZE_UNSCALED / ps->scale;
  footerx_charsize = ps->legx_charsize;
  
  fprintf(fp, "/%s findfont %.2f scalefont setfont\n", DEFAULT_FONT, footer_fontsize);
  if(! (esl_opt_GetBoolean(go, "--no-foot"))) { 
    /* draw alignment file name in lower left hand corner */
    if(ps->mask != NULL) { 
      if(esl_opt_GetString(go, "--mask-diff") != NULL) { 
	fprintf(fp, "(alifile: %s; mask 1 file: %s; mask 2 file: %s;) %.2f %.2f moveto show\n", esl_opt_GetArg(go, 1), esl_opt_GetString(go, "--mask"), esl_opt_GetString(go, "--mask-diff"), PAGE_SIDEBUF, PAGE_BOTBUF + (1.25 * footer_fontsize));
      }
      else { 
	fprintf(fp, "(alifile: %s; mask file: %s;) %.2f %.2f moveto show\n", esl_opt_GetArg(go, 1), esl_opt_GetString(go, "--mask"), PAGE_SIDEBUF, PAGE_BOTBUF + (1.25 * footer_fontsize));
      }
    }
    else { 
      fprintf(fp, "(alifile: %s) %.2f %.2f moveto show\n", esl_opt_GetArg(go, 1), PAGE_SIDEBUF, PAGE_BOTBUF + (1.25 * footer_fontsize));
    }

    /* put page number */
    /* determine ndigits */
    int tmp, ndigits;
    tmp = pageidx2print;
    ndigits = 1;
    while(tmp >= 10) { tmp /= 10; ndigits++; }
    x = ps->pagex_max - (PAGE_SIDEBUF) - (footerx_charsize * (5 + ndigits)); 
    fprintf(fp, "(page %d) %.2f %.2f moveto show\n", pageidx2print, x, PAGE_BOTBUF);
    
  }    
  fprintf(fp, "(structure diagram derived from CRW database: http://www.rna.ccbb.utexas.edu/) %.2f %.2f moveto show\n", PAGE_SIDEBUF , PAGE_BOTBUF);
  fprintf(fp, "%% end ignore\n");

  if(model2print != NULL) free(model2print);
  return eslOK;

 ERROR: ESL_FAIL(eslEINVAL, errbuf, "draw_header_and_footer(), memory error.");
}


/* Function: read_seq_list_file_bigmem
 * Date:     EPN, Thu Jun  5 13:21:36 2008
 * 
 * Read a file listing sequence names to draw individual
 * structure diagrams for, and find those names in an 
 * input MSA. This function is only called if --small
 * is not enabled, b/c if it is, then the msa does not
 * have any sequence info in it.
 *
 * <ret_useme> is an array specifying which sequences
 * were listed in the file. [0..i..msa->nseq-1]. It is
 * allocated here.
 * 
 * ret_useme[i] == TRUE if seq i was listed in the file,
 *                 FALSE otherwise.
 * <ret_nused> is filled with number of unique
 *             sequences read in the file that exist
 *             in the msa.
 * 
 * Returns eslOK on success.
 *
 * Dies with an error message if a sequence listed
 * in the file does not exist in the msa, or
 * some other problem is encountered.
 */
static int
read_seq_list_file_bigmem(char *filename, ESL_MSA *msa, int **ret_useme, int *ret_nused)
{
  int             status;
  ESL_FILEPARSER *efp;
  char           *seqname;
  int            *useme = NULL;
  int             nused = 0;
  int             seqidx;

  ESL_ALLOC(useme, sizeof(int) * msa->nseq);
  esl_vec_ISet(useme, msa->nseq, FALSE);

  if (esl_fileparser_Open(filename, NULL, &efp) != eslOK) esl_fatal("Error: failed to open list file %s\n", filename);
  
  while((status = esl_fileparser_GetToken(efp, &seqname, NULL)) != eslEOF) {
    status = esl_key_Lookup(msa->index, seqname, &seqidx);
    if(status == eslENOTFOUND) esl_fatal("Error while reading list file %s, sequence %s does not exist in the alignment.", filename, seqname);
    if(useme[seqidx] == FALSE) { 
      useme[seqidx] = TRUE;
      nused++;
    }
  }
  esl_fileparser_Close(efp);

  *ret_useme = useme;
  *ret_nused = nused;
  return eslOK;

 ERROR:
  if(useme != NULL) free(useme);
  esl_fatal("Memory allocation error while reading list file %s.", filename);
  return status; /* NEVERREACHED */
}

/* Function: read_seq_list_file_smallmem
 * Date:     EPN, Tue Jan 19 15:44:08 2010
 * 
 * Read a file listing sequence names to draw individual
 * structure diagrams for, and create a keyhash with
 * only those names in it. This function is only called if 
 * --small is enabled.
 *
 * Returns eslOK on success.
 *
 * Dies with an error message if a sequence listed
 * in the file does not exist in the msa, or
 * some other problem is encountered.
 */
static int
read_seq_list_file_smallmem(char *filename, ESL_KEYHASH **ret_useme_keyhash, int *ret_nused)
{
  int             status;
  ESL_FILEPARSER *efp;
  char           *seqname;
  int             nused = 0;
  ESL_KEYHASH    *useme_keyhash;

  useme_keyhash = esl_keyhash_Create();
  if(useme_keyhash == NULL) esl_fatal("Memory allocation error.");
  
  if (esl_fileparser_Open(filename, NULL, &efp) != eslOK) esl_fatal("Error: failed to open list file %s\n", filename);
  
  while((status = esl_fileparser_GetToken(efp, &seqname, NULL)) != eslEOF) {
    status = esl_key_Store(useme_keyhash, seqname, NULL);
    if(status == eslOK) nused++;
    else if(status != eslEDUP) esl_fatal("Error adding sequence %s to keyhash", seqname);
  }
  esl_fileparser_Close(efp);

  *ret_useme_keyhash = useme_keyhash;
  *ret_nused = nused;
  return eslOK;

}

/* Function: get_insert_info_from_msa
 * Date:     EPN, Fri Dec  4 13:52:53 2009
 * 
 * Read an MSA with #=GC RF annotation defining
 * consensus columns and count how many insertions
 * occur after each consensus column for each sequence.
 * 
 * msa         - the alignment
 * rflen       - expected nongap RF length (consensus length)
 * ret_nseq_with_ins_ct - [0..rflen]  number of sequences with >= 1 
 *               insert after each position. NULL if unwanted.
 * ret_ins_ct - [0..i..msa->nseq-1][0..rflen] number of inserts 
 *              insert after each position per sequence. NULL if unwanted.
 * 
 * Returns void. Dies with an informative error message upon an error.
 */
void
get_insert_info_from_msa(ESL_MSA *msa, int rflen, int **ret_nseq_with_ins_ct, int ***ret_ins_ct)
{
  int             status;
  int             i;
  int            *nseq_with_ins_ct = NULL;
  int           **ins_ct = NULL;
  int             rfpos, apos;

  /* contract check */
  if(msa->rf == NULL) esl_fatal("Error in get_insert_info_from_msa(), msa->rf is NULL.");

  /* allocate and initialize */
  ESL_ALLOC(nseq_with_ins_ct, sizeof(int) * (rflen+1));
  esl_vec_ISet(nseq_with_ins_ct, rflen+1, 0);
  ESL_ALLOC(ins_ct, sizeof(int *) * (msa->nseq));
  for(i = 0; i < msa->nseq; i++) { 
    ESL_ALLOC(ins_ct[i],  sizeof(int) * (rflen+1));
    esl_vec_ISet(ins_ct[i], (rflen+1), 0);
  }

  /* fill ins_ct, nseq_with_ins_ct */
  rfpos = 0;
  for(apos = 0; apos < msa->alen; apos++) { 
    if((! esl_abc_CIsGap(msa->abc, msa->rf[apos])) && 
       (! esl_abc_CIsMissing(msa->abc, msa->rf[apos])) && 
       (! esl_abc_CIsNonresidue(msa->abc, msa->rf[apos]))) { 
      rfpos++;
      if(rfpos > rflen) esl_fatal("Error in get_insert_info_from_msa(), expected consensus length (%d) is incorrect."); 
    }
    else { 
      for(i = 0; i < msa->nseq; i++) { 
	if(! esl_abc_CIsGap(msa->abc, msa->aseq[i][apos])) { 
	  ins_ct[i][rfpos]++;
	  if(ins_ct[i][rfpos] == 1) nseq_with_ins_ct[rfpos]++;
	}	  
      }
    }
  }

  if(ret_nseq_with_ins_ct != NULL) *ret_nseq_with_ins_ct = nseq_with_ins_ct;
  else free(nseq_with_ins_ct);
  if(ret_ins_ct != NULL) *ret_ins_ct = ins_ct;
  else esl_Free2D((void **) ins_ct, msa->nseq);

  return;

 ERROR:
  esl_fatal("Error in get_insert_info_from_msa(), memory allocation error.");
  return; /* NEVERREACHED */
}


/* Function: get_insert_info_from_ifile
 * Date:     EPN, Fri Dec  4 14:49:12 2009
 * 
 * Read a file output from Infernal's cmalign 
 * when run with '--matchonly --ifile <f>' and
 * fill *ret_ict[0..rflen][0..msa->nseq-1] number
 * of inserts after each consensus position for 
 * each sequence.
 *
 * Format of an insert file (from the commented header of an ifile):
 * This file includes 2+<nseq> non-'#' pre-fixed lines per model used for alignment,
 * where <nseq> is the number of sequences in the target file.
 * The first non-'#' prefixed line per model includes 2 tokens, separated by a single space (' '):
 * The first token is the model name and the second is the consensus length of the model (<clen>).
 * The following <nseq> lines include (1+3*<n>) whitespace delimited tokens per line.
 * The format for these <nseq> lines is:
 *   <seqname> <c_1> <u_1> <i_1> <c_2> <u_2> <i_2> .... <c_x> <u_x> <i_x> .... <c_n> <u_n> <i_n>
 *   indicating <seqname> has >= 1 inserted residues after <n> different consensus positions,
 *   <c_x> is a consensus position and
 *   <u_x> is the *unaligned* position in <seqname> of the first inserted residue after <c_x>.
 *   <i_x> is the number of inserted residues after position <c_x> for <seqname>.
 * Lines for sequences with 0 inserted residues will include only <seqname>.
 * The final non-'#' prefixed line per model includes only '//', indicating the end of info for a model.
 * 
 * Example: 
 * trna 72
 * trna-1
 * trna-2 19 20 1
 * trna-7 33 35 1 64 67 1
 * 
 * ifile         - name of ifile
 * rflen         - expected nongap RF length (consensus length)
 * msa_nseq      - expected number of sequences
 * useme_keyhash - keyhash with names of sequences for which we will store insert info
 *                 if NULL, we store insert info for all seqs.
 * ret_nseq_with_ins_ct  - [0..rflen] number of sequences with >= 1 inserts
 *                         after each RF position.
 * ret_ins_ct  - [0..i..msa->nseq-1][0..rflen number of inserts
 *               after each position for each sequence.
 *               Filled here.
 * 
 * Returns void. Dies with an informative error message on an error.
 */
void
get_insert_info_from_ifile(char *ifile, int rflen, int msa_nseq, ESL_KEYHASH *useme_keyhash, int **ret_nseq_with_ins_ct, int ***ret_ins_ct)
{
  int             status;
  ESL_FILEPARSER *efp;
  char           *tok;
  int             nseq_read = 0;
  int             nseq_stored = 0;
  int           **ins_ct = NULL;
  int            *nseq_with_ins_ct = NULL;
  int             nins;
  int             i;
  int             rfpos;
  int             seen_model_name_line = FALSE;
  int             seen_end_of_model_line = FALSE;
  int             nseq2store;     /* number of seqs we'll store insert info on */

  if (esl_fileparser_Open(ifile, NULL, &efp) != eslOK) esl_fatal("Error: failed to open list file %s\n", ifile);
  esl_fileparser_SetCommentChar(efp, '#');

  /* determine how many sequences we'll be storing info for */
  nseq2store = (useme_keyhash == NULL) ? msa_nseq : esl_keyhash_GetNumber(useme_keyhash);

  /* allocate and initialize */
  ESL_ALLOC(nseq_with_ins_ct, sizeof(int) * (rflen+1));
  esl_vec_ISet(nseq_with_ins_ct, rflen+1, 0);

  if(ret_ins_ct != NULL) { 
    ESL_ALLOC(ins_ct, sizeof(int *) * (nseq2store));
    for(i = 0; i < nseq2store; i++) { 
      ESL_ALLOC(ins_ct[i],  sizeof(int) * (rflen+1));
      esl_vec_ISet(ins_ct[i], (rflen+1), 0);
    }
  }
  /* Read the file, verify that it contains the correct number of sequences and the
   * consensus length(s) listed in the file agrees with expected rflen. 
   * Special care is taken to allow concatenated ifiles, so we may see more than 
   * one // lines, but the total number of seqs should match what we expect. 
   */
  while (esl_fileparser_NextLine(efp) == eslOK) { 
    if (esl_fileparser_GetTokenOnLine(efp, &tok, NULL) != eslOK) { 
      if(seen_model_name_line) esl_fatal("Error reading insert file, failed to read seq name on line %d of file %s\n", efp->linenumber, ifile);
      else                     esl_fatal("Error reading insert file, failed to read model name on line %d of file %s\n", efp->linenumber, ifile);
    }
    if(! seen_model_name_line) { /* this should be a special line, 2 tokens: <cmname> <rflen>, verify that <rflen> is what we expect  */
      seen_model_name_line   = TRUE;
      seen_end_of_model_line = FALSE;
      if (esl_fileparser_GetTokenOnLine(efp, &tok, NULL) != eslOK) {
	esl_fatal("Error reading insert file, failed to read consensus length on line %d of file %s\n", efp->linenumber, ifile); 
      }
      if(rflen != atoi(tok)) {
	esl_fatal("Error reading insert file, read consensus length of %d on line %d of file %s, but expected length %d\n", atoi(tok), rflen, efp->linenumber, ifile); 	
      } 
    }     
    else if (strncmp(tok, "//", 2) == 0) { /* end of data for an ifile, but we may have concatenated them, so we keep going */
      seen_model_name_line   = FALSE;
      seen_end_of_model_line = TRUE;
    }
    else { /* should be a seq line with 3*n tokens, <rfpos> <nins> ...., n can be 0 */
      /* determine if we're using this sequence */
      if(useme_keyhash == NULL || (esl_key_Lookup(useme_keyhash, tok, NULL) == eslOK)) { 
	while(esl_fileparser_GetTokenOnLine(efp, &tok, NULL) == eslOK) { 
	  rfpos = atoi(tok);
	  if(rfpos > rflen) { 
	    esl_fatal("Error reading insert file, read insert info for position %d that exceeds expected consensus length %don line %d of file %s.\n", rfpos, rflen, efp->linenumber, ifile);
	  }
	  if((status = esl_fileparser_GetTokenOnLine(efp, &tok, NULL)) != eslOK) { 
	    esl_fatal("Error reading insert file, didn't read unaligned sequence position for rfpos %d on line %d of file %s.\n", rfpos, efp->linenumber, ifile);
	  }
	  /* don't bother parsing this token */
	  
	  if((status = esl_fileparser_GetTokenOnLine(efp, &tok, NULL)) != eslOK) { 
	    esl_fatal("Error reading insert file, didn't read number of inserts for position %d on line %d of file %s.\n", rfpos, efp->linenumber, ifile);
	  }
	  nins = atoi(tok);
	  if(ins_ct != NULL) ins_ct[nseq_stored][rfpos] = nins;
	  if(nins > 0) nseq_with_ins_ct[rfpos]++; /* nins should always be > 0, but why not check?  */
	}
	nseq_stored++;
      }
      nseq_read++;
      if(nseq_read > msa_nseq) esl_fatal("Error reading insert file, read info for more sequences than expected (%d) at line %d of file %s.", msa_nseq, efp->linenumber, ifile);
    }
  }
  esl_fileparser_Close(efp);

  /* end of file, make sure we read a '//' at the end of it */
  if(! seen_end_of_model_line) esl_fatal("Error reading insert file, didn't read the special '//' line at the end of file %s.\n", rfpos, efp->linenumber, ifile);

  /* if useme_keyhash != NULL, make sure we read all the seqs we wanted to */
  if((useme_keyhash != NULL) && (nseq_stored != nseq2store)) { 
    esl_fatal("Error reading insert file, wanted to read insert info on %d seqs, but only found %d of them in the insert file %s\n", nseq2store, nseq_stored, ifile);      
  }
  if(nseq_read != msa_nseq)  esl_fatal("Error reading insert file, expected to read info on %d seqs, but only found %d in the insert file %s\n", msa_nseq, nseq_read, ifile);      

  if(ret_nseq_with_ins_ct != NULL) *ret_nseq_with_ins_ct = nseq_with_ins_ct;
  else free(nseq_with_ins_ct);
  if(ret_ins_ct != NULL) *ret_ins_ct = ins_ct;
  else if (ins_ct != NULL) esl_Free2D((void **) ins_ct, msa_nseq);

  return;

 ERROR:
  esl_fatal("Memory allocation error while reading insert file %s.", ifile);
  return;; /* NEVERREACHED */
}

/* Function: get_insert_info_from_abc_ct
 * Date:     EPN, Tue Jan 19 09:32:30 2010
 * 
 * Given an abc_ct array:
 * [0..apos..alen-1][0..abc->K]: per position count of 
 * each symbol in alphabet over all seqs. 
 * 
 * Determine the number of sequences with >= 1 inserted
 * residues after each RF position.
 *  
 * IMPORTANT NOTE: This is done based on the assumption that inserts
 * have been systematically placed in the MSA by a program like
 * Infernal or HMMER and that the maximum number of nongap residues in
 * any insert position after rfpos is the number of sequences with >=1
 * insert residues after rfpos. However, if the alignments were not
 * generated from HMMER nor Infernal, or have been modified, the
 * inserts may not follow this rule (the rule, in other words is that
 * there is always 1 insert column after each nongap RF position rfpos
 * that includes a residue from all sequences that have >= 1 inserted
 * residues after rfpos).
 * 
 * abc_ct      - [0..apos..msa_alen-1][0..abc->K] count of each residue in each position 
 * abc         - the alphabet
 * msa_rf      - msa->rf, the reference annotation 
 * msa->alen   - length of alignment
 * rflen       - expected nongap RF length (consensus length)
 * ret_nseq_with_ins_ct - [0..rflen]  number of sequences with >= 1 
 *               insert after each position. NULL if unwanted.
 * 
 * Returns void. Dies with an informative error message upon an error.
 */
void
get_insert_info_from_abc_ct(double **abc_ct, ESL_ALPHABET *abc, char *msa_rf, int64_t msa_alen, int rflen, int **ret_nseq_with_ins_ct)
{
  int             status;
  int            *nseq_with_ins_ct;
  int             nmaxins;
  int             rfpos, apos;

  /* allocate and initialize */
  ESL_ALLOC(nseq_with_ins_ct, sizeof(int) * (rflen+1));
  esl_vec_ISet(nseq_with_ins_ct, rflen+1, 0);

  nmaxins = 0;
  rfpos = 0;
  for(apos = 0; apos < msa_alen; apos++) { 
    if((! esl_abc_CIsGap(abc, msa_rf[apos])) && 
       (! esl_abc_CIsMissing(abc, msa_rf[apos])) && 
       (! esl_abc_CIsNonresidue(abc, msa_rf[apos]))) { 
      nseq_with_ins_ct[rfpos] = nmaxins;
      nmaxins = 0;
      rfpos++;
      if(rfpos > rflen) esl_fatal("Error in get_insert_info_from_abc_ct(), expected consensus length (%d) is incorrect."); 
    }
    else { 
      nmaxins = ESL_MAX(nmaxins, (int) esl_vec_DSum(abc_ct[apos], abc->K)); 
    }
  }
  /* get max_nseq with inserts after the final position */
  nseq_with_ins_ct[rfpos] = nmaxins;

  if(ret_nseq_with_ins_ct != NULL) *ret_nseq_with_ins_ct = nseq_with_ins_ct;
  else free(nseq_with_ins_ct);

  return;

 ERROR:
  esl_fatal("Error in get_insert_info_from_abc_ct(), memory allocation error.");
  return; /* NEVERREACHED */
}

/* get_pp_idx
 *                   
 * Given a #=GR PP or #=GC PP_cons character, return the appropriate index
 * in a pp_ct[] vector. 
 * '0' return 0;
 * '1' return 1;
 * '2' return 2;
 * '3' return 3;
 * '4' return 4;
 * '5' return 5;
 * '6' return 6;
 * '7' return 7;
 * '8' return 8;
 * '9' return 9;
 * '*' return 10;
 * gap return 11;
 * 
 * Anything else (including missing or nonresidue) return -1;
 */
int
get_pp_idx(ESL_ALPHABET *abc, char ppchar)
{
  if(esl_abc_CIsGap(abc, ppchar)) return 11;
  if(ppchar == '*')               return 10;
  if(ppchar == '9')               return 9;
  if(ppchar == '8')               return 8;
  if(ppchar == '7')               return 7;
  if(ppchar == '6')               return 6;
  if(ppchar == '5')               return 5;
  if(ppchar == '4')               return 4;
  if(ppchar == '3')               return 3;
  if(ppchar == '2')               return 2;
  if(ppchar == '1')               return 1;
  if(ppchar == '0')               return 0;
  return -1;
}

/* spos_and_epos2span_ct
 *                   
 * Given two arrays, srfpos_ct[0..rfpos..rflen-1] and erfpos_ct[0..rfpos..rflen-1], 
 * specifying the number of sequences for which the first and last non-gap RF (consensus)
 * position is rfpos, calculate the number of sequences that 'span' each RF
 * position. A sequence spans position rfpos if it has at least one consensus residue (non-gap)
 * in a RF position x <= rfpos, and at least one residue in a position y, y >= rfpos.
 * 
 * span_ct[0..rfpos..rflen-1] = x, x sequences 'span' position rfpos.
 */
int
get_span_ct(int rflen, int nseq, int *srfpos_ct, int *erfpos_ct, int **ret_span_ct)
{
  int status;
  int *nseq_start_after_rfpos = NULL;
  int *nseq_end_before_rfpos = NULL;
  int *span_ct = NULL;
  int rfpos;

  ESL_ALLOC(span_ct, sizeof(int) * rflen);
  ESL_ALLOC(nseq_start_after_rfpos, sizeof(int) * rflen);
  ESL_ALLOC(nseq_end_before_rfpos, sizeof(int) * rflen);

  esl_vec_ISet(nseq_start_after_rfpos, rflen, 0);
  esl_vec_ISet(nseq_end_before_rfpos, rflen, 0);

  /* first count number of seqs that start after each position */
  nseq_start_after_rfpos[rflen-1] = 0; /* initialize */
  for(rfpos = rflen-2; rfpos >= 0; rfpos--) { 
    nseq_start_after_rfpos[rfpos] = srfpos_ct[rfpos+1] + nseq_start_after_rfpos[rfpos+1];
  }
  /* count number of seqs that end before each position */
  nseq_end_before_rfpos[0] = 0; /* initialize */
  for(rfpos = 1; rfpos < rflen; rfpos++) { 
    nseq_end_before_rfpos[rfpos] += erfpos_ct[rfpos-1] + nseq_end_before_rfpos[rfpos-1];
  }

  /* we now know how many seqs start after each rfpos (a = nseq_start_after_rfpos[rfpos]) and
   *             how many seqs end  before each rfpos (b = nseq_end_before_rfpos[rfpos])
   * so c = [nseq-(a+b)] seqs span rfpos (b/c they don't start after it AND don't end before it).
   * 
   * Note: this is really confusing to me, but I've convinced myself it's true, and I think
   *       it is only true because we know each sequence's end position must be >= its start position
   */
  for(rfpos = 0; rfpos < rflen; rfpos++) { 
    span_ct[rfpos] = nseq - (nseq_start_after_rfpos[rfpos] + nseq_end_before_rfpos[rfpos]);
  }

  *ret_span_ct = span_ct;
  free(nseq_start_after_rfpos);
  free(nseq_end_before_rfpos);
  return eslOK;

 ERROR: 
  return eslEMEM;
}
