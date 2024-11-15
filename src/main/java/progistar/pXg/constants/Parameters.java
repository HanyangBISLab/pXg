package progistar.pXg.constants;

public class Parameters {
	// Number of SMA/BAM files
	public static int NUM_OF_SAM_FILES									=	1;
	public static int CURRENT_FILE_INDEX								=	0;


	// Input file paths
	public static final String CMD_GENOMIC_ANNOTATION_PATH				=	"--gtf_file".toLowerCase();
	public static String genomicAnnotationFilePath	=	null;

	public static final String CMD_GENOMIC_SEQUENCE_PATH				=	"--sam_file".toLowerCase();
	public static String[] sequenceFilePaths		=	null;

	public static final String CMD_PEPTIDE_ANNOTATION_PATH				=	"--psm_file".toLowerCase();
	public static String peptideFilePath			=	null;

	public static final String CMD_PROTEIN_SEQUENCE_PATH		=	"--fasta_file".toLowerCase();
	public static String proteinFastaPath			=	null;

	public static final String SEP_TYPE							=	"--sep".toLowerCase();
	public static String sepType					=	"tsv".toLowerCase();
	// Output file path
	public static final String CMD_OUTPUT_PATH				=	"--output".toLowerCase();
	public static String outputFilePath					=	null;
	public static String pinFilePath					=	null;
	public static String[] unmappedFilePaths				=	null;
	public static String[] exportSAMPaths					=	null;
	public static String[] tmpOutputFilePaths				=	null;


	// Peptide length
	public static final String CMD_LENGTH			=	"--lengths".toLowerCase();
	public static final String CMD_MAX_FLANK_SIZE	=	"--max_flank_size".toLowerCase();
	public static int minPeptLen					=	8;
	public static int maxPeptLen					=	15;
	public static int maxFlankNSize					=	1000;

	public static final String CMD_IL				=	"--ileq".toLowerCase();
	public static boolean leucineIsIsoleucine		=	true;

	// Those options cannot be accessed by users.
	// It is only used for test which method is more adaptable.
	public static byte	mocks						=	Constants.MOCK_REVERSE;
	public static byte	mockPolicy				=	Constants.MOCK_MAX_ONE;

	// Output format
	public static final String CMD_SAM_FORMAT		=	"--out_sam".toLowerCase();
	public static final String CMD_GTF_FORMAT		=	"--out_gtf".toLowerCase();
	public static final String CMD_NONCANONICAL		=	"--out_noncanonical".toLowerCase();
	public static final String CMD_CANONICAL		=	"--out_canonical".toLowerCase();
	public static final String CMD_UNMAPPED			=	"--out_unmapped".toLowerCase();

	// GTF partition size
	public static final String CMD_GENOMIC_ANNOTATION_PARTITION_SIZE		=	"--gtf_partition_size".toLowerCase();
	public static int partitionSize					=	10000000; // 10^7 * 10 * 8 = 0.8 G

	public static final String CMD_GENOMIC_SEQUENCE_PARTITION_SIZE			=	"--sam_partition_size".toLowerCase();
	public static int readSize						=	1000000; // 1 * 10^6

	// maxJunctionSize MUST not exceed partitionSize in GTF
	public static int maxJunctionSize				=	1000000; // 1 * 10^6

	// Peptide file
	// for user-friendly purpose, peptideColumnIndex is taken one-based and converted to zero-based.
	public static final String CMD_PEPTIDE_COLUMN_INDEX	=	"--pept_col".toLowerCase();
	public static int peptideColumnIndex			=	-1; // user-specific peptide index

	public static final String CMD_SCAN_COLUMN_INDEX	=	"--scan_col".toLowerCase();
	public static int scanColumnIndex			=	-1; // user-specific scan index

	public static final String CMD_FILE_COLUMN_INDEX	=	"--file_col".toLowerCase();
	public static int fileColumnIndex						=	-1;

	public static final String CMD_SCORE_COLUMN_INDEX	=	"--score_col".toLowerCase();
	public static int scoreColumnIndex				=	-1;

	public static final String CMD_CHARGE_COLUMN_INDEX	=	"--charge_col".toLowerCase();
	public static int chargeColumnIndex				=	-1;

	public static final String CMD_ADD_FEAT_COLUMN_INDICES	=	"--add_feat_cols".toLowerCase();
	public static int[] additionalFeatureIndices				=	null;

	public static final String CMD_CANDIDATE_RANK	=	"--rank".toLowerCase();
	public static int psmRank						=	100;

	public static final String pParserRegExr		=	"--aareg".toLowerCase();
	public static String peptideParserRegExr		=	"[A-Z]"; // read sequence matched to the RegExr.
	// Note that
	// comment is different from field.
	// In pXg definition, field indicates column names and comment shows some meta-data.
	// Comment line must be present on the top of records.
	// Field must follow "comment" when the comment presents in the file.
	public static final String cMarker				=	"cm".toLowerCase();
	public static String commentMarker				=	"#|@|%"; // if line starts with the pattern, the line will be skipped during parsing the file.

	public static int maxProteinOut					=	10;

	// Penalty
	public static final String CMD_PENALTY_MUTATION	=	"--penalty_mutation".toLowerCase();
	public static final String CMD_PENALTY_AS		=	"--penalty_AS".toLowerCase();
	public static final String CMD_PENALTY_5UTR		=	"--penalty_5UTR".toLowerCase();
	public static final String CMD_PENALTY_3UTR		=	"--penalty_3UTR".toLowerCase();
	public static final String CMD_PENALTY_FS		=	"--penalty_FS".toLowerCase();
	public static final String CMD_PENALTY_ncRNA	=	"--penalty_ncRNA".toLowerCase();
	public static final String CMD_PENALTY_IR		=	"--penalty_IR".toLowerCase();
	public static final String CMD_PENALTY_IGR		=	"--penalty_IGR".toLowerCase();
	public static final String CMD_PENALTY_asRNA	=	"--penalty_asRNA".toLowerCase();
	public static final String CMD_PENALTY_SOFTCLIP	=	"--penalty_softclip".toLowerCase();
	public static final String CMD_PENALTY_UNMAP	=	"--penalty_unknown".toLowerCase();

	public static double PENALTY_MUTATION			=	1;
	public static double PENALTY_AS					=	10;
	public static double PENALTY_5UTR				=	20;
	public static double PENALTY_3UTR				=	20;
	public static double PENALTY_FS					=	20;
	public static double PENALTY_ncRNA				=	20;
	public static double PENALTY_IR					=	30;
	public static double PENALTY_IGR				=	30;
	public static double PENALTY_asRNA				=	30;
	public static double PENALTY_SOFTCLIP			=	50;
	public static double PENALTY_UNMAP				=	100;

	// Export option
	public static boolean EXPORT_TSV				=	true;
	public static boolean EXPORT_UNMAPPED_SEQ		=	true;
	public static boolean EXPORT_SAM				=	false;
	public static boolean EXPORT_NONCANONICAL		=	true;
	public static boolean EXPORT_CANONICAL			=	true;

	// System Parameters
	public static final String CMD_THREADS			=	"--threads".toLowerCase();
	public static int nThreads						=	4;

	public static boolean isDecoyOut 				=	true;

	public static final String CMD_TRANSLATION		=	"--mode";

	// Three or Six frame translation
	public static int translationMethod				=	Constants.THREE_FRAME;

	public static final String CMD_REMOVE_QUOTES	=	"--rm_quotes";
	public static boolean rmQuotes					=	true;

	// Third-party path
	public static String netMHCpanPath				=	"/Users/gistar/tools/netMHCpan4.1/netMHCpan-4.1/netMHCpan";



	// hidden parameters for revision
	public static String CMD_PHRED_CAL					=	"--cal_phred";
	public static String PHRED_CAL						=	Constants.CAL_PHRED_AVG;
	
	
	// Primary count
	public static boolean COUNT_PRIMARY_ONLY			=	false;
	public static String strandedness					=	Constants.AUTO_STRANDED;



	///////////////////// Parameters for "Build Sequence Database" /////////////////////
	public static String sequencedbInputPath			=	null;
	public static String sequencedbOutputPath			=	null;
	public static String referenceSequencePath			=	null;
	public static boolean isIncludedCanonical			=	false;
	public static boolean isIncludedNoncanonical		=	false;
	public static boolean isIncludedFlankSequence		=	false;
	public static boolean isStringent					=	false;
	////////////////////////////////////////////////////////////////////////////////////

}
