package progistar.pXg.constants;

public class Parameters {
	// Number of SMA/BAM files
	public static int NUM_OF_SAM_FILES									=	1;
	public static int CURRENT_FILE_INDEX								=	0;


	// Input file paths
	public static String genomicAnnotationFilePath	=	null;
	public static String[] sequenceFilePaths		=	null;
	public static String peptideFilePath			=	null;
	public static String proteinFastaPath			=	null;
	public static String sepType					=	"tsv".toLowerCase();
	// Output file path
	public static String outputFilePath					=	null;
	public static String pinFilePath					=	null;
	public static String[] unmappedFilePaths				=	null;
	public static String[] exportSAMPaths					=	null;
	public static String[] tmpOutputFilePaths				=	null;


	// Peptide length
	public static int minPeptLen					=	8;
	public static int maxPeptLen					=	15;
	public static int maxFlankNSize					=	1000;
	public static boolean leucineIsIsoleucine		=	true;

	// Those options cannot be accessed by users.
	// It is only used for test which method is more adaptable.
	public static byte	mocks						=	Constants.MOCK_REVERSE;
	public static byte	mockPolicy				=	Constants.MOCK_MAX_ONE;

	// GTF partition size
	public static int partitionSize					=	10000000; // 10^7 * 10 * 8 = 0.8 G
	public static int readSize						=	1000000; // 1 * 10^6
	// maxJunctionSize MUST not exceed partitionSize in GTF
	public static int maxJunctionSize				=	1000000; // 1 * 10^6

	// Peptide file
	// for user-friendly purpose, peptideColumnIndex is taken one-based and converted to zero-based.
	public static int peptideColumnIndex			=	-1; // user-specific peptide index
	public static int[] identifierColumnIndicies = null;
	public static int scoreColumnIndex				=	-1;
	public static int chargeColumnIndex				=	-1;
	public static int[] additionalFeatureIndices				=	null;
	public static int psmRank						=	100;
	public static String peptideParserRegExr		=	"[A-Z]"; // read sequence matched to the RegExr.
	// Note that
	// comment is different from field.
	// In pXg definition, field indicates column names and comment shows some meta-data.
	// Comment line must be present on the top of records.
	// Field must follow "comment" when the comment presents in the file.
	public static String commentMarker				=	"#|@|%"; // if line starts with the pattern, the line will be skipped during parsing the file.

	public static int maxProteinOut					=	10;

	// Penalty
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
	public static int nThreads						=	4;
	public static boolean isDecoyOut 				=	true;
	public static boolean rmQuotes					=	true;

	// Third-party path
	public static String netMHCpanPath				=	"/Users/gistar/tools/netMHCpan4.1/netMHCpan-4.1/netMHCpan";



	// hidden parameters for revision
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
