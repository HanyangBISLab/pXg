package progistar.pXg.constants;

public class Constants {

	//
	public static final byte EXON = 30;
	
	// EXON is classified into three types such as CDS, UTR and NCDS (Non CoDing Sequence)
	public static final byte CDS = 100;
	public static final byte UTR5 = 25;
	public static final byte UTR3 = 23;
	public static final byte NCDS = 20;
	public static final byte INTRON = 0;
	public static final byte INTERGENIC = -4;
	
	// Regional Character
	public static final char MARK_CDS 			=	'C';
	public static final char MARK_5UTR			=	'F';
	public static final char MARK_3UTR			=	'T';
	public static final char MARK_NCDS			=	'N';
	public static final char MARK_INTRON		=	'I';
	public static final char MARK_INTERGENIC	=	'-';
	public static final char MARK_SOFTCLIP		=	'?';
	public static final char MARK_UNMAPPED		=	'*';
	
	// Alternative Splicing Character
	public static final char MARK_AS			=	'A'; // for alternative splicing form
	public static final char MARK_CA			=	'C'; // for canonical form
	
	// Transcript coding type
	public static final byte NON_CODING_TRANSCRIPT	=	0;
	public static final byte CODING_TRANSCRIPT		=	1;
	
	// Frame type
	public static final char NO_FRAME		=	'N';
	public static final char IN_FRAME		=	'I';
	public static final char OUT_OF_FRAME	=	'O';
	
	// Frame index
	// Note that FRAME_X denots NO_FRAME.
	public static final byte FRAME_0		=	0;
	public static final byte FARME_1		=	1;
	public static final byte FRAME_2		=	2;
	public static final byte FRAME_X		=	3;
	
	// Mutation type
	public static final byte SNP			=	0;
	public static final byte INS			=	1;
	public static final byte DEL			=	2;
	
	// RNA-Seq parameter
	public static final byte FORWARD_STRAND_READS	=	0;
	public static final byte REVERSE_STRAND_READS	=	1;
	
	// RNA-Seq Mock
	public static final byte MOCK_NONE			=	0;
	public static final byte MOCK_REVERSE		=	1;
	
	public static final byte MOCK_ALL			=	0; // use all multiply-assigned read count
	public static final byte MOCK_MAX_ONE		=	1; // use maximum one as a representative
	public static final byte MOCK_MEAN			=	2; // use mean as a representative
	
	// TASKS
	public static final int TASK_G_MAP				=	1;
	
	// Output Annotation
	public static final String OUTPUT_G_UNIQUE_ID	=	"[ID]";
	public static final String OUTPUT_G_SEQUENCE	=	"[SEQ]";
	public static final String OUTPUT_G_REGION		=	"[REGION]";
	public static final String OUTPUT_G_PEPTIDE		=	"[PEPTIDE]";
	
	// Events
	public static final String EVENT_ANTISENSE		=	"antisense";
	public static final String EVENT_SENSE			=	"sense";
	public static final String EVENT_5UTR			=	"5UTR";
	public static final String EVENT_3UTR			=	"3UTR";
	public static final String EVENT_NONCODING		=	"noncoding";
	public static final String EVENT_INTERGENIC		=	"intergenic";
	public static final String EVENT_INTRON			=	"intron";
	public static final String EVENT_FRAMESHIFT		=	"frameshift";
	public static final String EVENT_PROTEINCODING	=	"proteincoding";
	public static final String EVENT_AS				=	"alternativesplicing";
	public static final String EVENT_UNKNOWN		=	"unknown";
	
	// PSM Status
	public static final byte PSM_STATUS_RANDOM		=	0;
	public static final byte PSM_STATUS_DECOY		=	1;
	public static final byte PSM_STATUS_TARGET		=	2;

}
