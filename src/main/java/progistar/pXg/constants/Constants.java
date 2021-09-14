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
	public static final char MARK_UTR5			=	'F';
	public static final char MARK_UTR3			=	'T';
	public static final char MARK_NCDS			=	'N';
	public static final char MARK_INTRON		=	'I';
	public static final char MARK_INTERGENIC	=	'-';
	public static final char MARK_SOFTCLIP		=	'?';
	
	// Transcript coding type
	public static final byte NON_CODING_TRANSCRIPT	=	0;
	public static final byte CODING_TRANSCRIPT		=	1;
	
	// Frame annotation
	public static final byte IN_FRAME_0		=	0;
	public static final byte IN_FRAME_1		=	1;
	public static final byte IN_FRAME_2		=	2;
	
	// Frame type
	public static final byte IN_FRAME		=	0;
	public static final byte OUT_OF_FRAME	=	1;
	public static final byte NO_FRAME		=	2;
	
	// RNA-Seq parameter
	public static final byte FORWARD_STRAND_READS	=	0;
	public static final byte REVERSE_STRAND_READS	=	1;
	
	// RNA-Seq decoy
	public static final byte DECOY_NONE			=	0;
	public static final byte DECOY_REVERSE		=	1;
	public static final byte DECOY_SHUFFLE		=	2;
	
	// TASKS
	public static final int TASK_G_MAP				=	1;
}
