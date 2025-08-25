package progistar.tdc;

import java.io.File;
import java.io.IOException;

public class PercolatorTSV extends TDTSV {

	public static String PROTEIN_IDS	= "proteinids";
	public static String PSM_ID			= "psmid";
	public static String SCORE			= "score";
	public static String PEPTIDE		= "peptide";
	
	public PercolatorTSV(File file)
			throws IOException {
		super(file, PSM_ID, SCORE, PEPTIDE, PROTEIN_IDS);
		// TODO Auto-generated constructor stub
	}

	
}
