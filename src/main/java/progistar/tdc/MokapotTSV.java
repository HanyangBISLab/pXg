package progistar.tdc;

import java.io.File;
import java.io.IOException;

public class MokapotTSV extends TDTSV {

	public static String PROTEIN_IDS	= "protein_list";
	public static String PSM_ID			= "spectrum_id";
	public static String SCORE			= "mokapot score";
	public static String PEPTIDE		= "peptide";
	
	public MokapotTSV(File file) throws IOException {
		super(file, PSM_ID, SCORE, PEPTIDE, PROTEIN_IDS);
		// TODO Auto-generated constructor stub
	}
	
	@Override
	public String parseProteinId(String proteinId) {
		proteinId = proteinId.replace("['", "").replace("']", "");
		return super.parseProteinId(proteinId);
	}

}
