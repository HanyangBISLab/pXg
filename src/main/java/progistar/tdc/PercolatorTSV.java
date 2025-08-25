package progistar.tdc;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;

class PercolatorRecord {
	public String psmId = null;
	public String score = null;
	public String peptide = null;
	public String proteinIds = null;
	public boolean isDecoy = false;
	
	public PercolatorRecord (String psmId, String score, String peptide, String proteinIds) {
		this.psmId = psmId;
		this.score = score;
		this.peptide = peptide;
		this.proteinIds = proteinIds;
	}
	
	public String getKey () {
		return this.psmId + "@" + proteinIds + "@" + peptide;
	}
}

public class PercolatorTSV {

	public static final String PROTEIN_IDS	= "proteinids";
	public static final String PSM_ID		= "psmid";
	public static final String PEPTIDE		= "peptide";
	
	public ArrayList<PercolatorRecord> records = new ArrayList<PercolatorRecord>();
	public Hashtable<String, PercolatorRecord> keyToRecordMapper = new Hashtable<String, PercolatorRecord>();
	
	public PercolatorTSV (File file) throws IOException {
		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;
		
		// fields
		String[] fields = BR.readLine().split("\t");
		int proteinIdsIdx = -1;
		int psmIdIdx = -1;
		int scoreIdx = -1;
		int peptideIdx = -1;
		
		for(int i=0; i<fields.length; i++) {
			if(fields[i].equalsIgnoreCase(PROTEIN_IDS)) {
				proteinIdsIdx = i;
			} else if(fields[i].equalsIgnoreCase(PSM_ID)) {
				psmIdIdx = i;
			} else if(fields[i].equalsIgnoreCase(TDC.scoreName)) {
				scoreIdx = i;
			} else if(fields[i].equalsIgnoreCase(PEPTIDE)) {
				peptideIdx = i;
			}
		}
		
		while((line = BR.readLine()) != null) {
			fields = line.split("\t");
			String psmId = fields[psmIdIdx];
			String score = fields[scoreIdx];
			String peptide = fields[peptideIdx];
			String proteinIds = fields[proteinIdsIdx];
			
			PercolatorRecord pr = new PercolatorRecord(psmId, score, peptide, proteinIds);
			records.add(pr);
		}
		
		BR.close();
		
		System.out.println("Read "+records.size()+" from "+file.getName());
		
		// check and remove decoy prefix.
		checkAndRemoveDecoyPrefix();
		
		// generate mapper
		for(int i=0; i<records.size(); i++) {
			PercolatorRecord pr = records.get(i);
			keyToRecordMapper.put(pr.getKey(), pr);
		}
	}
	
	
	private void checkAndRemoveDecoyPrefix () {
		System.out.println("Check decoy prefix...");
		int endIdx = -1;
		boolean endOfPattern = false;
		
		while(!endOfPattern) {
			endIdx++;
			
			char previousChar = records.get(0).proteinIds.charAt(endIdx);
			for(int i=0; i<records.size(); i++) {
				PercolatorRecord pr = records.get(i);
				char curChar = pr.proteinIds.charAt(endIdx);
				
				if(previousChar != curChar) {
					endOfPattern = true;
					break;
				}
			}
		}
		
		if(endIdx > 0) {
			String decoyPrefix = records.get(0).proteinIds.substring(0, endIdx);
			System.out.println("Decoy prefix: "+decoyPrefix);
			
			for(int i=0; i<records.size(); i++) {
				PercolatorRecord pr = records.get(i);
				pr.isDecoy = true;
				pr.proteinIds = pr.proteinIds.replace(decoyPrefix, "");
			}
			
		} else {
			System.out.println("No decoy prefix found. It should be target PSMs.");
		}
	}
	
}
