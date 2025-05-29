package progistar.pXg.data;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import progistar.pXg.constants.Constants;

public class AAVariantTable {

	
	private static int endPos = 0;
	public static String[] aaRNAArray = new String[100];
	public static String[] aaPeptideArray = new String[100];
	
	public static void init () {
		// init
		aaRNAArray[0] = Constants.ID_NULL;
		aaPeptideArray[0] = Constants.ID_NULL;
		
	}
	
	public static String getAnnotation (int idx) {
		if(idx < endPos) {
			return aaRNAArray[idx] + ">" +aaPeptideArray[idx];
		} else {
			return null;
		}
	}
	
	public static void addAAVariant (String aaRNA, String aaPeptide) {
		assert endPos < aaRNAArray.length;
		
		aaRNAArray[endPos] = aaRNA;
		aaPeptideArray[endPos] = aaPeptide;
		
		// add end mark
		aaRNAArray[endPos+1] = Constants.ID_NULL;
		aaPeptideArray[endPos+1] = Constants.ID_NULL;
		
		endPos++;
	}
	
	public static int getSize () {
		return endPos;
	}
	
	public static ArrayList<AAVariant> getAAVariants (String peptide) {
		ArrayList<AAVariant> aaVars = new ArrayList<AAVariant>();
		
		for(int i=0; i<endPos; i++) {
			String aaPeptide = aaPeptideArray[i];
			int idx = peptide.indexOf(aaPeptide);
			
			while(idx != -1) {
				AAVariant aaVar = new AAVariant();
				aaVar.aaVariantIndex = i;
				aaVar.pos = idx;
				aaVars.add(aaVar);
				
				idx = peptide.indexOf(aaPeptide, idx + 1);
			}
		}
		
		return aaVars;
	}
	
	public static void parseTable (String tableFilePath) {
		// if no file provided?
		if(tableFilePath == null) {
			return;
		}
		
		try {
			AAVariantTable.init();
			
			BufferedReader BR = new BufferedReader(new FileReader(tableFilePath));
			int aaRNAIdx = -1;
			int aaPeptideIdx = -1;
			
			String header = BR.readLine();
			String[] fields = header.split("\t");
			
			for(int i=0; i<fields.length; i++) {
				if(fields[i].equalsIgnoreCase(Constants.AA_RNA_COLUMN_NAME)) {
					aaRNAIdx = i;
				} else if(fields[i].equalsIgnoreCase(Constants.AA_PEPTIDE_COLUMN_NAME)) {
					aaPeptideIdx = i;
				}
			}
			
			if(aaRNAIdx == -1 ||aaPeptideIdx == -1) {
				System.out.println(tableFilePath +" has wrong column names!");
				System.out.println(header);
				System.out.println("Header must have aaRNA and aaPeptide columns!");
				System.exit(1);
			}
			
			
			String line = null;
			while((line = BR.readLine()) != null) {
				fields = line.split("\t");
				String aaRNA = fields[aaRNAIdx];
				String aaPeptide = fields[aaPeptideIdx];
				
				if(aaRNA.equalsIgnoreCase(aaPeptide)) {
					System.out.println("Skip synounymous amino acid variation: "+aaRNA+">"+aaPeptide);
				} else {
					// add aaVariant
					AAVariantTable.addAAVariant(aaRNA, aaPeptide);
				}
				
				
			}
			
			BR.close();
			
		}catch(IOException ioe) {
			
		}
		if(AAVariantTable.getSize() > 0) {
			System.out.println("Load an amino acid variant table: "+AAVariantTable.getSize() +" were stored.");
		}
	}
}
