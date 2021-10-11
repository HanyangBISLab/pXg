package progistar.pXg.tset;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class SamSelector {
	
	private enum FieldIndex {
		QNAME(0), CHR(2), START_POS(3), CIGAR(5), SEQUENCE(9);

		private int value;

		FieldIndex(int value) {
			this.value = value;
		}
	}
	
	
	public static void main(String[] args) throws IOException {
		writeSam(args[0], args[1]);
	}
	
	/**
	 * Select some lines of records from a sam file. <br>
	 * Use for testing purpose. <br>
	 * 
	 * @param samFilePath
	 * @param numOfRecords
	 */
	public static void writeSam (String samFilePath, String targetChr) {
		try {
			File samFile = new File(samFilePath);
			File outSamFile = new File(samFile.getAbsolutePath().replace(".sam", ".sub."+targetChr+".sam"));
			
			BufferedReader BR = new BufferedReader(new FileReader(samFile));
			BufferedWriter BW = new BufferedWriter(new FileWriter(outSamFile));
			String line = null;
			
			while((line = BR.readLine()) != null) {

				// meta data is not counted.
				if(line.startsWith("@")) {
					BW.append(line);
					BW.newLine();
				}  else  {
					String[] fields = line.split("\t");
					if(fields[FieldIndex.CHR.value].equalsIgnoreCase(targetChr)) {
						BW.append(line);
						BW.newLine();
					}
				}
				
			}
			
			BR.close();
			BW.close();
			
		}catch (IOException ioe) {
			
		}
	}

}
