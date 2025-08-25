package progistar.tdc;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Hashtable;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import progistar.pXg.constants.Constants;
import progistar.pXg.data.pXgRecord;
import progistar.pXg.data.parser.pXgParser;

public class TDC {
	
	public static final String PERCOLATOR	= "percolator";
	public static final String MOKAPOT		= "mokapot";

	public static File targetTSVFile	= null;
	public static File decoyTSVFile 	= null;
	public static File pXgFile			= null;
	public static double fdr			= 0.01;
	public static boolean isLengthSpecific = false;
	public static boolean isPeptideLevel   = false;
	public static File prefixOfOutputFile		= null;
	public static String software		= PERCOLATOR;
	
	public static void main(String[] args) throws IOException {
		System.out.println(Constants.VERSION+" "+Constants.RELEASE);
		System.out.println(Constants.INTRODUCE);
		System.out.println();
		
		parseOptions(args);
		
		TDTSV targetTSV = null;
		TDTSV decoyTSV = null;
		
		if(software.equalsIgnoreCase(PERCOLATOR)) {
			targetTSV = new PercolatorTSV(targetTSVFile);
			decoyTSV = new PercolatorTSV(decoyTSVFile);
		} else if(software.equalsIgnoreCase(MOKAPOT)) {
			targetTSV = new MokapotTSV(targetTSVFile);
			decoyTSV = new MokapotTSV(decoyTSVFile);
		}
		
		// remove redundant psm_id
		TDTSV.refineTargetDecoyPSMs(targetTSV, decoyTSV);
		
		ArrayList<pXgRecord> records = pXgParser.parse(pXgFile, false);
		
		// connect between pXg and percolator results
		ArrayList<pXgRecord> referencePSMs = new ArrayList<pXgRecord>();
		ArrayList<pXgRecord> nonreferencePSMs = new ArrayList<pXgRecord>();
		
		int numReferenceTargetPSMs = 0;
		int numReferenceDecoyPSMs = 0;
		int numNonreferenceTargetPSMs = 0;
		int numNonreferenceDecoyPSMs = 0;
		
		for(int i=0; i<records.size(); i++) {
			pXgRecord record = records.get(i);
			String specId = record.getValueByFieldName("specid");
			String genomicId = record.getValueByFieldName("genomicid");
			String peptide = record.getValueByFieldName("InferredPeptide");
			
			String key = specId + "@" +genomicId + "@" +peptide;
			
			TDRecord pr = targetTSV.keyToRecordMapper.get(key);
			if(pr == null) {
				pr = decoyTSV.keyToRecordMapper.get(key);
			}
			
			if( pr != null ) {
				record.setValueByFieldName("Rescore", pr.score);
				if(record.isReference()) {
					referencePSMs.add(record);
					
					if(record.isTarget()) {
						numReferenceTargetPSMs++;
					} else {
						numReferenceDecoyPSMs++;
					}
					
				} else {
					nonreferencePSMs.add(record);
					
					if(record.isTarget()) {
						numNonreferenceTargetPSMs++;
					} else {
						numNonreferenceDecoyPSMs++;
					}
				}
			}
		}
		
		System.out.println("A total of PSMs: "+records.size());
		System.out.println(" -Top PSMs: "+ (referencePSMs.size() + nonreferencePSMs.size()));
		System.out.println("  -Reference target/decoy PSMs: "+numReferenceTargetPSMs+"/"+numReferenceDecoyPSMs);
		System.out.println("  -Non-reference target/decoy PSMs: "+numNonreferenceTargetPSMs+"/"+numNonreferenceDecoyPSMs);
		
		// sort by percolator score
		Collections.sort(referencePSMs, new Comparator<pXgRecord>() {
			@Override
			public int compare(pXgRecord o1, pXgRecord o2) {
				double s1 = Double.parseDouble(o1.getValueByFieldName("Rescore"));
				double s2 = Double.parseDouble(o2.getValueByFieldName("Rescore"));
				
				if(s1 > s2) {
					return -1;
				} else if(s1 < s2) {
					return 1;
				} else if(o1.isTarget() && (!o2.isTarget())) {
					return 1;
				}  else if((!o1.isTarget()) && o2.isTarget()) {
					return -1;
				}
				
				return 0;
			}
		});
		
		Collections.sort(nonreferencePSMs, new Comparator<pXgRecord>() {
			@Override
			public int compare(pXgRecord o1, pXgRecord o2) {
				double s1 = Double.parseDouble(o1.getValueByFieldName("Rescore"));
				double s2 = Double.parseDouble(o2.getValueByFieldName("Rescore"));
				
				if(s1 > s2) {
					return -1;
				} else if(s1 < s2) {
					return 1;
				} else if(o1.isTarget() && (!o2.isTarget())) {
					return 1;
				}  else if((!o1.isTarget()) && o2.isTarget()) {
					return -1;
				}
				
				return 0;
			}
		});
		
		if(isPeptideLevel) {
			
			numReferenceTargetPSMs = 0;
			numReferenceDecoyPSMs = 0;
			numNonreferenceTargetPSMs = 0;
			numNonreferenceDecoyPSMs = 0;
			
			ArrayList<pXgRecord> peptides = referencePSMs;
			referencePSMs = new ArrayList<pXgRecord>();
			Hashtable<String, Boolean> checks  = new Hashtable<String, Boolean>();
			for(pXgRecord record : peptides) {
				String peptide = record.getValueByFieldName("InferredPeptide");
				if(checks.get(peptide) == null) {
					checks.put(peptide, true);
					referencePSMs.add(record);
					
					if(record.isTarget()) {
						numReferenceTargetPSMs++;
					} else {
						numReferenceDecoyPSMs++;
					}
				}
			}
			
			checks.clear();
			peptides = nonreferencePSMs;
			nonreferencePSMs = new ArrayList<pXgRecord>();
			
			for(pXgRecord record : peptides) {
				String peptide = record.getValueByFieldName("InferredPeptide");
				if(checks.get(peptide) == null) {
					checks.put(peptide, true);
					nonreferencePSMs.add(record);
					
					if(record.isTarget()) {
						numNonreferenceTargetPSMs++;
					} else {
						numNonreferenceDecoyPSMs++;
					}
				}
			}
			
			System.out.println("  -Reference target/decoy peptides: "+numReferenceTargetPSMs+"/"+numReferenceDecoyPSMs);
			System.out.println("  -Non-reference target/decoy peptides: "+numNonreferenceTargetPSMs+"/"+numNonreferenceDecoyPSMs);
			
		}
		
		// do FDR
		ArrayList<pXgRecord> passedPSMs = new ArrayList<pXgRecord>();
		if(isLengthSpecific) {
			ArrayList<pXgRecord> referencePSMsLenShort = new ArrayList<pXgRecord>();
			ArrayList<pXgRecord> referencePSMsLenLong = new ArrayList<pXgRecord>();
			
			ArrayList<pXgRecord> nonreferencePSMsLenShort = new ArrayList<pXgRecord>();
			ArrayList<pXgRecord> nonreferencePSMsLenLong = new ArrayList<pXgRecord>();
			
			for(int i=0; i<referencePSMs.size(); i++) {
				pXgRecord record = referencePSMs.get(i);
				int length = Integer.parseInt(record.getValueByFieldName("PeptideLength"));
				if(length <= 8) {
					referencePSMsLenShort.add(record);
				} else {
					referencePSMsLenLong.add(record);
				}
			}
			
			for(int i=0; i<nonreferencePSMs.size(); i++) {
				pXgRecord record = nonreferencePSMs.get(i);
				int length = Integer.parseInt(record.getValueByFieldName("PeptideLength"));
				if(length <= 8) {
					nonreferencePSMsLenShort.add(record);
				} else {
					nonreferencePSMsLenLong.add(record);
				}
			}
			
			referencePSMsLenShort = getMarkedPSMFilter(referencePSMsLenShort, "Reference target PSMs (short peptides)", fdr);
			referencePSMsLenLong = getMarkedPSMFilter(referencePSMsLenLong, "Reference target PSMs (long peptides)", fdr);
			nonreferencePSMsLenShort = getMarkedPSMFilter(nonreferencePSMsLenShort, "Non-reference target PSMs (short peptides)", fdr);
			nonreferencePSMsLenLong = getMarkedPSMFilter(nonreferencePSMsLenLong, "Non-reference target PSMs (long peptides)", fdr);
			
			// add all
			passedPSMs.addAll(referencePSMsLenShort);
			passedPSMs.addAll(referencePSMsLenLong);
			passedPSMs.addAll(nonreferencePSMsLenShort);
			passedPSMs.addAll(nonreferencePSMsLenLong);
			
			// report
			
			System.out.println("Passed PSMs: "+passedPSMs.size());
		} else {
			referencePSMs = getMarkedPSMFilter(referencePSMs, "Reference target PSMs", fdr);
			nonreferencePSMs = getMarkedPSMFilter(nonreferencePSMs, "Non-reference target PSMs", fdr);
			
			passedPSMs.addAll(referencePSMs);
			passedPSMs.addAll(nonreferencePSMs);
			
			System.out.println("Passed PSMs: "+passedPSMs.size());
		}
		
		// write output
		BufferedWriter BW = new BufferedWriter(new FileWriter(prefixOfOutputFile.getAbsolutePath()+".fdr.tsv"));
		
		BW.append(pXgParser.header.get(0));
		for(int i=1; i<pXgParser.header.size(); i++) {
			BW.append("\t").append(pXgParser.header.get(i));
		}
		BW.newLine();
		
		for(int i=0; i<passedPSMs.size(); i++) {
			pXgRecord pr = passedPSMs.get(i);
			BW.append(pr.toString());
			BW.newLine();
		}
		
		BW.close();
	}
	
	public static ArrayList<pXgRecord> getMarkedPSMFilter (ArrayList<pXgRecord> records, String title, double fdr) {
		ArrayList<pXgRecord> passedPSMs = new ArrayList<pXgRecord>();
		
		double numTarget = 0;
		double numDecoy = 0;
		int passedTarget = 0;
		int passedDecoy = 0;
		double thisFDR = 1;
		int boundIdx = -1;
		
		for(int i=0; i<records.size(); i++) {
			pXgRecord record = records.get(i);
			
			if(record.isTarget()) {
				numTarget++;
				
			} else {
				numDecoy++;
			}
			
			if(numTarget != 0 && numDecoy/numTarget < fdr) {
				thisFDR = numDecoy/numTarget;
				boundIdx = i;
			}
		}
		
		for(int i=0; i<records.size(); i++) {
			pXgRecord record = records.get(i);
			
			if(i <= boundIdx) {
				record.setValueByFieldName("FDRFilter", "Pass");
				passedPSMs.add(record);
				if(record.isTarget()) {
					passedTarget++;
				} else {
					passedDecoy++;
				}
			} else {
				record.setValueByFieldName("FDRFilter", "Fail");
				passedPSMs.add(record);
			}
			
			
		}
		
		System.out.println(" "+title);
		System.out.println("  -target/decoy: "+passedTarget+"/"+passedDecoy);
		System.out.println("  -FDR: "+thisFDR);
		
		
		return passedPSMs;
	}
	

	public static void parseOptions (String[] args) {
		CommandLine cmd = null;
		Options options = new Options();
		
		// Mandatory
		Option optionInput = Option.builder("i")
				.longOpt("input").argName("pXg")
				.hasArg()
				.required(true)
				.desc("pXg result")
				.build();
		
		Option optionTargetTSV = Option.builder("t")
				.longOpt("target").argName("tsv")
				.hasArg()
				.required(true)
				.desc("A list of target PSMs reported by Percolator or Mokapot.")
				.build();
		
		Option optionDecoyTSV = Option.builder("d")
				.longOpt("decoy").argName("tsv")
				.hasArg()
				.required(true)
				.desc("A list of decoy PSMs reported by Percolator or Mokapot.")
				.build();
		
		Option optionOutput = Option.builder("o")
				.longOpt("output").argName("tsv")
				.hasArg()
				.required(true)
				.desc("specifies the output file prefix, generating the final result as [prefix].fdr.tsv.")
				.build();
		
		Option optionFDR = Option.builder("f")
				.longOpt("fdr").argName("float")
				.hasArg()
				.required(true)
				.desc("FDR value. Default value is 0.01.")
				.build();
		
		Option optionLengthSpecificFDR = Option.builder("l")
				.longOpt("length_specific")
				.required(false)
				.desc("Calculate a length-specific FDR")
				.build();
		
		Option optionPeptideFDR = Option.builder("p")
				.longOpt("peptide")
				.required(false)
				.desc("Calculate a peptide level FDR")
				.build();
		
		Option optionSoftware = Option.builder("s")
				.longOpt("software").argName("percolator|mokapot")
				.hasArg()
				.required(false)
				.desc("Specify the software tool used to generate target and decoy rescoring PSMs. The default is percolator.")
				.build();
		
		
		
		options.addOption(optionInput)
		.addOption(optionTargetTSV)
		.addOption(optionDecoyTSV)
		.addOption(optionPeptideFDR)
		.addOption(optionFDR)
		.addOption(optionOutput)
		.addOption(optionLengthSpecificFDR)
		.addOption(optionSoftware);
		
		CommandLineParser parser = new DefaultParser();
	    HelpFormatter helper = new HelpFormatter();
	    boolean isFail = false;
	    
		try {
		    cmd = parser.parse(options, args, false);
		    
		    if(cmd.hasOption("i")) {
		    	pXgFile = new File(cmd.getOptionValue("i"));
		    }
		    
		    if(cmd.hasOption("t")) {
		    	targetTSVFile = new File(cmd.getOptionValue("t"));
		    }
		    
		    if(cmd.hasOption("d")) {
		    	decoyTSVFile = new File(cmd.getOptionValue("d"));
		    }
		    
		    if(cmd.hasOption("p")) {
		    	isPeptideLevel = true;
		    }
		    
		    if(cmd.hasOption("f")) {
		    	fdr = Double.parseDouble(cmd.getOptionValue("f"));
		    }
		    
		    if(cmd.hasOption("s")) {
		    	software = cmd.getOptionValue("s");
		    }
		    
		    if(cmd.hasOption("l")) {
		    	isLengthSpecific = true;
		    }
		    if(cmd.hasOption("o")) {
		    	prefixOfOutputFile = new File(cmd.getOptionValue("o"));
		    }
		} catch (ParseException e) {
			System.out.println(e.getMessage());
			isFail = true;
		}
		
		if(isFail) {
		    helper.printHelp("Usage:", options);
		    System.exit(0);
		}
		
		System.out.println();
	}
}
