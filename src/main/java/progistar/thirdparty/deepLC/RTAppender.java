package progistar.thirdparty.deepLC;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class RTAppender {

	public static Pattern MODREG = Pattern.compile("([+0-9.]+)");
	public static String MET = "Oxidation";
	public static String DEAMI = "Deamidated";
	public static String CARBAM = "Carbamidomethyl";
	public static String PHOSPHO = "Phospho";
	public static String CYSTEINLY = "Cysteinyl";
	
	public static void main(String[] args) throws IOException {
		String pXgFileName = "/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/7.Unmodified/pXg/S4.UNMOD.PEAKS.pxg.BA.fdr";
		String deeplcFolder = "/Users/gistar/projects/pXg/Laumont_NatCommun2016/RT/S4_deeplc_predictions.csv";
		
		File file = new File(deeplcFolder);
		// key
		// fileName_InferredPeptide_PTM_RT
		// value
		// predicted RT
		Hashtable<String, String> RTMapper = new Hashtable<String, String>();
		
		BufferedReader BR = new BufferedReader(new FileReader(file));
		String line = null;
		
		BR.readLine(); // skip header
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\\,");
			String key = fields[1]+"_"+fields[2]+"_"+Double.parseDouble(fields[3]);
			
			RTMapper.put(key, fields[4]);
		}
		BR.close();
		
		// read PXG	
		Hashtable<String, String> massToModName = new Hashtable<String, String>();
		massToModName.put("+15.99", MET);
		massToModName.put("+.98", DEAMI);
		massToModName.put("+119.00", CYSTEINLY);
		massToModName.put("+79.97", PHOSPHO);
		massToModName.put("+58.01", CARBAM);
		
		int peaksPeptideIndex = 3;
		int rtIndex = 11;
		int peptideIndex = 20;
		
		BR = new BufferedReader(new FileReader(pXgFileName));
		
		// header
		System.out.println(BR.readLine()+"\tdeeplcRT");
		
		
		while((line = BR.readLine()) != null) {
			String[] fields = line.split("\t");
			String modifications = "";
			if(fields[peaksPeptideIndex].contains("+") || fields[peaksPeptideIndex].contains("-")) {
				Matcher matcher = MODREG.matcher(fields[peaksPeptideIndex]);
				while(matcher.find()) {
					String group = matcher.group();
					int pos = fields[peaksPeptideIndex].indexOf("("+group+")");
					fields[peaksPeptideIndex] = fields[peaksPeptideIndex].replaceFirst("\\([+0-9.]+\\)", "");
					modifications += pos+"|"+massToModName.get(group)+"|";
				}
				modifications = modifications.substring(0, modifications.length()-1);
			}
			//
			
			String key = fields[peptideIndex]+"_"+modifications+"_"+Double.parseDouble(fields[rtIndex]);
			String pRT = RTMapper.get(key);
			
			if(pRT == null) {
				// multiple genomic loci with ambiguous I/L
				System.out.println(line+"\t-");
			} else {
				System.out.println(line+"\t"+pRT);
			}
		}
		
		
		BR.close();
	}
}
