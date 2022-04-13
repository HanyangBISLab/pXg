package progistar.thirdparty.deepLC;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Hashtable;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

class DeepLCData implements Comparable<DeepLCData> {
	public String content = null;
	public double score = 0;
	public double rt = 0;
	@Override
	public int compareTo(DeepLCData o) {
		
		if(this.score > o.score) {
			return -1;
		} else if(this.score < o.score) {
			return 1;
		}
		
		return 0;
	}
	
	
}

public class pXg2DeepLCInput {

	public static Pattern MODREG = Pattern.compile("([+0-9.]+)");
	public static String MET = "Oxidation";
	public static String DEAMI = "Deamidated";
	public static String CARBAM = "Carbamidomethyl";
	public static String PHOSPHO = "Phospho";
	public static String CYSTEINLY = "Cysteinyl";
	
	public static void main(String[] args) throws IOException {
		
		Hashtable<String, String> massToModName = new Hashtable<String, String>();
		massToModName.put("+15.99", MET);
		massToModName.put("+.98", DEAMI);
		massToModName.put("+119.00", CYSTEINLY);
		massToModName.put("+79.97", PHOSPHO);
		massToModName.put("+58.01", CARBAM);
		
		int calSize = 25;
		int peaksPeptideIndex = 3;
		int rtIndex = 11;
		int fractionIndex = 1;
		int scanIndex = 4;
		int peptideIndex = 20;
		int peaksScoreIndex = 7;
		
		String fileName = "/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/6.CaliScanwoDeamiFix/pXg/PeptideAnnotationS3_5ppm_002_recal.scanNum.rep1.rank10.pXg.BA.fdr";
		
		BufferedReader BR = new BufferedReader(new FileReader(fileName));
		
		String line = null;
		
		
		BR.readLine(); // skip header
		Hashtable<String, String> rmDuplications = new Hashtable<String, String>();
		Hashtable<String, ArrayList<DeepLCData>> fractions = new Hashtable<String, ArrayList<DeepLCData>>();
		
		
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
			
			ArrayList<DeepLCData> fraction = fractions.get(fields[fractionIndex]);
			if(fraction == null) {
				fraction = new ArrayList<DeepLCData>();
				fractions.put(fields[fractionIndex], fraction);
			}
			
			if(rmDuplications.get(fields[fractionIndex]+"_"+fields[scanIndex]+"_"+fields[peptideIndex]) == null) {
				DeepLCData data = new DeepLCData();
				String key = fields[peptideIndex]+","+modifications+","+fields[rtIndex];
				
				data.content = key;
				data.score = Double.parseDouble(fields[peaksScoreIndex]);
				data.rt = Double.parseDouble(fields[rtIndex]);
				
				fraction.add(data);
				rmDuplications.put(fields[fractionIndex]+"_"+fields[scanIndex]+"_"+fields[peptideIndex], "");
			}
		}
		
		fractions.forEach( (fName, fraction) -> {
			try {
				BufferedWriter BW = new BufferedWriter(new FileWriter(fName+".csv"));
				BufferedWriter BWCal = new BufferedWriter(new FileWriter(fName+".cal.csv"));
				BW.append("seq,modifications,tr");
				BW.newLine();
				
				BWCal.append("seq,modifications,tr");
				BWCal.newLine();
				
				Collections.sort(fraction);
				
				int cnt = 0;
				while(cnt < 25) {
					for(DeepLCData data : fraction) {
						BW.append(data.content);
						BW.newLine();
						cnt++;
					}
				}
				
				// for calibration input file
				double maxRT = 0;
				for(int i=0; i<fraction.size(); i++) {
					maxRT = Math.max(fraction.get(i).rt, maxRT);
				}
				
				double intervalRT = maxRT/calSize;
				for(int i=0; i<calSize; i++) {
					double startRT = i * intervalRT;
					double endRT = (i+1) * intervalRT;
					
					ArrayList<DeepLCData> calData = new ArrayList<DeepLCData>();
					while(calData.size() < 10) {
						// select top 10 peptides
						for(int j=0; j<fraction.size(); j++) {
							DeepLCData data = fraction.get(j);
							if(data.rt >= startRT && data.rt <= endRT) {
								calData.add(data);
							}
						}
					}
					
					for(DeepLCData data : calData) {
						BWCal.append(data.content);
						BWCal.newLine();
					}
				}
				
				BW.close();
				BWCal.close();
				
				System.out.println(fName+":\t"+fraction.size()+"=> "+cnt);
			}catch(IOException ioe) {
				
			}
		});
		
		
		BR.close();
	}
}
