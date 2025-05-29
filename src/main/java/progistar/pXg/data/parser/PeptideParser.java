package progistar.pXg.data.parser;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import progistar.pXg.constants.Parameters;
import progistar.pXg.constants.RunInfo;
import progistar.pXg.data.PBlock;
import progistar.pXg.data.PeptideAnnotation;

public class PeptideParser {

	private static String[]	commentMarkers;

	private PeptideParser () {}

	/**
	 * The peptideFile must contain field (column names). <br>
	 * Header lines are optional and if exist, then the lines must be positioned on the top. <br>
	 *
	 *
	 * @param peptideFilePath
	 */
	public static void parseResult (String peptideFilePath) {
		
		PeptideAnnotation.pBlocks = new ArrayList<>();
		System.out.println("Parsing peptide file: "+peptideFilePath);
		long startTime = System.currentTimeMillis();

		// there are PTM patterns : May 1, 2025
		Pattern ptmPattern = Pattern.compile(Parameters.ptmParserRegExr);
		// set regular expressions
		Pattern peptideRegExr	= Pattern.compile(Parameters.peptideParserRegExr);
		commentMarkers = Parameters.commentMarker.split("\\|");

		StringBuilder pSeq = new StringBuilder();
		int fieldLength = 0;
		try {
			File file = new File(peptideFilePath);

			BufferedReader BR = new BufferedReader(new FileReader(file));
			String line = null;

			int recordCount = -1;
			while((line = BR.readLine()) != null) {
				// skip header marker
				// comment marker is not considered record.
				for(String headerMarker : commentMarkers) {
					if(line.startsWith(headerMarker)) {
						continue;
					}
				}

				String[] record = null;
				if(Parameters.sepType.equalsIgnoreCase("tsv")) {
					record = line.split("\t");
				} else if(Parameters.sepType.equalsIgnoreCase("csv")) {
					record = line.split(",");
				}

				if(Parameters.rmQuotes) {
					for(int i=0; i<record.length; i++) {
						record[i] = record[i].replace("\"", "");
					}
				}

				// the first line after headers must be field line.
				if(recordCount == -1) {
					PeptideAnnotation.setFields(record);
					fieldLength = PeptideAnnotation.getFieldLength();
				}
				// record
				else {
					// if the number of fields in the record is different from the number of fields in the header,
					// adjust
					if(record.length != fieldLength) {
						String[] adjustedRecord = new String[fieldLength];
						for(int i=0; i<fieldLength; i++) {
							if(record.length > i) {
								adjustedRecord[i] = record[i];
							} else {
								adjustedRecord[i] = "";
							}
						}
						
						// replace an original record to an adjust record.
						record = adjustedRecord;
					}
					
					String peptide = record[Parameters.peptideColumnIndex];
					
					
					// check PTM pattern
					Matcher matcher = ptmPattern.matcher(peptide);
					while(matcher.find()) {
						Parameters.detectedPTMTable.put(matcher.group(), true);
					}
					
					// remove unimod and mass relating patterns
					peptide = peptide.replaceAll(Parameters.ptmParserRegExr, "");
					
					// find peptide strip sequence
					matcher = peptideRegExr.matcher(peptide);

					while(matcher.find()) {
						pSeq.append(matcher.group());
					}

					PBlock pBlock = new PBlock(record, pSeq.toString(), recordCount);

					PeptideAnnotation.pBlocks.add(pBlock);
					pSeq.setLength(0);

				}

				recordCount ++;

			}

			BR.close();
			
			// print what kinds of PTM patterns in it.
			System.out.println("Found PTM patterns: "+Parameters.detectedPTMTable.size());
			if(Parameters.detectedPTMTable.size()>0) {
				Parameters.detectedPTMTable.forEach((ptm, nil)->{
					System.out.println(ptm);
				});
			}
		}catch (IOException ioe) {

		}


		// determine rank
		PeptideAnnotation.assignRank();

		RunInfo.initialPeptideNum = PeptideAnnotation.getPeptideSize();
		RunInfo.initialScanNum = PeptideAnnotation.getScanSize();

		// filter PSMs by delta-rank
		PeptideAnnotation.filter();

		RunInfo.rankFilterPeptideNum1 = PeptideAnnotation.getPeptideSize();
		RunInfo.rankFilterScanNum1 = PeptideAnnotation.getScanSize();

		// length filter
		PeptideAnnotation.peptideLengthFilter();

		RunInfo.lengthFilterPeptideNum2 = PeptideAnnotation.getPeptideSize();
		RunInfo.lengthFilterScanNum2 = PeptideAnnotation.getScanSize();

		
		// add AAVs
		PeptideAnnotation.setAAVs();
		
		// calculate total candidate peptides per rank
		for(PBlock pBlock : PeptideAnnotation.pBlocks) {
			RunInfo.totalRankPSMs[pBlock.rank]++;
		}

		long endTime = System.currentTimeMillis();
		System.out.println("\tElapsed time: "+((endTime-startTime)/1000) + " sec");


		// build keyword-trie
		PeptideAnnotation.buildKeywordTrie();
	}
}
