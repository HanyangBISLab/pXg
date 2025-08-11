package progistar.pXg.data.parser;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import progistar.pXg.constants.Constants;
import progistar.pXg.data.GenomicSequence;
import progistar.pXg.data.Global;
import progistar.pXg.data.PxGAnnotation;
import progistar.pXg.data.XBlock;

public class ResultParser {

	private ResultParser () {}

	public static PxGAnnotation parseResult (ArrayList<File> files) {
		PxGAnnotation annotation = new PxGAnnotation();

		try {
			for(File file : files) {
				System.out.println("parsing "+file.getName()+" ...");
				BufferedReader BR = new BufferedReader(new FileReader(file));
				String line = null;

				String uniqueID = null;
				boolean isDecoy = false;
				double meanQScore = 0;

				// for unmapped reads
				StringBuilder transcriptsWithOutBANlist = new StringBuilder();

				String[] field = null;
				while((line = BR.readLine()) != null) {
					// initialize transcript buffer
					transcriptsWithOutBANlist.setLength(0);

					field = line.split("\t");
					if(field[0].equalsIgnoreCase(Constants.OUTPUT_G_UNIQUE_ID)) {
						uniqueID = field[1];
						// decoy decision
						if(uniqueID.startsWith(Constants.DECOY_PREFIX)) {
							isDecoy = true;
						} else {
							isDecoy = false;
						}
					} else if(field[0].equalsIgnoreCase(Constants.OUTPUT_G_QSCORE)) {
						meanQScore = Double.parseDouble(field[1]);
					} 
					else if(field[0].equalsIgnoreCase(Constants.OUTPUT_G_PEPTIDE)) {

						try {
							String pSeq = field[1]; // peptide sequence without I/L consideration
							XBlock xBlock = new XBlock();
							xBlock.genomicLocus = field[2];
							xBlock.strand = field[3].charAt(0);

							xBlock.leftFlankSequenceIdx = Global.getSequenceValue(field[4]);
							xBlock.genomicSequenceIdx = Global.getSequenceValue(field[5]);
							xBlock.rightFlankSequenceIdx = Global.getSequenceValue(field[6]);

							xBlock.leftFlankRefSequenceIdx = Global.getSequenceValue(field[7]);
							xBlock.referenceSequenceIdx = Global.getSequenceValue(field[8]);
							xBlock.rightFlankRefSequenceIdx = Global.getSequenceValue(field[9]);

							xBlock.mutations = field[10];
							xBlock.mutationStatus = field[11];
							xBlock.tAnnotations = field[12];

							String[] exonLengths = field[13].split("\\|");
							String[] percentFullDistances = field[14].split("\\|");
							String[] percentExonDistances = field[15].split("\\|");
							String[] percentCDSDistances = field[16].split("\\|");
							String[] fromStartDistances = field[17].split("\\|");
							String[] fromStopDistances = field[18].split("\\|");

							// If unmapped reads, merging xBlocks and making a single contig xBlock.
							// antisense checker if three-frame
							/**
							 * Antisense event is not allowed.
							 * TODO: More events can be added in future.
							 */
							String[] transcripts = xBlock.tAnnotations.split("\\|");
							Boolean[] bans = new Boolean[transcripts.length];
							for(int i=0; i<transcripts.length; i++) {
								String transcript = transcripts[i];
								if(!transcript.contains(";antisense;")) {
									if(transcriptsWithOutBANlist.length() != 0) {
										transcriptsWithOutBANlist.append("|");
									}
									transcriptsWithOutBANlist.append(transcript);
									bans[i] = false;
								} else {
									bans[i] = true;
								}
							}

							// If all transcripts were discarded, then it means only possible interpretation is antisense
							// in this case, we accept antisense even three-frame translation.
							if(transcriptsWithOutBANlist.length() != 0) {
								xBlock.tAnnotations = transcriptsWithOutBANlist.toString();

								exonLengths = getWithoutBanList(exonLengths, bans);
								percentFullDistances = getWithoutBanList(percentFullDistances, bans);
								percentExonDistances = getWithoutBanList(percentExonDistances, bans);
								percentCDSDistances = getWithoutBanList(percentCDSDistances, bans);
								fromStartDistances = getWithoutBanList(fromStartDistances, bans);
								fromStopDistances = getWithoutBanList(fromStopDistances, bans);
							}

							if(xBlock.strand == '+') {
								xBlock.peptideSequence = GenomicSequence.translation(Global.SEQUENCE_ARRAYLIST.get(xBlock.genomicSequenceIdx), 0);
							} else {
								xBlock.peptideSequence = GenomicSequence.reverseComplementTranslation(Global.SEQUENCE_ARRAYLIST.get(xBlock.genomicSequenceIdx), 0);
							}

							if(xBlock.mockReadCount != 0 || xBlock.targetReadCount != 0 || xBlock.decoyQScore != 0 || xBlock.targetQScore != 0) {
								System.out.println("Line 100 at ResultParser.java: why they are not zero?");
							}

							if(isDecoy) {
								xBlock.mockReadCount++;
								xBlock.decoyQScore += meanQScore;
							} else {
								xBlock.targetReadCount++;
								xBlock.targetQScore += meanQScore;
							}

							xBlock.sequenceID = uniqueID;
							if(xBlock.tAnnotations.contains(Constants.MARK_UNMAPPED+"")) {
								// those xBlocks are needed to merge into contig.
								xBlock.genomicLocus = Constants.ID_NULL;
								// TODO: assembly?
								// Assembler.addXBlock(pSeq, xBlock);
							}

							// If two records share the same genomic loci, their relative genomic information should be the same.
							Global.putExonLengths(xBlock.getRelInfoKey(), exonLengths, false);
							Global.putPercentFullDist(xBlock.getRelInfoKey(), percentFullDistances, false);
							Global.putPercentCDSDist(xBlock.getRelInfoKey(), percentCDSDistances, false);
							Global.putPercentExonDist(xBlock.getRelInfoKey(), percentExonDistances, false);
							Global.putStartDist(xBlock.getRelInfoKey(), fromStartDistances, false);
							Global.putStopDist(xBlock.getRelInfoKey(), fromStopDistances, false);

							annotation.putXBlock(pSeq, xBlock);
						}catch(Exception e) {
							e.printStackTrace();
							System.out.println(line);
						}
					}
				}

				BR.close();
			}
		}catch(IOException ioe ) {

		}

		// just write down unmapped reads (which mapped to any peptides)
		//Assembler.write();

		/*
		// assemble xBlocks
		ArrayList<XBlock> xBlocks = Assembler.assemble();

		xBlocks.forEach(xBlock -> {
			String pSeq = Assembler.getPSeq(xBlock);
			annotation.putXBlock(pSeq, xBlock);
		});
		*/

		return annotation;
	}

	public static String[] getWithoutBanList (String[] list, Boolean[] isBan) {
		// if the list has been already updated by other records,
		if(list.length != isBan.length) {
			return list;
		}

		StringBuilder str = new StringBuilder();

		for(int i=0; i<list.length; i++) {
			if(!isBan[i]) {
				if(str.length() != 0) {
					str.append("|");
				}
				str.append(list[i]);
			}
		}

		return str.toString().split("\\|");
	}
}
