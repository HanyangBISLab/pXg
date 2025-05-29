package progistar.pXg.data;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Hashtable;
import java.util.Iterator;

import progistar.pXg.constants.Constants;
import progistar.pXg.constants.Parameters;
import progistar.pXg.data.parser.SAMExportor;

public class PxGAnnotation {

	// the first key: peptide sequence without I/L consideration
	// the first value: xBlocks corresponding to the key
	// and the second key: peptide sequence from nucleotides + "_" + genomic locus
	private Hashtable<String, Hashtable<String, XBlock>> xBlockMapper = new Hashtable<>();
	// this is a subset of xBlockMapper to map target xBlocks (pass the rna cutoff)
	private Hashtable<String, Hashtable<String, XBlock>> targetXBlockMapper = new Hashtable<>();
	// this is a subset of xBlockMapper to map decoy xBlocks (pass the rna cutoff)
	private Hashtable<String, Hashtable<String, XBlock>> decoyXBlockMapper = new Hashtable<>();
	private int maxNGSreadCount = 0;

	/**
	 * Do not modify by this return object.<br>
	 * Strongly prevent modifying object! Look-up is okay. <br>
	 *
	 * @return
	 */
	public Hashtable<String, Hashtable<String, XBlock>> getXBlockMapper () {
		return this.xBlockMapper;
	}

	public void putXBlock (String pSeq, XBlock xBlock) {
		Hashtable<String, XBlock> xBlocks = this.xBlockMapper.get(pSeq);

		if(xBlocks == null) {
			xBlocks = new Hashtable<>();
			this.xBlockMapper.put(pSeq, xBlocks);
		}

		String key = xBlock.getKey();

		XBlock thisXBlock = xBlocks.get(key);
		if(thisXBlock == null) {
			thisXBlock = xBlock;
			xBlocks.put(key, thisXBlock);
		} else {
			thisXBlock.mockReadCount += xBlock.mockReadCount;
			thisXBlock.targetReadCount += xBlock.targetReadCount;

			if(Parameters.PHRED_CAL.equalsIgnoreCase(Constants.CAL_PHRED_AVG)) {
				thisXBlock.decoyQScore += xBlock.decoyQScore;
				thisXBlock.targetQScore += xBlock.targetQScore;
			} else if(Parameters.PHRED_CAL.equalsIgnoreCase(Constants.CAL_PHRED_MAX)) {
				thisXBlock.decoyQScore = Math.max(xBlock.decoyQScore, thisXBlock.decoyQScore);
				thisXBlock.targetQScore = Math.max(xBlock.targetQScore, thisXBlock.targetQScore);
			} else if(Parameters.PHRED_CAL.equalsIgnoreCase(Constants.CAL_PHRED_MIN)) {
				thisXBlock.decoyQScore = Math.min(xBlock.decoyQScore, thisXBlock.decoyQScore);
				thisXBlock.targetQScore = Math.min(xBlock.targetQScore, thisXBlock.targetQScore);
			}

			thisXBlock.siblingXBlocks.add(xBlock); //
		}

		maxNGSreadCount = Math.max(thisXBlock.mockReadCount, maxNGSreadCount);
		maxNGSreadCount = Math.max(thisXBlock.targetReadCount, maxNGSreadCount);
	}

	/**
	 *
	 *
	 */
	public void markFasta () {

		System.out.print("Marking fasta ids... ");
		long startTime = System.currentTimeMillis();

		Fasta fasta = new Fasta(Parameters.proteinFastaPath);

		ArrayList<String> sequences = PeptideAnnotation.enumerateSequence();
		Hashtable<String, ArrayList<String>> matchedList = fasta.findAll(sequences);

		// assign fasta IDs to pBlock
		for(PBlock pBlock : PeptideAnnotation.pBlocks) {
			ArrayList<String> ids = matchedList.get(pBlock.getPeptideSequence(Parameters.leucineIsIsoleucine));
			if(ids != null) {
				pBlock.fastaIDs = new String[ids.size()];
				for(int i=0; i<pBlock.fastaIDs.length; i++) {
					pBlock.fastaIDs[i] = ids.get(i);
				}
				// 22.01.24
				// if the pBlock has fasta!
				// 23.12.18
				// it was deprecated.
				// pBlock.isCannonical = true;
			}
		}

		long endTime = System.currentTimeMillis();
		System.out.println("\tElapsed time: "+((endTime-startTime)/1000) + " sec");
	}

	public void write (String fileName) {
		try {
			// TSV output
			File tsvFile = new File(fileName);
			BufferedWriter BW = new BufferedWriter(new FileWriter(tsvFile));
			File samFile = null;
			BufferedWriter BWSAM = null;

			if(Parameters.EXPORT_SAM) {
				samFile = new File(Parameters.exportSAMPaths[Parameters.CURRENT_FILE_INDEX]);
				BWSAM = new BufferedWriter(new FileWriter(samFile));
			}

			ArrayList<PBlock> pBlocks = PeptideAnnotation.pBlocks;
			BW.append("SpecID").append("\t");
			BW.append("GenomicID").append("\t");
			BW.append("Label").append("\t");
			BW.append(PeptideAnnotation.toFields()).append("\t");
			BW.append("DeltaScore").append("\t");
			BW.append("Rank").append("\t");
			BW.append("GenomicLociCount").append("\t");
			BW.append(Constants.AA_VARIANT_COLUMN_NAME).append("\t");
			BW.append(Constants.INFERRED_PEPTIDE_COLUMN_NAME).append("\t");
			
			BW.append("GenomicLoci").append("\t");
			BW.append("Strand").append("\t");

			BW.append("ObservedLeftFlankNucleotide").append("\t");
			BW.append("ObservedNucleotide").append("\t");
			BW.append("ObservedRightFlankNucleotide").append("\t");

			BW.append("ReferenceLeftFlankNucleotide").append("\t");
			BW.append("ReferenceNucleotide").append("\t");
			BW.append("ReferenceRightFlankNucleotide").append("\t");

			BW.append("Mutations").append("\t");
			BW.append("MutationStatus").append("\t");
			BW.append("TranscriptIDs").append("\t");
			BW.append("TranscriptIDCount").append("\t");
			BW.append("GeneIDs").append("\t");
			BW.append("GeneIDCount").append("\t");
			BW.append("GeneNames").append("\t");
			BW.append("GeneNameCount").append("\t");

			BW.append("PercentFullDistance").append("\t");
			BW.append("PercentExonDistance").append("\t");
			BW.append("PercentCDSDistance").append("\t");
			BW.append("FromCDSStartSite").append("\t");
			BW.append("FromCDSStopSite").append("\t");

			BW.append("Events").append("\t");
			BW.append("EventCount").append("\t");
			BW.append("FastaIDs").append("\t");
			BW.append("FastaIDCount").append("\t");
			BW.append("Reads").append("\t");
			BW.append("MeanQScore").append("\t");
			BW.append(Constants.CLASS_COLUMN_NAME);
			BW.newLine();

			File outFile = new File(Parameters.unmappedFilePaths[Parameters.CURRENT_FILE_INDEX]);
			BufferedWriter BWUnmapped = new BufferedWriter(new FileWriter(outFile));

			// calculate genomic ID
			// The combination of genomic loci + strand + nucleotide sequence is mapping to unique gneomic ID.
			Hashtable<String, Integer> genomicIDMapper = new Hashtable<>();
			boolean[] isTarget = new boolean[1];
			StringBuilder aaVariants = new StringBuilder();
			for(PBlock pBlock : pBlocks) {
				// peptide sequence without I/L consideration
				String key = pBlock.getPeptideSequence(Parameters.leucineIsIsoleucine);
				
				// generate aaVariant annotation
				aaVariants.setLength(0);
				if(pBlock.aaVariant == null) {
					aaVariants.append(Constants.ID_NULL);
				} else {
					AAVariant aaVar = pBlock.aaVariant;
					String aaRNA = AAVariantTable.aaRNAArray[aaVar.aaVariantIndex];
					String aaPeptide = AAVariantTable.aaPeptideArray[aaVar.aaVariantIndex];
					// position: zero to one base
					aaVariants.append(aaVar.pos+1).append(":").append(aaRNA).append(">").append(aaPeptide);
				}

				// for target
				isTarget[0] = true;
				Hashtable<String, XBlock> xBlocksTmp = this.targetXBlockMapper.get(key);
				if(xBlocksTmp == null) {
					xBlocksTmp = this.decoyXBlockMapper.get(key);
					isTarget[0] = false;
				}

				// for final state
				Hashtable<String, XBlock> xBlocks = xBlocksTmp;

				// there is a mapping.
				if(xBlocks != null) {
					xBlocks.forEach((pSeq, xBlock) -> {
						try {
							// assign fastaIDs.
							xBlock.fastaIDs = pBlock.fastaIDs;
							// assigne genomic ID
							String genomicIDKey = xBlock.genomicLocus+"|"+xBlock.strand+"|"+Global.SEQUENCE_ARRAYLIST.get(xBlock.genomicSequenceIdx);
							Integer genomicID = genomicIDMapper.get(genomicIDKey);
							if(genomicID == null) {
								genomicID = genomicIDMapper.size()+1;
								genomicIDMapper.put(genomicIDKey, genomicID);
							}

							// we treated unmapped reads as '0' genomic loci count.
							String gLociCount = "0";
							if(!xBlock.isMappedAmbiguous()) {
								gLociCount = xBlocks.size()+"";
							}

							// TSV writer
							BW.append(pBlock.toString(genomicID)).append("\t")
							.append(pBlock.deltaScore+"\t")
							.append(pBlock.rank+"\t")
							.append(gLociCount+"\t")
							.append(aaVariants.toString()+"\t")
							.append(xBlock.toString(pBlock.psmStatus))
							.append("\t"+pBlock.isReference);
							BW.newLine();

							// Jul 23, 2024
							// For additional information, 
							// decoy information is neglected.
							
							// if this is unmapped, then store.
							if(xBlock.isMappedAmbiguous() && isTarget[0]) {
								ArrayList<XBlock> unmappedXBlocks = new ArrayList<>();
								unmappedXBlocks.add(xBlock);
								unmappedXBlocks.addAll(xBlock.siblingXBlocks);

								BWUnmapped.append(">"+xBlock.peptideSequence);
								BWUnmapped.newLine();
								for(XBlock thisXBlock : unmappedXBlocks) {
									BWUnmapped.append(thisXBlock.sequenceID).append("\t")
									.append(thisXBlock.fullReadSequence).append("\t")
									.append(Global.SEQUENCE_ARRAYLIST.get(thisXBlock.genomicSequenceIdx));
									BWUnmapped.newLine();
								}
 							} else if(isTarget[0]){
 								// SAM ID Mapper
 								if(Parameters.EXPORT_SAM) {
 									if(		(Parameters.EXPORT_REFERENCE && pBlock.isReference) ||
 											(Parameters.EXPORT_NON_REFERENCE && !pBlock.isReference)) {
 										SAMExportor.putSequenceID(xBlock);
 									}
 								}
 							}
						}catch(IOException ioe) {

						}
					});
				}
			}

			// write exportSAM
			BW.close();
			BWUnmapped.close();
			if(Parameters.EXPORT_SAM) {
				SAMExportor.writeSAM(BWSAM);
				BWSAM.close();
			}
		}catch(IOException ioe) {

		}
	}
	
	/**
	 * Marking target PSM <br>
	 *
	 */
	public void assignXBlocks () {
		long startTime = System.currentTimeMillis();

		ArrayList<PBlock> pBlocks = PeptideAnnotation.pBlocks;
		PeptideAnnotation.pBlocks = new ArrayList<>();
		// update pBlocks!
		int size = pBlocks.size();
		
		for(int i=0; i<size; i++) {
			PBlock pBlock = pBlocks.get(i);
			// peptide sequence with I/L consideration
			String key = pBlock.getPeptideSequence(Parameters.leucineIsIsoleucine);

			// check error!
			boolean[] expAndMocks = new boolean[2];

			Hashtable<String, XBlock> xBlocks = this.targetXBlockMapper.get(key);
			if(xBlocks != null) {
				xBlocks.forEach((key_, xBlock) -> {
					if(xBlock.targetReadCount > 0) {
						pBlock.psmStatus = Constants.PSM_STATUS_TARGET;
						pBlock.isReference |= xBlock.isCannonical();
						pBlock.targetXBlocks.put(xBlock.getKey(), xBlock);
						expAndMocks[0] = true;
					}
				});
			}

			xBlocks = this.decoyXBlockMapper.get(key);
			if(xBlocks != null) {
				xBlocks.forEach((key_, xBlock) -> {
					if(xBlock.mockReadCount > 0) {
						if(expAndMocks[0]) {
							pBlock.psmStatus = Constants.PSM_STATUS_BOTH;
						} else {
							pBlock.psmStatus = Constants.PSM_STATUS_DECOY;
							pBlock.isReference |= xBlock.isCannonical();
						}
						pBlock.decoyXBlocks.put(xBlock.getKey(), xBlock);
						expAndMocks[1] = true;
					}
				});
			}
			
			// if the pBlock has aa variant?
			// it should be non-canonical
			if(pBlock.aaVariant != null) {
				pBlock.isReference = false;
			}

			// remove unassigned pBlock
			if(expAndMocks[0] || expAndMocks[1]) {
				PeptideAnnotation.pBlocks.add(pBlock);
			}
		}
		pBlocks.clear();
		long endTime = System.currentTimeMillis();
		
		System.out.println("\tElapsed time: "+((endTime-startTime)/1000) + " sec");

	}
	
	/**
	 * This method expects that: <br>
	 * 1) following estimatePvalueTreshold.<br>
	 * 2) so, pSeq has exp or mock read only.<br>
	 *
	 *
	 */
	public void regionScoreFilter () {
		// prevent zeor-size
		if(this.xBlockMapper.size() == 0) {
			return;
		}
		// filter regions in the same locus and nucleotides.
		Iterator<String> pSeqs = (Iterator<String>) this.xBlockMapper.keys();

		while(pSeqs.hasNext()) {
			String pSeq = pSeqs.next();
			Hashtable<String, XBlock> xBlocks = this.xBlockMapper.get(pSeq);
			if(xBlocks == null || xBlocks.size() == 0) {
				continue;
			}

			Iterator<String> keys = (Iterator<String>) xBlocks.keys();

			double minExpXBlockPenalty = Double.MAX_VALUE;
			double minMockXBlockPenalty = Double.MAX_VALUE;
			while(keys.hasNext()) {
				String key = keys.next();
				XBlock xBlock = xBlocks.get(key);

				// filter-out
				xBlock.filterRegions();

				// find min penalty
				if(xBlock.targetReadCount >=1) {
					minExpXBlockPenalty = Math.min(xBlock.bestRegionPriority, minExpXBlockPenalty);
				}
				if(xBlock.mockReadCount >= 1) {
					minMockXBlockPenalty = Math.min(xBlock.bestRegionPriority, minMockXBlockPenalty);
				}
			}

			keys = (Iterator<String>) xBlocks.keys();

			// target and decoy xBlocks are assigned.
			while(keys.hasNext()) {
				String key = keys.next();
				XBlock xBlock = xBlocks.get(key);

				if(xBlock.targetReadCount >= 1) {
					if(xBlock.bestRegionPriority <= minExpXBlockPenalty) {
						// for target xBlock mapper
						Hashtable<String, XBlock> xBlocks_ = this.targetXBlockMapper.get(pSeq);
						if(xBlocks_ == null) {
							xBlocks_ = new Hashtable<>();
							this.targetXBlockMapper.put(pSeq, xBlocks_);
						}
						xBlocks_.put(key, xBlock);
					}
				}
				if(xBlock.mockReadCount >= 1) {
					if(xBlock.bestRegionPriority <= minMockXBlockPenalty) {
						// for decoy xBlock mapper
						Hashtable<String, XBlock> xBlocks_ = this.decoyXBlockMapper.get(pSeq);
						if(xBlocks_ == null) {
							xBlocks_ = new Hashtable<>();
							this.decoyXBlockMapper.put(pSeq, xBlocks_);
						}
						xBlocks_.put(key, xBlock);
					}
				}
			}
		}
	}
	
	/**
	 * A peptide with/without AAVariants is aggregated by a record id.<br>
	 * If there is a peptide without AAVariants, all the peptides with AAVariants are ignored. 
	 * 
	 * Priority<br>
	 * Target > Decoy > Non-AAVariant > AAVariant
	 * 
	 * SAAV filter + TD filter
	 * 
	 */
	public void prioritizePeptidesInEachCandidate () {
		ArrayList<PBlock> newPBlocks = new ArrayList<PBlock>();
		
		int size = PeptideAnnotation.pBlocks.size();
		Hashtable<Integer, ArrayList<PBlock>> aggregates = new Hashtable<Integer, ArrayList<PBlock>>();
		for(int i=0; i<size; i++) {
			PBlock pBlock = PeptideAnnotation.pBlocks.get(i);
			ArrayList<PBlock> pBlocks = aggregates.get(pBlock.recordId);
			if(pBlocks == null) {
				pBlocks = new ArrayList<PBlock>();
			}
			
			pBlocks.add(pBlock);
			aggregates.put(pBlock.recordId, pBlocks);
		}
		
		aggregates.forEach((id, pBlocks)->{
			byte maxScore = 0;
			// check AAVariants
			for(PBlock pBlock : pBlocks) {
				byte score = 0;
				// target = 20. Decoy = 10.
				if(pBlock.psmStatus == Constants.PSM_STATUS_DECOY) {
					score = 10;
				} else {
					score = 20;
				}
				
				// Non-AAvariant = + 2. AAVariant = + 1.
				if(pBlock.isAAVariant()) {
					score += 1;
				} else {
					score += 2;
				}
				
				maxScore = score > maxScore ? score : maxScore;
			}
			
			for(PBlock pBlock : pBlocks) {
				byte score = 0;
				// target = 20. Decoy = 10.
				if(pBlock.psmStatus == Constants.PSM_STATUS_DECOY) {
					score = 10;
				} else {
					score = 20;
				}
				
				// Non-AAvariant = + 2. AAVariant = + 1.
				if(pBlock.isAAVariant()) {
					score += 1;
				} else {
					score += 2;
				}
				
				if(score < maxScore) {
					pBlock.recordId = -1; // mark as delete
				}
			}
		});
		
		
		// remove lower priority peptides in each candidate.
		for(int i=0; i<size; i++) {
			PBlock pBlock = PeptideAnnotation.pBlocks.get(i);
			if(pBlock.recordId != -1) {
				newPBlocks.add(pBlock);
			}
		}
		
		PeptideAnnotation.pBlocks.clear();
		PeptideAnnotation.pBlocks = newPBlocks;
		
		Collections.sort(PeptideAnnotation.pBlocks, new Comparator<PBlock>() {
			@Override
			public int compare(PBlock o1, PBlock o2) {
				if(o1.recordId < o2.recordId) {
					return -1;
				} else if(o1.recordId > o2.recordId) {
					return 1;
				}
				return 0;
			}
        });
		
	}
}
