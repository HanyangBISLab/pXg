package progistar.pXg.data;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Iterator;

import progistar.pXg.constants.Constants;
import progistar.pXg.constants.Parameters;
import progistar.pXg.constants.RunInfo;
import progistar.pXg.data.parser.GTFExportor;
import progistar.pXg.data.parser.SAMExportor;
import progistar.pXg.utils.Logger;

public class PxGAnnotation {

	
	
	// the first key: peptide sequence without I/L consideration
	// the first value: xBlocks corresponding to the key 
	// and the second key: peptide sequence from nucleotides + "_" + genomic locus
	private Hashtable<String, Hashtable<String, XBlock>> xBlockMapper = new Hashtable<String, Hashtable<String, XBlock>>();
	// this is a subset of xBlockMapper to map target xBlocks (pass the rna cutoff)
	private Hashtable<String, Hashtable<String, XBlock>> targetXBlockMapper = new Hashtable<String, Hashtable<String, XBlock>>();
	// this is a subset of xBlockMapper to map decoy xBlocks (pass the rna cutoff)
	private Hashtable<String, Hashtable<String, XBlock>> decoyXBlockMapper = new Hashtable<String, Hashtable<String, XBlock>>();
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
			xBlocks = new Hashtable<String, XBlock>();
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
			thisXBlock.decoyMeanQScore += xBlock.decoyMeanQScore;
			thisXBlock.targetMeanQScore += xBlock.targetMeanQScore;
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
			ArrayList<String> ids = matchedList.get(pBlock.getPeptideSequence());
			if(ids != null) {
				pBlock.fastaIDs = new String[ids.size()];
				for(int i=0; i<pBlock.fastaIDs.length; i++) {
					pBlock.fastaIDs[i] = ids.get(i);
				}
				// 22.01.24
				// if the pBlock has fasta!
				pBlock.isCannonical = true;
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
			
			File gtfFile = null;
			BufferedWriter BWGTF = null;
			
			if(Parameters.EXPORT_GTF) {
				gtfFile = new File(Parameters.exportGTFPath);
				BWGTF = new BufferedWriter(new FileWriter(gtfFile));
			}
			final BufferedWriter BWGTF_ = BWGTF;
			
			File samFile = null;
			BufferedWriter BWSAM = null;
			
			if(Parameters.EXPORT_SAM) {
				samFile = new File(Parameters.exportSAMPath);
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
			BW.append("InferredPeptide").append("\t");
			BW.append("GenomicLoci").append("\t");
			BW.append("Strand").append("\t");
			BW.append("Nucleotide").append("\t");
			BW.append("Mutations").append("\t");
			BW.append("TranscriptIDs").append("\t");
			BW.append("TranscriptIDCount").append("\t");
			BW.append("GeneIDs").append("\t");
			BW.append("GeneIDCount").append("\t");
			BW.append("GeneNames").append("\t");
			BW.append("GeneNameCount").append("\t");
			BW.append("Events").append("\t");
			BW.append("EventCount").append("\t");
			BW.append("FastaIDs").append("\t");
			BW.append("FastaIDCount").append("\t");
			BW.append("Reads").append("\t");
			BW.append("MeanQScore").append("\t");
			BW.append("IsCanonical");
			BW.newLine();
			
			File outFile = new File(Parameters.unmappedFilePath);
			BufferedWriter BWUnmapped = new BufferedWriter(new FileWriter(outFile));
			
			// calculate genomic ID
			// The combination of genomic loci + strand + nucleotide sequence is mapping to unique gneomic ID.
			Hashtable<String, Integer> genomicIDMapper = new Hashtable<String, Integer>();
			for(PBlock pBlock : pBlocks) {
				// peptide sequence without I/L consideration
				String key = pBlock.getPeptideSequence();
				
				// for target
				Hashtable<String, XBlock> xBlocksTmp = this.targetXBlockMapper.get(key);
				if(xBlocksTmp == null && Parameters.isDecoyOut) {
					xBlocksTmp = this.decoyXBlockMapper.get(key);
				}
				
				// for final state
				Hashtable<String, XBlock> xBlocks = xBlocksTmp;
				
				// there is no available mapping.
				if(xBlocks != null) {
					xBlocks.forEach((pSeq, xBlock) -> {
						try {
							// assign fastaIDs.
							xBlock.fastaIDs = pBlock.fastaIDs;
							// assigne genomic ID
							String genomicIDKey = xBlock.genomicLocus+"|"+xBlock.strand+"|"+xBlock.genomicSequence;
							Integer genomicID = genomicIDMapper.get(genomicIDKey);
							if(genomicID == null) {
								genomicID = genomicIDMapper.size()+1;
								genomicIDMapper.put(genomicIDKey, genomicID);
							}
							
							// we treated unmapped reads as '0' genomic loci count.
							String gLociCount = "0";
							if(xBlock.isMapped()) {
								gLociCount = xBlocks.size()+"";
							}
							
							// TSV writer
							BW.append(pBlock.toString(genomicID)).append("\t").append(pBlock.deltaScore+"\t").append(pBlock.rank+"\t").append(gLociCount+"\t").append(xBlock.toString(pBlock.psmStatus)).append("\t"+pBlock.isCannonical);
							BW.newLine();
							
							// if this is unmapped, then store.
							if(!xBlock.isMapped()) {
								ArrayList<XBlock> unmappedXBlocks = new ArrayList<XBlock>();
								unmappedXBlocks.add(xBlock);
								unmappedXBlocks.addAll(xBlock.siblingXBlocks);
								
								BWUnmapped.append(">"+xBlock.peptideSequence);
								BWUnmapped.newLine();
								for(XBlock thisXBlock : unmappedXBlocks) {
									BWUnmapped.append(thisXBlock.sequenceID).append("\t").append(thisXBlock.fullReadSequence).append("\t").append(thisXBlock.genomicSequence);
									BWUnmapped.newLine();
								}
 							} else {
								// GTF writer
 								if(Parameters.EXPORT_GTF) {
 									if(			(Parameters.EXPORT_CANONICAL && pBlock.isCannonical) ||
 											(Parameters.EXPORT_NONCANONICAL && !pBlock.isCannonical)) {
 										GTFExportor.writeGTF(pBlock, xBlock, BWGTF_);
 									}
 								}
 								
 								// SAM ID Mapper
 								if(Parameters.EXPORT_SAM) {
 									if(			(Parameters.EXPORT_CANONICAL && pBlock.isCannonical) ||
 											(Parameters.EXPORT_NONCANONICAL && !pBlock.isCannonical)) {
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
			
			if(Parameters.EXPORT_GTF) {
				BWGTF.close();
			}
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
		
		// update pBlocks!
		for(int i=0; i<pBlocks.size(); i++) {
			PBlock pBlock = pBlocks.get(i);
			// peptide sequence without I/L consideration
			String key = pBlock.getPeptideSequence();
			
			// check error!
			boolean[] expAndMocks = new boolean[2];
			
			Hashtable<String, XBlock> xBlocks = this.targetXBlockMapper.get(key);
			if(xBlocks != null) {
				// only select significantly mapped PSMs
				// this is because we are interested in P( decoy (score > X) | significantly mapped PSMs)
				xBlocks.forEach((key_, xBlock) -> {
					if(xBlock.targetReadCount >= 1) {
						pBlock.psmStatus = Constants.PSM_STATUS_TARGET;
						pBlock.isCannonical |= xBlock.isCannonical();
						pBlock.targetXBlocks.put(xBlock.getKey(), xBlock);
						expAndMocks[0] = true;
					} 
				});
			}
			
			xBlocks = this.decoyXBlockMapper.get(key);
			if(xBlocks != null) {
				// only select significantly mapped PSMs
				// this is because we are interested in P( decoy (score > X) | significantly mapped PSMs)
				xBlocks.forEach((key_, xBlock) -> {
					if(xBlock.mockReadCount >= 1) {
						if(expAndMocks[0]) {
							pBlock.psmStatus = Constants.PSM_STATUS_BOTH;
						} else {
							pBlock.psmStatus = Constants.PSM_STATUS_DECOY;
							pBlock.isCannonical |= xBlock.isCannonical();
						}
						pBlock.decoyXBlocks.put(xBlock.getKey(), xBlock);
						expAndMocks[1] = true;
					}
				});
			}
			
			// remove unassigned pBlock
			if(!expAndMocks[0] && !expAndMocks[1]) {
				pBlocks.remove(i--);
			}
			
		}
		
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
		if(this.xBlockMapper.size() == 0) return;
		// filter regions in the same locus and nucleotides.
		Iterator<String> pSeqs = (Iterator<String>) this.xBlockMapper.keys();

		while(pSeqs.hasNext()) {
			String pSeq = pSeqs.next();
			Hashtable<String, XBlock> xBlocks = this.xBlockMapper.get(pSeq);
			if(xBlocks == null || xBlocks.size() == 0) continue;
			
			Iterator<String> keys = (Iterator<String>) xBlocks.keys();
			
			double minExpXBlockPenalty = Double.MAX_VALUE;
			double minMockXBlockPenalty = Double.MAX_VALUE;
			int peptLen = pSeq.length();
			
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
							xBlocks_ = new Hashtable<String, XBlock>();
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
							xBlocks_ = new Hashtable<String, XBlock>();
							this.decoyXBlockMapper.put(pSeq, xBlocks_);
						}
						xBlocks_.put(key, xBlock);
					}
				}
			}
		}
	}
}
