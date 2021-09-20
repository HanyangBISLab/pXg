package progistar.pXg.data;

import progistar.pXg.constants.Constants;
import progistar.pXg.utils.Codon;
import progistar.pXg.utils.IndexConvertor;

public class Output {

	// positions in NGS-read
	// zero-based
	// [startPos, endPos] both are inclusive.
	// forward-oriented position.
	// it means that it is not care about strand.
	public int startPosInNGS;
	public int endPosInNGS;
	
	public int startGenomicPosition;
	public int endGenomicPosition;
	
	public int peptideIndex;
	public boolean strand;

	public GenomicSequence gSeq;
	
	public Output (GenomicSequence gSeq, int peptideIndex, int startPos, int endPos, boolean strand) {
		this.gSeq = gSeq;
		this.peptideIndex = peptideIndex;
		this.startPosInNGS = startPos;
		this.endPosInNGS = endPos;
		this.strand = strand;
	}
	
	public String getPeptide () {
		return PeptideAnnotation.indexedPeptide.get(this.peptideIndex);
	}
	
	public String getMatchedNucleotide () {
		String nucleotide = this.gSeq.getNucleotideString();
		return nucleotide.substring(this.startPosInNGS, this.endPosInNGS+1);
	}
	
	public byte getFrame (int transcriptNum) {
		TBlock tBlock = gSeq.tBlocks[transcriptNum];
		// this is intergenic
		// or non-coding
		if(tBlock == null || tBlock.transcriptCodingType == Constants.NON_CODING_TRANSCRIPT) return Constants.NO_FRAME;
		// with soft-clip
		if(this.startGenomicPosition == -1 || this.endGenomicPosition == -1) return Constants.NO_FRAME;
		
		
		return Constants.NO_FRAME;
		
	}
	
	public String getLocus () {
		// unknown locus
		if(this.startGenomicPosition == -1 || this.endGenomicPosition == -1) {
			return IndexConvertor.indexToChr(gSeq.chrIndex)+":?";
		}
		
		return IndexConvertor.indexToChr(gSeq.chrIndex) +":" +startGenomicPosition+"-"+endGenomicPosition;
	}
	
	/**
	 * 
	 * Genomic region matched to peptide.<br>
	 * 
	 * @param transcriptNum
	 * @return
	 */
	public String getGenomicRegion (int transcriptNum) {
		StringBuilder genomicRegion = new StringBuilder();
		int relPos = 0;
		for(Cigar cigar : gSeq.cigars) {
			// append sequence
			if(cigar.operation == 'S' || cigar.operation == 'M' || cigar.operation == 'I') {
				for(int i=0; i<cigar.annotations.length; i++) {
					if(this.startPosInNGS <= relPos && relPos <= this.endPosInNGS) {
						genomicRegion.append(cigar.annotations[i][transcriptNum]);
					}
					relPos++;
				}
			}
		}
		
		return genomicRegion.toString();
	}
	
	/**
	 * 
	 * get genomic regional information at AA level. <br>
	 * 
	 * @param transcriptNum
	 * @return
	 */
	public String getAARegionAnnotation (int transcriptNum) {
		String genomicRegion = getGenomicRegion(transcriptNum);
		StringBuilder aaRegionAnnotation = new StringBuilder();
		
		for(int i=0; i<genomicRegion.length(); i+=3) {
			char ntMark1 = genomicRegion.charAt(i);
			char ntMark2 = genomicRegion.charAt(i+1);
			char ntMark3 = genomicRegion.charAt(i+2);
			
			aaRegionAnnotation.append(Codon.getAARegion(ntMark1, ntMark2, ntMark3));
		}
		
		
		if(!strand) {
			aaRegionAnnotation = aaRegionAnnotation.reverse();
		}
		// region shrinkage
		StringBuilder shortAnnotation = new StringBuilder();
		int count = 1;
		char mark = aaRegionAnnotation.charAt(0);
		for(int i=1; i<aaRegionAnnotation.length(); i++) {
			char nextMark = aaRegionAnnotation.charAt(i);
			if(mark == nextMark) count++;
			else {
				shortAnnotation.append(count).append(mark);
				mark = nextMark;
				count = 1;
			}
		}
		shortAnnotation.append(count).append(mark);
		
		return shortAnnotation.toString();
		
	}
	/**
	 * Genomic positions of mapped peptides <br>
	 * Plus, frame-decision <br>
	 * 
	 */
	public void mapGenomicAnnotation () {
		// genomic position setting
		this.startGenomicPosition = -1;
		this.endGenomicPosition = -1;
		
		int relPos = 0;
		for(Cigar cigar : gSeq.cigars) {
			// append sequence
			if(cigar.operation == 'S' || cigar.operation == 'M' || cigar.operation == 'I') {
				for(int i=0; i<cigar.annotations.length; i++) {
					if(this.startPosInNGS <= relPos && relPos <= this.endPosInNGS) {
						// there is no relativePosition for soft-clip!
						if(cigar.operation == 'S') return;
						
						if(startGenomicPosition == -1) startGenomicPosition = cigar.relativePositions[i] + gSeq.startPosition;
						if(endGenomicPosition == -1) endGenomicPosition = cigar.relativePositions[i] + gSeq.startPosition;
						
						startGenomicPosition = Math.min(startGenomicPosition, cigar.relativePositions[i] + gSeq.startPosition);
						endGenomicPosition = Math.max(endGenomicPosition, cigar.relativePositions[i] + gSeq.startPosition);
					}
					relPos++;
				}
			}
		}
	}
}

