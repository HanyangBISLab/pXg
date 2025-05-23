package progistar.pXg.data;

import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMTag;
import progistar.pXg.constants.Constants;
import progistar.pXg.constants.Parameters;
import progistar.pXg.utils.Codon;
import progistar.pXg.utils.IndexConvertor;

/**
 * Maintain sequence read information: <br>
 * UniqueID: qName <br>
 * chrIndex: index of chr string, which is an auto-increment byte generated from IndexConvertor <br>
 * startPosition: start position of mapped read. 1-based index <br>
 * ArrayList<Cigar>: list of cigars. <br>
 *
 * @author gistar
 *
 */
public class GenomicSequence {

	/**
	 * Note that about MD tag.
	 * MD tag only presents about SNP and DEL reference sequences.
	 * So, basically Cigar 'M' string is only considered (in case of DEL, we can infer from 'M').
	 * Therefore, MD tag must be resolved by Cigar 'M' only! (do not use other Cigars to resolve MD tag).
	 *
	 *
	 */
	private static final Pattern EACH_MD_REGEX = Pattern.compile("(([0-9]+)|([A-Z]+|\\^[A-Z]+))");
	private static final Pattern EACH_CIGAR_REGEX = Pattern.compile("([0-9]+)([MINDSHPX=])");

	public String uniqueID;
	public int chrIndex;
	public int startPosition;
	public int endPosition;
	public String mdString;
	public ArrayList<Cigar> cigars;
	public double meanQScore;
	public int flags;

	// annotation
	// if there is no model on the mapped region, then it means this read maps on intergenic region only.
	// => in this case, tBlocks has only one element with null.
	// otherwise, there is at least one annotated model.
	// => in this case, tBlocks has annotated transcript models and fill with their object.
	// The assignment, see Mapper Class
	public TBlock[]	tBlocks; // we can get strand from tBlock!
	public int 		matchedTxds = 1; // in the case of only mapping to intergenic, the number is 1.

	// default is three frame translation.
	// it has a specific translation frame as if it can infer from annotation.

	//public byte[] 	transFrames; // byte[sizeOfTranscripts]
	// In case of intergenic (it implies that this sequence cannot be explained by annotations), transFrames = new byte[1].

	public GenomicSequence (SAMRecord read) {
		
		// TODO: unmapped reads cannot have genomic position information.
		// In case of unmapped read,  must consider it!
		String phred33 = read.getBaseQualityString();
		// average of phred33 QScore
		int length = phred33.length();
		this.meanQScore = 0;
		for(int i=0; i<length; i++) {
			char qChar = phred33.charAt(i);
			if(qChar != '*') {
				this.meanQScore += (qChar-33);
			}
		}
		this.meanQScore /= length;

		// Note that
		// Chr of unmapped reads are marked as *
		// From this, we can recognize unmapped reads

		// Cigar has nucleotides and relative positions to the start position.
		this.cigars = parseCigarString(read.getCigarString(), read.getReadString());
		
		// find MD string
		Object mdStr = read.getAttribute(SAMTag.MD);
		if(mdStr == null) {
			this.mdString = null;
		} else {
			this.mdString = (String) mdStr;
		}
		
		this.chrIndex = IndexConvertor.chrToIndex(read.getReferenceName());
		this.uniqueID = read.getReadName();
		this.startPosition = read.getAlignmentStart();
		this.endPosition = startPosition;
		this.flags = read.getFlags();

		for(Cigar cigar : this.cigars) {
			char op = cigar.operation;

			switch (op) {
	    	case 'M': case 'I':// match or mismatch or Insertion
	    		this.endPosition = Math.max(this.endPosition, this.startPosition+cigar.relativePositions[cigar.relativePositions.length-1]);
	    		break;
	    	case '*': // for unmapped
	    		this.startPosition = 1;
	    		this.endPosition = this.startPosition + cigar.markerSize - 1;
	    		break;
    		default :
    			break;
	    	}
		}
	}

	public void setNonTranscripts () {
		this.matchedTxds = 1;
		this.tBlocks = new TBlock[1];
		this.tBlocks[0] = null;

		for(Cigar cigar : this.cigars) {
			char op = cigar.operation;

			switch (op) {
	    	case 'M': // match or mismatch
	    		cigar.annotations = new char[cigar.relativePositions.length][1];

	    		for(int i=0; i<cigar.annotations.length; i++) {
	    			cigar.annotations[i][0] = Constants.MARK_INTERGENIC;
	    		}

	    		break;

	    	case 'I': // insertion
	    		cigar.annotations = new char[cigar.relativePositions.length][1];

	    		for(int i=0; i<cigar.annotations.length; i++) {
	    			cigar.annotations[i][0] = Constants.MARK_INTERGENIC;
	    		}

	    		break;

	    	case 'D': // deletion
	    		break;

	    	case 'N': // skip (ex> exon junction)
	    		break;

	    	case '*': // unmapped
	    		cigar.annotations = new char[cigar.relativePositions.length][1];

	    		for(int i=0; i<cigar.annotations.length; i++) {
	    			cigar.annotations[i][0] = Constants.MARK_UNMAPPED;
	    		}

	    		break;

	    	case 'S': // soft-clip
	    		cigar.annotations = new char[cigar.relativePositions.length][1];

	    		for(int i=0; i<cigar.annotations.length; i++) {
	    			cigar.annotations[i][0] = Constants.MARK_SOFTCLIP;
	    		}

	    		break;
    		default :
    			break;
	    	}
		}
	}

	/**
	 * Check unmapped/softclip status by cigar string.<br>
	 *
	 *
	 * If cigars have '*', then it returns 'MARK_UNMAPPED' <br>
	 * If cigars have 'S', then it returns 'MARK_SOFTCLIP' <br>
	 * Otherwise, it returns 'MARK_MAPPED'
	 *
	 * @return
	 */
	public char getMappingStatus () {
		assert this.cigars.size() != 0;

		if(this.cigars.get(0).operation == '*') {
			return Constants.MARK_UNMAPPED;
		}

		for(Cigar cigar : this.cigars) {
			if(cigar.operation == 'S') {
				return Constants.MARK_SOFTCLIP;
			}
		}

		return Constants.MARK_MAPPED;
	}

	public String getLocus () {
		return IndexConvertor.indexToChr(chrIndex) +":" +this.startPosition+"-"+this.endPosition;
	}

	public String getNucleotideString () {
		StringBuilder nucleotides = new StringBuilder();
		for(Cigar cigar : cigars) {
			// append sequence
			nucleotides.append(cigar.nucleotides);
		}
		return nucleotides.toString();
	}
	/**
	 * [start, end] zero-based. <br>
	 *
	 * @param start
	 * @param end
	 * @return
	 */
	public ArrayList<Mutation> getMutationsByPositionInNGS (int start, int end) {
		ArrayList<Mutation> allMutations = new ArrayList<>();
		ArrayList<Mutation> inMutations = new ArrayList<>();

		// decoy... has no md string.
		// OR there is no available MD tag: because of alignment option.
		// In the case of STAR2, add --outSAMattributes MD
		boolean isMDtag = mdString == null ? false : true;
		if(!isMDtag) {
			return inMutations;
		}

		// MD parsing
		Matcher mdMatcher = EACH_MD_REGEX.matcher(mdString);
		int mRelPos = 0;

		while(mdMatcher.find()) {
			String md = mdMatcher.group();
			char sign = md.charAt(0);

			// match size
			if(Character.isDigit(sign)) {
				mRelPos += Integer.parseInt(md);
			}
			// nt change
			else if(Character.isAlphabetic(sign)) {
				for(int i=0; i<md.length(); i++) {
					mRelPos++;
					Mutation mutation = new Mutation();
					mutation.relPos = mRelPos - 1; // to zero-based
					mutation.refSeq = md.charAt(i)+"";
					mutation.type = Constants.SNP;
					allMutations.add(mutation);
				}
			}
			// deletion sequence
			else if(sign == '^') {
				Mutation mutation = new Mutation();
				mutation.relPos = mRelPos; // it means that the next position of the current position as a zero-based
				mutation.refSeq = md.substring(1);
				mutation.type = Constants.DEL;
				allMutations.add(mutation);
			}
		}

		int relPos = 0;
		mRelPos = 0;
		for(Cigar cigar : this.cigars) {
			if(cigar.operation == 'M') {
				for(int i=0; i<cigar.relativePositions.length; i++) {
					if(start <= relPos && relPos <= end) {

						for(int j=0; j<allMutations.size(); j++) {
							if(allMutations.get(j).relPos == mRelPos && allMutations.get(j).type == Constants.SNP) {
								allMutations.get(j).relPos = relPos; // because of softclip
								allMutations.get(j).altSeq = cigar.nucleotides.charAt(i) +"";
								allMutations.get(j).chrIndex = this.chrIndex;
								allMutations.get(j).genomicPosition = this.startPosition + cigar.relativePositions[i];
								inMutations.add(allMutations.get(j));
								allMutations.remove(j);
							}
						}

					}
					mRelPos++;
					relPos++;
				}
			} else if(cigar.operation == 'I') {
				boolean isIncluded = false;
				for(int i=0; i<cigar.relativePositions.length; i++) {
					if(!isIncluded) {
						if(start <= relPos && relPos <= end) {
							Mutation mutation = new Mutation();
							mutation.altSeq = cigar.nucleotides;
							mutation.chrIndex = this.chrIndex;
							mutation.relPos = relPos;
							// the relative position of insertion is shifted by + 1 when parsing Cigar.
							// this is a tiny issue, so just shift by -1.
							mutation.genomicPosition = this.startPosition + cigar.relativePositions[0] -1;
							mutation.type = Constants.INS;
							inMutations.add(mutation);
							isIncluded = true;
						}
					}
					relPos++;
				}
			} else if(cigar.operation == 'D') {
				for(int i=0; i<cigar.relativePositions.length; i++) {
					if(start < relPos && relPos <= end) {
						for(int j=0; j<allMutations.size(); j++) {
							if(allMutations.get(j).relPos == mRelPos && allMutations.get(j).type == Constants.DEL) {
								allMutations.get(j).relPos = relPos; // because of softclip
								allMutations.get(j).chrIndex = this.chrIndex;
								allMutations.get(j).genomicPosition = this.startPosition + cigar.relativePositions[0];
								inMutations.add(allMutations.get(j));
								allMutations.remove(j);
							}
						}
					}
				}
			} else if(cigar.operation == 'S') {
				relPos += cigar.markerSize;
			}
		}

		return inMutations;
	}

	public String getGenomieRegion (int transcriptNum) {
		StringBuilder genomicRegion = new StringBuilder();
		for(Cigar cigar : cigars) {
			// append sequence
			if(cigar.operation == 'M' || cigar.operation == 'I') {
				for (char[] annotation : cigar.annotations) {
					genomicRegion.append(annotation[transcriptNum]);
				}
			}
		}
		return genomicRegion.toString();
	}

	public String getForwardStrandTranslation (int frame) {
		if(Parameters.leucineIsIsoleucine) {
			return GenomicSequence.translation(getNucleotideString(), frame).replace("I","L");
		} else {
			return GenomicSequence.translation(getNucleotideString(), frame);
		}
	}

	public String getReverseStrandTranslation (int frame) {
		if(Parameters.leucineIsIsoleucine) {
			return GenomicSequence.reverseComplementTranslation(this.getNucleotideString(), frame).replace("I", "L");
		} else {
			return GenomicSequence.reverseComplementTranslation(this.getNucleotideString(), frame);
		}
	}
	/**
	 * frame is a start position. This is zero-base.
	 *
	 * @param nucleotides
	 * @param frame
	 * @return
	 */
	public static String translation (String nucleotides, int frame) {
		StringBuilder peptides = new StringBuilder();
		int length = nucleotides.length();
		for(int position=frame; position<length-2; position+=3) {
			char aa = Codon.nuclToAmino(nucleotides.substring(position,position+3));
			peptides.append(aa);
		}
		return peptides.toString();
	}

	public static String reverseComplementTranslation (String nucleotides, int frame) {
		StringBuilder peptides = new StringBuilder();
		StringBuilder reverseComplementNTs = new StringBuilder(nucleotides);
		int length = nucleotides.length();
		for(int i=0; i<length; i++) {
			switch(reverseComplementNTs.charAt(i)) {
				case 'A': reverseComplementNTs.setCharAt(i, 'T'); break;
				case 'C': reverseComplementNTs.setCharAt(i, 'G'); break;
				case 'T': reverseComplementNTs.setCharAt(i, 'A'); break;
				case 'G': reverseComplementNTs.setCharAt(i, 'C'); break;
				default : break;
			}
		}

		nucleotides = reverseComplementNTs.reverse().toString();
		for(int position=frame; position<nucleotides.length()-2; position+=3) {
			char aa = Codon.nuclToAmino(nucleotides.substring(position,position+3));
			peptides.append(aa);
		}
		return peptides.toString();
	}
	
	private static ArrayList<Cigar> parseCigarString (String cigarString, String nucleotides) {
		Matcher matcher = EACH_CIGAR_REGEX.matcher(cigarString);
		ArrayList<Cigar> results = new ArrayList<>();
		ArrayList<Cigar> filterResults = new ArrayList<>();

		// unmapped reads have no matcher...!
	    while (matcher.find()) {
	      int markerSize = Integer.parseInt(matcher.group(1));
	      char operation = matcher.group(2).charAt(0);

	      results.add(new Cigar(markerSize, operation));
	    }

	    // Unmapped read checker
	    if(results.size() == 0 && cigarString.equalsIgnoreCase("*")) {
	    	// add unmapped read cigar
	    	results.add(new Cigar(nucleotides.length(), '*'));
	    }

	    int ntIndex = 0;
	    int relPos = 0;
	    int[] relativePositions = null;
	    // drop all cigars without MIND
	    for(int i=0; i<results.size(); i++) {
	    	Cigar cigar = results.get(i);
	    	char op = cigar.operation;

	    	switch (op) {
	    	case 'M': // match or mismatch
	    		cigar.nucleotides = nucleotides.substring(ntIndex, ntIndex + cigar.markerSize);
	    		ntIndex += cigar.markerSize;

	    		relativePositions = new int[cigar.markerSize];
	    		for(int j=0; j<relativePositions.length; j++) {
	    			relativePositions[j] = relPos++;
	    		}

	    		cigar.relativePositions = relativePositions;
	    		filterResults.add(cigar);
	    		break;

	    	case 'S': // soft clip
	    		cigar.nucleotides = nucleotides.substring(ntIndex, ntIndex + cigar.markerSize);
	    		ntIndex += cigar.markerSize;

	    		relativePositions = new int[cigar.markerSize];
	    		for(int j=0; j<relativePositions.length; j++) {
	    			relativePositions[j] = relPos; // relPos is not changed... consistent!
	    		}

	    		cigar.relativePositions = relativePositions;
	    		filterResults.add(cigar);
	    		break;

	    	case 'I': // insertion
	    		cigar.nucleotides = nucleotides.substring(ntIndex, ntIndex + cigar.markerSize);
	    		ntIndex += cigar.markerSize;

	    		relativePositions = new int[cigar.markerSize];
	    		for(int j=0; j<relativePositions.length; j++) {
	    			relativePositions[j] = relPos; // relPos is not changed... consistent!
	    		}

	    		cigar.relativePositions = relativePositions;
	    		filterResults.add(cigar);
	    		break;

	    	case 'D': // deletion
	    		relativePositions = new int[cigar.markerSize];
	    		for(int j=0; j<relativePositions.length; j++) {
	    			relativePositions[j] = relPos++;
	    		}

	    		cigar.relativePositions = relativePositions;
	    		filterResults.add(cigar);
	    		break;

	    	case 'N': // skip (ex> exon junction)
	    		relPos += cigar.markerSize;
	    		filterResults.add(cigar);
	    		break;

	    	case '*': // unmapped
	    		cigar.nucleotides = nucleotides;
	    		relativePositions = new int[cigar.markerSize];
	    		for(int j=0; j<relativePositions.length; j++) {
	    			relativePositions[j] = relPos++;
	    		}
	    		cigar.relativePositions = relativePositions;

	    		filterResults.add(cigar);
	    		break;
	    	}
	    }

	    return filterResults;
	}
	

	public ArrayList<Character> getStrandedness () {
		ArrayList<Character> strands = new ArrayList<Character>();
		boolean isFirstSegment = (0x40 & this.flags) == 0x40 ? true : false;
		boolean isForward = (0x10 & this.flags) == 0x10 ? false : true;
		
		// non-stranded
		if(Parameters.strandedness.equalsIgnoreCase(Constants.NON_STRANDED)) {
			strands.add('+');
			strands.add('-');
		}  
		// Single-end
		else if(Parameters.strandedness.equalsIgnoreCase(Constants.F_STRANDED)) {
			if(isForward) {
				strands.add('+');
			} else {
				strands.add('-');
			}
		} else if(Parameters.strandedness.equalsIgnoreCase(Constants.R_STRANDED)) {
			if(isForward) {
				strands.add('-');
			} else {
				strands.add('+');
			}
		} 
		// Paired-end
		else {
			// R1
			if(isFirstSegment) {
				if(Parameters.strandedness.equalsIgnoreCase(Constants.FR_STRANDED)) {
					if(isForward) {
						strands.add('+');
					} else {
						strands.add('-');
					}
				} else if(Parameters.strandedness.equalsIgnoreCase(Constants.RF_STRANDED)) {
					if(isForward) {
						strands.add('-');
					} else {
						strands.add('+');
					}
				}
			} 
			// R2
			else {
				if(Parameters.strandedness.equalsIgnoreCase(Constants.FR_STRANDED)) {
					if(isForward) {
						strands.add('-');
					} else {
						strands.add('+');
					}
				} else if(Parameters.strandedness.equalsIgnoreCase(Constants.RF_STRANDED)) {
					if(isForward) {
						strands.add('+');
					} else {
						strands.add('-');
					}
				} 
			}
		}
		
		return strands;
	}
	
}
