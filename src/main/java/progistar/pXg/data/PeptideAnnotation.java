package progistar.pXg.data;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Hashtable;
import java.util.Iterator;

import org.ahocorasick.trie.Emit;
import org.ahocorasick.trie.Trie;

import progistar.pXg.constants.Parameters;
import progistar.pXg.constants.RunInfo;

public class PeptideAnnotation {


	private static String[]			fields;
	private static Trie				trie;
	public static ArrayList<PBlock>	pBlocks = new ArrayList<>();
	public static Hashtable<String, Integer>	peptideIndexer	=	new Hashtable<>();
	public static Hashtable<Integer, String>	indexedPeptide	=	new Hashtable<>();

	private PeptideAnnotation() {}

	public static void setFields (String[] fields) {
		PeptideAnnotation.fields = fields;
	}

	public static String toFields () {
		StringBuilder fieldStr = new StringBuilder();

		for(int i=0; i<fields.length; i++) {
			if(i!=0) {
				fieldStr.append("\t");
			}
			fieldStr.append(fields[i]);
		}

		return fieldStr.toString();
	}

	public static void buildKeywordTrie () {
		System.out.println("Enumerate peptide sequences...");
		ArrayList<String> sequences = enumerateSequence();

		RunInfo.totalProcessedPeptides = sequences.size();

		System.out.println("A total of "+sequences.size()+" peptides were detected without duplications.");
		System.out.println("Build keyword trie...");
		trie = Trie.builder().addKeywords(sequences).build();
		System.out.println("Done!");
	}



	public static ArrayList<Output> find (GenomicSequence gSeq) {
		ArrayList<Output> outputs = new ArrayList<>();
		int ntLen = gSeq.getNucleotideString().length();
		ArrayList<Character> strandedness = gSeq.getStrandedness();
		
		StringBuilder decoySequence = new StringBuilder();
		// forward translation
		for(char strand : strandedness) {
			if(strand == '+') {
				for(int i=0; i<3; i++) {
					decoySequence.setLength(0);
					String target = gSeq.getForwardStrandTranslation(i);
					String decoy = decoySequence.append(target).reverse().toString();
					int decoyLen = decoy.length();


					// for target
					Collection<Emit> emits = trie.parseText(target);
					for(Emit emit : emits) {
						// convert peptide index to nucleotide index
						int start = emit.getStart() * 3 + i;
						int end = (emit.getEnd()+1) * 3 + i - 1;
						Output output = new Output(gSeq, peptideIndexer.get(emit.getKeyword()), start, end, true, true);
						outputs.add(output);
					}

					// for decoy
					if(!Parameters.printTargetOnly) {
						emits = trie.parseText(decoy);
						for(Emit emit : emits) {
							// convert peptide index to nucleotide index
							int tmpStart = ((decoyLen - 1) - emit.getEnd());
							int tmpEnd = ((decoyLen-1) - emit.getStart());
							
							int start = tmpStart * 3 + i;
							int end = (tmpEnd+1) * 3 + i - 1;
							
							Output output = new Output(gSeq, peptideIndexer.get(emit.getKeyword()), start, end, true, false);
							outputs.add(output);
						}
					}
				}
			} else if(strand == '-') {
				// reverse complement translation
				for(int i=0; i<3; i++) {
					decoySequence.setLength(0);
					String target = gSeq.getReverseStrandTranslation(i);
					String decoy = decoySequence.append(target).reverse().toString();
					int decoyLen = decoy.length();

					// for target
					Collection<Emit> emits = trie.parseText(target);
					for(Emit emit : emits) {
						// convert peptide index to nucleotide index
						int start = emit.getStart() * 3 + i;
						int end = (emit.getEnd()+1) * 3 + i - 1;

						// convert reverse index to forward index
						int tmp = start;
						start = ntLen - end - 1;
						end = ntLen - tmp - 1;

						Output output = new Output(gSeq, peptideIndexer.get(emit.getKeyword()), start, end, false, true);
						outputs.add(output);
					}

					// for decoy
					if(!Parameters.printTargetOnly) {
						emits = trie.parseText(decoy);
						for(Emit emit : emits) {
							// convert peptide index to nucleotide index
							int tmpStart = ((decoyLen - 1) - emit.getEnd());
							int tmpEnd = ((decoyLen-1) - emit.getStart());
							
							int start = tmpStart * 3 + i;
							int end = (tmpEnd+1) * 3 + i - 1;
							
							// convert reverse index to forward index
							int tmp = start;
							start = ntLen - end - 1;
							end = ntLen - tmp - 1;
							
							Output output = new Output(gSeq, peptideIndexer.get(emit.getKeyword()), start, end, false, false);
							outputs.add(output);
						}
					}
				}

			}
		}
		

		

		return outputs;
	}

	/**
	 * Get non-duplicated peptide sequences from PeptideAnnotation.<br>
	 *
	 *
	 * @return
	 */
	public static ArrayList<String> enumerateSequence () {
		ArrayList<String> sequences = new ArrayList<>();

		indexedPeptide.clear();
		peptideIndexer.clear();

		// put peptide sequences into checks
		pBlocks.forEach(pBlock ->
			{
				if(peptideIndexer.get(pBlock.getPeptideSequence()) == null) {
					indexedPeptide.put(peptideIndexer.size()+1, pBlock.getPeptideSequence());
					peptideIndexer.put(pBlock.getPeptideSequence(), peptideIndexer.size()+1);

				}
			}
		);

		// add to ArrayList without duplications
		peptideIndexer.forEach((k, v) -> { sequences.add(k); });

		return sequences;
	}
	/**
	 * assign rank of candidates. <br>
	 * If the scores are tied, than they are assigned the same rank.<br>
	 *
	 */
	public static void assignRank () {
		Hashtable<String, ArrayList<PBlock>> pBlocksByScan = aggregatePBlocksByScan();

		pBlocksByScan.forEach((scanID, scanPBlocks) -> {
			// sort pBlocks by scores, decreasing order.
			Collections.sort(scanPBlocks);

			// cutoff
			int size = scanPBlocks.size();
			double topScore = scanPBlocks.get(0).score;
			int rank = 0;
			for(int i=0; i<size; i++) {
				PBlock pBlock = scanPBlocks.get(i);
				pBlock.rank = ++rank;
				pBlock.deltaScore = topScore - pBlock.score;
			}
		});
	}

	/**
	 * filter by rank <br>
	 *
	 */
	public static void filter () {
		Hashtable<String, ArrayList<PBlock>> pBlocksByScan = aggregatePBlocksByScan();
		// remove original pBlocks
		pBlocks.clear();

		pBlocksByScan.forEach((scanID, scanPBlocks) -> {
			// sort pBlocks by scores, decreasing order.
			Collections.sort(scanPBlocks);

			// cutoff
			int size = scanPBlocks.size();
			for(int i=0; i<size; i++) {
				// add pBlocks which delta-score is equal or less than "score-threshold"

				// include the candidates satisfying the rank creterion
				if(scanPBlocks.get(i).rank <= Parameters.psmRank) {
					pBlocks.add(scanPBlocks.get(i));
				}
			}
		});
	}

	public static void peptideLengthFilter () {
		ArrayList<PBlock> tmpPBlocks = pBlocks;
		pBlocks = new ArrayList<>();

		int size = tmpPBlocks.size();
		for(int i=0; i<size; i++) {
			int peptLength = tmpPBlocks.get(i).getPeptideSequence().length();

			// filter-out non-interesting
			if(peptLength >= Parameters.minPeptLen && peptLength <= Parameters.maxPeptLen) {
				pBlocks.add(tmpPBlocks.get(i));
			}
		}
		tmpPBlocks.clear();
	}

	/**
	 * aggregate pblocks by scan. <br>
	 *
	 * @return
	 */
	private static Hashtable<String, ArrayList<PBlock>> aggregatePBlocksByScan () {
		Hashtable<String, ArrayList<PBlock>> pBlocksByScan = new Hashtable<>();

		// aggregate pBlocks by scanID
		pBlocks.forEach(pBlock -> {
			String scanID = pBlock.getSpecID();
			ArrayList<PBlock> scanPBlocks = pBlocksByScan.get(scanID);
			if(scanPBlocks == null) {
				scanPBlocks = new ArrayList<>();
				pBlocksByScan.put(scanID, scanPBlocks);
			}
			scanPBlocks.add(pBlock);
		});

		return pBlocksByScan;
	}

	/**
	 * return size of scans <br>
	 *
	 * @return
	 */
	public static int getScanSize () {
		return aggregatePBlocksByScan().size();
	}

	/**
	 * Count peptides with containing at least one matched experiment reads.<br>
	 *
	 *
	 * @param xBlockMapper
	 * @return
	 */
	public static int getPeptideSizeWithXBlocks (Hashtable<String, Hashtable<String, XBlock>> xBlockMapper) {
		Hashtable<String, String> peptides = new Hashtable<>();
		for (PBlock pBlock : pBlocks) {
			peptides.put(pBlock.getPeptideSequence(), "");
		}

		// if there is at least one xBlock with exp reads?
		ArrayList<String> passedPeptideIL = new ArrayList<>();

		peptides.forEach((peptideIL, NA) -> {
			boolean is = false;
			Hashtable<String, XBlock> xBlocks = xBlockMapper.get(peptideIL);

			if(xBlocks != null && xBlocks.size() != 0) {
				Iterator<String> keys = (Iterator<String>) xBlocks.keys();
				while(keys.hasNext()) {
					String key = keys.next();
					XBlock xBlock = xBlocks.get(key);
					if(xBlock != null && xBlock.targetReadCount > 0) {
						is = true;
					}
				}
			}

			// if at least one xBlock with exp read?
			if(is) {
				passedPeptideIL.add(peptideIL);
			}
		});

		return passedPeptideIL.size();
	}


	/**
	 * Count scans with containing at least one pBlock matched to experiment reads.<br>
	 *
	 *
	 * @param xBlockMapper
	 * @return
	 */
	public static int getScanSizeWithXBlocks (Hashtable<String, Hashtable<String, XBlock>> xBlockMapper) {
		Hashtable<String, ArrayList<PBlock>> pBlocksByScan = aggregatePBlocksByScan();

		// if there is at least one xBlock with exp. reads
		ArrayList<String> passedScans = new ArrayList<>();

		pBlocksByScan.forEach((scanID, pBlocks) -> {
			boolean is = false;
			for(PBlock pBlock : pBlocks) {
				String peptideIL = pBlock.getPeptideSequence();
				Hashtable<String, XBlock> xBlocks = xBlockMapper.get(peptideIL);

				if(xBlocks != null && xBlocks.size() != 0) {
					Iterator<String> keys = (Iterator<String>) xBlocks.keys();
					while(keys.hasNext()) {
						String key = keys.next();
						XBlock xBlock = xBlocks.get(key);
						if(xBlock != null && xBlock.targetReadCount > 0) {
							is = true;
						}
					}
				}
			}

			// if at least one pBlock in the scan has xBlock with exp read?
			if(is) {
				passedScans.add(scanID);
			}
		});

		return passedScans.size();
	}

	/**
	 * return size of peptides <br>
	 *
	 * @return
	 */
	public static int getPeptideSize () {
		Hashtable<String, String> peptides = new Hashtable<>();

		for (PBlock pBlock : pBlocks) {
			peptides.put(pBlock.getPeptideSequence(), "");
		}

		return peptides.size();
	}
	
	public static int getFieldLength () {
		return fields.length;
	}
}
