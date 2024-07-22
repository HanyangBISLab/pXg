package progistar.pXg.data;

import java.util.ArrayList;

import progistar.pXg.constants.Constants;
import progistar.pXg.constants.Parameters;

public class ConsensusSequence {
	
	public String leftFlankConsensus = Constants.ID_NULL;
	public String rightFlankConsensus = Constants.ID_NULL;
	public String leftFlankRefConsensus = Constants.ID_NULL;
	public String rightFlankRefConsensus = Constants.ID_NULL;
	
	public int getCharToIdx (char n) {
		if(n == 'A') {
			return 0;
		} else if(n == 'C') {
			return 1;
		} else if(n == 'T') {
			return 2;
		} else if(n == 'G') {
			return 3;
		} else if(n == 'a') {
			return 4;
		} else if(n == 'c') {
			return 5;
		} else if(n == 't') {
			return 6;
		} else if(n == 'g') {
			return 7;
		} else if(n == Constants.ID_NULL.charAt(0)) {
			return 8;
		} else {
			return 9;
		}
	}
	
	public char getIdxToChar (int idx) {
		if(idx == 0) {
			return 'A';
		} else if(idx == 1) {
			return 'C';
		} else if(idx == 2) {
			return 'T';
		} else if(idx == 3) {
			return 'G';
		} else if(idx == 4) {
			return 'a';
		} else if(idx == 5) {
			return 'c';
		} else if(idx == 6) {
			return 't';
		} else if(idx == 7) {
			return 'g';
		} else if(idx == 8) {
			return Constants.ID_NULL.charAt(0);
		} else {
			return 'N';
		}
	}
	
	private String[] getConsensusSequence(ArrayList<XBlock> list, char type) {

		StringBuilder leftConsensusSequence = new StringBuilder();
		StringBuilder rightConsensusSequence = new StringBuilder();
		
		int maxLeftSize = 0;
		int maxRightSize = 0;
		for(XBlock xBlock_ : list) {
			String leftFlankSequence = null;
			String rightFlankSequence = null;

			if(type == 'R') {
				leftFlankSequence = Global.SEQUENCE_ARRAYLIST.get(xBlock_.leftFlankRefSequenceIdx);
				rightFlankSequence = Global.SEQUENCE_ARRAYLIST.get(xBlock_.rightFlankRefSequenceIdx);
			} else {
				leftFlankSequence = Global.SEQUENCE_ARRAYLIST.get(xBlock_.leftFlankSequenceIdx);
				rightFlankSequence = Global.SEQUENCE_ARRAYLIST.get(xBlock_.rightFlankSequenceIdx);
			}
			
			maxLeftSize = Math.max(leftFlankSequence.length(), maxLeftSize);
			maxRightSize = Math.max(rightFlankSequence.length(), maxRightSize);
		}

		if(maxLeftSize > Parameters.maxFlankNSize) {
			maxLeftSize = Parameters.maxFlankNSize;
		}
		if(maxRightSize > Parameters.maxFlankNSize) {
			maxRightSize = Parameters.maxFlankNSize;
		}

		int[][] leftConsensusScore = new int[maxLeftSize][10];
		int[][] rightConsensusScore = new int[maxRightSize][10];


		// calculate score matrix
		for(XBlock xBlock_ : list) {
			// left
			String sequence = null;
			if(type == 'R') {
				sequence = Global.SEQUENCE_ARRAYLIST.get(xBlock_.leftFlankRefSequenceIdx);
			} else {
				sequence = Global.SEQUENCE_ARRAYLIST.get(xBlock_.leftFlankSequenceIdx);
			}
			
			int idx = maxLeftSize-1;
			for(int i=sequence.length()-1; i>=0; i--) {
				char nt = sequence.charAt(i);
				int ntIdx = getCharToIdx(nt);
				leftConsensusScore[idx][ntIdx] ++;
				idx--;
			}

			// right
			sequence = null;
			if(type == 'R') {
				sequence = Global.SEQUENCE_ARRAYLIST.get(xBlock_.rightFlankRefSequenceIdx);
			} else {
				sequence = Global.SEQUENCE_ARRAYLIST.get(xBlock_.rightFlankSequenceIdx);
			}
			
			for(int i=0; i<sequence.length(); i++) {
				char nt = sequence.charAt(i);
				int ntIdx = getCharToIdx(nt);
				rightConsensusScore[i][ntIdx] ++;
			}
		}
		
		// assign a consensus sequence (left flank)
		for(int i=0; i<maxLeftSize; i++) {
			int bestIdx = 9;
			for(int j=0; j<9; j++) {
				if(leftConsensusScore[i][j] > leftConsensusScore[i][bestIdx]) {
					bestIdx = j;
				}
			}
			leftConsensusSequence.append(getIdxToChar(bestIdx));
		}
		
		for(int i=0; i<maxRightSize; i++) {
			int bestIdx = 9;
			for(int j=0; j<9; j++) {
				if(rightConsensusScore[i][j] > rightConsensusScore[i][bestIdx]) {
					bestIdx = j;
				}
			}
			rightConsensusSequence.append(getIdxToChar(bestIdx));
		}
		
		if(leftConsensusSequence.length() == 0) {
			leftConsensusSequence.append(Constants.ID_NULL);
		}
		if(rightConsensusSequence.length() == 0) {
			rightConsensusSequence.append(Constants.ID_NULL);
		}
		
		return new String[] {leftConsensusSequence.toString(), rightConsensusSequence.toString()};
	}
	
	public ConsensusSequence (XBlock xBlock) {
		ArrayList<XBlock> list = new ArrayList<>();
		list.add(xBlock);
		list.addAll(xBlock.siblingXBlocks);
		
		// for NGS-based consensus sequence
		String[] lr = getConsensusSequence(list, 'N');
		this.leftFlankConsensus = lr[0];
		this.rightFlankConsensus = lr[1];
		
		// for reference-based consensus sequence
		lr = getConsensusSequence(list, 'R');
		this.leftFlankRefConsensus = lr[0];
		this.rightFlankRefConsensus = lr[1];
	}
}
