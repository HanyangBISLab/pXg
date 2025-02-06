package progistar.pXg.utils;

import progistar.pXg.constants.Parameters;
import progistar.pXg.constants.RunInfo;

public class ReadNormalizer {

	// RPHM = read per hundred million.
	private static double factor = Math.pow(10, 8);
	
	public static double normalize (int count) {
		double nCount = count;
		
		if(Parameters.COUNT_NORMALIZATION) {
			nCount = (nCount/RunInfo.totalProcessedReads)*factor;
		}
		
		return nCount;
	}
}
