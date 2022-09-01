package progistar.thirdparty.sSim;

import java.io.IOException;
import java.util.ArrayList;

public class FindPeakMatch {

	public static void main(String[] args) throws IOException {
		PeptideLoader.loadPeptideList("/Users/gistar/projects/pXg/Laumont_NatCommun2016/Results/9.Unmodified_10ppm/ProteomToolsHLA/Matched_NCPeptideList.tsv");
		ProteomTools.loadProteomeTools("/Volumes/One Touch/pxgSelectedRes", 
				"/Volumes/One Touch/pxgSelected", 
				PeptideLoader.matchedPeptides);
		
		ArrayList<Spectra> list = ProteomTools.spectraList;
		System.out.println(list.size()+" were retrived");
		int[] msLevel = {2};
		for(Spectra spectra : list) {
			System.out.println(spectra.getFileName()+": " +spectra.sizeOfSpectra());
			spectra.writeToMGF(spectra.getFileName().replace(".mgf", "selected.mgf"), msLevel);
		}
	}
}
