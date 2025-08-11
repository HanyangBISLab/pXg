package progistar.pXg;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import progistar.pXg.data.PxGAnnotation;
import progistar.pXg.data.parser.ResultParser;

public class Test {

	public static void main(String[] args) throws IOException {
		ArrayList<File>  tmpOutputFiles = new ArrayList<File>();
		
		tmpOutputFiles.add(new File("/Users/schoi/Documents/project/neoflow2/Error_reports/20250804_pXgMemory/C3L-06174.T.markdup.sorted.bam.1754889873214.worker1.tmp"));
		
		PxGAnnotation pXgA = ResultParser.parseResult(tmpOutputFiles);
	}
}
