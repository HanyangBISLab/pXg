package progistar.pXg;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import progistar.pXg.data.PxGAnnotation;
import progistar.pXg.data.parser.ResultParser;

public class Test {

	public static void main(String[] args) throws IOException {
		ArrayList<File>  tmpOutputFiles = new ArrayList<File>();
		
//		extract(new File("/Users/schoi/Documents/project/neoflow2/Error_reports/20250804_pXgMemory/C3L-06174.T.markdup.sorted.bam"));
		
		tmpOutputFiles.add(new File("/Users/schoi/Documents/project/neoflow2/Error_reports/20250804_pXgMemory/C3L-06174.T.markdup.sorted.bam.1754901466035.worker7.tmp"));
		
		
		PxGAnnotation pXgA = ResultParser.parseResult(tmpOutputFiles);
	}
	
	public static void extract(File bamFile) throws IOException {
		String[] ids = {"A00434:678:HCJ3CDSX7:4:1114:20003:30342"};
        
        // Open BAM file
        SamReaderFactory factory = SamReaderFactory.makeDefault();
        try (SamReader reader = factory.open(bamFile)) {
        	
        	SAMRecordIterator it = reader.queryUnmapped();
        	
        	while(it.hasNext()) {
        		SAMRecord record = it.next();
        		String readName = record.getReadName();
            	
            	for(String id : ids) {
            		if(id.equalsIgnoreCase(readName)) {
            			System.out.println(record.getSAMString());
            		}
            	}
        	}
        } catch (Exception e) {
            e.printStackTrace();
        }
		
	}
}
