package progistar.pXg.utils;

import java.io.File;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import progistar.pXg.constants.Constants;
import progistar.pXg.constants.Parameters;

public class CheckStrandedness {
	
	private static int R1F = 0;
	private static int R2F = 0;
	private static int R1R = 0;
	private static int R2R = 0;
	
	public static byte detect (File file) {
		long startTime = System.currentTimeMillis();
		// to prevent racing
		try (SamReader samReader = SamReaderFactory.makeDefault().open(file)) {
			SAMRecordIterator iterator = samReader.iterator();
			int size = 0;
			
			while((size++) < 1000000 && iterator.hasNext()) {
				SAMRecord samRecord = iterator.next();
				
				boolean isPass = false;
				if(Parameters.COUNT_PRIMARY_ONLY && samRecord.isSecondaryAlignment()) {
            		isPass = true;
            	}
				
				Object xsTag = samRecord.getAttribute("XS");
				if(xsTag == null) {
					isPass = true;
				}
				
				if(isPass) {
					continue;
				}
				
				
				int flags = samRecord.getFlags();
				boolean isFirstSegment = (0x40 & flags) == 0x40 ? true : false;
				boolean isForward = (0x10 & flags) == 0x10 ? false : true;
				boolean strand = ((Character) xsTag) == '+' ? true : false;
				
				// first segment
				if(isFirstSegment) {
					if(isForward == strand) {
						R1F++;
					} else {
						R1R++;
					}
				} 
				// second segment
				else {
					if(isForward == strand) {
						R2F++;
					} else {
						R2R++;
					}
				}
				
			}
			
			iterator.close();
		} catch(Exception e) {
			e.printStackTrace();
			System.exit(1);
		}
		
		if(R1F*10 < R1R && R2F > R2R*10) {
			Parameters.strandedness = Constants.RF_STRANDED;
		} else if(R1F > R1R*10 && R2F*10 < R2R) {
			Parameters.strandedness = Constants.FR_STRANDED;
		} else {
			Parameters.strandedness = Constants.NON_STRANDED;
		}
		System.out.println("Estimate strandedness");
		System.out.println("1F\t1R\t2F\t2R");
		System.out.println(R1F+"\t"+R1R+"\t"+R2F+"\t"+R2R);
		
		if(R1F+R1R+R2F+R2R == 0) {
			System.out.println("Fail to estimate stradedness!");
			System.out.println("It looks single-end RNA-seq experiement. Please specify strandedness.");
			System.exit(1);
		} else {
			System.out.println("Strandedness: "+Parameters.strandedness+"-stranded");
		}
		
		long endTime = System.currentTimeMillis();
		System.out.println("Check strandedness: "+(endTime-startTime)/1000+" sec");
		
		return 0;
	}
}
