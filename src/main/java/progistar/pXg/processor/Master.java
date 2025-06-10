package progistar.pXg.processor;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.Vector;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import progistar.pXg.constants.Constants;
import progistar.pXg.constants.Parameters;
import progistar.pXg.constants.RunInfo;
import progistar.pXg.data.GenomicAnnotation;
import progistar.pXg.data.PeptideAnnotation;
import progistar.pXg.data.PxGAnnotation;
import progistar.pXg.data.parser.GTFParser;
import progistar.pXg.data.parser.PeptideParser;
import progistar.pXg.data.parser.ResultParser;
import progistar.pXg.utils.BAMIndex;
import progistar.pXg.utils.CheckStrandedness;
import progistar.pXg.utils.Codon;
import progistar.pXg.utils.IndexConvertor;

public class Master {

	private static GenomicAnnotation genomicAnnotation = null;
	private static int taskCount = 0;
	private static Hashtable<String, BufferedWriter> tmpOutputFilePaths = null;
	private static File SAM_FILE = null;
	private static int[] chrIndices = null;
	private static int[] startPositions = null;
	private static boolean[] assignedArray = null;
	private static SAMRecord[] reads = null;

	private Master() {}


	/**
	 * Load GTF and peptide file and ready to read SAM file <br>
	 *
	 * @param genomicAnnotationFilePath
	 * @param sequenceFilePath
	 */
	public static void ready (String genomicAnnotationFilePath, 
			String sequenceFilePath, 
			String peptideFilePath) {
		// GTF parser and Peptide parser
		genomicAnnotation = GTFParser.parseGTF(genomicAnnotationFilePath);
		PeptideParser.parseResult(peptideFilePath); // static..!
		
		tmpOutputFilePaths = new Hashtable<>();

		SAM_FILE = new File(sequenceFilePath);

		try {
			BAMIndex.index(SAM_FILE);
		} catch(IOException ioe) {
		}
		
		// for SAM - GTF associated task assignment
		// TASK-related variables
		chrIndices = new int[Parameters.readSize];
		startPositions = new int[Parameters.readSize];
		assignedArray = new boolean[Parameters.readSize];
		reads = new SAMRecord[Parameters.readSize];

		// loading Codon.
		Codon.mapping();
		
		if(Parameters.strandedness.equalsIgnoreCase(Constants.AUTO_STRANDED)) {
			CheckStrandedness.detect(SAM_FILE);
		}
	}


	/**
	 * Start to map peptides to NGS-reads <br>
	 *
	 */
	public static void run () {
		assert genomicAnnotation != null;
		// TODO:
		// Auto detection of already made index files.

		try {
			RunInfo.totalProcessedReads = 0;
			RunInfo.workerProcessedReads = new long[Parameters.nThreads+1];
			Worker[] workers = new Worker[Parameters.nThreads];

			try (SamReader samReader = SamReaderFactory.makeDefault().open(SAM_FILE)) {

				// initialize
				Vector<Task> taskQueue = new Vector<>(); // for synchronized
				int readCount = 0;
				long readPartitionSize = Parameters.readSize;

	        	// get records
	            for (SAMRecord samRecord : samReader) {
	            	// Is it primary only?
	            	if(Parameters.COUNT_PRIMARY_ONLY && samRecord.isSecondaryAlignment()) {
	            		continue;
	            	}
	            	
	                // Process each SAMRecord as needed
	            	reads[readCount] = samRecord;

	            	String chr = samRecord.getReferenceName();
	            	int startPosition = samRecord.getAlignmentStart();

					// the index for that chr is automatically assigned by auto-increment key.
					IndexConvertor.putChrIndexer(chr);
					int chrIndex_ = IndexConvertor.chrToIndex(chr);
					// check all chromosomes are well processed.
					RunInfo.processedChromosomes.put(chr, chrIndex_);

					// store chr and start positions
					chrIndices[readCount] = chrIndex_;
					startPositions[readCount] = startPosition;
					readCount ++;
					if(readCount == readPartitionSize) {
						// the array is initialized as "false"
						Arrays.fill(assignedArray, false);
						// assign tasks to workers
						while(!assignTasks(taskQueue, workers, readCount)) {
						}
						// once give the tasks, remove reads
						RunInfo.totalProcessedReads += readCount;
						readCount = 0;
					}
	            }

	            // do last tasks
	            if(readCount > 0) {
	            	Arrays.fill(assignedArray, false);
					// assign tasks to workers
					while(!assignTasks(taskQueue, workers, readCount)) {
						
					}
					// once give the tasks, remove reads
					RunInfo.totalProcessedReads += readCount;
					readCount = 0;
	            }


	        } catch (IOException e) {
	            e.printStackTrace();
	        }

			// wait for finishing all tasks from workers
			waitUntilAllWorkersDone(workers);

			// read tmp output files
			ArrayList<File> tmpOutputFiles = new ArrayList<>();
			tmpOutputFilePaths.forEach((path, tmpBW) ->{
				tmpOutputFiles.add(new File(path));
				Master.closeOutputBW(path);
			});

			PxGAnnotation pXgA = ResultParser.parseResult(tmpOutputFiles);

			// removing tmpOutputFiles
			tmpOutputFiles.forEach(file -> {file.delete();});

			// count peptides and scans matching to exp.reads
			RunInfo.mappingFilterPeptideNum3 = PeptideAnnotation.getPeptideSizeWithXBlocks(pXgA.getXBlockMapper());
			RunInfo.mappingFilterScanNum3 = PeptideAnnotation.getScanSizeWithXBlocks(pXgA.getXBlockMapper());

			// count peptides and scans after p-value
			RunInfo.pvalueFilterPeptideNum4 = PeptideAnnotation.getPeptideSizeWithXBlocks(pXgA.getXBlockMapper());
			RunInfo.pvalueFilterScanNum4 = PeptideAnnotation.getScanSizeWithXBlocks(pXgA.getXBlockMapper());

			// filter regions
			pXgA.regionScoreFilter();
			// mark fasta result
			// to distinguish ambiguous interpretation
			pXgA.markFasta();
			
			/**
			 * Filtering priority:
			 * Target > Decoy > Non-AAVariant peptide > AAVariant peptide.
			 * It means that 
			 */
			
			// marking target PSMs
			pXgA.assignXBlocks();
			// SAAV filter + TD filter
			pXgA.prioritizePeptidesInEachCandidate();

			pXgA.write(Parameters.tmpOutputFilePaths[Parameters.CURRENT_FILE_INDEX]);
		}catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static boolean assignTasks (Vector<Task> taskQueue, Worker[] workers, int readCount) {
		Task[] tasks = getTasks(readCount);
		boolean isDone = true;
		for (Task task : tasks) {
			if(task.isAssigned) {
				// add task into taskQueue
				taskQueue.add(task);
				isDone = false;
			}
		}

		while(!taskQueue.isEmpty()) {
			Task task = taskQueue.firstElement();
			taskQueue.remove(0);

			boolean isAssigned = false;

			while(!isAssigned) {
				for(int i=0; i<workers.length; i++) {
					if(workers[i] == null || !workers[i].isAlive()) {
						workers[i] = new Worker(i+1, task);
						workers[i].start();
						isAssigned = true;
						break;
					}
				}

				Thread.yield();
			}
		}

		return isDone;
	}

	/**
	 * waiting for all workers are done with their tasks.<br>
	 *
	 * @param workers
	 */
	private static void waitUntilAllWorkersDone (Worker[] workers) {
		boolean isProcessing = true;
		while(isProcessing) {
			isProcessing = false;
			for (Worker worker : workers) {
				if(worker != null) {
					isProcessing |= worker.isAlive();
				}
			}
			Thread.yield();
		}
	}

	/**
	 * Partitioning tasks. <br>
	 * @param gSeqs
	 * @param assignedArray
	 * @return
	 */
	private static Task[] getTasks (int readCount) {

		assert readCount != 0;


		Task[] tasks = new Task[Parameters.nThreads];

		// Task class contains "isAssigned" feature.
		// The default value of the feature is "false"
		// The value becomes "true" when task has something to do.
		for(int i=0; i<tasks.length; i++) {
			tasks[i] = new Task();
		}

		int gSeqSize	=	readCount;
		int chrIndex	=	0;
		int start		=	0;
		int end			=	0;

		// Setting the pivot start position information
		for(int i=0; i<gSeqSize; i++) {
			// start position of NOT treated sequence
			if(!assignedArray[i]) {
				chrIndex	=	chrIndices[i];
				start		=	startPositions[i];
				end			=	start + Parameters.partitionSize - 1;

				break;
			}
		}

		ArrayList<SAMRecord> readPartition = new ArrayList<>();

		for(int i=0; i<gSeqSize; i++) {
			// already treated sequence
			if(assignedArray[i]) {
				continue;
			}

			if(chrIndices[i] == chrIndex) {
				if(startPositions[i] >= start && startPositions[i] + Parameters.maxJunctionSize <= end) {
					assignedArray[i] = true;
					readPartition.add(reads[i]);
				}
			}

		}

		// have a task at least one.
		int partitionInSize = readPartition.size();
		if(partitionInSize != 0) {
			// Annotation index
			int[][] gIndex = genomicAnnotation.getIndexingBlocks(chrIndex, start, end);
			int taskIndex = 0;
			for(int i=0; i<partitionInSize; i++) {
				// target NGS-read
				tasks[taskIndex].samReads.add(readPartition.get(i));
				taskIndex++;
				if(taskIndex == tasks.length) {
					taskIndex = 0;
				}
			}

			for(int i=0; i<tasks.length; i++) {
				if(tasks[i].samReads.size() != 0) {
					tasks[i].isAssigned = true;
					tasks[i].genomicAnnotationIndex = gIndex;
					tasks[i].genomicAnnotation = genomicAnnotation;
					tasks[i].gIndexStart = start;
					tasks[i].taskID = ++Master.taskCount;
					tasks[i].taskType = Constants.TASK_G_MAP;

//					System.out.println(tasks[i].description());
				}
			}
		}

		return tasks;
	}

	/**
	 * Enroll temporary output file path. <br>
	 * The enrolled file paths will be processed when making the final output file. <br>
	 *
	 * @param outputFilePath
	 */
	public static void enrollTmpOutputFilePath (String outputFilePath) {
		if(tmpOutputFilePaths.get(outputFilePath) == null) {
			try {
				tmpOutputFilePaths.put(outputFilePath, new BufferedWriter(new FileWriter(outputFilePath)));
			}catch(IOException ioe) {

			}
		}
	}

	public static BufferedWriter getOutputBW (String outputFilePath) {
		BufferedWriter BW = tmpOutputFilePaths.get(outputFilePath);
		if(BW == null) {
			enrollTmpOutputFilePath(outputFilePath);
			BW = tmpOutputFilePaths.get(outputFilePath);
		}

		return BW;
	}

	public static void closeOutputBW (String outputFilePath) {
		BufferedWriter BW = tmpOutputFilePaths.get(outputFilePath);
		if(BW != null) {
			try {
				BW.close();
			}catch(IOException ioe) {

			}
		}
	}
}
