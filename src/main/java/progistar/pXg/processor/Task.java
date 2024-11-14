package progistar.pXg.processor;

import java.util.ArrayList;

import htsjdk.samtools.SAMRecord;
import progistar.pXg.data.GenomicAnnotation;

public class Task {
	protected int[][] genomicAnnotationIndex;
	protected GenomicAnnotation genomicAnnotation;
	protected ArrayList<SAMRecord> samReads = new ArrayList<>();
	protected int gIndexStart;
	protected boolean isAssigned = false;
	protected int taskID;
	protected int taskType;

	/**
	 * default values: <br>
	 * isAssigned = false <br>
	 * genomicSequences are allocated <br>
	 *
	 */
	public Task() {}

	public Task(int[][] genomicAnnotationIndex, GenomicAnnotation genomicAnnotation,
			ArrayList<SAMRecord> samReads, int start, int taskID, int taskType) {
		super();
		this.genomicAnnotationIndex = genomicAnnotationIndex;
		this.genomicAnnotation = genomicAnnotation;
		this.samReads = samReads;
		this.gIndexStart = start;
		this.taskID = taskID;
		this.taskType = taskType;
	}

	public int getTaskType () {
		return this.taskType;
	}

	public int getTaskID () {
		return taskID;
	}

	public String description () {
		return "task " + this.taskID +": " + this.samReads.size() +" sequences";
	}
}
