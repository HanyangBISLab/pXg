package progistar.pXg.data.parser;

import java.io.File;
import java.util.Comparator;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.io.FilenameUtils;

import progistar.pXg.constants.Constants;
import progistar.pXg.constants.Parameters;
import progistar.pXg.utils.Logger;

public class ParameterParser {

	/**
	 * If parsing successfully, return 0. else -1.
	 *
	 * @param args
	 * @return
	 */
	public static int parseParams (String[] args) {
		System.out.println(Constants.VERSION+" "+Constants.RELEASE);
		System.out.println(Constants.INTRODUCE);
		System.out.println();

		
		CommandLine cmd = null;
		
		
		// Mandatory
		Option optionPSM = Option.builder("p")
				.longOpt("psm").argName("madatory, tsv|csv")
				.hasArg()
				.required(true)
				.desc("PSM file path. It is expected that the psm file is derived from proteomics search by de novo or database search engine.")
				.build();
		Option optionBAM = Option.builder("b")
				.longOpt("bam").argName("mandatory, bam|sam")
				.hasArg()
				.required(true)
				.desc("SAM/BAM file path. The sam/bam file must be sorted by coordinate. Multiple SAM/BAM files should be separated by comma (,).")
				.build();
		Option optionGTF = Option.builder("g")
				.longOpt("gtf").argName("mandatory, gtf")
				.hasArg()
				.required(true)
				.desc("GTF file path. We recommand to use the same gtf corresponding to alignment.")
				.build();
		Option optionOutput = Option.builder("o")
				.longOpt("output").argName("mandatory, string")
				.hasArg()
				.required(true)
				.desc("Output file path of pXg.")
				.build();
		
		Option optionPeptideIdx = Option.builder("pi")
				.longOpt("peptide_index").argName("mandatory, integer")
				.hasArg()
				.required(true)
				.desc("Peptide column index in the psm file (one-based).")
				.build();
		Option optionIdIdx = Option.builder("idi")
				.longOpt("identifier_index").argName("madatory, integer")
				.hasArg()
				.required(true)
				.desc("PSM identifier indicies (one-based). One or more indicies can be specified by comma separated. ex> 3,5,7")
				.build();
		Option optionChargeIdx = Option.builder("ci")
				.longOpt("charge_index").argName("madatory, integer")
				.hasArg()
				.required(true)
				.desc("Charge state index (one-based).")
				.build();
		Option optionScoreIdx = Option.builder("si")
				.longOpt("score_index").argName("madatory, integer")
				.hasArg()
				.required(true)
				.desc("Main search score index (one-based).")
				.build();
		
		// optional
		Option optionCountStrategy = Option.builder("c")
				.longOpt("count").argName("optional, primary|all")
				.hasArg()
				.required(false)
				.desc("Specify which reads to consider for counting. The default is primary.")
				.build();
		Option optionNormalization = Option.builder("n")
				.longOpt("normalize").argName("optional, rphm|raw")
				.hasArg()
				.required(false)
				.desc("Specify the moethod of normalization to be used. The default is RPHM (read per a hundred million).")
				.build();
		Option optionAdditionalFeatures = Option.builder("ai")
				.longOpt("add_index").argName("optional, integer")
				.hasArg()
				.required(false)
				.desc("Specify the indices for additional features to generate PIN file (one-based). Several features can be added by comma separator. ex> 8,10,12")
				.build();
		Option optionMode = Option.builder("m")
				.longOpt("mode").argName("optional, auto|fr|rf|r|f|none")
				.hasArg()
				.required(false)
				.desc("Specify strandedness (default is auto). "
						+ "auto: auto-detection. only available in paired-ends. "
						+ "fr: first-forward second-reverse. "
						+ "rf: first-reverse second-forward. "
						+ "r: reverse single end. "
						+ "f: forward single end. "
						+ "none: non-strandedness.")
				.build();
		Option optionSeparator = Option.builder("s")
				.longOpt("sep").argName("optional, tsv|csv")
				.hasArg()
				.required(false)
				.desc("Specify the column separator. Possible values are csv or tsv. Default is tsv.")
				.build();
		Option optionAAVariant = Option.builder("a")
				.longOpt("aa_variant").argName("optional, string")
				.hasArg()
				.required(false)
				.desc("File path of amino acid variant table.")
				.build();
		Option optionPrintTargetOnly = Option.builder("to")
				.longOpt("print_target_only").argName("optional")
				.required(false)
				.desc("Print target only.")
				.build();
		Option optionIL = Option.builder("il")
				.longOpt("il_equivalent").argName("optional, true|false")
				.hasArg()
				.required(false)
				.desc("Controls whether pXg treats isoleucine (I) and leucine (L) as the same/equivalent with respect to a peptide identification. Default is true.")
				.build();
		Option optionLengths = Option.builder("l")
				.longOpt("lengths").argName("optional, integer")
				.hasArg()
				.required(false)
				.desc("Range of peptide length to consider. Default is 8-15. You can write in this way (min-max, both inclusive) : 8-13.")
				.build();
		Option optionMaxFlankSize = Option.builder("fs")
				.longOpt("flank_size").argName("optional, integer")
				.hasArg()
				.required(false)
				.desc("Specify to print maximum flank nuleotides from the matched sequence. Default is 1,000.")
				.build();
		Option optionFasta = Option.builder("f")
				.longOpt("fasta").argName("optional, fasta|fa")
				.hasArg()
				.required(false)
				.desc("Reference sequence database to avoid ambiguous assignment of non-reference peptides.")
				.build();
		Option optionRank = Option.builder("r")
				.longOpt("rank").argName("optional, integer")
				.hasArg()
				.required(false)
				.desc("How many candidates will be considered per scan. Default is 100 (in other words, use all candidates).")
				.build();
		Option optionOutputSAM = Option.builder("os")
				.longOpt("output_sam").argName("optional, true|false")
				.hasArg()
				.required(false)
				.desc("Report matched reads as SAM format (true or false). Default is false.")
				.build();
		Option optionOutputNonreference = Option.builder("on")
				.longOpt("output_non_reference").argName("optional, true|false")
				.hasArg()
				.required(false)
				.desc("Report non-reference peptides for SAM and/or GTF formats (true or false). Default is true.")
				.build();
		Option optionOutputReference = Option.builder("or")
				.longOpt("output_reference").argName("optional, true|false")
				.hasArg()
				.required(false)
				.desc("Report reference peptides for SAM and/or GTF formats (true or false). Default is true.")
				.build();
		Option optionPenaltyMutation = Option.builder("pm")
				.longOpt("penalty_mutation").argName("optional, integer")
				.hasArg()
				.required(false)
				.desc("Penalty per a mutation. Default is "+Parameters.PENALTY_MUTATION+".")
				.build();
		Option optionPenaltyAlternativeSplicing = Option.builder("pas")
				.longOpt("penalty_alternative_splicing").argName("optional, integer")
				.hasArg()
				.required(false)
				.desc("Penalty for alternative splicing. Default is "+Parameters.PENALTY_AS+".")
				.build();
		Option optionPenalty5UTR = Option.builder("p5")
				.longOpt("penalty_5utr").argName("optional, integer")
				.hasArg()
				.required(false)
				.desc("Penalty for 5`-UTR. Default is "+Parameters.PENALTY_5UTR+".")
				.build();
		Option optionPenalty3UTR = Option.builder("p3")
				.longOpt("penalty_3utr").argName("optional, integer")
				.hasArg()
				.required(false)
				.desc("Penalty for 3`-UTR. Default is "+Parameters.PENALTY_3UTR+".")
				.build();
		Option optionPenaltyncRNA = Option.builder("pn")
				.longOpt("penalty_ncrna").argName("optional, integer")
				.hasArg()
				.required(false)
				.desc("Penalty for noncoding RNA. Default is "+Parameters.PENALTY_ncRNA+".")
				.build();
		Option optionPenaltyncFS = Option.builder("pf")
				.longOpt("penalty_frameshift").argName("optional, integer")
				.hasArg()
				.required(false)
				.desc("Penalty for frame shift. Default is "+Parameters.PENALTY_FS+".")
				.build();
		Option optionPenaltyncIR = Option.builder("pir")
				.longOpt("penalty_intron_retention").argName("optional, integer")
				.hasArg()
				.required(false)
				.desc("Penalty for intron region. Default is "+Parameters.PENALTY_IR+".")
				.build();
		Option optionPenaltyncIGR = Option.builder("pigr")
				.longOpt("penalty_intergenic_region").argName("optional, integer")
				.hasArg()
				.required(false)
				.desc("Penalty for intergenic region. Default is "+Parameters.PENALTY_IGR+".")
				.build();
		Option optionPenaltyncasRNA = Option.builder("pa")
				.longOpt("penalty_antisense_region").argName("optional, integer")
				.hasArg()
				.required(false)
				.desc("Penalty for antisense RNA. Default is "+Parameters.PENALTY_asRNA+".")
				.build();
		Option optionPenaltyncSoftClip = Option.builder("ps")
				.longOpt("penalty_softclip").argName("optional, integer")
				.hasArg()
				.required(false)
				.desc("Penalty for softclip reads. Default is "+Parameters.PENALTY_SOFTCLIP+".")
				.build();
		Option optionPenaltyncUnknown = Option.builder("pu")
				.longOpt("penalty_unknown").argName("optional, integer")
				.hasArg()
				.required(false)
				.desc("Penalty for unmapped reads. Default is "+Parameters.PENALTY_UNMAP+".")
				.build();
		Option optionGTFPartitionSize = Option.builder("gps")
				.longOpt("gtf_partition_size").argName("optional, integer")
				.hasArg()
				.required(false)
				.desc("The block size of loading genomic region at once. Default is "+Parameters.partitionSize+".")
				.build();
		Option optionBAMPartitionSize = Option.builder("bps")
				.longOpt("bam_partition_size").argName("optional, integer")
				.hasArg()
				.required(false)
				.desc("The number of loading reads at once. Default is "+Parameters.readSize+".")
				.build();
		Option optionThreads = Option.builder("t")
				.longOpt("threads").argName("optional, integer")
				.hasArg()
				.required(false)
				.desc("The number of threads. Default is "+Parameters.nThreads+".")
				.build();
		
		
		Options options = new Options();
		options
		.addOption(optionPSM)
		.addOption(optionBAM)
		.addOption(optionGTF)
		.addOption(optionOutput)
		.addOption(optionPeptideIdx)
		.addOption(optionIdIdx)
		.addOption(optionChargeIdx)
		.addOption(optionScoreIdx)
		
		.addOption(optionCountStrategy)
		.addOption(optionNormalization)
		.addOption(optionAdditionalFeatures)
		.addOption(optionFasta)
		.addOption(optionMode)
		.addOption(optionSeparator)
		.addOption(optionPrintTargetOnly)
		.addOption(optionIL)
		.addOption(optionLengths)
		.addOption(optionMaxFlankSize)
		.addOption(optionRank)
		.addOption(optionAAVariant)
		.addOption(optionOutputSAM)
		.addOption(optionOutputNonreference)
		.addOption(optionOutputReference)
		.addOption(optionPenaltyMutation)
		.addOption(optionPenaltyAlternativeSplicing)
		.addOption(optionPenalty5UTR)
		.addOption(optionPenalty3UTR)
		.addOption(optionPenaltyncRNA)
		.addOption(optionPenaltyncFS)
		.addOption(optionPenaltyncIR)
		.addOption(optionPenaltyncIGR)
		.addOption(optionPenaltyncasRNA)
		.addOption(optionPenaltyncSoftClip)
		.addOption(optionPenaltyncUnknown)
		.addOption(optionGTFPartitionSize)
		.addOption(optionBAMPartitionSize)
		.addOption(optionThreads);
		
		
		CommandLineParser parser = new DefaultParser();
	    boolean isFail = false;

		try {
		    cmd = parser.parse(options, args);
		    
		    // --bam
		    if(cmd.hasOption("b")) {
		    	String[] paths = cmd.getOptionValue("b").split(",");
				Parameters.NUM_OF_SAM_FILES = paths.length;
				Parameters.sequenceFilePaths = new String[Parameters.NUM_OF_SAM_FILES];
				Parameters.unmappedFilePaths = new String[Parameters.NUM_OF_SAM_FILES];
				Parameters.exportSAMPaths	 = new String[Parameters.NUM_OF_SAM_FILES];
				Parameters.tmpOutputFilePaths= new String[Parameters.NUM_OF_SAM_FILES];

				for(int idx=0; idx < Parameters.NUM_OF_SAM_FILES; idx++) {
					Parameters.sequenceFilePaths[idx] = paths[idx];
					if(!(Parameters.sequenceFilePaths[idx].toLowerCase().endsWith(".bam") ||
							Parameters.sequenceFilePaths[idx].toLowerCase().endsWith(".sam"))) {
						System.out.println(Parameters.sequenceFilePaths[idx] +" must be .sam or .bam file.");
						isFail = true;
					}

					if(!isExist(Parameters.sequenceFilePaths[idx])) {
						printNoSuchFileOrDirectory(Parameters.sequenceFilePaths[idx]);
						isFail = true;
					}
					
					int lastIdx = Parameters.sequenceFilePaths[idx].lastIndexOf(".");
					
					Parameters.unmappedFilePaths[idx] = Parameters.sequenceFilePaths[idx].substring(0, lastIdx) + ".unknown.seq";
					Parameters.exportSAMPaths[idx]	  = Parameters.sequenceFilePaths[idx].substring(0, lastIdx) + ".ided.sam";
					Parameters.tmpOutputFilePaths[idx]= Parameters.sequenceFilePaths[idx].substring(0, lastIdx) + "."+Constants.UNIQUE_RUN_ID;
					
				}
		    }
		    // --gtf
		    if(cmd.hasOption("g")) {
		    	Parameters.genomicAnnotationFilePath = cmd.getOptionValue("g");
				if(!isExist(Parameters.genomicAnnotationFilePath)) {
					printNoSuchFileOrDirectory(Parameters.genomicAnnotationFilePath);
					return -1;
				}
		    }
		    
		    // --psm
		    if(cmd.hasOption("p")) {
		    	Parameters.peptideFilePath = cmd.getOptionValue("p");
				if(!isExist(Parameters.peptideFilePath)) {
					printNoSuchFileOrDirectory(Parameters.peptideFilePath);
					return -1;
				}
		    }
		    
		    // --aa_variant
		    if(cmd.hasOption("a")) {
		    	Parameters.aaVariantTableFilePath = cmd.getOptionValue("a");
		    	if(!isExist(Parameters.aaVariantTableFilePath)) {
					printNoSuchFileOrDirectory(Parameters.aaVariantTableFilePath);
					return -1;
				}
		    }
		    
		    // --peptide_index
		    if(cmd.hasOption("pi")) {
		    	Parameters.peptideColumnIndex = Integer.parseInt(cmd.getOptionValue("pi"));
		    }
		    
		    // --identifier_index
		    if(cmd.hasOption("idi")) {
		    	String[] indicies = cmd.getOptionValue("idi").split(",");
		    	
		    	Parameters.identifierColumnIndicies = new int[indicies.length];
		    	for(int i=0; i<indicies.length; i++) {
		    		Parameters.identifierColumnIndicies[i] = Integer.parseInt(indicies[i]);
		    	}
		    }
		    
		    // --charge_index
		    if(cmd.hasOption("ci")) {
		    	Parameters.chargeColumnIndex = Integer.parseInt(cmd.getOptionValue("ci"));
		    }
		    
		    // --score_index
		    if(cmd.hasOption("si")) {
		    	Parameters.scoreColumnIndex = Integer.parseInt(cmd.getOptionValue("si"));
		    }
		    
		    
		    // --output
		    if(cmd.hasOption("o")) {
		    	Parameters.outputFilePath = cmd.getOptionValue("o") +".pXg";
				Parameters.pinFilePath = Parameters.outputFilePath +".pin";
				if(isExist(Parameters.outputFilePath)) {
					printAlreadyExistFileOrDirectory(Parameters.outputFilePath);
					isFail = true;
				}
				
				File baseFile = new File(Parameters.outputFilePath);
				String basePath = FilenameUtils.getFullPathNoEndSeparator(baseFile.getAbsolutePath());
				if(basePath == null) {
					basePath = "";
				}
				
				
				File outputFileWOExtension = new File(cmd.getOptionValue("o"));
				// rename and relocation of sam/bam-related outputs
				for(int idx=0; idx < Parameters.NUM_OF_SAM_FILES; idx++) {
					
					File file = new File(Parameters.unmappedFilePaths[idx]);
					Parameters.unmappedFilePaths[idx] = basePath +"/"+outputFileWOExtension.getName()+"."+file.getName();
					
					file = new File(Parameters.exportSAMPaths[idx]);
					Parameters.exportSAMPaths[idx] = basePath +"/"+outputFileWOExtension.getName()+"."+file.getName();
					
					file = new File(Parameters.tmpOutputFilePaths[idx]);
					Parameters.tmpOutputFilePaths[idx] = basePath +"/"+outputFileWOExtension.getName()+"."+file.getName();
				}
		    }
		    
		    // --count
		    if(cmd.hasOption("c")) {
		    	String count = cmd.getOptionValue("c");
		    	if(count.equalsIgnoreCase("primary")) {
		    		Parameters.COUNT_PRIMARY_ONLY = true;
		    	} else if(count.equalsIgnoreCase("all")) {
		    		Parameters.COUNT_PRIMARY_ONLY = false;
		    	} else {
		    		System.out.println(count +" is not supported value.");
					isFail = true;
		    	}
		    }
		    
		 // --count
		    if(cmd.hasOption("n")) {
		    	String normalization = cmd.getOptionValue("n");
		    	if(normalization.equalsIgnoreCase("rphm")) {
		    		Parameters.COUNT_NORMALIZATION = true;
		    	} else if(normalization.equalsIgnoreCase("raw")) {
		    		Parameters.COUNT_NORMALIZATION = false;
		    	} else {
		    		System.out.println(normalization +" is not supported value.");
					isFail = true;
		    	}
		    }
		    
		    // --add_index
		    if(cmd.hasOption("ai")) {
		    	String[] indicies = cmd.getOptionValue("ai").split(",");
				Parameters.additionalFeatureIndices = new int[indicies.length];
				for(int idx=0; idx<indicies.length; idx++) {
					Parameters.additionalFeatureIndices[idx] = Integer.parseInt(indicies[idx]);
				}
		    }
		    
		    // --sep
		    if(cmd.hasOption("s")) {
		    	Parameters.sepType = cmd.getOptionValue("s");
		    }
		    
		 // --sep
		    if(cmd.hasOption("to")) {
		    	Parameters.printTargetOnly = true;
		    }
		    
		    // --mode
		    if(cmd.hasOption("m")) {
		    	String mode = cmd.getOptionValue("m");
		    	if(!mode.equalsIgnoreCase(Constants.NON_STRANDED) && 
				   !mode.equalsIgnoreCase(Constants.F_STRANDED) &&
				   !mode.equalsIgnoreCase(Constants.R_STRANDED) &&
				   !mode.equalsIgnoreCase(Constants.RF_STRANDED) &&
				   !mode.equalsIgnoreCase(Constants.FR_STRANDED) &&
				   !mode.equalsIgnoreCase(Constants.AUTO_STRANDED)) {
					System.out.println(mode +" is not supported value.");
					isFail = true;
				} else {
					Parameters.strandedness = mode;
				}
		    }
		    
		    // --il_equivalence
		    if(cmd.hasOption("il")) {
		    	if(cmd.getOptionValue("il").equalsIgnoreCase("false")) {
					Parameters.leucineIsIsoleucine = false;
				}
		    }
		    
		    // --lengths
		    if(cmd.hasOption("l")) {
		    	String[] range = cmd.getOptionValue("l").split("\\-");
				Parameters.minPeptLen = Integer.parseInt(range[0]);
				Parameters.maxPeptLen = Integer.parseInt(range[1]);
		    }
		    
		    // --flank_size
		    if(cmd.hasOption("fs")) {
		    	Integer mFlankSize = Integer.parseInt(cmd.getOptionValue("fs"));
				Parameters.maxFlankNSize= mFlankSize;
		    }
		    
		    // --rank
		    if(cmd.hasOption("r")) {
		    	Parameters.psmRank = Integer.parseInt(cmd.getOptionValue("r"));
		    }
		    
		    // --output_sam
		    if(cmd.hasOption("os")) {
		    	if(cmd.getOptionValue("os").equalsIgnoreCase("true")) {
					Parameters.EXPORT_SAM = true;
				}
		    }
		    
		    // --output_non_reference
		    if(cmd.hasOption("on")) {
		    	if(cmd.getOptionValue("on").equalsIgnoreCase("false")) {
					Parameters.EXPORT_NON_REFERENCE = false;
				}
		    }
		    
		    // --output_reference
		    if(cmd.hasOption("or")) {
		    	if(cmd.getOptionValue("oc").equalsIgnoreCase("false")) {
					Parameters.EXPORT_REFERENCE= false;
				}
		    }
		    
		    // --penalty_mutation
		    if(cmd.hasOption("pm")) {
	    		Double penalty = Double.parseDouble(cmd.getOptionValue("pm"));
				Parameters.PENALTY_MUTATION = penalty;
		    }
		    
		    // --penalty_alternative_splicing
		    if(cmd.hasOption("pas")) {
	    		Double penalty = Double.parseDouble(cmd.getOptionValue("pas"));
				Parameters.PENALTY_AS = penalty;
		    }
		    
		    // --penalty_5utr
		    if(cmd.hasOption("p5")) {
	    		Double penalty = Double.parseDouble(cmd.getOptionValue("p5"));
				Parameters.PENALTY_5UTR = penalty;
		    }
		    
		    // --penalty_3utr
		    if(cmd.hasOption("p3")) {
	    		Double penalty = Double.parseDouble(cmd.getOptionValue("p3"));
				Parameters.PENALTY_3UTR = penalty;
		    }
		    
		    // --penalty_ncrna
		    if(cmd.hasOption("pn")) {
	    		Double penalty = Double.parseDouble(cmd.getOptionValue("pn"));
				Parameters.PENALTY_ncRNA = penalty;
		    }
		    
		    // --penalty_frameshift
		    if(cmd.hasOption("pf")) {
	    		Double penalty = Double.parseDouble(cmd.getOptionValue("pf"));
				Parameters.PENALTY_FS = penalty;
		    }
		    
		    // --penalty_intron_retention
		    if(cmd.hasOption("pir")) {
	    		Double penalty = Double.parseDouble(cmd.getOptionValue("pir"));
				Parameters.PENALTY_IR = penalty;
		    }
		    
		    // --penalty_intergenic_region
		    if(cmd.hasOption("pigr")) {
	    		Double penalty = Double.parseDouble(cmd.getOptionValue("pigr"));
				Parameters.PENALTY_IGR = penalty;
		    }
		    
		    // --penalty_asRNA
		    if(cmd.hasOption("pa")) {
	    		Double penalty = Double.parseDouble(cmd.getOptionValue("pa"));
				Parameters.PENALTY_asRNA = penalty;
		    }
		    
		    // --penalty_softclip
		    if(cmd.hasOption("ps")) {
	    		Double penalty = Double.parseDouble(cmd.getOptionValue("ps"));
				Parameters.PENALTY_SOFTCLIP = penalty;
		    }
		    
		    // --penalty_unknown
		    if(cmd.hasOption("pu")) {
	    		Double penalty = Double.parseDouble(cmd.getOptionValue("pu"));
				Parameters.PENALTY_UNMAP = penalty;
		    }
		    
		    // --gtf_partition_size
		    if(cmd.hasOption("gps")) {
		    	Parameters.partitionSize = Integer.parseInt(cmd.getOptionValue("gps"));
		    }
		    
		    // --bam_partition_size
		    if(cmd.hasOption("bps")) {
		    	Parameters.readSize = Integer.parseInt(cmd.getOptionValue("bps"));
		    }
		    
		    // --fasta
		    if(cmd.hasOption("f")) {
		    	Parameters.proteinFastaPath = cmd.getOptionValue("f");
				if(!isExist(Parameters.proteinFastaPath)) {
					printNoSuchFileOrDirectory(Parameters.proteinFastaPath);
					isFail = true;
				}
		    }
		    
		    
		    // --threads
		    if(cmd.hasOption("t")) {
		    	Parameters.nThreads = Integer.parseInt(cmd.getOptionValue("t"));
		    }
		    
		} catch (ParseException e) {
			System.out.println(e.getMessage());
			isFail = true;
		}
		
		if(isFail) {
			System.out.println("Usage:");
		    for(Option option : options.getOptions()) {
		    	System.out.println("-"+option.getOpt()+", --"+option.getLongOpt()+" <"+option.getArgName()+">");
		    	System.out.println("\t"+option.getDescription());
		    	System.out.println();
		    }
		    
		    System.out.println("Example1");
			System.out.println("java -Xmx30G -jar pXg.jar --gtf gencode.gtf --bam aligned_1.sorted.bam,aligned_2.sorted.bam --psm peaks.result --identifier_index 2,5 --peptide_index 4 --charge_index 11 --score_index 8 --add_index 14,15  --output test");
			System.out.println("Example2");
			System.out.println("java -Xmx30G -jar pXg.jar --gtf gencode.gtf --bam aligned.sorted.sam -psm peaks.result --identifier_index 2,5 --peptide_index 4 --charge_index 11 --score_index 8 --add_index 14,15 --lengths 8-11 --output test");
			
		    System.exit(0);
		}
	

		// open logger
    	Logger.create(Parameters.outputFilePath+".log");

		printSetting();

		// change one-based to zero-based
		Parameters.peptideColumnIndex--;
		Parameters.scoreColumnIndex--;
		Parameters.chargeColumnIndex--;
		
		for(int i=0; i<Parameters.identifierColumnIndicies.length; i++) {
			Parameters.identifierColumnIndicies[i]--;
		}
		if(Parameters.additionalFeatureIndices != null) {
			for(int i=0; i<Parameters.additionalFeatureIndices.length; i++) {
				Parameters.additionalFeatureIndices[i]--;
			}
		}

		return 0;
	}
	
	/**
	 * GTF
	 * SAM
	 * PSM
	 * FASTA
	 * OUT
	 * P-value
	 * FDR
	 * SCAN_COLS
	 * PEPT_COL
	 * SCORE_COL
	 * RANKS TO CONSIDER
	 * LENGTH
	 * GTF_PARTITION_SIZE
	 * SAM_PARTITION_SIZE
	 * THREADS
	 */
	private static void printSetting () {
		System.out.println("Running info");
		System.out.println(" GTF: "+Parameters.genomicAnnotationFilePath);
		System.out.println("  GTF_PARTITION_SIZE: "+Parameters.partitionSize);
		String samPaths = Parameters.sequenceFilePaths[0];
		for(int i=1; i<Parameters.sequenceFilePaths.length; i++) {
			samPaths += "," + Parameters.sequenceFilePaths[i];
		}
		System.out.println(" SAM: "+samPaths);
		System.out.println("  SAM_PARTITION_SIZE: "+Parameters.readSize);

		System.out.println("  STRANDEDNESS_MODE: "+Parameters.strandedness);
		if(Parameters.proteinFastaPath != null) {
			System.out.println(" PROTEIN_DB: "+Parameters.proteinFastaPath);
		}
		System.out.println(" PSM: "+Parameters.peptideFilePath);
		String identifierIndex = "";
		for(int identifier : Parameters.identifierColumnIndicies) {
			identifierIndex += "," + identifier;
		}
		identifierIndex = identifierIndex.substring(1);
		System.out.println("  IDENTIFIER_COLS: "+identifierIndex);
		System.out.println("  PEPT_COL: "+Parameters.peptideColumnIndex);
		System.out.println("  SCORE_COL: "+Parameters.scoreColumnIndex);
		
		System.out.println("  CHARGE_COL: "+Parameters.chargeColumnIndex);
		// to display array
		String addFeatCols = "NA";
		if(Parameters.additionalFeatureIndices != null) {
			addFeatCols = "";
			for (int additionalFeatureIndex : Parameters.additionalFeatureIndices) {
				addFeatCols += "," + additionalFeatureIndex;
			}
			addFeatCols = addFeatCols.substring(1);
		}
		
		String countReads = "all";
		if(Parameters.COUNT_PRIMARY_ONLY) {
			countReads = "primary";
		}
		
		String normalization = "raw";
		if(Parameters.COUNT_NORMALIZATION) {
			normalization = "RPHM";
		}
		
		String printTargetOnly = "target and decoy";
		if(Parameters.printTargetOnly) {
			printTargetOnly = "target only";
		}
		
		String aaVariant = "NA";
		if(Parameters.aaVariantTableFilePath != null) {
			aaVariant = Parameters.aaVariantTableFilePath;
		}
		
		System.out.println("  COUNT_READS: "+countReads);
		System.out.println("  NORMALIZATED_READS: "+normalization);
		System.out.println("  ADDITIONAL_FEATURE_COLS: "+addFeatCols);
		System.out.println("  RANK TO CONSIDER: "+Parameters.psmRank);
		System.out.println("  PEPTIDE_LENGTHS: "+Parameters.minPeptLen+"-"+Parameters.maxPeptLen);
		System.out.println("  MAX_FLANK_SIZE: "+Parameters.maxFlankNSize);
		System.out.println(" OUT_RESULT: "+Parameters.outputFilePath);
		System.out.println("  OUT_PIN.: "+Parameters.pinFilePath);
		System.out.println("  OUT_UNKNOWN: "+Parameters.EXPORT_UNMAPPED_SEQ);
		System.out.println("  OUT_SAM: "+Parameters.EXPORT_SAM);
		System.out.println("  OUT_REFERENCE: "+Parameters.EXPORT_REFERENCE);
		System.out.println("  OUT_NON_REFERENCE: "+Parameters.EXPORT_NON_REFERENCE);
		System.out.println("  PRINT_TYPE: "+printTargetOnly);
		System.out.println(" AA_VARIANT: "+aaVariant);
		System.out.println(" penalty_mutation: "+Parameters.PENALTY_MUTATION);
		System.out.println(" penalty_AS: "+Parameters.PENALTY_AS);
		System.out.println(" penalty_5UTR: "+Parameters.PENALTY_5UTR);
		System.out.println(" penalty_3UTR: "+Parameters.PENALTY_3UTR);
		System.out.println(" penalty_ncRNA: "+Parameters.PENALTY_ncRNA);
		System.out.println(" penalty_FS: "+Parameters.PENALTY_FS);
		System.out.println(" penalty_IR: "+Parameters.PENALTY_IR);
		System.out.println(" penalty_IGR: "+Parameters.PENALTY_IGR);
		System.out.println(" penalty_asRNA: "+Parameters.PENALTY_asRNA);
		System.out.println(" penalty_softclip: "+Parameters.PENALTY_SOFTCLIP);
		System.out.println(" penalty_unknown: "+Parameters.PENALTY_UNMAP);
		System.out.println(" THREADS: "+Parameters.nThreads);

		// append to logger
		Logger.append("Running info");
		Logger.newLine();
		Logger.append(" GTF: "+Parameters.genomicAnnotationFilePath);
		Logger.newLine();
		Logger.append("  GTF_PARTITION_SIZE: "+Parameters.partitionSize);
		Logger.newLine();
		Logger.append(" SAM: "+samPaths);
		Logger.newLine();
		Logger.append("  SAM_PARTITION_SIZE: "+Parameters.readSize);
		Logger.newLine();
		Logger.append("  STRANDEDNESS_MODE: "+Parameters.strandedness);
		Logger.newLine();
		if(Parameters.proteinFastaPath != null) {
			Logger.append(" PROTEIN_DB: "+Parameters.proteinFastaPath);
			Logger.newLine();
		}
		Logger.append(" PSM: "+Parameters.peptideFilePath);
		Logger.newLine();
		Logger.append("  IDENTIFIER_COLS: "+identifierIndex);
		Logger.newLine();
		Logger.append("  PEPT_COL: "+Parameters.peptideColumnIndex);
		Logger.newLine();
		Logger.append("  SCORE_COL: "+Parameters.scoreColumnIndex);
		Logger.newLine();
		
		Logger.append("  COUNT_READS: "+countReads);
		Logger.newLine();
		Logger.append("  NORMALIZATED_READS: "+normalization);
		Logger.newLine();
		Logger.append("  ADDITIONAL_FEATURE_COLS: "+addFeatCols);
		Logger.newLine();
		Logger.append("  RANK TO CONSIDER: "+Parameters.psmRank);
		Logger.newLine();
		Logger.append("  PEPTIDE_LENGTHS: "+Parameters.minPeptLen+"-"+Parameters.maxPeptLen);
		Logger.newLine();
		Logger.append("  MAX_FLANK_SIZE: "+Parameters.maxFlankNSize);
		Logger.newLine();
		Logger.append(" OUT_RESULT: "+Parameters.outputFilePath);
		Logger.newLine();
		Logger.append("  OUT_PIN: "+Parameters.pinFilePath);
		Logger.newLine();
		Logger.append("  OUT_UNMAPPED: "+Parameters.EXPORT_UNMAPPED_SEQ);
		Logger.newLine();
		Logger.append("  OUT_SAM: "+Parameters.EXPORT_SAM);
		Logger.newLine();
		Logger.append("  OUT_REFERENCE: "+Parameters.EXPORT_REFERENCE);
		Logger.newLine();
		Logger.append("  OUT_NON_REFERENCE: "+Parameters.EXPORT_NON_REFERENCE);
		Logger.newLine();
		Logger.append("  PRINT_TYPE: "+printTargetOnly);
		Logger.newLine();
		Logger.append(" AA_VARIANT: "+aaVariant);
		Logger.newLine();
		Logger.append(" penalty_mutation: "+Parameters.PENALTY_MUTATION);
		Logger.newLine();
		Logger.append(" penalty_AS: "+Parameters.PENALTY_AS);
		Logger.newLine();
		Logger.append(" penalty_5UTR: "+Parameters.PENALTY_5UTR);
		Logger.newLine();
		Logger.append(" penalty_3UTR: "+Parameters.PENALTY_3UTR);
		Logger.newLine();
		Logger.append(" penalty_ncRNA: "+Parameters.PENALTY_ncRNA);
		Logger.newLine();
		Logger.append(" penalty_FS: "+Parameters.PENALTY_FS);
		Logger.newLine();
		Logger.append(" penalty_IR: "+Parameters.PENALTY_IR);
		Logger.newLine();
		Logger.append(" penalty_IGR: "+Parameters.PENALTY_IGR);
		Logger.newLine();
		Logger.append(" penalty_asRNA: "+Parameters.PENALTY_asRNA);
		Logger.newLine();
		Logger.append(" penalty_unknown: "+Parameters.PENALTY_UNMAP);
		Logger.newLine();
		Logger.append(" THREADS: "+Parameters.nThreads);
		Logger.newLine();
	}

	private static void printNoSuchFileOrDirectory (String fileName) {
		System.out.println("No such file or directory: "+fileName);
	}

	private static void printAlreadyExistFileOrDirectory (String fileName) {
		System.out.println("The file already exists: "+fileName);
	}

	private static boolean isExist (String fileName) {
		File file = new File(fileName);
		return file.exists();
	}
}
