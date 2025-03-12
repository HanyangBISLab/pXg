package progistar.pXg;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import progistar.pXg.constants.Constants;
import progistar.pXg.constants.Parameters;
import progistar.pXg.data.pXgRecord;
import progistar.pXg.data.parser.pXgParser;

public class BuildSequenceDB {

	public static void main(String[] args) throws IOException {
		long startTime = System.currentTimeMillis();
		System.out.println(Constants.VERSION+" "+Constants.RELEASE);
		System.out.println(Constants.INTRODUCE);
		System.out.println();

		// parse the options
		parseOptions(args);
		
		File output = new File(Parameters.sequencedbOutputPath);
		System.out.println("write "+output.getName());
		BufferedWriter BW = new BufferedWriter(new FileWriter(output));

		// pXg records
		for(int i=0; i<5; i++) {
			String filePath = Parameters.sequencedbpXgPath[i];
			if(filePath != null) {
				ArrayList<pXgRecord> records = pXgParser.parse(new File(filePath), true);

				int writeCnt = 0;
				for (pXgRecord record : records) {
					
					boolean isWrite = false;
					if(Parameters.isIncludedCanonical && record.isCanonical()) {
						isWrite = true;
					}
					else if(Parameters.isIncludedNoncanonical && !record.isCanonical()) {
						isWrite = true;
					}

					if(isWrite) {
						writeCnt++;
						BW.append(record.getHeader());
						BW.newLine();
						BW.append(record.getTranslatedSequence());
						BW.newLine();
					}
				}
				

				// write the records
				System.out.println("A total of "+writeCnt+" were written in the sequence database");
				System.out.println();

			}
		}
		
		// fasta records
		for(int i=0; i<5; i++) {
			int PE = i+1;
			String filePath = Parameters.sequencedbFastaPath[i];
			if(filePath != null) {
				File referenceFile = new File(filePath);
				System.out.println("Read and write fasta file: "+referenceFile.getName());
				System.out.println("All entries in the file enforce to have PE="+PE+" at header section");

				Pattern PEPattern = Pattern.compile("PE=[0-9]+");
				BufferedReader BR = new BufferedReader(new FileReader(referenceFile));
				String line = null;
				int refCnt = 0;
				while((line = BR.readLine()) != null) {
					if(line.startsWith(">")) {
						refCnt++;
						Matcher matcher = PEPattern.matcher(line);
						if(matcher.find()) {
							line = line.replace(matcher.group(), "PE="+PE);
						} else {
							line += " PE="+PE;
						}
					}
					BW.append(line);
					BW.newLine();
				}

				BR.close();
				System.out.println("A total of "+refCnt +" were written in the sequence database");
				System.out.println();
			}
		}

		long endTime = System.currentTimeMillis();
		System.out.println("Total elapsed time: "+(endTime - startTime)/1000 +" sec");
		BW.close();
	}

	/**
	 * Parse and apply arguments
	 *
	 * @param args
	 */
	public static void parseOptions (String[] args) {
		CommandLine cmd = null;
		Options options = new Options();

		// Mandatory
		Option optionP1 = Option.builder("p1")
				.longOpt("pxg_pe1").argName("pXg input")
				.hasArg()
				.required(false)
				.desc("pXg file path")
				.build();
		
		Option optionP2 = Option.builder("p2")
				.longOpt("pxg_pe2").argName("pXg input")
				.hasArg()
				.required(false)
				.desc("pXg file path")
				.build();
		
		Option optionP3 = Option.builder("p3")
				.longOpt("pxg_pe3").argName("pXg input")
				.hasArg()
				.required(false)
				.desc("pXg file path")
				.build();
		
		Option optionP4 = Option.builder("p4")
				.longOpt("pxg_pe4").argName("pXg input")
				.hasArg()
				.required(false)
				.desc("pXg file path")
				.build();
		
		Option optionP5 = Option.builder("p5")
				.longOpt("pxg_pe5").argName("pXg input")
				.hasArg()
				.required(false)
				.desc("pXg file path")
				.build();
		
		Option optionF1 = Option.builder("f1")
				.longOpt("fasta_pe1").argName("fasta file")
				.hasArg()
				.required(false)
				.desc("fasta file path")
				.build();
		
		Option optionF2 = Option.builder("f2")
				.longOpt("fasta_pe2").argName("fasta file")
				.hasArg()
				.required(false)
				.desc("fasta file path")
				.build();
		
		Option optionF3 = Option.builder("f3")
				.longOpt("fasta_pe3").argName("fasta file")
				.hasArg()
				.required(false)
				.desc("fasta file path")
				.build();
		
		Option optionF4 = Option.builder("f4")
				.longOpt("fasta_pe4").argName("fasta file")
				.hasArg()
				.required(false)
				.desc("fasta file path")
				.build();
		
		Option optionF5 = Option.builder("f5")
				.longOpt("fasta_pe5").argName("fasta file")
				.hasArg()
				.required(false)
				.desc("fasta file path")
				.build();
		
		Option optionOutput = Option.builder("o")
				.longOpt("output").argName("Database output")
				.hasArg()
				.required(true)
				.desc("Sequence database output path")
				.build();

		// Optional
		Option optionCanonical = Option.builder("c")
				.longOpt("canonical").argName("Canonical peptides")
				.required(false)
				.desc("Include canonical peptides")
				.build();
		Option optionNoncanonical = Option.builder("n")
				.longOpt("noncanonical").argName("Non-canonical peptides")
				.required(false)
				.desc("Include non-canonical peptides")
				.build();
		Option optionFlank = Option.builder("f")
				.longOpt("flank").argName("Flank sequences")
				.required(false)
				.desc("Include flank sequences")
				.build();
		Option optionStringent = Option.builder("s")
				.longOpt("stringent").argName("Unambiguous non-canonical peptides")
				.required(false)
				.desc("Exclude non-canonical peptides with FastaIDs")
				.build();


		options
		.addOption(optionP1)
		.addOption(optionP2)
		.addOption(optionP3)
		.addOption(optionP4)
		.addOption(optionP5)
		.addOption(optionF1)
		.addOption(optionF2)
		.addOption(optionF3)
		.addOption(optionF4)
		.addOption(optionF5)
		.addOption(optionOutput)
		.addOption(optionCanonical)
		.addOption(optionNoncanonical)
		.addOption(optionFlank)
		.addOption(optionStringent);

		CommandLineParser parser = new DefaultParser();
	    HelpFormatter helper = new HelpFormatter();
	    boolean isFail = false;

		try {
		    cmd = parser.parse(options, args);

		    if(!cmd.hasOption("c") && !cmd.hasOption("n")) {
		    	isFail = true;
		    	System.out.println("-c or -n must be included.");
		    }
		    if(cmd.hasOption("c")) {
		    	Parameters.isIncludedCanonical = true;
		    	System.out.println("Canonical peptides are included in the sequence database.");
		    }
		    if(cmd.hasOption("n")) {
		    	Parameters.isIncludedNoncanonical = true;
		    	System.out.println("Non-canonical peptides are included in the sequence database.");
		    }
		    if(cmd.hasOption("f")) {
		    	Parameters.isIncludedFlankSequence = true;
		    	System.out.println("Flank sequences are included in the sequence database.");
		    }
		    if(cmd.hasOption("s")) {
		    	Parameters.isStringent = true;
		    	System.out.println("Excluding non-canonical peptides with FastaIDs.");
		    }
		    
		    // pXg inputs
		    boolean isTakenInput = false;
		    if(cmd.hasOption("p1")) {
		    	Parameters.sequencedbpXgPath[0] = cmd.getOptionValue("p1");
		    	isTakenInput = true;
		    }
		    if(cmd.hasOption("p2")) {
		    	Parameters.sequencedbpXgPath[1] = cmd.getOptionValue("p2");
		    	isTakenInput = true;
		    }
		    if(cmd.hasOption("p3")) {
		    	Parameters.sequencedbpXgPath[2] = cmd.getOptionValue("p3");
		    	isTakenInput = true;
		    }
		    if(cmd.hasOption("p4")) {
		    	Parameters.sequencedbpXgPath[3] = cmd.getOptionValue("p4");
		    	isTakenInput = true;
		    }
		    if(cmd.hasOption("p5")) {
		    	Parameters.sequencedbpXgPath[4] = cmd.getOptionValue("p5");
		    	isTakenInput = true;
		    }
		    // fasta inputs
		    if(cmd.hasOption("f1")) {
		    	Parameters.sequencedbFastaPath[0] = cmd.getOptionValue("f1");
		    	isTakenInput = true;
		    }
		    if(cmd.hasOption("f2")) {
		    	Parameters.sequencedbFastaPath[1] = cmd.getOptionValue("f2");
		    	isTakenInput = true;
		    }
		    if(cmd.hasOption("f3")) {
		    	Parameters.sequencedbFastaPath[2] = cmd.getOptionValue("f3");
		    	isTakenInput = true;
		    }
		    if(cmd.hasOption("f4")) {
		    	Parameters.sequencedbFastaPath[3] = cmd.getOptionValue("f4");
		    	isTakenInput = true;
		    }
		    if(cmd.hasOption("f5")) {
		    	Parameters.sequencedbFastaPath[4] = cmd.getOptionValue("f5");
		    	isTakenInput = true;
		    }
		    // check input
		    if(!isTakenInput) {
		    	isFail = true;
		    	System.out.println("pXg or fasta files must be provided.");
		    }
		    
		    if(cmd.hasOption("o")) {
		    	Parameters.sequencedbOutputPath = cmd.getOptionValue("o");
		    	System.out.println("Output the generated sequence database: "+Parameters.sequencedbOutputPath);
		    }
		} catch (ParseException e) {
			System.out.println(e.getMessage());
			isFail = true;
		}

		if(isFail) {
		    helper.printHelp("Usage:", options);
		    System.exit(0);
		}

		System.out.println();
	}
}
