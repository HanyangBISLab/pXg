package progistar.pXg.data;

import progistar.pXg.constants.Parameters;
import progistar.pXg.data.parser.pXgParser;

public class pXgRecord {
	private String[] fields = null;

	public pXgRecord (String[] fields) {
		this.fields = fields;
	}
	
	public String getHeader (int pe) {
		StringBuilder header = new StringBuilder(">pXg");

		String id = getID();
		String gn = getValueByFieldName("GeneNames").replace("|", ",");
		String ev = getValueByFieldName("Events").replace("|", ",");
		String gCount = getValueByFieldName("GenomicLociCount");
		String exp = getValueByFieldName("Reads");
		String var = getValueByFieldName("Mutations").replace("|", ",");
		String alt = getValueByFieldName("MutationStatus");
		String gId = getValueByFieldName("GeneIDs").replace("|", ",");
		String rna = getNucleotideSequence();
		String genomicLoci = getValueByFieldName("GenomicLoci").replace("|", " ");
		String strand = getValueByFieldName("Strand");
		String td = isTarget() ? "Target" : "Decoy";
		
		header.append("|").append(id)
		.append("|").append(gId)
		.append(" ").append("TD="+td)
		.append(" ").append("GN="+gn)
		.append(" ").append("NUM="+gCount)
		.append(" ").append("FR=0")
		.append(" ").append("EV="+ev)
		.append(" ").append("EXP="+exp)
		.append(" ").append("VAR="+var)
		.append(" ").append("ALT="+alt)
		.append(" ").append("SR="+strand)
		.append(" ").append("gene_site="+genomicLoci)
		.append(" ").append("PE="+pe);

		if(Parameters.isIncludedFlankSequence) {
			header.append(" ").append("RNA=").append(rna.replace("|", ","));
		}
		
		return header.toString();
	}
	


	public String getID () {
		String id = null;

		String genomicLoci = getValueByFieldName("GenomicLoci").replace("|", ",");
		String centerSeuqnece = getValueByFieldName("ObservedNucleotide");
		String strand = getValueByFieldName("Strand");

		id = genomicLoci+"_"+strand+"_"+centerSeuqnece;

		return id;
	}

	public boolean hasFastaID () {
		return getValueByFieldName("FastaIDs").equalsIgnoreCase("-") ? false : true;
	}

	public String getTranslatedSequence () {
		String strand = getValueByFieldName("Strand");
		String nucleotide = getNucleotideSequence();
		if(Parameters.isIncludedFlankSequence) {
			nucleotide = nucleotide.replaceAll("[-\\|]", "");
		} else {
			nucleotide = nucleotide.split("\\|")[1];
		}

		String peptide = null;
		if(strand.equalsIgnoreCase("+")) {
			peptide = GenomicSequence.translation(nucleotide, 0);
		} else {
			peptide = GenomicSequence.reverseComplementTranslation(nucleotide, 0);
		}
		
		// reverse decoy
		if(!isTarget()) {
			peptide = new StringBuilder(peptide).reverse().toString();
		}
		
		return peptide;
	}

	public boolean isTarget () {
		return getValueByFieldName("Label").equalsIgnoreCase("1") ? true : false;
	}
	
	public boolean isCanonical () {
		return getValueByFieldName("isCanonical").equalsIgnoreCase("true") ? true : false;
	}

	/**
	 * return
	 * left-flank|center|right-flank
	 *
	 * @return
	 */
	private String getNucleotideSequence () {
		StringBuilder sequence = new StringBuilder();
		String centerSeuqnece = getValueByFieldName("ObservedNucleotide").toUpperCase();
		String leftFlank = getValueByFieldName("ObservedLeftFlankNucleotide").toUpperCase();
		String rightFlank = getValueByFieldName("ObservedRightFlankNucleotide").toUpperCase();
		sequence.append(leftFlank);
		sequence.append("|");
		sequence.append(centerSeuqnece);
		sequence.append("|");
		sequence.append(rightFlank);

		return sequence.toString();
	}


	public String getValueByFieldName (String fieldName) {
		String[] header = pXgParser.header;
		String value = null;
		for(int i=0; i<header.length; i++) {
			if(header[i].equalsIgnoreCase(fieldName)) {
				if(value != null) {
					System.out.println(fieldName+" is duplciated");
				} else {
					value = fields[i];
				}
			}
		}
		return value;
	}

	public void setValueByFieldName (String fieldName, String value) {
		String[] header = pXgParser.header;
		for(int i=0; i<header.length; i++) {
			if(header[i].equalsIgnoreCase(fieldName)) {
				fields[i] = value;
			}
		}
	}

	@Override
	public String toString() {
		StringBuilder str = new StringBuilder();

		for(int i=0; i<fields.length; i++) {
			if(i != 0) {
				str.append("\t");
			}
			str.append(fields[i]);
		}

		return str.toString();
	}
}
