package progistar.pXg.data;

public class AAVariant {

	public int pos = -1; // position in the peptide (zero-based)
	public int aaVariantIndex = -1; // this value indicates the index in the AAVriantTable
	
	public String getAnnotation () {
		// position: zero to one base
		return (this.pos+1) + ":" + AAVariantTable.getAnnotation(this.aaVariantIndex);
	}
}
