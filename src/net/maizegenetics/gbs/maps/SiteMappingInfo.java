
package net.maizegenetics.gbs.maps;

/**
 * Container class for passing information on where a tags maps.
 * @deprecated Should use {@link TagMappingInfo}
 * @author edbuckler
 */
@Deprecated
public class SiteMappingInfo {
    public int chromosome=Integer.MIN_VALUE;  // 4 bytes
    public byte strand=Byte.MIN_VALUE; // 1 = same sense as reference FASTA file.  -1 = opposite sense.  unknown = Byte.MIN_VALUE  // 1 byte
    public int position=Integer.MIN_VALUE;  // chromosomal position of the barcoded end of the tag  // 4 bytes
    public float r2=Float.NaN;
    public float mapP=Float.NaN;  //Round(Log2(P)), unknown Byte.MIN_VALUE
    public int site=Integer.MIN_VALUE;
    
    public SiteMappingInfo() {   
    }
    
    public SiteMappingInfo(int chromosome, byte strand, int position, 
                float r2, float mapP, int site) {
        this.chromosome=chromosome;
        this.strand=strand;
        this.position=position;
        this.r2=r2;
        this.mapP=mapP;
        this.site=site;

    } 
    
        

}
