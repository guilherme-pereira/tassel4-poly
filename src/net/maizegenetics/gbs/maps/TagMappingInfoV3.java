
package net.maizegenetics.gbs.maps;

/**
 * Container class for information on mapping position.  
 * <p>
 * TODO: Needs to be generalized for describing all mapping positions.
 * <p>
 * The mappingSource can be an identifier. Byte.MIN_VALUE means the tag doesn't exist in the mapping result, other values means it exist either mapped or unmapped.
 * @author edbuckler Fei Lu
 */
public class TagMappingInfoV3 {
    /** Chromosome as an integer */
    public int chromosome=Integer.MIN_VALUE;  // 4 bytes
    /** Strand relative to reference genome. 1 = same sense as reference FASTA file.  -1 = opposite sense.  unknown = Byte.MIN_VALUE */
    public byte strand=Byte.MIN_VALUE;   // 1 byte
    /**Chromosomal position of the barcoded end of the tag */
    public int startPosition=Integer.MIN_VALUE;  //   // 4 bytes
    /**Chromosomal position of the common adapter end of the tag (smaller than startPosition if tag matches minus strand), inclusive */
    public int endPosition=Integer.MIN_VALUE;  //   // 4 bytes
    /**Number of diverging bp (edit distance) from reference, unknown = Byte.MIN_VALUE*/
    public byte divergence=Byte.MIN_VALUE;  
    /**Code of mappingSource.0: Bowtie2; 1: BWA; 2: BLAST; 3: BWAMEM; 4: PE one end; 5: PE the other end; 6: Genetic Mapping*/
    public byte mappingSource = Byte.MIN_VALUE;
    /**The rank of this mapping based on the scores from one aligner*/
    public byte mappingRank = Byte.MIN_VALUE;
    /**The mapping score of this mapping*/
    public short mappingScore = Short.MIN_VALUE;
    /**Double cross-over probability Round(Log2(P)), unknown Byte.MIN_VALUE */
    public byte dcoP=Byte.MIN_VALUE;
    /**Genetic mapping probability Round(Log2(P)), unknown Byte.MIN_VALUE */
    public byte mapP=Byte.MIN_VALUE;  //Round(Log2(P)), unknown Byte.MIN_VALUE
    
    /**
     * MappingSource definition
     */
    public static enum Aligner {
        Bowtie2((byte)0, "Bowtie2"),
        BWA((byte)1, "BWA"),
        Blast((byte)2, "Blast"),
        BWAMEM((byte)3,"BWAMEM"),
        PEEnd1((byte)4, "PEEnd1"),
        PEEnd2((byte)5, "PEEnd2");
        private byte value;
        private String name;
        Aligner (byte mappingSource, String alignerName) {
            value = mappingSource;
            name = alignerName;
        }
        public byte getValue () {
            return value;
        }
        public String getName () {
            return name;
        }
        public static byte getValueFromName (String name) {
            if (name.equalsIgnoreCase(Bowtie2.name)) return Bowtie2.getValue();
            if (name.equalsIgnoreCase(BWA.name)) return BWA.getValue();
            if (name.equalsIgnoreCase(Blast.name)) return Blast.getValue();
            if (name.equalsIgnoreCase(BWAMEM.name)) return BWAMEM.getValue();
            if (name.equalsIgnoreCase(PEEnd1.name)) return PEEnd1.getValue();
            if (name.equalsIgnoreCase(PEEnd2.name)) return PEEnd2.getValue();
            return Byte.MIN_VALUE;
        }
        public static Aligner getAlignerFromName (String name) {
            if (name.equalsIgnoreCase(Bowtie2.name)) return Bowtie2;
            if (name.equalsIgnoreCase(BWA.name)) return BWA;
            if (name.equalsIgnoreCase(Blast.name)) return Blast;
            if (name.equalsIgnoreCase(BWAMEM.name)) return BWAMEM;
            if (name.equalsIgnoreCase(PEEnd1.name)) return PEEnd1;
            if (name.equalsIgnoreCase(PEEnd2.name)) return PEEnd2;
            return null;
        }
    }

    public TagMappingInfoV3() {   
    }
    
    public TagMappingInfoV3(int chromosome, byte strand, int startPosition, int endPosition, byte divergence, byte mappingSource, short mappingScore) {
        this.chromosome=chromosome;
        this.strand=strand;
        this.startPosition=startPosition;
        this.endPosition=endPosition;
        this.divergence=divergence;
        this.mappingSource = mappingSource;
        this.mappingScore = mappingScore;
    } 
    
    public TagMappingInfoV3(int chromosome, byte strand, int startPosition, 
                int endPosition, byte divergence, byte mappingSource, short mappingScore, byte[] variantPosOff, byte[] variantDef,
                byte dcoP, byte mapP) {
        this(chromosome, strand, startPosition, endPosition, divergence, mappingSource, mappingScore);
        this.dcoP=dcoP;
        this.mapP=mapP;
    }
    
    public void setMappingRank (byte mappingRank) {
        this.mappingRank = mappingRank;
    }
    
    public void setMappingSource (byte mappingSource) {
        this.mappingSource = mappingSource;
    }
}
