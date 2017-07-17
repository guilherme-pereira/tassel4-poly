/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.pal.site;

import com.google.common.collect.ComparisonChain;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;
import net.maizegenetics.pal.util.GeneralAnnotationUtils;

import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;
//import java.util.Objects;

/**
 * Provide information on a site and its annotations.  This includes information
 * on position, MAF, coverage.  This class is immutable.
 * <p></p>
 * The annotations are all set using the builder.
 * <p></p>
 * This class has been optimized for memory size compared to the other version.  Every annotated position takes about 84
 * bytes, even with a names and declared variants.  Positions with extended annotation take roughly 117 bytes, and as
 * long as the annotation are repeated extensively the take an additional 8 bytes each.  So a site with 100 annotation would
 * take <1000 bytes.
 *
 * @author Ed Buckler
 */
public final class GeneralPosition implements Position {
       /**Locus of the site (required)*/
    private final Chromosome myChromosome;
    /**Physical position of the site (unknown = Float.NaN)*/
    private final int myPosition;
    /**Strand of the site (unknown = Byte.MIN_VALUE)*/
    private final byte myStrand;
    /**Genetic position in centiMorgans (unknown = Float.NaN)*/
    private final float myCM;
    /**Is type Nucleotide or Text*/
    private final boolean isNucleotide;
    /**Whether the variant define the nature of the indel*/
    private final boolean isIndel;
    private final float myMAF;
    private final float mySiteCoverage;
    private final long myAlleleValue;
    /**Name of the site (default = SLocus_Position)*/
    private final byte[] mySNPIDAsBytes;
//    /**Define the nature of the polymorphism {"ACTAT","-"} or {"A","C","G"} or {"100","103","106"}
//     */
   private final Map.Entry<String, String>[] myVariantsAndAnno;
    private final int hashCode;


    //since there are numerous redundant annotations and variant descriptions, this class use a annotation hash, so that
    //only the pointers are stored.  It takes a little longer to instantiate, but save 3-fold in memory.
    private static final ConcurrentMap<Map.Entry<String, String>,Map.Entry<String, String>> ANNO_HASH = new ConcurrentHashMap<>(500_000);

    public static Map.Entry<String, String> getCanonicalAnnotation(String key, String value) {
        if (ANNO_HASH.size() > 100000) {
            ANNO_HASH.clear();
        }
        Map.Entry<String, String> str= new AbstractMap.SimpleImmutableEntry(key,value);
        Map.Entry<String, String> canon = ANNO_HASH.putIfAbsent(str, str);
        return (canon == null) ? str : canon;
    }

    /**
     * A builder for creating immutable CoreAnnotatedPosition instances. AnnotatedPositions are
     * built off a base of a CorePosition, so build it first.
     *<p> Example:
     * <pre>   {@code
     * Position cp= new CorePosition.Builder(new Chromosome("1"),1232).build();
     * CoreAnnotatedPosition ap= new CoreAnnotatedPosition.Builder(cp)
     *    .maf(0.05f)
     *    .ancAllele(NucleotideAlignmentConstants.C_ALLELE)
     *    .build();}</pre>
     * <p>This would create nucleotide position on chromosome 1 at position 1232.  The MAF is 0.05 and the ancestral allele
     * is C.
     */
    public static class Builder {
        // Required parameters
        private final Chromosome myChromosome;
        private final int myPosition;
        // Optional parameters - initialized to default values
        private byte myStrand=1;
        private float myCM=Float.NaN;
        private String mySNPID=null;
        private boolean isNucleotide=true;
        private boolean isIndel=false;
        private Map.Entry<String, String> myKnownVariants=null;

        //in an allele annotation objects
        private float myMAF = Float.NaN;
        private float mySiteCoverage = Float.NaN;
        private byte[] myAlleles=new byte[Allele.COUNT];
        private long myAllelesAsLong;
        //in an general annotation object
        private ArrayList<Map.Entry<String, String>> myAnnotations=new ArrayList<>(0);

        /**Constructor requires a Position before annotation of the position*/
        public Builder(Chromosome chr, int position) {
            myChromosome = chr;
            myPosition = position;
            Arrays.fill(myAlleles,Alignment.UNKNOWN_ALLELE);
        }

        /**Constructor requires a Position before annotation of the position*/
        public Builder(Position aCorePosition) {
           this(aCorePosition.getChromosome(),aCorePosition.getPosition());
           myStrand=aCorePosition.getStrand();
           myCM=aCorePosition.getCM();
           mySNPID=aCorePosition.getSNPID();
           isNucleotide=aCorePosition.isNucleotide();
           isIndel=aCorePosition.isIndel();
           //myKnownVariants=aCorePosition.getKnownVariants(); //todo Fix
        }

        /**Set strand (default=1)*/
        public Builder strand(byte val) {myStrand = val; return this;}
        /**Set strand (default=Float.NaN)*/
        public Builder cM(float val) {myCM = val; return this;}
        /**Set SNP name (default="S"+Chromosome+"_"+position)*/
        public Builder snpName(String val) {mySNPID = val; return this;}
        /**Set whether position is nucleotide (default=true)*/
        public Builder nucleotide(boolean val) {isNucleotide = val; return this; }
        /**Set whether position is indel (default=false)*/
        public Builder indel(boolean val) {isIndel = val; return this;}
        /**Set text definition of variants (default=null)*/
        public Builder knownVariants(String[] val) {
            Map.Entry<String, String> ent=getCanonicalAnnotation("VARIANT",val[0]+"/"+val[1]);
            myKnownVariants=ent;
            return this;
        }
        /**Set text definition of variants (default=null)*/
        public Builder knownVariants(String val) {
            Map.Entry<String, String> ent=getCanonicalAnnotation("VARIANT",val);
            myKnownVariants=ent;
            return this;
        }

        /**Set Minor Allele Frequency annotation (default=Float.NaN)*/
        public Builder maf(float val) {myMAF = val; return this;}
        /**Set site coverage annotation (default=Float.NaN)*/
        public Builder siteCoverage(float val) {mySiteCoverage = val; return this;}
        /**Set allele annotation by Allele type (default=Alignment.UNKNOWN_ALLELE)*/
        public Builder allele(Allele aT, byte val) {myAlleles[aT.index()] = val; return this;}

        /**Add non-standard annotation*/
        public Builder addAnno(String key, String value) {
            Map.Entry<String, String> ent=getCanonicalAnnotation(key,value);
            myAnnotations.add(ent);
            return this;
        }
        /**Add non-standard annotation*/
        public Builder addAnno(String key, Number value) {
            Map.Entry<String, String> ent=getCanonicalAnnotation(key,value.toString());
            myAnnotations.add(ent);
            return this;
        }

        public GeneralPosition build() {
            for (int i = myAlleles.length-1; i >=0 ; i--) {
                myAllelesAsLong=(myAllelesAsLong<<8)|myAlleles[i];
            }
            if (mySNPID != null) {
                String defaultS=(new StringBuilder("S").append(myChromosome.getName()).append("_").append(myPosition)).toString();
                if(defaultS.equals(mySNPID)) mySNPID=null;
            }
            return new GeneralPosition(this);
        }
    }
    private GeneralPosition(Builder builder) {
        myChromosome = builder.myChromosome;
        myPosition = builder.myPosition;
        myStrand = builder.myStrand;
        myCM = builder.myCM;
        if(builder.mySNPID==null) {
            mySNPIDAsBytes=null;
        } else {
            mySNPIDAsBytes=builder.mySNPID.getBytes();
        }
        isNucleotide = builder.isNucleotide;
        isIndel = builder.isIndel;
        //this looks crazy because it java doesn't support generic arrays
        myVariantsAndAnno=(Map.Entry<String, String>[])new Map.Entry<?,?>[1+builder.myAnnotations.size()];
        myVariantsAndAnno[0]=builder.myKnownVariants;
        for (int i = 0; i < builder.myAnnotations.size(); i++) {
            myVariantsAndAnno[i+1]=builder.myAnnotations.get(i);
        }
        hashCode=calcHashCode();

        myMAF = builder.myMAF;
        mySiteCoverage = builder.mySiteCoverage;
        myAlleleValue=builder.myAllelesAsLong;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == this) {return true;}
        if (!(obj instanceof Position)) {return false;}
        Position o=(Position)obj;
        int result= ComparisonChain.start()
                .compare(myPosition, o.getPosition())  //position is most discriminating for speed
                .compare(myChromosome, o.getChromosome())
                .compare(myCM, o.getCM())
                .compare(myStrand,o.getStrand())
                .result();
        if(result!=0) return false;
        return getSNPID().equals(o.getSNPID()); //This is done last as the string comparison is expensive
    }

    @Override
    public int compareTo(Position o) {
        int result=ComparisonChain.start()
                .compare(myChromosome,o.getChromosome())
                .compare(myPosition,o.getPosition())
                .compare(myCM, o.getCM())
                .compare(myStrand,o.getStrand())
                .result();
        if(result!=0) return result;
        return getSNPID().compareTo(o.getSNPID()); //This is done last as the string comparison is expensive
    }

    @Override
    public String toString() {
        StringBuilder sb=new StringBuilder("Position");
        sb.append("\tChr:").append(getChromosome().getName());
        sb.append("\tPos:").append(getPosition());
        sb.append("\tName:").append(getSNPID());
        sb.append("\tMAF:").append(getGlobalMAF());
        sb.append("\tRef:").append(NucleotideAlignmentConstants.getHaplotypeNucleotide(getAllele(Allele.REF)));
        return sb.toString();
    }

    private int calcHashCode() {
        int hash = 7;
        hash = 37 * hash + this.myChromosome.hashCode();
        hash = 37 * hash + this.myPosition;
        hash = 37 * hash + this.myStrand;
        hash = 37 * hash + Float.floatToIntBits(this.myCM);
        if(mySNPIDAsBytes!=null) hash = 37 * hash + this.mySNPIDAsBytes.hashCode();
        return hash;
    }

    @Override
    public float getGlobalMAF() {
        return myMAF;
    }

    @Override
    public float getGlobalSiteCoverage() {
        return mySiteCoverage;
    }

    @Override
    public byte getAllele(Allele alleleType) {
        return (byte)((myAlleleValue>>(alleleType.index()*8))&0xF);
    }

    @Override
    public int hashCode() {
        return hashCode;
    }

    @Override
    public Chromosome getChromosome() {
        return myChromosome;
    }

    @Override
    public int getPosition() {
        return myPosition;
    }

    @Override
    public byte getStrand() {
        return myStrand;
    }

    @Override
    public float getCM() {
        return myCM;
    }

    @Override
    public String getSNPID() {
         if (mySNPIDAsBytes == null) {
            return (new StringBuilder("S").append(getChromosome().getName()).append("_").append(myPosition)).toString();
        } else {
            return new String(mySNPIDAsBytes);
        }
    }

    @Override
    public boolean isNucleotide() {
        return isNucleotide;
    }

    @Override
    public boolean isIndel() {
        return isIndel;
    }

    @Override
    public String[] getKnownVariants() {
        if((myVariantsAndAnno==null)||(myVariantsAndAnno[0]==null)) return new String[0];
        return myVariantsAndAnno[0].getValue().replace("[", "").replace("]", "").split(", ");
    }

    @Override
    public Object[] getAnnotation(String annoName) {
        return GeneralAnnotationUtils.getAnnotation(myVariantsAndAnno, annoName);
        //return myGA.getAnnotation(annoName);
//        switch (annoName) {  //TODO: uncomment once in Java 7
//            case "locus":return myLocus;
//            case "position":return new Integer[]{myPosition};
//            case "myCM":return myCM;
//            case "strand":return myStrand;
//            case "snpID":return mySNPID;
//        }
//        return myGA.getAnnotation(annoName);
    }

    @Override
    public String[] getTextAnnotation(String annoName) {
        return GeneralAnnotationUtils.getTextAnnotation(myVariantsAndAnno,annoName);
  //      return myGA.getTextAnnotation(annoName);
    }

    @Override
    public double[] getQuantAnnotation(String annoName) {
        return GeneralAnnotationUtils.getQuantAnnotation(myVariantsAndAnno,annoName);
    //    return myGA.getQuantAnnotation(annoName);
    }


    @Override
    public String getConsensusAnnotation(String annoName) {
        return GeneralAnnotationUtils.getConsensusAnnotation(myVariantsAndAnno,annoName);
    }

    @Override
    public double getAverageAnnotation(String annoName) {
        return GeneralAnnotationUtils.getAverageAnnotation(myVariantsAndAnno,annoName);
    }
    
}
