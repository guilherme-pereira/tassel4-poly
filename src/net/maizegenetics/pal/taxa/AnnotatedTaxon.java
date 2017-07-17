/*
 *  AnnotatedTaxon
 */
package net.maizegenetics.pal.taxa;

import com.google.common.collect.ImmutableMultimap;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.util.AbstractAnnotation;
import net.maizegenetics.pal.util.GeneralAnnotation;

/**
 * The generally used class for defining a taxon. Contains its name, plus a
 * series of optional annotations. Use the builder to create immutable
 * instances. //TODO change Identifer name to Taxon and move the class to this
 * package
 *
 * @author Ed Buckler
 */
public final class AnnotatedTaxon extends Identifier implements GeneralAnnotation {

    private final GeneralAnnotation myGA;
    private final String myParent1;  // generally female
    private final String myParent2;  // generally male
    private final float myInbreedF; // inbreeding coefficient
    private final byte mySex;  // 0=both, 1=female, 2=male
    private final String myPedigree;
    //Custom annotation are stored in the map
    //private final Multimap<String, Object> myAnnoMap;

    public String getParent1() {
        return myParent1;
    }

    public String getParent2() {
        return myParent2;
    }

    public float getInbreedF() {
        return myInbreedF;
    }

    public byte getSex() {
        return mySex;
    }

    public String getPedigree() {
        return myPedigree;
    }

    @Override
    public Object[] getAnnotation(String annoName) {
        return myGA.getAnnotation(annoName);
//        switch (annoName) {  //TODO: uncomment once in Java 7
//            case "myParent1":return myLocus;
//            case "myParent2":return myPosition;
//            case "myInbreedF":return myCM;
//            case "mySex":return myStrand;
//            case "pedigree":return mySNPID;
//        }
//       }
    }

    @Override
    public String[] getTextAnnotation(String annoName) {
        return myGA.getTextAnnotation(annoName);
    }

    @Override
    public double[] getQuantAnnotation(String annoName) {
        return myGA.getQuantAnnotation(annoName);
    }

    @Override
    public String getConsensusAnnotation(String annoName) {
        return myGA.getConsensusAnnotation(annoName);
    }

    @Override
    public double getAverageAnnotation(String annoName) {
        return myGA.getAverageAnnotation(annoName);
    }

    /**
     * A builder for creating immutable AnnotatedTaxon instances.
     * <p> Example:
     * <pre>   {@code
     * AnnotatedTaxon cp= new AnnotatedTaxon.Builder("Z001E0001:Line:mays:Zea")
     *   .inbreedF(0.99)
     *   .parents("B73","B97")
     *   .pedigree("(B73xB97)S6I1")
     *   .build();}</pre>
     * <p>This would create an annotatedTaxon.
     */
    public static class Builder {

        // Required parameters
        private final String myTaxonFullName;
        // Optional parameters - initialized to default values
        private String parent1 = null;  //generally female
        private String parent2 = null;  //generally male
        private float inbreedF = Float.NaN;
        private byte sex = 0;  //0=both, 1=female, 2=male
        private String pedigree = null;
        private ImmutableMultimap.Builder<String, Object> myAnnoMapBld = null;
        private ImmutableMultimap<String, Object> myAnnoMap = null;

        /**
         * Constructor for Builder, requires a Taxon object
         *
         * @param aTaxon taxon object
         */
        public Builder(Identifier aTaxon) {
            myTaxonFullName = aTaxon.getFullName();
        }

        /**
         * Constructor for Builder, requires a Taxon name
         *
         * @param aTaxonName name of the taxon
         */
        public Builder(String aTaxonName) {
            myTaxonFullName = aTaxonName;
        }

        /**
         * Set sex: 0=both, 1=female, 2=male (default=0 Both)
         */
        public Builder sex(byte val) {
            sex = val;
            return this;
        }

        /**
         * Set inbreeding coefficient (default=Float.NaN)
         */
        public Builder inbreedF(float val) {
            inbreedF = val;
            return this;
        }

        /**
         * Set text definition of parents (default=null)
         */
        public Builder parents(String mom, String dad) {
            parent1 = mom;
            parent2 = dad;
            return this;
        }

        /**
         * Set text definition of pedigree (default=null)
         */
        public Builder pedigree(String val) {
            pedigree = val;
            return this;
        }

        /**
         * Add non-standard annotation
         */
        public Builder addAnno(String key, String value) {
            if (myAnnoMapBld == null) {
                myAnnoMapBld = new ImmutableMultimap.Builder();
            }
            myAnnoMapBld.put(key, value);
            return this;
        }

        /**
         * Add non-standard annotation
         */
        public Builder addAnno(String key, Number value) {
            if (myAnnoMapBld == null) {
                myAnnoMapBld = new ImmutableMultimap.Builder();
            }
            myAnnoMapBld.put(key, value);
            return this;
        }

        public AnnotatedTaxon build() {
            if (myAnnoMapBld != null) {
                myAnnoMap = myAnnoMapBld.build();
            }
            return new AnnotatedTaxon(this);
        }
    }

    private AnnotatedTaxon(Builder builder) {
        super(builder.myTaxonFullName);
        this.myParent1 = builder.parent1;
        this.myParent2 = builder.parent2;
        this.mySex = builder.sex;
        this.myPedigree = builder.pedigree;
        this.myInbreedF = builder.inbreedF;
        this.myGA = new AbstractAnnotation(builder.myAnnoMap);
    }
}
