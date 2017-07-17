package net.maizegenetics.pal.site;

import net.maizegenetics.pal.alignment.Locus;
import net.maizegenetics.pal.util.GeneralAnnotation;

/**
 * Defines the chromosome structure and length. The name and length recorded for
 * each chromosome.
 *
 * @author Terry Casstevens and Ed Buckler
 */
public class Chromosome implements Comparable<Chromosome>, GeneralAnnotation {

    private static final long serialVersionUID = -5197800047652332969L;
    public static Locus UNKNOWN = new Locus("Unknown");
    private final String myName;
    private final int myChromosomeNumber;
    private final int myLength;
    private final GeneralAnnotation myGA;
    private final int hashCode;

    /**
     *
     * @param name Name of the chromosome
     * @param length Length of chromosome in base pairs
     * @param features Map of features about the chromosome
     */
    public Chromosome(String name, int length, GeneralAnnotation features) {
        myName = name;
        myLength = length;
        int convChr = Integer.MAX_VALUE;
        try {
            convChr = Integer.parseInt(name);
        } catch (NumberFormatException ne) {
        }
        myChromosomeNumber = convChr;
        myGA = features;
        hashCode = calcHashCode();
    }

    public Chromosome(String name) {
        this(name, -1, null);
    }

    public String getName() {
        return myName;
    }

    /**
     * Returns the interger value of the chromosome (if name is not a number
     * then Integer.MAX_VALUE is returned)
     */
    public int getChromosomeNumber() {
        return myChromosomeNumber;
    }

    public int getLength() {
        return myLength;
    }

    @Override
    public Object[] getAnnotation(String annoName) {
        return myGA.getAnnotation(annoName);
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

    @Override
    public String toString() {
        return getName();
    }

    @Override
    public int hashCode() {
        return hashCode;
    }

    private int calcHashCode() {
        int hash = 7;
        hash = 79 * hash + (this.myName != null ? this.myName.hashCode() : 0);
        return hash;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == this) {
            return true;
        }
        if (!(obj instanceof Chromosome)) {
            return false;
        }
        return (compareTo((Chromosome) obj) == 0);
    }

    @Override
    public int compareTo(Chromosome o) {
        //int result=Integer.compare(myChromosomeNumber,o.getChromosomeNumber());
        int result = Integer.valueOf(myChromosomeNumber).compareTo(Integer.valueOf(o.getChromosomeNumber()));
        if (result != 0) {
            return result;
        }
        return myName.compareTo(o.getName());
    }
}
