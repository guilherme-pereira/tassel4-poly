/*
 * Locus
 */
package net.maizegenetics.pal.alignment;

import com.google.common.collect.ComparisonChain;
import java.io.Serializable;

import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

/**
 *
 * @author terry
 */
public class Locus implements Serializable, Comparable<Locus>{

    private static final long serialVersionUID = -5197800047652332969L;
    public static Locus UNKNOWN = new Locus("Unknown");
    private final String myName;
    private final String myChromosome;
    private final int myChromosomeNumber;
    private final int myStart;
    private final int myEnd;
    private final Map<String, Integer> myFeatures = new HashMap<String, Integer>();

    public Locus(String name, String chromosome, int start, int end, String[] featureNames, int[] featurePositions) {

        myName = name;
        myChromosome = chromosome;
        myStart = start;
        myEnd = end;
        int convChr=Integer.MAX_VALUE;
        try{convChr=Integer.parseInt(myChromosome);
        }
        catch(NumberFormatException ne) {}
        myChromosomeNumber=convChr;

        if ((featureNames == null) || (featurePositions == null)) {
            return;
        }
        if ((featureNames != null) && (featurePositions != null) && (featureNames.length != featurePositions.length)) {
            throw new IllegalArgumentException("Locus: init: number of feature names should equals number of feature positions.");
        }

        for (int i = 0, n = featureNames.length; i < n; i++) {
            myFeatures.put(featureNames[i], featurePositions[i]);
        }

    }

    public static Locus getMergedInstance(Locus locus1, Locus locus2) {
        String chromosome = locus1.getChromosomeName();
        if (!chromosome.equals(locus2.getChromosomeName())) {
            throw new IllegalArgumentException("Locus: getInstance: Chromosome Names must be the same.  locus 1: " + chromosome + "  locus 2: " + locus2.getChromosomeName());
        }
        int start = Math.min(locus1.getStart(), locus2.getStart());
        int end = Math.max(locus1.getEnd(), locus2.getEnd());
        Map features1 = locus1.getFeatures();
        Map features2 = locus2.getFeatures();
        String[] featureNames = new String[features1.size() + features2.size()];
        int[] featurePositions = new int[features1.size() + features2.size()];
        int count = 0;
        Iterator itr = locus1.getFeatures().keySet().iterator();
        while (itr.hasNext()) {
            featureNames[count] = (String) itr.next();
            featurePositions[count] = (Integer) features1.get(featureNames[count]);
            count++;
        }
        itr = locus2.getFeatures().keySet().iterator();
        while (itr.hasNext()) {
            featureNames[count] = (String) itr.next();
            featurePositions[count] = (Integer) features2.get(featureNames[count]);
            count++;
        }
        return new Locus(chromosome, chromosome, start, end, featureNames, featurePositions);
    }

    public Locus(String name) {
        this(name, name, -1, -1, null, null);
    }

    public String getName() {
        return myName;
    }

    public int getStart() {
        return myStart;
    }

    public int getEnd() {
        return myEnd;
    }

    public String getChromosomeName() {
        return myChromosome;
    }
    
    /**
     * Returns the interger value of the chromosome (if name is not a number then
     * Integer.MAX_VALUE is returned)
     */
    public int getChromosomeNumber() {
        return myChromosomeNumber;
    }

    public boolean isChromosome() {
        if (myChromosome == null) {
            return false;
        }
        return myChromosome.equals(myName);
    }

    public Map getFeatures() {
        return Collections.unmodifiableMap(myFeatures);
    }

    @Override
    public String toString() {
        return getName();
    }

    public boolean equalName(Object obj) {
        if (obj == this) {
            return true;
        }

        if (!(obj instanceof Locus)) {
            return false;
        }
        Locus other = (Locus) obj;

        if (!myName.equals(other.getName())) {
            return false;
        }

        if (!myChromosome.equals(other.getChromosomeName())) {
            return false;
        }
        return true;
    }
    
    
    @Override
    public boolean equals(Object obj) {

        if (obj == this) {
            return true;
        }

        if (!(obj instanceof Locus)) {
            return false;
        }
        Locus other = (Locus) obj;
        if (myChromosomeNumber != other.getChromosomeNumber()) {
            return false;
        }
        if (!myName.equals(other.getName())) {
            return false;
        }

        if (!myChromosome.equals(other.getChromosomeName())) {
            return false;
        }

        if (myStart != other.getStart()) {
            return false;
        }

        if (myEnd != other.getEnd()) {
            return false;
        }

        Map<String, Integer> otherFeatures = other.getFeatures();
        if (myFeatures.size() != otherFeatures.size()) {
            return false;
        }

        Iterator itr = myFeatures.keySet().iterator();
        while (itr.hasNext()) {
            String key = (String) itr.next();
            Integer value = myFeatures.get(key);
            Integer otherValue = otherFeatures.get(key);
            if ((otherValue == null) && (value != otherValue)) {
                return false;
            }
        }

        return true;
    }

    @Override
    public int compareTo(Locus o) {
        return ComparisonChain.start()
                .compare(myChromosomeNumber,o.getChromosomeNumber())
                .compare(myName,o.getName())
                .compare(myChromosome,o.getChromosomeName())
                .compare(myStart,o.getStart())
                .compare(myEnd, o.getEnd())
                .result();
    }
}
