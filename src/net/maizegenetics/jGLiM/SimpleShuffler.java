package net.maizegenetics.jGLiM;
//package net.maizegenetics.jGLiM;

import cern.jet.random.engine.MersenneTwister;

import java.util.Date;
import java.util.Arrays;
import java.util.Comparator;

import net.maizegenetics.jGLiM.Shuffler;

/**
 * Created by IntelliJ IDEA.
 * User: Peter Bradbury
 * Date: Mar 31, 2006
 * Time: 1:40:00 PM
 * returns the integers 0 to n in random order
 */
public class SimpleShuffler implements Shuffler {
    Date seed;
    MersenneTwister randomizer;
    int numberOfElements;
    int[][] pairs;

    public SimpleShuffler(int n) {
        seed = new Date();
        resetRandomGenerator();
        numberOfElements = n;
        pairs = new int[numberOfElements][];
        for (int i=0;i<numberOfElements;i++) {
            pairs[i] = new int[]{i, 0};
        }

        for (int i=0;i<1000;i++) randomizer.nextInt();

    }

    public SimpleShuffler() {
        seed = new Date();
        resetRandomGenerator();
        for (int i=0;i<1000;i++) randomizer.nextInt();
    }
    
    public void setRandomSeed(Date seedDate) {
        seed = seedDate;
    }

    public void resetRandomGenerator() {
        randomizer = new MersenneTwister(seed);
    }

    public int[] getNextPermutationIndex() {
        return nextSequence();
    }

    private int[] nextSequence() {
        int[] result = new int[numberOfElements];

        for (int i = 0; i < numberOfElements; i++) {
            ((int[]) pairs[i])[1] = randomizer.nextInt();
        }
        Arrays.sort(pairs, new Comparator<int[]>() {
            public int compare(int[] a1, int[] a2) {
                if (a1[1] > a2[1]) return 1;
                if (a1[1] < a2[1]) return -1;
                return 0;
            }
        });
        for (int i = 0; i < numberOfElements; i++) result[i] = ((int[]) pairs[i])[0];

        return result;
    }

    public void setNumberOfElements(int numberOfElements) {
        this.numberOfElements = numberOfElements;
        pairs = new int[numberOfElements][];
        for (int i=0;i<numberOfElements;i++) {
            pairs[i] = new int[]{i, 0};
        }
    }


}
