package net.maizegenetics.jGLiM;
//package net.maizegenetics.jGLiM;

import java.util.Date;

/**
 * Created by IntelliJ IDEA.
 * User: Peter Bradbury
 * Date: Feb 1, 2006
 * Time: 10:14:32 AM
 */
public interface Shuffler {
    void setRandomSeed(Date seedDate);
    public void resetRandomGenerator();
    int[] getNextPermutationIndex();
}
