package net.maizegenetics.jGLiM;
//package net.maizegenetics.jGLiM;

/**
 * Created by IntelliJ IDEA.
 * User: Peter Bradbury
 * Date: Apr 24, 2006
 * Time: 9:50:40 AM
 */
public interface ProgressIndicator {
    void setProgress(int percent);
    void incrementProgress(int numberOfSteps);
    void setTotalSteps(int totalSteps);
}
