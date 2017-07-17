/*
 * ProgressListener
 */

package net.maizegenetics.util;

/**
 *
 * @author terry
 */
public interface ProgressListener {

    /**
     * Returns progress of execution.
     *
     * @param percent percent complete
     * @param meta meta data
     */
    public void progress (int percent, Object meta);


}