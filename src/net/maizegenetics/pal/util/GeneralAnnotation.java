package net.maizegenetics.pal.util;

/**
 * Provide generalized annotations (descriptors) for a taxon or site.
 *
 */
public interface GeneralAnnotation {

    /**
     * Returns all annotation value for a given annotation key
     * @param annoName annotation key
     * @return array of annotation values (if not present new String[0])
     */
    public Object[] getAnnotation(String annoName);

    /**
     * Returns all annotation value for a given annotation key
     * @param annoName annotation key
     * @return array of annotation values (if not present new String[0])
     */
    public String[] getTextAnnotation(String annoName);

    /**
     * Returns all annotation value for a given annotation key
     * @param annoName annotation key
     * @return array of annotation values (if not present new String[0])
     */
    public String getConsensusAnnotation(String annoName);

    /**
     * Returns all annotation value for a given annotation key
     * @param annoName annotation key
     * @return array of annotation values (if not present new double[0])
     */
    public double[] getQuantAnnotation(String annoName);

    /**
     * Returns average annotation for a given annotation key
     * @param annoName annotation key
     * @return average value (if not present - return Double.NaN)
     */
    public double getAverageAnnotation(String annoName);

    //should we provide methods, to average the quantitative annotations, the first annotation
    //


}
