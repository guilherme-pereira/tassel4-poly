/*
 *  SiteScore
 */
package net.maizegenetics.pal.alignment.score;

import net.maizegenetics.pal.alignment.AlignmentNew;

/**
 *
 * @author Terry Casstevens
 */
public interface SiteScore {

    /**
     * Returns the site score of the given sequence and site.
     *
     * @param seq sequence index
     * @param site site
     *
     * @return site score.
     */
    public float getSiteScore(int seq, int site);

    /**
     * Returns the site scores.
     *
     * @return site scores.
     */
    public float[][] getSiteScores();

    /**
     * Returns true if this alignment has site scores.
     *
     * @return true if this alignment has site scores.
     */
    public boolean hasSiteScores();

    /**
     * Return what type of site scores this alignment has.
     *
     * @return site score type.
     */
    public AlignmentNew.SITE_SCORE_TYPE getSiteScoreType();
}
