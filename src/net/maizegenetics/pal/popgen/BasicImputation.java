/*
 * BasicImputation
 */
package net.maizegenetics.pal.popgen;

import net.maizegenetics.pal.alignment.Alignment;

/**
 *
 * @author ed
 */
public class BasicImputation {

    /**
     * Finds the maximum shared haplotype length for a pair of taxa centered around the
     * initialSite.  Note
     * @param align Alignment
     * @param initialSite the site at which the search begins
     * @param taxa1
     * @param taxa2
     * @param mismatch the number mismatching sites in the best match length
     * @return array of three ints with maxLength, leftStartSite, rightStartSite
     */
    public static int[] maxHaplotypeLength(Alignment align, int initialSite, int taxa1, int taxa2, int mismatch) {
        int[] hapDescription = new int[3];  //coded length, left start site, right start site
        //hapDescription[0]=-1;
        int b, rightMissCnt = 0, leftMissCnt = 0;
        byte s1b, s2b;
        boolean ibd;
        int[][] leftMis = new int[2][mismatch + 1];
        int[][] rightMis = new int[2][mismatch + 1];
        ibd = true;
        for (b = initialSite; (b >= 0) && ibd; b--) {
            s1b = align.getBase(taxa1, b);
            s2b = align.getBase(taxa2, b);
            if ((s1b != Alignment.UNKNOWN_DIPLOID_ALLELE) && (s2b != Alignment.UNKNOWN_DIPLOID_ALLELE)) {
                if (s1b == s2b) {
                    leftMis[0][leftMissCnt]++;
                    //hapDescription[1]=b;
                } else {
                    leftMis[1][leftMissCnt] = b;
                    leftMissCnt++;
                    if (leftMissCnt > mismatch) {
                        ibd = false;
                    } else {
                        leftMis[0][leftMissCnt] = leftMis[0][leftMissCnt - 1];
                    }
                }
            }
        }
        ibd = true;
        for (b = initialSite; (b < align.getSiteCount()) && ibd; b++) {
            s1b = align.getBase(taxa1, b);
            s2b = align.getBase(taxa2, b);
            if ((s1b != Alignment.UNKNOWN_DIPLOID_ALLELE) && (s2b != Alignment.UNKNOWN_DIPLOID_ALLELE)) {
                if (s1b == s2b) {
                    if (b != initialSite) {
                        rightMis[0][rightMissCnt]++;
                    }
                    //hapDescription[2]=b;
                } else {
                    rightMis[1][rightMissCnt] = b;
                    rightMissCnt++;
                    if (rightMissCnt > mismatch) {
                        ibd = false;
                    } else {
                        rightMis[0][rightMissCnt] = rightMis[0][rightMissCnt - 1];
                    }
                }
            }
        }
        for (int i = 0; i <= mismatch; i++) {
            int length = leftMis[0][mismatch - i] + rightMis[0][i];
            //   if(leftMis[1][0]==initialSite) length=0;  //mismatch at the initial site
            if (hapDescription[0] < length) {
                hapDescription[0] = length;
                hapDescription[1] = leftMis[1][mismatch - i];
                hapDescription[2] = rightMis[1][mismatch - i];
            }
        }
        return hapDescription;
    }

    /**
     *
     * NOTE in the return matrix I extend the length of the haplotype by 1 for all, so
     *there is a difference between missing for zero.
     * @param align
     * @param initialSite
     * @param maxMismatch
     * @param maskKnown if true, pair of sites where all is known are ignored and set
     * @return matrix of distances between
     */
    public static int[][] maxHaplotypeLengthMatrix(Alignment align, int initialSite, int maxMismatch, boolean maskKnown) {
        int[][] m = new int[align.getSequenceCount()][align.getSequenceCount()];
        for (int i = 0; i < align.getSequenceCount(); i++) {
            for (int j = 0; j < i; j++) {
                byte bj = align.getBase(j, initialSite);
                byte bi = align.getBase(i, initialSite);
                if (maskKnown && (bj != Alignment.UNKNOWN_DIPLOID_ALLELE) && (bi != Alignment.UNKNOWN_DIPLOID_ALLELE)) {
                    m[j][i] = m[i][j] = Integer.MIN_VALUE;
                } else {
                    m[j][i] = m[i][j] = maxHaplotypeLength(align, initialSite, i, j, maxMismatch)[0] + 1;
                    if (bj == Alignment.UNKNOWN_DIPLOID_ALLELE) {
                        m[i][j] *= -1;
                    }
                    if (bi == Alignment.UNKNOWN_DIPLOID_ALLELE) {
                        m[j][i] *= -1;
                    }
                }
            }
        }
        return m;
    }

    /**
     * This imputation approach looks for the longest matching haplotype surrounding
     * every missing SNP. It finds the greatest length by only counting sites where the two haplotypes
     * have non-missing data, and then will find the longest matching stretch in either direction with
     * the number of mismatches.
     *
     * The a taxa with a missing data has no match of the minLength, then the major allele is used to replace the SNP
     *
     * Currently this requires a Pack1Alignment as there not a duplication method for simple alignment.
     * @param align Alignment used
     * @param minLength threshold below which majority base is used.
     * @param mismatchNum the number of misMatching bases in the length search
     * @return Imputed Pack1Alignment
     */
    public static Alignment imputeBySite(Alignment align, int minLength, int mismatchNum) {
        throw new UnsupportedOperationException("BasicImputation: imputeBySite: Not Supported.");
    }
}
