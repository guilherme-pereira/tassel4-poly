package net.maizegenetics.gbs.homology;

/**
 * Created by IntelliJ IDEA.
 * User: ed
 * Date: May 31, 2008
 * Time: 9:35:16 AM
 * To change this template use File | Settings | File Templates.
 */

/*
 * SmithWaterman.java
 *
 * Copyright 2003 Sergio Anibal de Carvalho Junior
 *
 * This file is part of NeoBio.
 *
 * NeoBio is free software; you can redistribute it and/or modify it under the terms of
 * the GNU General Public License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * NeoBio is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE. See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with NeoBio;
 * if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330,
 * Boston, MA 02111-1307, USA.
 *
 * Proper attribution of the author as the source of the software would be appreciated.
 *
 * Sergio Anibal de Carvalho Junior		mailto:sergioanibaljr@users.sourceforge.net
 * Department of Computer Science		http://www.dcs.kcl.ac.uk
 * King's College London, UK			http://www.kcl.ac.uk
 *
 * Please visit http://neobio.sourceforge.net
 *
 * This project was supervised by Professor Maxime Crochemore.
 *
 */
import java.util.Arrays;
import org.biojava3.alignment.NeedlemanWunsch;

/**
 * This class implement the classic local alignment algorithm (with linear gap penalty
 * function) due to T.F.Smith and M.S.Waterman (1981).
 *
 * <P>This algorithm is very similar to the {@linkplain NeedlemanWunsch} algorithm for
 * global alignment. The idea here also consists of building an (n+1 x m+1) matrix M given
 * two sequences A and B of sizes n and m, respectively. However, unlike in the global
 * alignment case, every position M[i,j] in the matrix contains the similarity score of
 * <B>suffixes</B> of A[1..i] and B[1..j].</P>
 *
 * <P>Starting from row 0, column 0, the {@link #computeMatrix computeMatrix} method
 * computes each position M[i,j] with the following recurrence:</P>
 *
 * <CODE><BLOCKQUOTE><PRE>
 * M[0,0] = <B>M[0,j]</B> = <B>M[i,0]</B> = 0
 * M[i,j] = max { M[i,j-1]   + scoreInsertion (B[j]),
 *                M[i-1,j-1] + scoreSubstitution (A[i], B[j]),
 *                M[i-1,j]   + scoreDeletion(A[i])             }
 * </PRE></BLOCKQUOTE></CODE>
 *
 * <P>Note that, here, all cells in the first row and column are set to zero. The best
 * local alignment score is the highest value found anywhere in the matrix.</P>
 *
 * <P>Just like in global alignment case, this algorithm has quadratic space complexity
 * because it needs to keep an (n+1 x m+1) matrix in memory. And since the work of
 * computing each cell is constant, it also has quadratic time complexity.</P>
 *
 * <P>After the matrix has been computed, the alignment can be retrieved by tracing a path
 * back in the matrix from the position of the highest score until a cell of value zero is
 * reached. This step is performed by the {
 * buildOptimalAlignment} method, and its time complexity is linear on the size of the
 * alignment.
 *
 * <P>If the similarity value only is needed (and not the alignment itself), it is easy to
 * reduce the space requirement to O(n) by keeping just the last row or column in memory.
 * This is precisely what is done by the {@link #computeScore computeScore} method. Note
 * that it still requires O(n<SUP>2</SUP>) time.</P>
 *
 * <P>For a more efficient approach to the local alignment problem, see the
 * {@linkplain CrochemoreLandauZivUkelson} algorithm. For global alignment, see the
 * {@linkplain NeedlemanWunsch} algorithm.</P>
 *
 * @author Sergio A. de Carvalho Jr.

 */
public class SmithWaterman {

    /**
     * The first sequence of an alignment.
     */
    protected String seq1 = "CGGGTGTGACAGTCGTGCAGTCGACCGTTGGG";
    /**
     * The second sequence of an alignment.
     */
    protected String seq2 = "XXXXXCGGGTGTGACAGTCGTGCAGTCGACCGTTGGGXXXXXXX";
//   protected String seq2="CGGGTGTGACAGTCGTGCAGTCGACCGTTGGG";
    /**
     * The dynamic programming matrix. Each position (i, j) represents the best score
     * between a suffic of the firsts i characters of <CODE>seq1</CODE> and a suffix of
     * the first j characters of <CODE>seq2</CODE>.
     */
    protected int[][] matrix;
    /**
     * Indicate the row of where an optimal local alignment can be found in the matrix..
     */
    protected int max_row;
    /**
     * Indicate the column of where an optimal local alignment can be found in the matrix.
     */
    protected int max_col;
    int rows, cols;
    byte[] bseq1, bseq2;
    int[] array;

    public SmithWaterman() {
        long time = System.currentTimeMillis();
        int maxScore = 0;
        rows = seq1.length() + 1;
        cols = seq2.length() + 1;
        if (rows <= cols) {
            array = new int[rows];
        } else {
            array = new int[cols];
        }
        bseq1 = seq1.getBytes();
        bseq2 = seq2.getBytes();
        //        matrix = new int[rows][cols];
        for (int i = 0; i < 30000; i++) {
            //      maxScore=computeMatrix();
            maxScore = computeScore();
        }
        System.out.println(Arrays.deepToString(matrix));
        System.out.println("MaxScore:" + maxScore);
        System.out.println("time:" + (System.currentTimeMillis() - time));
    }

    public SmithWaterman(int rows, int cols) {
        this.rows = rows + 1;
        this.cols = cols + 1;
        int max = (rows > cols) ? rows : cols;
        array = new int[max];
        //       if (rows <= cols) {array = new int [this.rows];}
        //          else{array = new int [this.cols];}
    }

    int computeScore(byte[] b1, byte[] b2) {
        bseq1 = b1;
        bseq2 = b2;
        this.rows = b1.length + 1;
        this.cols = b2.length + 1;
        return computeScore();
    }

    /**
     * Computes the dynamic programming matrix.
     *
     *  If the scoring scheme is not compatible
     * with the loaded sequences.
     */
    protected int computeMatrix() {
        int r, c, ins, sub, del, max_score;

        // initiate first row
        for (c = 0; c < cols; c++) {
            matrix[0][c] = 0;
        }

        // keep track of the maximum score
        this.max_row = this.max_col = max_score = 0;

        // calculates the similarity matrix (row-wise)
        for (r = 1; r < rows; r++) {
            // initiate first column
            matrix[r][0] = 0;

            for (c = 1; c < cols; c++) {
//				ins = matrix[r][c-1] + scoreInsertion(seq2.charAt(c));
//				sub = matrix[r-1][c-1] + scoreSubstitution(seq1.charAt(r),seq2.charAt(c));
//				del = matrix[r-1][c] + scoreDeletion(seq1.charAt(r));

                ins = matrix[r][c - 1] - 1;
//				sub = matrix[r-1][c-1] + scoreSubstitution(seq1.charAt(r-1),seq2.charAt(c-1));
                sub = matrix[r - 1][c - 1] + (bseq1[r - 1] == bseq2[c - 1] ? 2 : 0);
                del = matrix[r - 1][c] - 1;

                // choose the greatest
                matrix[r][c] = max(ins, sub, del, 0);

                if (matrix[r][c] > max_score) {
                    // keep track of the maximum score
                    max_score = matrix[r][c];
                    this.max_row = r;
                    this.max_col = c;
                }
            }
        }
//        System.out.println("MaxSore:"+max_score);
        return max_score;
    }

    /**
     * Computes the score of the best local alignment between the two sequences using the
     * scoring scheme previously set. This method calculates the similarity value only
     * (doesn't build the whole matrix so the alignment cannot be recovered, however it
     * has the advantage of requiring O(n) space only).
     *
     * @return the score of the best local alignment between the loaded sequences
     * with the loaded sequences.
     */
    protected int computeScore() {
//		int[]	array;
//		int 	rows = seq1.length()+1, cols = seq2.length()+1;
        int r, c, tmp, ins, del, sub, max_score;

        // keep track of the maximum score
        max_score = 0;

        if (rows <= cols) {
//			// goes columnwise
//			array = new int [rows];

            // initiate first column
            for (r = 0; r < rows; r++) {
                array[r] = 0;
            }

            // calculate the similarity matrix (keep current column only)
            for (c = 1; c < cols; c++) {
                // set first position to zero (tmp hold values
                // that will be later moved to the array)
                tmp = 0;

                for (r = 1; r < rows; r++) {
                    ins = array[r] - 1;
//                 if(r>30) {
//                     System.out.println();
//                 }
                    sub = array[r - 1] + (bseq1[r - 1] == bseq2[c - 1] ? 2 : 0);
                    del = tmp - 1;

                    // move the temp value to the array
                    array[r - 1] = tmp;

                    // choose the greatest (or zero if all negative)
                    tmp = max(ins, sub, del, 0);

                    // keep track of the maximum score
                    if (tmp > max_score) {
                        max_score = tmp;
                    }
                }

                // move the temp value to the array
                array[rows - 1] = tmp;
            }
        } else {
            // goes rowwise
//			array = new int [cols];

            // initiate first row
            for (c = 0; c < cols; c++) {
                array[c] = 0;
            }

            // calculate the similarity matrix (keep current row only)
            for (r = 1; r < rows; r++) {
                // set first position to zero (tmp hold values
                // that will be later moved to the array)
                tmp = 0;

                for (c = 1; c < cols; c++) {
                    ins = tmp - 1;
                    sub = array[c - 1] + (bseq1[r - 1] == bseq2[c - 1] ? 2 : 0);
                    del = array[c] - 1;

                    // move the temp value to the array
                    array[c - 1] = tmp;

                    // choose the greatest (or zero if all negative)
                    tmp = max(ins, sub, del, 0);

                    // keep track of the maximum score
                    if (tmp > max_score) {
                        max_score = tmp;
                    }
                }

                // move the temp value to the array
                array[cols - 1] = tmp;
            }
        }

        return max_score;
    }

    protected final int scoreInsertion(char a) {
        return -1;

    }

    protected final int scoreSubstitution(char a, char b) {
        if (a == b) {
            return 2;
        }
        return 0;

    }

    protected final int scoreDeletion(char a) {
        return -1;

    }

    /**
     * Helper method to compute the the greater of two values.
     *
     * @param v1 first value
     * @param v2 second value
     * @return the larger of <CODE>v1</CODE> and <CODE>v2</CODE>
     */
    protected final int max(int v1, int v2) {
        return (v1 >= v2) ? v1 : v2;
    }

    /**
     * Helper method to compute the the greater of three values.
     *
     * @param v1 first value
     * @param v2 second value
     * @param v3 third value
     * @return the larger of <CODE>v1</CODE>, <CODE>v2</CODE> and <CODE>v3</CODE>
     */
    protected final int max(int v1, int v2, int v3) {
        return (v1 >= v2) ? ((v1 >= v3) ? v1 : v3) : ((v2 >= v3) ? v2 : v3);
    }

    /**
     * Helper method to compute the the greater of four values.
     *
     * @param v1 first value
     * @param v2 second value
     * @param v3 third value
     * @param v4 fourth value
     * @return the larger of <CODE>v1</CODE>, <CODE>v2</CODE> <CODE>v3</CODE> and
     * <CODE>v4</CODE>
     */
    protected final int max(int v1, int v2, int v3, int v4) {
        int m1 = ((v1 >= v2) ? v1 : v2);
        int m2 = ((v3 >= v4) ? v3 : v4);

        return (m1 >= m2) ? m1 : m2;
    }
}
