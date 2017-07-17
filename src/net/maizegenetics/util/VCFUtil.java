/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.util;

import java.util.HashMap;
import net.maizegenetics.pal.alignment.Alignment;

/**
 *
 * @author qs24, Gabriel Rodrigues Alves Margarido
 */
public class VCFUtil {
    // variables for calculating OS and PL for VCF, might not be in the correct class
    private static double error;
    private static double v1;
    private static double v2;
    private static double v3;
    private static int[][][] myGenoScoreMap;
    public static final int VCF_DEFAULT_MAX_NUM_ALLELES = 3;
        
    private VCFUtil ()
    {
        
    }
    
    private static void initVCFScoreMap() {
        error = 0.001;
        v1 = Math.log10(1.0 - error * 3.0 /4.0);
        v2 = Math.log10(error/4);
        v3 = Math.log10(0.5 - (error/4.0));
        myGenoScoreMap = new int[501][501][];
        for (int i = 0; i < 501; i++) {
            for (int j = 0; j < 501; j++) {
                myGenoScoreMap[i][j]= calcScore(i, j);
            }
        }
    }
    
    public static int[] getScore(int i, int j) {
        if (myGenoScoreMap == null) {
            initVCFScoreMap();
        }
        
        if (i < 501 && j < 501) {
            return myGenoScoreMap[i][j];
        } else {
            return calcScore(i, j);
        }
    }
    
    // Calculate QS and PL for VCF might not be in the correct class
    private static int[] calcScore (int a, int b)
    {   
        int[] results= new int[4];
        int n = a + b;
        int m = a;
        if (b > m) {
            m = b;
        }

        double fact = 0;
        if (n > m) {
            for (int i = n; i > m; i--) {
               fact += Math.log10(i);
            }
            for (int i = 1; i <= (n - m); i++){
               fact -= Math.log10(i);
            }
        }
        double aad = Math.pow(10, fact + (double)a * v1 + (double)b * v2);
        double abd = Math.pow(10, fact + (double)n * v3);
        double bbd = Math.pow(10, fact + (double)b * v1 + (double)a * v2);
        double md = aad;
        if (md < abd) {
            md = abd;
        }
        if (md < bbd) {
            md = bbd;
        }
        int gq = 0;
        if ((aad + abd + bbd) > 0) {
            gq = (int)(md / (aad + abd + bbd) * 100);
        }
        
        int aa =(int) (-10 * (fact + (double)a * v1 + (double)b * v2));
        int ab =(int) (-10 * (fact + (double)n * v3));
        int bb =(int) (-10 * (fact + (double)b * v1 + (double)a * v2));
        
        m = aa;
        if (m > ab) {
            m = ab;
        }
        if (m>bb) {
            m = bb;
        }
        aa -= m;
        ab -= m;
        bb -= m;
        results[0] = aa > 255 ? 255 : aa;
        results[1] = ab > 255 ? 255 : ab;
        results[2] = bb > 255 ? 255 : bb;
        results[3] = gq;
        
        return results;
    }
    
     public static byte resolveVCFGeno(byte[] alleles, int[][] allelesInTaxa, int tx) {
        int[] alleleDepth = new int[allelesInTaxa.length];
        for (int i=0; i<allelesInTaxa.length; i++)
        {
            alleleDepth[i] = allelesInTaxa[i][tx];
        }
        return resolveVCFGeno(alleles, alleleDepth);
    }
     
     public static byte resolveVCFGeno(byte[] alleles, int[] alleleDepth) { 
        int depth = 0;
        for (int i = 0; i < alleleDepth.length; i++) {
            depth += alleleDepth[i];
        }
        if (depth == 0) {
            return (byte)((Alignment.UNKNOWN_ALLELE << 4) | Alignment.UNKNOWN_ALLELE);
        }
        int max = 0;
        byte maxAllele = Alignment.UNKNOWN_ALLELE;
        int nextMax = 0;
        byte nextMaxAllele = Alignment.UNKNOWN_ALLELE;
        for (int i = 0; i < alleles.length; i++) {
            if (alleleDepth[i] > max) {
                nextMax = max;
                nextMaxAllele = maxAllele;
                max = alleleDepth[i];
                maxAllele = alleles[i];
            } else if (alleleDepth[i] > nextMax) {
                nextMax = alleleDepth[i];
                nextMaxAllele = alleles[i];
            }
        }
        if (alleles.length == 1) {
            return (byte)((alleles[0] << 4) | alleles[0]);
        } else {
            max = (max > 32767) ? 32767 : max;
            nextMax = (nextMax > 32767) ? 32767 : nextMax;
            int[] scores = getScore(max, nextMax);
            if ((scores[1] <= scores[0]) && (scores[1] <= scores[2])) {
                return (byte)((maxAllele << 4) | nextMaxAllele);
            } else if ((scores[0] <= scores[1]) && (scores[0] <= scores[2])) {
                return (byte)((maxAllele << 4) | maxAllele);
            } else {
                return (byte)((nextMaxAllele << 4) | nextMaxAllele);
            }
        }
     }
     

}
