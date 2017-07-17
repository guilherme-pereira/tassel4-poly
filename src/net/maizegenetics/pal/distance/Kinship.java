package net.maizegenetics.pal.distance;

import net.maizegenetics.pal.distance.DistanceMatrix;
import net.maizegenetics.pal.distance.IBSDistanceMatrix;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.AlignmentUtils;
import net.maizegenetics.pal.alignment.SimplePhenotype;

import java.awt.*;

/**
 * Created by IntelliJ IDEA.
 * User: Zhiwu
 * Date: Apr 29, 2007
 * Time: 4:21:44 PM
 * To change this template use File | Settings | File Templates.
 */
public class Kinship extends DistanceMatrix {

    Frame parentFrame;
    DistanceMatrix dm;
    Alignment mar;
    SimplePhenotype ped;
    int[][] parents;
    private double kMin = 99999;
    private double kMax = -99999;
    private double kAvg = 0;
    private double kSD = 0;
    private double cutOff = 2;
    private int numSeqs;

    public Kinship(Alignment mar) {
        this(mar, false, true);
    }

    public Kinship(Alignment mar, boolean areHetsRelated, boolean rescaleKinship) {
        this.mar = mar;
        numSeqs = this.mar.getIdGroup().getIdCount();
        buildFromMarker();
    }

    public Kinship(SimplePhenotype ped) {
        this.ped = ped;
        buildFromPed();
    }

    public Kinship(DistanceMatrix dm) {
        this.dm = dm;
    }

    public void buildFromMarker() {

    	IBSDistanceMatrix adm = new IBSDistanceMatrix(mar, 0, true, null);
    	dm = new DistanceMatrix(adm.getDistances(), mar.getIdGroup());
    	toSimilarity();
    	getKStatistics();
//    	pullBackExtrem();
    	//cutOff();
    	rescale();
    	System.out.println("Kinship was built from markers");
    }

    public void buildFromPed() {
        // get data from ped (SimplePhenotype) to parents (int[][]);
        //to do;
        System.out.println("Building Kinship From pedigree");
        parents = new int[ped.getNumberOfTaxa()][ped.getNumberOfTraits()];
        try {
            for (int row = 0; row < ped.getNumberOfTaxa(); row++) {
                for (int col = 0; col < ped.getNumberOfTraits(); col++) {
                    parents[row][col] = (int) ped.getData(row, col);
                }
            }
        } catch (NumberFormatException e) {
            e.printStackTrace();
        }

        dm = new DistanceMatrix(kinshipRelation(parents), ped.getTaxa());

        System.out.println("Kinship was build from pedigree");
    }

    public static double[][] kinshipRelation(int[][] ped) {
        int n = ped.length;
        int maleParent;
        int femaleParent;
        double[][] aMatrix = new double[n][n];
        //System.out.println("size of ped: "+n);



        //initial: diagonal 1, 0 otherwise;
        for (int i = 0; i < n; i++) {
            aMatrix[i][i] = 1;
            for (int j = i + 1; j < n; j++) {
                aMatrix[i][j] = 0;
                aMatrix[j][i] = 0;
            }
        }

        System.out.println("initial: diagonal 1, 0 otherwise");
        for (int i = 0; i < n; i++) {
            //diagonals
            femaleParent = ped[i][1];
            maleParent = ped[i][2];
            if ((femaleParent > 0) && (maleParent > 0)) {
                aMatrix[i][i] = aMatrix[i][i] + .5 * aMatrix[maleParent - 1][femaleParent - 1];
            }
            //Off diagonals
            for (int j = i + 1; j < n; j++) {
                femaleParent = ped[j][1];
                maleParent = ped[j][2];

                if ((femaleParent > 0) && (maleParent > 0)) {
                    aMatrix[i][j] = .5 * (aMatrix[i][femaleParent - 1] + aMatrix[i][maleParent - 1]);
                } else if (maleParent > 0) {
                    aMatrix[i][j] = .5 * aMatrix[i][maleParent - 1];
                } else if (femaleParent > 0) {
                    aMatrix[i][j] = .5 * aMatrix[i][femaleParent - 1];
                } else {
                    //do nothing
                }
                aMatrix[j][i] = aMatrix[i][j];
            }
        }
        System.out.println("A matrix finished");


        return aMatrix;
    }

    //Convert distance to similarity
    //By Zhiwu Zhang
    public void toSimilarity() {
        double s;
        System.out.println("toSimilarity " + numSeqs);

        for (int i = 0; i < numSeqs; i++) {
            for (int j = i; j < numSeqs; j++) {
                s = cutOff - dm.getDistance(i, j);
                dm.setDistance(i, j, s);
                dm.setDistance(j, i, s);
            }
        }
        System.out.println("toSimilarity finish" + numSeqs);
    }

    public void getKStatistics() {
        //get average
        double total = 0;
        double totalsq = 0;
        double nk = numSeqs * (numSeqs - 1) / 2;
        for (int i = 0; i < numSeqs - 1; i++) {
            for (int j = i + 1; j < numSeqs; j++) {
                total += dm.getDistance(i, j);
                totalsq += (dm.getDistance(i, j) * dm.getDistance(i, j));
                if (dm.getDistance(i, j) < kMin) {
                    kMin = dm.getDistance(i, j);
                }
                if (dm.getDistance(i, j) > kMax) {
                    kMax = dm.getDistance(i, j);
                }
            }
        }
        kAvg = total / nk;
        kSD = Math.sqrt((totalsq - nk * kAvg * kAvg) / (nk - 1));
        System.out.println(kAvg);
    }

    public void pullBackExtrem() {
        //take values beyond 3 sd from mean back
        //By Zhiwu Zhang
        for (int i = 0; i < numSeqs - 1; i++) {
            for (int j = i + 1; j < numSeqs; j++) {
                if (dm.getDistance(i, j) < kAvg - cutOff * kSD) {
                    dm.setDistance(i, j, kAvg - cutOff * kSD);
                    kMin = dm.getDistance(i, j);
                }
            }
        }
        System.out.println("values beyond 3 sd from mean were pulled back");
    }

    public void cutOff() {
        //Set vale to 0 if below than avg
        //By Zhiwu Zhang
        double s;
        for (int i = 0; i < numSeqs; i++) {
            for (int j = i + 0; j < numSeqs; j++) {
                if (dm.getDistance(i, j) < kAvg) {
                    dm.setDistance(i, j, kAvg);
                }
            }
        }
        kMin = kAvg;
    }

    public void rescale() {
        //rescale from theMin~2 to 0~2
        //By Zhiwu Zhang
        double s;
        for (int i = 0; i < numSeqs; i++) {
            for (int j = i; j < numSeqs; j++) {
                s = (dm.getDistance(i, j) - kMin) * cutOff / (cutOff - kMin);
                dm.setDistance(i, j, s);
                dm.setDistance(j, i, s);
            }
        }
        System.out.println("K rescaled");

    }

    public DistanceMatrix getDm() {
        return dm;
    }


}
