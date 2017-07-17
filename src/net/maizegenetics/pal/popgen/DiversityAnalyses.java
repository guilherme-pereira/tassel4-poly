/*
 * TASSEL - Trait Analysis by a aSSociation Evolution & Linkage
 * Copyright (C) 2003 Ed Buckler
 *
 * This software evaluates linkage disequilibrium nucletide diversity and 
 * associations. For more information visit http://www.maizegenetics.net
 *
 * This software is distributed under GNU general public license and without
 * any warranty ot technical support.
 *
 * You can redistribute and/or modify it under the terms of GNU General 
 * public license. 
 *
 */
package net.maizegenetics.pal.popgen;

import net.maizegenetics.pal.alignment.*;
import net.maizegenetics.pal.report.TableReport;
import net.maizegenetics.pal.report.AbstractTableReport;

import java.io.Serializable;
import java.util.Vector;
import net.maizegenetics.pal.distance.IBSDistanceMatrix;

/**
 *This method calculated estimates of nucleotide diversity (pi, theta, etc).
 *
 * In order to provide high performance, it establishes one nucleotide SitePattern, and then
 * manipulates the weighting patterns for sliding windows and different types of sites.
 *
 *Total segregating sites needs to be adjusted for missing and gap data, as this will modify Theta and Tajima estimates
 *
 * @author Ed Buckler
 * @version 1.0
 */
public class DiversityAnalyses extends AbstractTableReport implements TableReport, Serializable {

    /** the limits of analysis */
    int startSite, endSite;
    /** values for the sliding window*/
    int step, window;
    /** number of site in the subsetted SitePattern*/
    int currNumSites;
    /** whether to use a sliding window*/
    boolean slideWindow = false;
    /** the base annotated alignment, it should be unprocessed (raw alignment)*/
    Alignment theAAlignment;
    /** the nature of the sites being evaluated, currently set to ALL POLYMORPHISMs, need to add specifics indels and SNPs.*/
    byte siteGroup = Alignment.POSITION_TYPE_ALL_GROUP;
    Vector diversityResultsVector = new Vector();
    PolymorphismDistribution thePolymorphismDistribution = null;

    /**
     * Constructor that uses a sliding window
     * @param aa an annotated alignment
     * @param sitePositionGroups a vector that contains the various site groups found in PositionTypeSubsetter
     * @param slidingWindow if true a sliding window will be run
     *
     */
    public DiversityAnalyses(Alignment aa, Vector sitePositionGroups, boolean slidingWindow, int start, int end, int window, int step,
            PolymorphismDistribution thePolymorphismDistribution) {
        this.startSite = start;
        this.endSite = end;
        this.step = step;
        this.window = window;
        this.slideWindow = slidingWindow;
        this.theAAlignment = aa;

        this.thePolymorphismDistribution = thePolymorphismDistribution;
        runAnalyses();
    }

    public DiversityAnalyses(Alignment aa, Vector sitePositionGroups, boolean slidingWindow, int start, int end, int window, int step) {
        this(aa, sitePositionGroups, slidingWindow, start, end, window, step, null);
    }

    /**
     * This will run the basic analyses and controls the sliding window analyses
     */
    private void runAnalyses() {
        if (!slideWindow) {
            runAnalysisForRegion(startSite, endSite);
            return;
        }
        //slide the window
        for (int i = startSite; i < (endSite - window); i += step) {
            runAnalysisForRegion(i, i + window - 1);
        }
    }

    /**
     * This will determine what analyses are to be run and run them
     */
    private void runAnalysisForRegion(int start, int end) {
        Locus locus = theAAlignment.getLocus(start);
        int chromosome = -1;
        try {
            chromosome = Integer.parseInt(locus.getName());
        } catch (Exception e) {
            //do nothing
        }
        double startChrPosition = theAAlignment.getPositionInLocus(start);
        double endChrPosition = theAAlignment.getPositionInLocus(end);
        Alignment theFilteredAlignment = FilterAlignment.getInstance(theAAlignment, start, end);
        IBSDistanceMatrix adm = new IBSDistanceMatrix(theFilteredAlignment);
        diversityResultsVector.add(evaluate(siteGroup, theFilteredAlignment, adm, start, end, chromosome, startChrPosition, endChrPosition));
        if (thePolymorphismDistribution != null) {
            thePolymorphismDistribution.addDistribution(Alignment.POSITION_TYPE_GROUP_TEXT[siteGroup] + "s" + start + "-e" + end, theFilteredAlignment, true);
        }
    }

    DiversityResults evaluate(byte siteGroup, Alignment theAlignment, IBSDistanceMatrix dm,
            int start, int end, int chromosome, double startChrPosition, double endChrPosition) {
        int sites = end - start + 1;
        DiversityResults theDiversityResults = new DiversityResults(start, end, siteGroup, chromosome, startChrPosition, endChrPosition);
        if (dm == null) {
            theDiversityResults.totalSites = 0;
            theDiversityResults.pipbp = Double.NaN;
            theDiversityResults.thetapbp = Double.NaN;
            theDiversityResults.segregatingSites = 0;
            return theDiversityResults;
        }
        double pipbp = dm.meanDistance();
        int segSites = countSegregatingSites(theAlignment);
        int taxa = theAlignment.getSequenceCount();
        theDiversityResults.pipbp = pipbp;
        theDiversityResults.avgSiteCoverage = dm.getAverageTotalSites();
        theDiversityResults.totalSites = sites;
        theDiversityResults.segregatingSites = segSites;
        theDiversityResults.thetapbp = estimateThetaPerbp(segSites, sites, theDiversityResults.avgSiteCoverage, taxa);

        // theDiversityResults.theta=estimateTheta(segSites,sites,theDiversityResults.avgSiteCoverage, theAlignment.getSequenceCount());
        theDiversityResults.tajimaD = estimateTajimaD(segSites, theDiversityResults.totalSites, theDiversityResults.avgSiteCoverage,
                theAlignment.getSequenceCount(), theDiversityResults.pipbp, theDiversityResults.thetapbp);
        return theDiversityResults;
    }

    public static double estimateTheta(int segSites, int totalSites, double averageSiteCoverage, int taxa) {
        double a = 0.0;
        double n = taxa * averageSiteCoverage / totalSites;
        for (double i = 1; i < n; i += 1.0) {
            a += 1 / i;
        }
        return segSites / a;
    }

    public static double estimateThetaPerbp(int segSites, int totalSites, double averageSiteCoverage, int taxa) {
        return estimateTheta(segSites, totalSites, averageSiteCoverage, taxa) / totalSites;
    }

    public static double estimatePi(int totalSites, double avgPairwiseDivergence, double averageSiteCoverage) {
        //this function is only needed to clearly show what we are doing with missing sites
        return totalSites * avgPairwiseDivergence / averageSiteCoverage;
    }

    public static double estimatePiPerbp(int totalSites, double avgPairwiseDivergence, double averageSiteCoverage) {
        //this function is only needed to clearly show what we are doing with missing sites
        return estimatePi(totalSites, avgPairwiseDivergence, averageSiteCoverage) / averageSiteCoverage;
    }

    public static double estimateTajimaD(int segSites, double totalSites, double averageSiteCoverage, double taxa, double pipbp, double thetapbp) {
        //this is the Tajima 1989 formula, but I don't know how it will deal with the
        //differential missing data
        double a1 = 0.0, a2 = 0.0;
        double n = taxa * averageSiteCoverage / totalSites;
        for (double i = 1; i < n; i += 1.0) {
            a1 += (1 / i);
            a2 += (1 / (i * i));
        }
        double b1 = (n + 1) / (3 * (n - 1));
        double b2 = (2 * (n * n + n + 3)) / (9 * n * (n - 1));
        double c1 = b1 - (1 / a1);
        double c2 = b2 - ((n + 2) / (a1 * n)) + (a2 / (a1 * a1));
        double e1 = c1 / a1;
        double e2 = c2 / (a1 * a1 + a2);
        double D = (pipbp - thetapbp) / (Math.sqrt((e1 * segSites) + (e2 * segSites * (segSites - 1))) / totalSites);
        return D;
    }

    int countSegregatingSites(Alignment theAlignment) {
        int total = 0;
        if (theAlignment.isAllPolymorphic()) {
            return theAlignment.getSiteCount();
        }
        for (int i = 0; i < theAlignment.getSiteCount(); i++) {
            if (theAlignment.isPolymorphic(i)) {
                total++;
            }
        }
        return total;
    }

    //Implementation of TableReport Interface
    public Object[] getTableColumnNames() {
        String[] basicLabels = {"Site_Type", "Chromosome", "StartChrPosition", "EndChrPosition", "StartSite", "EndSite", "MidSite",
            "SiteCount", "AvgSiteCount", "SegSites", "PiPerBP",
            "ThetaPerBP", "Haplotypes", "TajimaD"};
        return basicLabels;
    }

    /**
     * Returns specified row.
     *
     * @param row row number
     *
     * @return row
     */
    public Object[] getRow(int row) {

        Object[] data;
        java.text.NumberFormat nf = new java.text.DecimalFormat();
        nf.setMaximumFractionDigits(5);
        java.text.NumberFormat nf2 = new java.text.DecimalFormat();
        nf2.setMaximumFractionDigits(1);
        int basicCols = 14, labelOffset;
        data = new String[basicCols];
        DiversityResults theDiversityResults;
        theDiversityResults = (DiversityResults) diversityResultsVector.get(row);
        labelOffset = 0;
        data[labelOffset++] = theDiversityResults.getSiteGroupString();
        data[labelOffset++] = "" + theDiversityResults.chromosome;
        data[labelOffset++] = "" + nf2.format(theDiversityResults.startChrPosition);
        data[labelOffset++] = "" + nf2.format(theDiversityResults.endChrPosition);
        data[labelOffset++] = "" + theDiversityResults.startSite;
        data[labelOffset++] = "" + theDiversityResults.endSite;
        data[labelOffset++] = "" + (theDiversityResults.startSite + theDiversityResults.endSite) / 2;
        data[labelOffset++] = "" + theDiversityResults.totalSites;
        data[labelOffset++] = "" + nf.format(theDiversityResults.avgSiteCoverage);
        data[labelOffset++] = "" + theDiversityResults.segregatingSites;
        data[labelOffset++] = "" + nf.format(theDiversityResults.pipbp);
        data[labelOffset++] = "" + nf.format(theDiversityResults.thetapbp);
        data[labelOffset++] = "NotAvail";//+theDiversityResults.haplotypes;
        data[labelOffset++] = "" + nf.format(theDiversityResults.tajimaD);
        return data;

    }

    public String getTableTitle() {
        return "Diversity estimates";
    }

    public String toString() {
        if (diversityResultsVector.size() == 0) {
            return "Needs to be run";
        }
        StringBuffer cs = new StringBuffer();
        Object[] header = this.getTableColumnNames();
        for (int i = 0; i < header.length; i++) {
            cs.append(header[i]);
        }
        cs.append("\n");
        Object[][] data = this.getTableData();
        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[i].length; j++) {
                cs.append(data[i][j]);
            }
            cs.append("\n");
        }
        return cs.toString();
    }

    public int getRowCount() {
        return diversityResultsVector.size();
    }

    public int getElementCount() {
        throw new UnsupportedOperationException();
    }

    public int getColumnCount() {
        throw new UnsupportedOperationException();
    }

    void testCode() {
        int n = 10, sites = 41, segSites = 16;
        double pi = 3.88;
        double theta = DiversityAnalyses.estimateTheta(segSites, sites, sites, n);
        double pipbp = pi / sites;
        double td = DiversityAnalyses.estimateTajimaD((int) segSites, sites, sites, n, pipbp, theta);

    }
}

class DiversityResults implements Serializable {

    protected double pipbp, thetapbp, totalSites, avgSiteCoverage, tajimaD, startChrPosition, endChrPosition;
    protected int startSite, endSite, haplotypes, segregatingSites, chromosome;
    protected int siteGroup;
    private int index;

    public DiversityResults(int start, int end, byte siteGroup, int chromosome,
            double startChrPosition, double endChrPosition) {
        this.startSite = start;
        this.endSite = end;
        this.siteGroup = siteGroup;
        this.chromosome = chromosome;
        this.startChrPosition = startChrPosition;
        this.endChrPosition = endChrPosition;
    }

    public boolean equals(Object anObject) {
        DiversityResults x = (DiversityResults) anObject;
        if (x.index == this.index) {
            return true;
        } else {
            return false;
        }
    }

    public void setIndex(int theIndex) {
        this.index = theIndex;
    }

    public String getSiteGroupString() {
        return Alignment.POSITION_TYPE_GROUP_TEXT[siteGroup];
    }

    public String toString() {
        String result = "Diversity Results for " + getSiteGroupString() + "\n";
        result += "Pi =" + pipbp + "\n";
        result += "Theta =" + thetapbp + "\n";
        result += "Segregrating Sites =" + segregatingSites + "\n";
        result += "Total Sites =" + totalSites + "\n";
        result += "Start Site =" + startSite + "\n";
        result += "End Site =" + endSite + "\n";
        return result;
    }
}
