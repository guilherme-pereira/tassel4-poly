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
package net.maizegenetics.stats.logisticReg;

import net.maizegenetics.pal.report.TableReport;
import net.maizegenetics.pal.report.AbstractTableReport;
import net.maizegenetics.util.RightJustifiedFormat;

import java.io.Serializable;
import java.util.Vector;


/**
 * Title:        TASSEL
 * Description:  calculate permutations across environments, across locusjnmhhhhhhhhhhhhhhhhhhh
 * @author Ed Buckler
 * @version 1.0
 */
public class LocusEnvironmentPermutation extends AbstractTableReport implements TableReport, Serializable {

    int traitCol, environCol, siteCol, locusCol;
    Vector origStats, permStats, traits, environments, locusSites;
    Vector permutationResultsVector = new Vector();
    double[][][] origTestStats;


    //calculate permutations over locs
    public LocusEnvironmentPermutation(Vector sr) {
        origStats = sr;
    }

    public void setCol(int tc, int ec, int lc, int sc) {
        traitCol = tc;
        environCol = ec;
        locusCol = lc;
        siteCol = sc;
        traits = new Vector();
        locusSites = new Vector();
        environments = new Vector();
        String ct, en, ls;
        StatisticResult sr;
        for (int i = 0; i < origStats.size(); i++) {
            sr = (StatisticResult) origStats.get(i);
            ct = sr.getAttribute(traitCol);
            if (!traits.contains(ct)) {
                traits.add(ct);
            }
            en = sr.getAttribute(environCol);
            if (!environments.contains(en)) {
                environments.add(en);
            }
            ls = sr.getAttribute(locusCol) + "#" + sr.getAttribute(siteCol);
            if (!locusSites.contains(ls)) {
                locusSites.add(ls);
            }
        }
        origTestStats = new double[traits.size()][environments.size()][locusSites.size()];
        double[][][] origTestPValues = new double[traits.size()][environments.size()][locusSites.size()];
        for (int i = 0; i < origStats.size(); i++) {
            sr = (StatisticResult) origStats.get(i);
//      if(!Double.isNaN(sr.testLikelihood()))
            //{
            origTestStats[traits.indexOf(sr.getAttribute(traitCol))][environments.indexOf(sr.getAttribute(environCol))][locusSites.indexOf(sr.getAttribute(locusCol) + "#" + sr.getAttribute(siteCol))] = sr.testLikelihood();
            origTestPValues[traits.indexOf(sr.getAttribute(traitCol))][environments.indexOf(sr.getAttribute(environCol))][locusSites.indexOf(sr.getAttribute(locusCol) + "#" + sr.getAttribute(siteCol))] = sr.testPValue();
//        }
        }
//    TestStatAndSite[] ssoebt=siteSumOverEnvironmentByTrait(origTestStats);
        createPermutationResults("Site_Statistic_Summed_Over_Environment_By_Trait", siteSumOverEnvironmentByTrait(origTestStats));
//    TestStatAndSite ssoeat=siteSumOverEnvironmentAndTrait(origTestStats);
        createPermutationResults("Site_Statistic_Summed_Over_Environment_And_Trait", siteSumOverEnvironmentAndTrait(origTestStats));
//    TestStatAndSite[] soebt=siteOverEnvironmentByTrait(origTestStats);
        createPermutationResults("Site_Statistic_Over_Environment_By_Trait", siteOverEnvironmentByTrait(origTestStats));
//    TestStatAndSite soeat=siteOverEnvironmentAndTrait(origTestStats);
        createPermutationResults("Site_Statistic_Over_Environment_And_Trait", siteOverEnvironmentAndTrait(origTestStats));
        createPermutationResults("SumSite_Statistic_Over_Environment_And_Trait", siteOverEnvironmentAndTraitAddAllLoci(origTestStats));
        createPermutationResults("Proportion<0.100", proportionSignificant(origTestPValues, 0.1));
        createPermutationResults("Proportion<0.010", proportionSignificant(origTestPValues, 0.01));
        createPermutationResults("Proportion<0.001", proportionSignificant(origTestPValues, 0.001));
    }

    private void createPermutationResults(String category, TestStatAndSite[] tsas) {
        for (int i = 0; i < tsas.length; i++) {
            String[] annotation = {category, tsas[i].trait, tsas[i].locus, tsas[i].site};
            PermutationResults pr = new PermutationResults(annotation, tsas[i].testStat);
            permutationResultsVector.add(pr);
        }
    }

    private void createPermutationResults(String category, TestStatAndSite tsas) {
        TestStatAndSite[] ntsas = {tsas};
        createPermutationResults(category, ntsas);
    }
    /*
    public void permutation(Vector pr) {
    double[][][] permTestStats=new double[traits.size()][environments.size()][locusSites.size()];
    for(int i=0; i<pr.size(); i++)
    {sr=(StatisticResult)pr.get(i);
    permTestStats[traits.indexOf(sr.getAttribute(traitCol))]
    [environments.indexOf(sr.getAttribute(environCol))]
    [locusSites.indexOf(sr.getAttribute(locusCol)+"#"+sr.getAttribute(siteCol))]=sr.testLikelihood();
    }
    updatePermutationResults("Site_Statistic_Summed_Over_Environment_By_Trait", siteSumOverEnvironmentByTrait(permTestStats));
    updatePermutationResults("Site_Statistic_Summed_Over_Environment_And_Trait", siteSumOverEnvironmentAndTrait(permTestStats));
    updatePermutationResults("Site_Statistic_Over_Environment_By_Trait", siteOverEnvironmentByTrait(permTestStats));
    updatePermutationResults("Site_Statistic_Over_Environment_And_Trait", siteOverEnvironmentAndTrait(permTestStats));
    }
     */

    public void addPermutation(Vector psr) {
        double[][][] permTestStats = new double[traits.size()][environments.size()][locusSites.size()];
        double[][][] permTestPValues = new double[traits.size()][environments.size()][locusSites.size()];
        StatisticResult sr;
        for (int i = 0; i < psr.size(); i++) {
            sr = (StatisticResult) psr.get(i);
            permTestStats[traits.indexOf(sr.getAttribute(traitCol))][environments.indexOf(sr.getAttribute(environCol))][locusSites.indexOf(sr.getAttribute(locusCol) + "#" + sr.getAttribute(siteCol))] = sr.testLikelihood();
            permTestPValues[traits.indexOf(sr.getAttribute(traitCol))][environments.indexOf(sr.getAttribute(environCol))][locusSites.indexOf(sr.getAttribute(locusCol) + "#" + sr.getAttribute(siteCol))] = sr.testPValue();
        }
        updatePermutationResults("Site_Statistic_Summed_Over_Environment_By_Trait", siteSumOverEnvironmentByTrait(permTestStats));
        updatePermutationResults("Site_Statistic_Summed_Over_Environment_And_Trait", siteSumOverEnvironmentAndTrait(permTestStats));
        updatePermutationResults("Site_Statistic_Over_Environment_By_Trait", siteOverEnvironmentByTrait(permTestStats));
        updatePermutationResults("Site_Statistic_Over_Environment_And_Trait", siteOverEnvironmentAndTrait(permTestStats));
        updatePermutationResults("SumSite_Statistic_Over_Environment_And_Trait", siteOverEnvironmentAndTraitAddAllLoci(permTestStats));
        updatePermutationResults("Proportion<0.100", proportionSignificant(permTestPValues, 0.1));
        updatePermutationResults("Proportion<0.010", proportionSignificant(permTestPValues, 0.01));
        updatePermutationResults("Proportion<0.001", proportionSignificant(permTestPValues, 0.001));
    }

    private void updatePermutationResults(String category, TestStatAndSite[] tsas) {
        for (int i = 0; i < tsas.length; i++) {
            //the nulls means the seach is only on category and trait
            String[] annotation = {category, tsas[i].trait, null, null};
            PermutationResults ppr = new PermutationResults(annotation, tsas[i].testStat);
            int index = permutationResultsVector.indexOf(ppr);
            PermutationResults pr = (PermutationResults) permutationResultsVector.get(index);
            pr.addPermutation(ppr.testStat);
        }
    }

    private void updatePermutationResults(String category, TestStatAndSite tsas) {
        TestStatAndSite[] ntsas = {tsas};
        updatePermutationResults(category, ntsas);
    }


    //Create a test statistic by combining all environments, and then returns the best site value for a given trait
    private TestStatAndSite[] siteSumOverEnvironmentByTrait(double[][][] test) {
        TestStatAndSite[] traitSums = new TestStatAndSite[test.length];
        double siteSum;
        for (int t = 0; t < test.length; t++) {
            traitSums[t] = new TestStatAndSite(-999, (String) traits.get(t), "", "");
            for (int s = 0; s < test[0][0].length; s++) {
                siteSum = 0;
                for (int e = 0; e < test[0].length; e++) {
                    if (!Double.isNaN(test[t][e][s])) {
                        siteSum += test[t][e][s];
                    }
                }//end of environments
                if (siteSum > traitSums[t].testStat) {//traitSums[t]=siteSum;
                    traitSums[t].testStat = siteSum;
                    traitSums[t].setLocusSite((String) locusSites.get(s));
                }
            }//end of sites
        }//end traits
        return traitSums;
    }

    //Create a test statistic by combining all traits and environments, and then returns the best site value
    private TestStatAndSite siteSumOverEnvironmentAndTrait(double[][][] test) {
        TestStatAndSite bestSiteSum = new TestStatAndSite(-999, "All", "", "");
        double siteSum;
        for (int s = 0; s < test[0][0].length; s++) {
            siteSum = 0;
            for (int e = 0; e < test[0].length; e++) {
                for (int t = 0; t < test.length; t++) {
                    if (!Double.isNaN(test[t][e][s])) {
                        siteSum += test[t][e][s];
                    }
                }//end traits
            }//end of environments
            if (siteSum > bestSiteSum.testStat) {//bestSiteSum=siteSum;
                bestSiteSum.testStat = siteSum;
                bestSiteSum.setLocusSite((String) locusSites.get(s));
            }
        }//end of sites
        return bestSiteSum;
    }

    //Returns the best site value for a given trait
    private TestStatAndSite[] siteOverEnvironmentByTrait(double[][][] test) {
        TestStatAndSite[] traitBest = new TestStatAndSite[test.length];
        for (int t = 0; t < test.length; t++) {
            traitBest[t] = new TestStatAndSite(-999, (String) traits.get(t), "", "");
            for (int s = 0; s < test[0][0].length; s++) {
                for (int e = 0; e < test[0].length; e++) {
                    if ((!Double.isNaN(test[t][e][s])) && (test[t][e][s] > traitBest[t].testStat)) {
                        traitBest[t].testStat = test[t][e][s];
                        traitBest[t].setLocusSite((String) locusSites.get(s));
                    }
                }//end of environments
            }//end of sites
        }//end traits
        return traitBest;
    }

    //returns the best site value
    private TestStatAndSite siteOverEnvironmentAndTrait(double[][][] test) {
        TestStatAndSite bestSite = new TestStatAndSite(-999, "All", "", "");
        for (int s = 0; s < test[0][0].length; s++) {
            for (int e = 0; e < test[0].length; e++) {
                for (int t = 0; t < test.length; t++) {
                    if ((!Double.isNaN(test[t][e][s])) && (test[t][e][s] > bestSite.testStat)) {//bestSite=test[t][e][s];
                        bestSite.testStat = test[t][e][s];
                        bestSite.setLocusSite((String) locusSites.get(s));
                    }
                }//end traits
            }//end of environments
        }//end of sites
        return bestSite;
    }

    //returns a test statistic that combines all the loci into one statistic
    private TestStatAndSite siteOverEnvironmentAndTraitAddAllLoci(double[][][] test) {
        TestStatAndSite bestSite = new TestStatAndSite(-999, "All", "", "");
        bestSite.testStat = 0;
        for (int s = 0; s < test[0][0].length; s++) {
            for (int e = 0; e < test[0].length; e++) {
                for (int t = 0; t < test.length; t++) {
                    //if(test[t][e][s]>bestSite.testStat)
                    //  {//bestSite=test[t][e][s];
                    if (!Double.isNaN(test[t][e][s])) {
                        bestSite.testStat = test[t][e][s] + bestSite.testStat;
                    }
                    bestSite.setLocusSite("All");
                //   }
                }//end traits
            }//end of environments
        }//end of sites
        return bestSite;
    }

    //returns a test statistic that is the proportion of sites more significant than the cutoff
    //this is still dealing if NaN wrong
    private TestStatAndSite proportionSignificant(double[][][] test, double cutoff) {
        TestStatAndSite bestSite = new TestStatAndSite(-999, "All", "", "");
        bestSite.testStat = 0;
        int count = 0;
        for (int s = 0; s < test[0][0].length; s++) {
            for (int e = 0; e < test[0].length; e++) {
                for (int t = 0; t < test.length; t++) {
                    if (!Double.isNaN(test[t][e][s])) {
                        count++;
                        if (test[t][e][s] < cutoff) {
                            bestSite.testStat += 1;
                        }
                    }
                    bestSite.setLocusSite("All");
                //   }
                }//end traits
            }//end of environments
        }//end of sites
        bestSite.testStat /= (double) count;
//    return Double.NaN;
        return bestSite;
    }

    private double[] locusSumOverEnvironmentByTrait(Vector psr) {
        return null;
    }

    //returns the best site value,
    private double locusSumOverEnvironmentAndTrait(Vector psr) {
        return -999;
    }

    private double[] locusOverEnvironmentByTrait(Vector psr) {
        return null;
    }

    //returns the best site value,
    private double locusOverEnvironmentAndTrait(Vector psr) {
        return -999;
    }

    //Implementation of TableReport Interface
    public Object[] getTableColumnNames() {
        String[] theLabels = {"Perm_Test", "Trait", "Locus", "Site", "Test_Stat", "P_value", "Permutations", "#Higher"};
        return theLabels;
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
        nf.setMaximumFractionDigits(8);
        int labelOffset, basicCols = 8;
        data = new String[basicCols];
        PermutationResults thePermResults;
        thePermResults = (PermutationResults) permutationResultsVector.get(row);
        labelOffset = 0;
        for (int j = 0; j < thePermResults.annotation.length; j++) {
            data[labelOffset++] = thePermResults.annotation[j];
        }
        data[labelOffset++] = "" + thePermResults.testStat;
        data[labelOffset++] = "" + thePermResults.permutedP;
        data[labelOffset++] = "" + thePermResults.permutations;
        data[labelOffset++] = "" + thePermResults.permutationsGreaterThenTestStat;
        return data;

    }

    public String getTableTitle() {
        //return "Permutations Of Single Factor ANOVA - No pop structure";
        return "Permutations of Logistic regression with pop structure";
    }

    public String toString() {
        if (permutationResultsVector.size() == 0) {
            return "Needs to be run";
        }
        StringBuffer cs = new StringBuffer();
        Object[] header = this.getTableColumnNames();
        for (int i = 0; i < header.length; i++) {
            cs.append(RightJustifiedFormat.monoString(header[i]));
        }
        cs.append("\n");
        Object[][] data = this.getTableData();
        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[i].length; j++) {
                cs.append(RightJustifiedFormat.monoString(data[i][j]));
            }
            cs.append("\n");
        }
        return cs.toString();
    }

    public int getRowCount() {
        return permutationResultsVector.size();
    }

    public int getElementCount() {
        throw new UnsupportedOperationException();
    }

    public int getColumnCount() {
        throw new UnsupportedOperationException();
    }
}

class TestStatAndSite {

    double testStat;
    String locus, site, trait;

    public TestStatAndSite(double ts, String t, String l, String s) {
        testStat = ts;
        trait = t;
        locus = l;
        site = s;
    }

    public TestStatAndSite(double ts, String t, String ls) {
        testStat = ts;
        trait = t;
        setLocusSite(ls);
    }

    void setLocusSite(String s) {
        int p = s.indexOf("#");
        if (p > 0) {
            locus = s.substring(0, p);
            site = s.substring(p + 1, s.length());
        } else {
            locus = site = s;
        }
    }
}

class PermutationResults implements Serializable, Cloneable, StatisticResult {

    String[] annotation;
    int index = -1;
    double testStat = 0;							// SST = total sum of squares
    int permutations = 0;
    int permutationsGreaterThenTestStat;
    double permutedP;

    public PermutationResults(String[] annotationM, double ts) {
        testStat = ts;
        annotation = new String[annotationM.length];
        for (int i = 0; i < annotationM.length; i++) {
            annotation[i] = annotationM[i];
        }
    }

    public Object clone() {
        try {
            PermutationResults pr = (PermutationResults) super.clone();
            pr.annotation = (String[]) this.annotation.clone();
            return pr;
        } catch (CloneNotSupportedException e) {
            throw new InternalError(e.toString());
        }
    }

    public void addPermutation(double pTS) {
        if (Double.isInfinite(pTS)) {
            return;  //this handles infinite likelihoods that get caused by missing data, resulting in variation in the dependent state
        }
        permutations++;
        if (pTS >= testStat) {
            permutationsGreaterThenTestStat++;
        }
        permutedP = (double) permutationsGreaterThenTestStat / (double) permutations;
    }

    public boolean equals(Object anObject) {
        PermutationResults x = (PermutationResults) anObject;
        for (int i = 0; i < this.annotation.length; i++) {
            if ((x.annotation[i] != null) && (this.annotation[i] != null) && (!x.annotation[i].equals(this.annotation[i]))) {
                return false;
            }
        }
        return true;
    }

    public void setIndex(int theIndex) {
        this.index = theIndex;
    }

    public String toString() {
        String result = "Permutation test\n";
        result += "Test Stat =" + testStat + "   \n";
        for (int i = 0; i < annotation.length; i++) {
            result += annotation[i] + "\t";
        }
        result += "\n";
        result += "Permutations=" + permutations + " P-value=" + permutedP + " \n";
        return result;
    }

    //StatisticResult interface methods
    public double testStatistic() {
        return testStat;
    }

    public double testPValue() {
        return permutedP;
    }

    public double testLikelihood() {
        return Math.abs(Math.log(permutedP));
    }

    public String[] getAttributes() {
        return annotation;
    }

    public String getAttribute(int a) {
        return annotation[a];
    }
}