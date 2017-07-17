/*
 * AlignmentFilterByGBSUtils
 */
package net.maizegenetics.gbs.pipeline;

import cern.colt.list.IntArrayList;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;
import java.util.TreeMap;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.AlignmentUtils;
import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.IdGroupUtils;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.pal.popgen.LinkageDisequilibrium;

import org.apache.poi.util.IntList;

/**
 * Filter methods for reference versus non-reference map. These are not really
 * appropriate for a general solution. GBSHapMapFilters are working towards the
 * more general filters.
 *
 * If there are useful methods they should be move to pal.
 * 
 * @author  Ed Buckler
 */
@Deprecated
public class AlignmentFilterByGBSUtils {

    public static final byte refAllele = (byte) 0x00; //A
    public static final byte altAllele = (byte) 0x11; //C
    public static final byte hetAllele = (byte) 0x01; //M

    private AlignmentFilterByGBSUtils() {
    }

    public static IdGroup getFilteredIdGroupByName(IdGroup orig, String[] filter, boolean keepTaxaMatchFilter) {
        ArrayList<Identifier> keepTaxa = new ArrayList<Identifier>();
        for (int i = 0; i < orig.getIdCount(); i++) {
            boolean matchFlag = false;
            for (String s : filter) {
                if (orig.getIdentifier(i).getFullName().contains(s)) {
                    matchFlag = true;
                }
            }
            if (matchFlag == keepTaxaMatchFilter) {
                keepTaxa.add(orig.getIdentifier(i));
            }
        }
        Identifier[] ids = new Identifier[keepTaxa.size()];
        keepTaxa.toArray(ids);
        IdGroup newGroup = new SimpleIdGroup(ids);
        return newGroup;
    }

    public static int[][] hetsByLine(Alignment a, boolean isRefAltCoded, boolean printToScreen) {
        int[][] counts = new int[2][a.getSequenceCount()];
        int totalScored = 0, totalHets = 0;
        for (int j = 0; j < a.getSequenceCount(); j++) {
            for (int i = 0; i < a.getSiteCount(); i++) {
                if (a.getBase(j, i) != Alignment.UNKNOWN_DIPLOID_ALLELE) {
                    counts[0][j]++;
                    if (isRefAltCoded) {
                        if (a.getBase(j, i) == hetAllele) {
                            counts[1][j]++;
                        }
                    } else {
                        if (AlignmentUtils.isHeterozygous(a.getBase(j, i))) {
                            counts[1][j]++;
                        }
                    }
                }
            }
            //            if(printToScreen) System.out.println(a.getIdGroup().getIdentifier(j).getName()+"\t"+
            //                    counts[0][j]+"\t"+counts[1][j]);
        }
        for (int c : counts[0]) {
            totalScored += c;
        }
        for (int c : counts[1]) {
            totalHets += c;
        }
        if (printToScreen) {
            System.out.println("Total Alleles:" + totalScored + " TotalHets:" + totalHets);
        }
        return counts;
    }

    public static IdGroup getLowHetIdGroup(Alignment a, boolean isRefAltCoded, double maxHets, int minCount) {
        int[][] hetCnt = hetsByLine(a, isRefAltCoded, false);
        boolean[] include = new boolean[a.getSequenceCount()];
        for (int i = 0; i < hetCnt[0].length; i++) {
            if (((double) hetCnt[1][i] / (double) hetCnt[0][i] > maxHets) || (hetCnt[0][i] < minCount)) {
                include[i] = false;
            } else {
                include[i] = true;
            }
        }
        return IdGroupUtils.idGroupSubset(a.getIdGroup(), include);
    }

    public static int[][] genotypicCountsBySite(Alignment a, boolean isRefAltCoded, boolean printToScreen) {
        int[][] counts = new int[5][a.getSiteCount()];  //total count [0], hets count Aa [1], major homozygous AA [2], minor homo aa [3], gap [4]
        if (printToScreen) {
            System.out.println("Locus\tMAF\tTaxaWKnown\tHetNum\tHetRate\tMajorAllele\tMajorCnt\tMinorAllele\tMinorCnt\tGapCnt");
        }
        for (int i = 0; i < a.getSiteCount(); i++) {
            byte majorAllele = a.getMajorAllele(i);
            majorAllele = (byte) (majorAllele << 4 | majorAllele);
            byte minorAllele = a.getMinorAllele(i);
            minorAllele = (byte) (minorAllele << 4 | minorAllele);
            for (int j = 0; j < a.getSequenceCount(); j++) {
                byte currentBase = a.getBase(j, i);
                if (currentBase != Alignment.UNKNOWN_DIPLOID_ALLELE) {
                    counts[0][i]++;
                    if (isRefAltCoded) {
                        if (AlignmentUtils.isEqual(currentBase, hetAllele)) {
                            counts[1][i]++;
                        }
                    } else {
                        if (AlignmentUtils.isHeterozygous(currentBase)) {
                            counts[1][i]++;
                        }
                    }
                    if (AlignmentUtils.isEqual(currentBase, majorAllele)) {
                        counts[2][i]++;
                    }
                    if (AlignmentUtils.isEqual(currentBase, minorAllele)) {
                        counts[3][i]++;
                    }
                    if (AlignmentUtils.isEqual(currentBase, NucleotideAlignmentConstants.GAP_DIPLOID_ALLELE)) {
                        counts[4][i]++;
                    }
                }
            }
            if (printToScreen) {
                System.out.println(a.getSNPID(i) + "\t" + a.getMinorAlleleFrequency(i) + "\t"
                        + counts[0][i] + "\t" + counts[1][i] + "\t" + (double) counts[1][i] / counts[0][i] + "\t"
                        + (char) majorAllele + "\t" + counts[2][i] + "\t" + (char) minorAllele + "\t" + counts[3][i] + "\t" + counts[4][i]);
            }
        }
        // System.out.println("counts"+Arrays.deepToString(counts));
        return counts;
    }

    public static int[] getLowHetSNPs(Alignment a, boolean isRefAltCoded, double minF, int minCount, double minMAF, double maxMAF) {
        return getLowHetSNPs(a, isRefAltCoded, minF, minCount, minMAF, maxMAF, null, null);
    }

    public static int[] getLowHetSNPs(Alignment a, boolean isRefAltCoded, double minF, int minCount, double minMAF, double maxMAF, SNPLogging snpLogging, String logTest) {
        String REMOVED_STATUS = "Removed";
        int[][] hetCnt = genotypicCountsBySite(a, isRefAltCoded, false);
        IntArrayList goodSites = new IntArrayList();
        //        System.out.println("Site Genotypes Hets propHets theMAF expHets relHets obsF");
        for (int i = 0; i < hetCnt[0].length; i++) {
            double propHets = hetCnt[0][i] > 0 ? (double) hetCnt[1][i] / (double) hetCnt[0][i] : 0.0;
            double theMAF = hetCnt[0][i] > 0 ? (double) (hetCnt[3][i] + ((double) hetCnt[1][i] / 2.0)) / (double) hetCnt[0][i] : 0.0;
            double expHets = 2.0 * theMAF * (1 - theMAF);
            double obsF = expHets > 0 ? 1.0 - (propHets / expHets) : 1.0;  // if all missing data or the SNP is monomorphic, don't reject merely on the basis of minF
            //            double relHets=(double)hetCnt[1][i]/expHets;
            //            System.out.printf("%d %d %d %g %g %g %g %g %n",i, hetCnt[0][i], hetCnt[1][i], propHets,theMAF, expHets, relHets, obsF);
            if ((obsF > minF) && (hetCnt[0][i] >= minCount) && (theMAF >= minMAF) && (theMAF <= maxMAF)) {
                goodSites.add(i);
            } else if (snpLogging != null) {
                String cutoff = "isRefAltCoded: " + isRefAltCoded + " minF: " + minF + " minCount: " + minCount + " minMAF: " + minMAF + " maxMAF: " + maxMAF;
                String value = "obsF: " + obsF + " obsMinCount: " + hetCnt[0][i] + " MAF: " + theMAF;
                snpLogging.writeEntry(a, i, null, null, logTest, REMOVED_STATUS, value, cutoff);
            }
        }
        goodSites.trimToSize();
        return goodSites.elements();
    }

    public static void getCoverage_MAF_F_Dist(Alignment a, boolean isRefAltCoded) {
        if (a.getSiteCount() == 0) {
            return;
        }
        int[][] hetCnt = genotypicCountsBySite(a, isRefAltCoded, false);
        double[] coverageD = new double[a.getSiteCount()];
        double[] mafD = new double[a.getSiteCount()];
        double[] fD = new double[a.getSiteCount()];
        double[] gapD = new double[a.getSiteCount()];
        //       System.out.println("Site Genotypes Hets propHets theMAF expHets relHets obsF");
        for (int i = 0; i < hetCnt[0].length; i++) {
            double propHets = (double) hetCnt[1][i] / (double) hetCnt[0][i];
            double theMAF = mafD[i] = (double) (hetCnt[3][i] + ((double) hetCnt[1][i] / 2.0)) / (double) hetCnt[0][i];
            //        if(theMAF<0.0001) System.out.printf("%d %s %d %d %d %g %n",i,a.getSNPID(i),hetCnt[3][i],hetCnt[1][i],hetCnt[0][i], theMAF);
            double expHets = 2.0 * theMAF * (1 - theMAF);
            fD[i] = 1.0 - (propHets / expHets);
            coverageD[i] = (double) hetCnt[0][i] / (double) a.getSequenceCount();
            gapD[i] = (double) hetCnt[4][i] / (double) a.getSequenceCount();
        }
        Arrays.sort(mafD);
        Arrays.sort(fD);
        Arrays.sort(coverageD);
        Arrays.sort(gapD);
        System.out.println("Percentile\tCoverage\tMAF\tF(inbreeding)\tGapProp");
        if (coverageD.length >= 20) {
            for (int i = 0; i < coverageD.length; i += (coverageD.length / 20)) {
                System.out.printf("%.2g\t%.3g\t%.3g\t%.3g\t%.3g%n", ((double) i / (double) coverageD.length), coverageD[i], mafD[i], fD[i], gapD[i]);
            }
        } else {
            System.out.printf("%.2g\t%.3g\t%.3g\t%.3g\t%.3g%n", 1.0, coverageD[coverageD.length - 1], mafD[coverageD.length - 1], fD[coverageD.length - 1], gapD[coverageD.length - 1]);
        }

    }

    public static double getErrorRateForDuplicatedTaxa(Alignment a, boolean ignoreHets, boolean random, boolean printToScreen) {
        IdGroup idg = a.getIdGroup();
        TreeMap<String, Integer> sortedIds = new TreeMap<String, Integer>();
        for (int i = 0; i < idg.getIdCount(); i++) {
            sortedIds.put(idg.getIdentifier(i).getFullName().toUpperCase(), i);
        }
        long cntDiff = 0, cntTotal = 0;
        Map.Entry<String, Integer> priorEntry = sortedIds.lastEntry();
        if (printToScreen) {
            System.out.println("Entry1 Entry2 SNPDiff SNPsCompared PropDiff");
        }
        for (Map.Entry<String, Integer> entry : sortedIds.entrySet()) {
            if (priorEntry.getKey().split(":")[0].equals(entry.getKey().split(":")[0])) {
                int cntDiffTaxa = 0, cntTotalTaxa = 0;
                for (int i = 0; i < a.getSiteCount(); i++) {
                    byte pB = Alignment.UNKNOWN_DIPLOID_ALLELE;
                    if (random) {
                        int t = (int) Math.round(Math.random() * (a.getSequenceCount() - 1));
                        pB = a.getBase(t, i);
                    } else {
                        pB = a.getBase(priorEntry.getValue(), i);
                    }
                    byte cB = a.getBase(entry.getValue(), i);
                    if (ignoreHets && AlignmentUtils.isHeterozygous(pB)) {
                        continue;
                    }
                    if (ignoreHets && AlignmentUtils.isHeterozygous(cB)) {
                        continue;
                    }
                    if ((pB != Alignment.UNKNOWN_DIPLOID_ALLELE) && (cB != Alignment.UNKNOWN_DIPLOID_ALLELE)) {
                        cntTotalTaxa++;
                        if (!AlignmentUtils.isEqual(pB, cB)) {
                            cntDiffTaxa++;
                        }
                    }
                }
                if (printToScreen) {
                    System.out.printf("%s %s %d %d %g %n", priorEntry.getKey(), entry.getKey(), cntDiffTaxa, cntTotalTaxa, (double) cntDiffTaxa / (double) cntTotalTaxa);
                }
                cntDiff += cntDiffTaxa;
                cntTotal += cntTotalTaxa;
            }
            priorEntry = entry;
        }
        double totalER = (double) cntDiff / (double) cntTotal;
        if (random) {
            System.out.print("RandomDups ");
        } else {
            System.out.print("RealDups ");
        }
        if (ignoreHets) {
            System.out.print("ignoreHets ");
        } else {
            System.out.print("testHets ");
        }
        System.out.printf("ALL %d %d %g %n", cntDiff, cntTotal, (double) cntDiff / (double) cntTotal);
        return totalER;
    }

    public static int[] getGoodSitesByLD(Alignment a, double minR2, double minBonferroniP, int minPosDist, int windowSize,
            int minCnt, boolean keepUnproven) {
        IntArrayList goodSites = new IntArrayList();
        //LinkageDisequilibrium theLD = new LinkageDisequilibrium(a, true, 100,
        //        minCnt, windowSize, LinkageDisequilibrium.testDesign.SlidingWindow);

        LinkageDisequilibrium theLD = new LinkageDisequilibrium(a, windowSize,
                LinkageDisequilibrium.testDesign.SlidingWindow, -1, null, false, -1, null, LinkageDisequilibrium.HetTreatment.Homozygous);
        theLD.run();
        for (int i = 0; i < a.getSiteCount(); i++) {
            int cntInformative = 0;
            double obsMaxR2 = -1;
            double obsMinP = 1;
            for (int j = i - (windowSize / 2); j < i + (windowSize / 2); j++) {
                double rValue = theLD.getRSqr(i, j);
                if (Double.isNaN(rValue)) {
                    continue;
                }
                if (Math.abs(a.getPositionInLocus(i) - a.getPositionInLocus(j)) < minPosDist) {
                    continue;
                }
                cntInformative++;
                if (rValue > obsMaxR2) {
                    obsMaxR2 = rValue;
                }
                if (rValue < obsMinP) {
                    obsMinP = rValue;
                }

            }
            double bonP = (cntInformative > 0) ? obsMinP * cntInformative : 1;
            if (i % 1000 == 0) {
                System.out.printf(a.getSNPID(i) + " %d %g %g %g %n", cntInformative, obsMaxR2, obsMinP, bonP);
            }
            //           if(obsMaxR2<minR2) System.out.println(i+" "+cntInformative+" "+obsMaxR2+" "+obsMinP);

            if (((obsMaxR2 >= minR2) && (bonP <= minBonferroniP)) || (keepUnproven && (cntInformative == 0))) {
                goodSites.add(i);
            }
        }
        goodSites.trimToSize();
        return goodSites.elements();
    }

    public static int[][] countCrossoversByLine(Alignment a) {
        int[][] cos = new int[2][a.getSequenceCount()];  //total count of useful markers in
        //row 0, crossovers in row 1
        int sumGood = 0, sumCO = 0;
        for (int t = 0; t < a.getSequenceCount(); t++) {
            byte lastHomozygous = Alignment.UNKNOWN_DIPLOID_ALLELE;
            for (int i = 0; i < a.getSiteCount(); i++) {
                byte currBase = a.getBase(t, i);
                if ((currBase != refAllele) && (currBase != altAllele)) {
                    continue;  //not useful
                }
                cos[0][t]++;
                if (currBase != lastHomozygous) {
                    cos[1][t]++;
                }
                lastHomozygous = currBase;
            }
            sumCO += cos[1][t];
            sumGood += cos[0][t];
            //  System.out.println(a.getIdGroup().getIdentifier(t).getName()+" "+cos[0][t]+" "+cos[1][t]);
        }
        System.out.println("TotalHomoMarkers:" + sumGood + " TotalCrossover:" + sumCO);
        //        System.out.println("By Line:"+Arrays.toString(cos[0]));
        //        System.out.println("By Line:"+Arrays.toString(cos[1]));
        return cos;
    }

    public static int[][] countDCO(Alignment a, boolean byTaxa) {
        int[][] cos;
        if (byTaxa) {
            cos = new int[2][a.getSequenceCount()];
        } else {
            cos = new int[2][a.getSiteCount()];
        }
        int sumGood = 0, sumDCO = 0;
        for (int t = 0; t < a.getSequenceCount(); t++) {
            byte[] base = {-1, -1, -1};
            int[] site = {-1, -1, -1};
            for (int s = 0; s < a.getSiteCount(); s++) {
                if ((a.getBase(t, s) == refAllele) || (a.getBase(t, s) == altAllele)) {
                    base[0] = base[1];
                    base[1] = base[2];
                    base[2] = a.getBase(t, s);
                    site[0] = site[1];
                    site[1] = site[2];
                    site[2] = s;
                    if (base[0] == -1) {
                        continue;  //need to load the arrays first before checking
                    }
                    sumGood++;
                    if (byTaxa) {
                        cos[0][t]++;
                    } else {
                        cos[0][site[1]]++;
                    }
                    if ((base[0] == base[2]) && (base[0] != base[1])) {
                        sumDCO++;
                        if (byTaxa) {
                            cos[1][t]++;
                        } else {
                            cos[1][site[1]]++;
                        }
                    }
                }
            }
            //           System.out.println(s+" "+cos[0][s]+" "+cos[1][s]);
        }
        if (byTaxa) {
            System.out.print("ByTaxa:");
        } else {
            System.out.print("BySite:");
        }
        System.out.println("TotalHomoMarkers:" + sumGood + " TotalDoubleCrossover:" + sumDCO);
        //  System.out.println("By Line:"+Arrays.toString(cos[1]));
        return cos;
    }

    public static int[] getLowDCOSNPs(Alignment a, double maxDCOrate, int minCount) {
        int[][] dcoCnt = countDCO(a, false);
        IntList goodSites = new IntList();
        for (int i = 0; i < dcoCnt[0].length; i++) {
            if (((double) dcoCnt[1][i] / (double) dcoCnt[0][i] < maxDCOrate)
                    && (dcoCnt[0][i] >= minCount)) {
                goodSites.add(i);
            }
        }
        return goodSites.toArray();
    }

    public static IdGroup getLowDCOIdGroup(Alignment a, boolean isRefAltCoded, double maxDCO, int minCount) {
        int[][] dcoCnt = hetsByLine(a, isRefAltCoded, false);
        boolean[] include = new boolean[a.getSequenceCount()];
        for (int i = 0; i < dcoCnt[0].length; i++) {
            if (((double) dcoCnt[1][i] / (double) dcoCnt[0][i] > maxDCO) || (dcoCnt[0][i] < minCount)) {
                include[i] = false;
            } else {
                include[i] = true;
            }
        }
        return IdGroupUtils.idGroupSubset(a.getIdGroup(), include);
    }
}
