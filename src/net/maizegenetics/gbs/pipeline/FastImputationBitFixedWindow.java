/*
 * FastImputationBitFixedWindow
 */
package net.maizegenetics.gbs.pipeline;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.TreeMap;
import net.maizegenetics.baseplugins.ConvertSBitTBitPlugin;

import net.maizegenetics.util.OpenBitSet;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.AlignmentUtils;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.alignment.MutableNucleotideAlignment;
import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;
import net.maizegenetics.pal.statistics.ChiSquareTest;

/**
 * An extremely fast imputation approach that uses KNN, but calculates the matrix
 * with a 64 bit window, and it does all operations with bits.  This has a fixed window size
 * unlike the other BDI imputation approaches.
 *
 * All heterozygous sets are to missing.
 *
 * Definitely a work in progress.
 * Perhaps we should do ratio to major/minor
 * Perhaps we should scan for massive nearly identical regions
 * Need to decide whether P or Length/Identity approach is better
 *
 * @deprecated Replaced by better methods {@link MinorWindowViterbiImputationPlugin}
 * @author ed
 */
@Deprecated  
public class FastImputationBitFixedWindow {

    short[][] matchInWin, diffInWin;
    float[] presentProp;  //Present proportion
    int[] presentCntForSites;
    int[] presentCntForTaxa;
    int[] hetCntForTaxa;
    boolean[] highHet = null;
    long totalPresent = 0;
    Alignment anchorAlignment = null;
    MutableNucleotideAlignment impAlign = null;
    static int windowSize = 64 * 64;  //left and right distance of 64bp window, so window 2 is 2+1+2=5 or 320bp
    static int minLengthOfMatch = 50;
    static double minIdentity = 0.95;  //std 0.99, 0.95 gave some pretty good results also
    static int minCountNNperLine = 4;  //std 4
    static int minCountNNperSite = 2;  //std 2
    static int segments = 1;
    static double majorityRule = 0.76;
    static double maxLineErrorRate = 0.025; //std 0.01
    static boolean imputeGaps = false;
    static boolean maskHighHets = true;
    static double maxHetStatic = 0.01; //HetStat = hetSiteCnt/covSiteCnt :  empirically derived

    // static int minimumMajorityRule=4;  //NN within  is used, and tie are not imputed
    public FastImputationBitFixedWindow(Alignment a, boolean[] highHet) {
        this.highHet = highHet;
        this.anchorAlignment = ConvertSBitTBitPlugin.convertAlignment(a, ConvertSBitTBitPlugin.CONVERT_TYPE.tbit, null);
        double avgMissing = calcPropMissingByLine();
        System.out.println("Average Missing:" + avgMissing);
        avgMissing = calcPropMissing();
        System.out.println("Average Missing (site and taxa):" + avgMissing + " totalPresent" + totalPresent);
        System.out.println("Average MAF:" + avgMinorAlleleFrequency());
        reportTaxaStats();
        AlignmentFilterByGBSUtils.getCoverage_MAF_F_Dist(anchorAlignment, false);
        double realDist = AlignmentFilterByGBSUtils.getErrorRateForDuplicatedTaxa(anchorAlignment, true, false, true);
        double randomDist = AlignmentFilterByGBSUtils.getErrorRateForDuplicatedTaxa(anchorAlignment, true, true, false);
        System.out.println("Ratio of RandomToReal:" + randomDist / realDist);
        System.out.println("Creating mutable alignment");
        impAlign = MutableNucleotideAlignment.getInstance(a);
        //impAlign = new MutableSimpleAlignment(a);
        for (int i = 4096; i >= 1024; i /= 2) {
            //            segments=i;  
            System.out.println("Starting imputation");
            //           windowSize=((a.getSiteCount()/64)/segments)*64;
            windowSize = i;
            int offset = 0;
            System.out.println("Window size:" + windowSize + " offset:" + offset);
            imputeBySiteJump(windowSize, 0, minLengthOfMatch, minIdentity);
            //            realDist=AlignmentFilterByGBSUtils.getErrorRateForDuplicatedTaxa(impAlign, true,false,false);
            //            randomDist=AlignmentFilterByGBSUtils.getErrorRateForDuplicatedTaxa(impAlign, true,true,false);
            //            System.out.println("Ratio of RandomToReal:"+randomDist/realDist);
            System.out.println("Window size:" + windowSize + " offset:" + offset);
            AlignmentFilterByGBSUtils.getCoverage_MAF_F_Dist(impAlign, false);
            offset = windowSize / 2;
            System.out.println("Window size:" + windowSize + " offset:" + offset);
            imputeBySiteJump(windowSize, offset, minLengthOfMatch, minIdentity);
            System.out.println("Window size:" + windowSize + " offset:" + offset);
            AlignmentFilterByGBSUtils.getCoverage_MAF_F_Dist(impAlign, false);
            System.out.println("Xratio: " + countXRateInImputation());
        }
    }

    public FastImputationBitFixedWindow(Alignment a) {
        this(a, null);
    }

    private double avgMinorAlleleFrequency() {
        long totmj = 0, totmn = 0;
        for (int i = 0; i < anchorAlignment.getSequenceCount(); i++) {
            totmj += anchorAlignment.getAllelePresenceForAllSites(i, 0).cardinality();
            totmn += anchorAlignment.getAllelePresenceForAllSites(i, 1).cardinality();
            //totmj += anchorAlignment.getTaxaBitsNoClone(i, 0).cardinality();
            //totmn += anchorAlignment.getTaxaBitsNoClone(i, 1).cardinality();
        }
        double theMAF = (double) totmn / (double) (totmj + totmn);
        return theMAF;
    }

    private double countXRateInImputation() {
        long xCnt = 0, knownCnt = 0;
        for (int i = 0; i < impAlign.getSequenceCount(); i++) {
            for (int j = 0; j < impAlign.getSiteCount(); j++) {
                if (impAlign.getBase(i, j) != Alignment.UNKNOWN_DIPLOID_ALLELE) {
                    knownCnt++;
                    if (impAlign.getBase(i, j) == (byte) 'X') {
                        xCnt++;
                    }
                }

            }
        }
        double theXRate = (double) xCnt / (double) knownCnt;
        return theXRate;
    }

    public void reportTaxaStats() {
        double sites = anchorAlignment.getSiteCount();
        for (int t = 0; t < anchorAlignment.getSequenceCount(); t++) {
            System.out.printf("%s %d %.5g %d %.5g %s %n", anchorAlignment.getIdGroup().getIdentifier(t).getFullName(), presentCntForTaxa[t],
                    (double) presentCntForTaxa[t] / sites, hetCntForTaxa[t], (double) hetCntForTaxa[t] / sites, highHet[t]);
        }
    }

    private double calcPropMissing() {
        boolean calculateHighHet = false;
        presentCntForTaxa = new int[anchorAlignment.getSequenceCount()];
        presentCntForSites = new int[anchorAlignment.getSiteCount()];
        hetCntForTaxa = new int[anchorAlignment.getSequenceCount()];
        if (highHet == null) {
            highHet = new boolean[anchorAlignment.getSequenceCount()];
            calculateHighHet = true;
        }
        totalPresent = 0;
        double sd = anchorAlignment.getSiteCount();
        for (int t = 0; t < anchorAlignment.getSequenceCount(); t++) {
            for (int s = 0; s < anchorAlignment.getSiteCount(); s++) {
                byte cb = anchorAlignment.getBase(t, s);
                if ((cb != Alignment.UNKNOWN_DIPLOID_ALLELE) && (cb != NucleotideAlignmentConstants.GAP_DIPLOID_ALLELE)) {
                    presentCntForTaxa[t]++;
                    presentCntForSites[s]++;
                    if (AlignmentUtils.isHeterozygous(cb)) {
                        hetCntForTaxa[t]++;
                    }
                    totalPresent++;
                }
            }
            //            double covProp=(double)presentCntForTaxa[t]/sd;
            double hetProp = (double) hetCntForTaxa[t] / presentCntForTaxa[t];
            if (calculateHighHet) {
                highHet[t] = (hetProp > maxHetStatic) ? true : false;
            }

        }
        return (double) totalPresent / ((double) anchorAlignment.getSequenceCount() * (double) anchorAlignment.getSiteCount());
    }

    private double calcPropMissingByLine() {
        presentProp = new float[anchorAlignment.getSequenceCount()];
        double avgMissing = 0;
        for (int i = 0; i < anchorAlignment.getSequenceCount(); i++) {
            long present = OpenBitSet.unionCount(anchorAlignment.getAllelePresenceForAllSites(i, 0), anchorAlignment.getAllelePresenceForAllSites(i, 1));
            //long present = OpenBitSet.unionCount(anchorAlignment.getTaxaBitsNoClone(i, 0), anchorAlignment.getTaxaBitsNoClone(i, 1));
            presentProp[i] = (float) present / (float) anchorAlignment.getSiteCount();
            avgMissing += presentProp[i];
        }
        return avgMissing / (double) anchorAlignment.getSequenceCount();
    }

    private void imputeBySiteJump(int window, int offset, int minLength, double minIdentity) {
        long time = System.currentTimeMillis();
        int knownSNPs = 0, unknownSNPs = 0, imputedSNPs = 0;
        int numSeqs = anchorAlignment.getSequenceCount();
        System.out.println("Initial matrix created in " + (System.currentTimeMillis() - time));
        time = System.currentTimeMillis();
        int corrCnt = 0, wrongCnt = 0, gapCnt = 0;
        int imputableLines = 0, imputedWithHighError = 0;
        TreeMap<Double, Integer> lenBestLine = new TreeMap<Double, Integer>(Collections.reverseOrder());
        for (int b = 0 + offset; b < anchorAlignment.getSiteCount() - window; b += window) {
            //       currWord=b>>6;
            int endBase = b + window;
            initHapLengths(b, endBase);
            double rate = (double) unknownSNPs / (double) (System.currentTimeMillis() - time);
            System.out.println("Imputed base:" + b + " known:" + knownSNPs + " unknownSNPs:" + unknownSNPs
                    + " imputed:" + imputedSNPs + " Rate:" + rate);
            for (int i = 0; i < numSeqs; i++) {
                if (maskHighHets && highHet[i]) {
                    continue;
                }
                lenBestLine.clear();
                for (int j = 0; j < numSeqs; j++) {
                    if (i == j) {
                        continue;
                    }
                    if (maskHighHets && highHet[j]) {
                        continue;
                    }
                    int sum = matchInWin[i][j] + diffInWin[i][j];
                    double identity = (double) matchInWin[i][j] / (double) sum;
                    if ((sum > minLength) && (identity > minIdentity)) {
                        lenBestLine.put(identity, j);
                    }
                }
                if (lenBestLine.size() < minCountNNperLine) {
                    continue;
                }
                imputableLines++;
                if (i == 0) {
                    System.out.println(" cnt" + lenBestLine.size());
                }
                //                byte[] calls=consensusCalls(anchorAlignment, i, b, endBase, lenBestLine, false, majorityRule, minCountNNperSite);
                //                System.out.println(Arrays.toString(calls));
                byte[] calls = consensusCallBit(anchorAlignment, i, b, endBase, lenBestLine,
                        false, majorityRule, minCountNNperSite, true, imputeGaps);  //determine error while ignoring known
                int[] callError = compareCallsWithActual(anchorAlignment, b, i, calls);
                double lineError = (double) callError[1] / (double) (callError[1] + callError[0]);
                if (lineError > maxLineErrorRate) {
                    imputedWithHighError++;
                    //                    System.out.println("high error line:"+anchorAlignment.getTaxaName(i)+" Error:"+lineError);
                    continue;
                }
                calls = consensusCallBit(anchorAlignment, i, b, endBase, lenBestLine,
                        false, majorityRule, minCountNNperSite, false, imputeGaps);  //recall with known; currently conflict are set to unknown
                setBaseInImputedAlignment(impAlign, b, i, calls);
                //              if(b>3000) System.out.println(Arrays.toString(callError));
                corrCnt += callError[0];
                wrongCnt += callError[1];
                gapCnt += callError[3];
                double errorRate = (double) wrongCnt / (double) (corrCnt + wrongCnt);

                if (i % 100 == 0) {
                    System.out.printf("%d-%d %d R: %d W: %d G: %d ErrorRate: %g Imputable:%d ImputedHighError:%d %n", b, endBase, i,
                            corrCnt, wrongCnt, gapCnt, errorRate, imputableLines, imputedWithHighError);
                }
            }
            double errorRate = (double) wrongCnt / (double) (corrCnt + wrongCnt);
            System.out.printf("R: %d W: %d ErrorRate: %g %n", corrCnt, wrongCnt, errorRate);
        }
    }

    private int[] compareCallsWithActual(Alignment a, int startBase, int taxon, byte[] calls) {
        int[] result = new int[4];  //agree[0], disagree[1], noncomparison[2], gaps[3]
        byte alignB;
        for (int s = 0; s < calls.length; s++) {
            alignB = a.getBase(taxon, s + startBase);
            if (calls[s] == NucleotideAlignmentConstants.GAP_DIPLOID_ALLELE) {
                result[3]++;
            }
            if (calls[s] == Alignment.UNKNOWN_DIPLOID_ALLELE) {
                result[2]++;
            } else if (alignB == Alignment.UNKNOWN_DIPLOID_ALLELE) {
                result[2]++;
            } else if (alignB == calls[s]) {
                result[0]++;
            } else {
                result[1]++;
            }
        }
        return result;
    }

    private void setBaseInImputedAlignment(MutableNucleotideAlignment a, int startBase, int taxon, byte[] calls) {
        for (int s = 0; s < calls.length; s++) {
            if (calls[s] != Alignment.UNKNOWN_DIPLOID_ALLELE) {
                byte cb = a.getBase(taxon, s + startBase);
                if (cb == Alignment.UNKNOWN_DIPLOID_ALLELE) {
                    a.setBase(taxon, s + startBase, calls[s]);
                } else if (cb != calls[s]) {
                    a.setBase(taxon, s + startBase, Alignment.UNKNOWN_DIPLOID_ALLELE);
                }
            }
        }
    }

    private byte[] consensusCalls(Alignment a, int taxon, int startBase, int endBase, TreeMap<Double, Integer> taxa,
            boolean callhets, double majority) {
        short[][] siteCnt = new short[2][a.getSiteCount()];
        int[] taxaIndex = new int[taxa.size()];
        ArrayList<Integer> taxaList = new ArrayList(taxa.values());
        for (int t = 0; t < taxaIndex.length; t++) {
            taxaIndex[t] = taxaList.get(t);
        }
        byte[] calls = new byte[endBase - startBase];
        Arrays.fill(calls, Alignment.UNKNOWN_DIPLOID_ALLELE);
        for (int s = startBase; s < endBase; s++) {
            calls[s - startBase] = a.getBase(taxon, s);
            byte mj = a.getMajorAllele(s);
            byte mn = a.getMinorAllele(s);
            byte het = AlignmentUtils.getDiploidValue(mj, mn);
            //byte[] snpValue = {mj, mn};
            //byte het = IUPACNucleotides.getDegerateSNPByteFromTwoSNPs(snpValue);
            for (int t = 0; t < taxaIndex.length; t++) {
                byte ob = a.getBase(taxaIndex[t], s);
                if (ob == Alignment.UNKNOWN_DIPLOID_ALLELE) {
                    continue;
                }
                if (ob == mj) {
                    siteCnt[0][s]++;
                } else if (ob == mn) {
                    siteCnt[1][s]++;
                } else if (ob == het) {
                    siteCnt[0][s]++;
                    siteCnt[1][s]++;
                }
            }
            int totalCnt = siteCnt[0][s] + siteCnt[1][s];
            if (totalCnt == 0) {
                continue;  //no data leave missing
            }
            if ((double) siteCnt[0][s] / (double) totalCnt > majority) {
                calls[s - startBase] = mj;
            } else if ((double) siteCnt[1][s] / (double) totalCnt > majority) {
                calls[s - startBase] = mn;
            } else if (callhets) {
                calls[s - startBase] = het;
            }
        }
        //       System.out.println("Byt:"+Arrays.toString(siteCnt[0]));
        return calls;
    }

    private byte[] consensusCallBit(Alignment a, int taxon, int startBase, int endBase, TreeMap<Double, Integer> taxa,
            boolean callhets, double majority, int minCount, boolean ignoreKnownBases, boolean imputeGaps) {
        int[] taxaIndex = new int[taxa.size()];
        ArrayList<Integer> taxaList = new ArrayList(taxa.values());
        for (int t = 0; t < taxaIndex.length; t++) {
            taxaIndex[t] = taxaList.get(t);
        }
        short[][] siteCnt = new short[2][endBase - startBase];
        double[] sumExpPresent = new double[endBase - startBase];
        int[] sumNNxSitePresent = new int[endBase - startBase];
        int sumNNPresent = 0;
        for (int t = 0; t < taxaIndex.length; t++) {
            sumNNPresent += presentCntForTaxa[taxaIndex[t]];
        }
        byte[] calls = new byte[endBase - startBase];
        Arrays.fill(calls, Alignment.UNKNOWN_DIPLOID_ALLELE);
        for (int alignS = startBase; alignS < endBase; alignS += 64) {
            int currWord = alignS / 64;
            int callSite = alignS - startBase;
            for (int t = 0; t < taxaIndex.length; t++) {
                long bmj = anchorAlignment.getAllelePresenceForAllSites(taxaIndex[t], 0).getBits()[currWord];
                //long bmj = anchorAlignment.getTaxaBitsNoClone(taxaIndex[t], 0).getBits()[currWord];
                long bmn = anchorAlignment.getAllelePresenceForAllSites(taxaIndex[t], 1).getBits()[currWord];
                //long bmn = anchorAlignment.getTaxaBitsNoClone(taxaIndex[t], 1).getBits()[currWord];
                int cs = callSite;
                for (int j = 0; j < 64; j++) {
                    boolean presentFlag = false;
                    if ((bmj & 0x01) != 0) {
                        siteCnt[0][cs]++;
                        presentFlag = true;
                    }
                    bmj = bmj >> 1;
                    if ((bmn & 0x01) != 0) {
                        siteCnt[1][cs]++;
                        presentFlag = true;
                    }
                    bmn = bmn >> 1;
                    sumExpPresent[cs] += presentProp[taxaIndex[t]];
                    if (presentFlag) {
                        sumNNxSitePresent[cs]++;
                    }
                    cs++;
                }
            }
        }
        //       System.out.println("Bit:"+Arrays.toString(siteCnt[0]));
        for (int alignS = startBase; alignS < endBase; alignS++) {
            int callS = alignS - startBase;
            byte ob = a.getBase(taxon, alignS);
            if (ignoreKnownBases) {
                ob = Alignment.UNKNOWN_DIPLOID_ALLELE;
            }
            calls[callS] = ob;
            byte mj = a.getMajorAllele(alignS);
            byte mn = a.getMinorAllele(alignS);
            int totalCnt = siteCnt[0][callS] + siteCnt[1][callS];
//            double expPres=sumExpPresent[callS]/(double)taxaIndex.length;

//            System.out.println(Arrays.toString(exp)+Arrays.toString(obs));
            if (imputeGaps && ((double) totalCnt < ((double) sumExpPresent[callS] / 2.0))) {
//                double[] exp= {expPres,1-expPres};
//                int[] obs={totalCnt,taxaIndex.length-totalCnt};
                int nonNNPresentCnt = (int) totalPresent - sumNNPresent;  //marginals
                int nonSitePresentCnt = (int) totalPresent - presentCntForSites[callS];  //marginals
                int[] obs2 = {sumNNxSitePresent[callS], presentCntForSites[callS] - sumNNxSitePresent[callS],
                    sumNNPresent - sumNNxSitePresent[callS], (int) totalPresent - presentCntForSites[callS] - sumNNPresent + sumNNxSitePresent[callS]};
                double totalPresentSqr = (double) totalPresent * (double) totalPresent;
                double[] exp2 = {(double) sumNNPresent * (double) presentCntForSites[callS] / totalPresentSqr, (double) nonNNPresentCnt * (double) presentCntForSites[callS] / totalPresentSqr,
                    (double) sumNNPresent * (double) nonSitePresentCnt / totalPresentSqr, (double) nonNNPresentCnt * (double) nonSitePresentCnt / totalPresentSqr};
//                double presProbX2=ChiSquareTest.compare(exp, obs);
                double presProbX2v2 = ChiSquareTest.compare(exp2, obs2);
//                System.out.printf("%d %d %d %d %g %g %n",taxaIndex.length, siteCnt[0][callS],siteCnt[1][callS],totalCnt,
//                     sumExpPresent[callS], presProbX2);
//                System.out.println(Arrays.toString(obs2)+Arrays.toString(exp2)+" P:"+presProbX2v2);
//                System.out.printf("%d %d %d %d %g %g %n",taxaIndex.length, siteCnt[0][callS],siteCnt[1][callS],totalCnt,
//                     sumExpPresent[callS], presProbX2);
                //       if((presProbX2>=0)&&(presProbX2<0.01)) {calls[callS]=NucleotideAlignmentConstants.GAP_DIPLOID_ALLELE; continue;}
                if ((presProbX2v2 >= 0) && (presProbX2v2 < 0.01)) {
                    calls[callS] = NucleotideAlignmentConstants.GAP_DIPLOID_ALLELE;
                    continue;
                }
            }

            if (totalCnt < minCount) {
                continue;  //no data leave missing
            }
            if ((double) siteCnt[0][callS] / (double) totalCnt > majority) {
                if ((ob != Alignment.UNKNOWN_DIPLOID_ALLELE) && (ob != mj)) {
                    calls[callS] = Alignment.UNKNOWN_DIPLOID_ALLELE;
                } else {
                   // calls[callS] = mj;
                    calls[callS] = AlignmentUtils.getDiploidValue(mj, mj);
                }
            } else if ((double) siteCnt[1][callS] / (double) totalCnt > majority) {
                if ((ob != Alignment.UNKNOWN_DIPLOID_ALLELE) && (ob != mn)) {
                    calls[callS] = Alignment.UNKNOWN_DIPLOID_ALLELE;
                } else {
                   // calls[callS] = mn;
                    calls[callS] = AlignmentUtils.getDiploidValue(mn, mn);
                }
            } else if (callhets) {
                //byte[] snpValue = {mj, mn};
                byte het = AlignmentUtils.getDiploidValue(mj, mn);
                //byte het = IUPACNucleotides.getDegerateSNPByteFromTwoSNPs(snpValue);
                calls[callS] = het;
            }
        }
        return calls;
    }

    private void initHapLengths(int startSite, int endSite) {
        matchInWin = new short[anchorAlignment.getSequenceCount()][anchorAlignment.getSequenceCount()];
        diffInWin = new short[anchorAlignment.getSequenceCount()][anchorAlignment.getSequenceCount()];
        //       pOfMatch=new float[anchorAlignment.getSequenceCount()][anchorAlignment.getSequenceCount()];
//        int currWord=initialSite>>6;
        int maxWord = anchorAlignment.getAllelePresenceForAllSites(0, 0).getNumWords() - 1;
        //int maxWord = anchorAlignment.getTaxaBitsNoClone(0, 0).getNumWords() - 1;
        int startWord = startSite >> 6;
        int endWord = endSite >> 6;
        if (endWord > maxWord) {
            endWord = maxWord;
        }
//        int startWord=(currWord-windowSize<0)?0:currWord-windowSize;
//        int endWord=(currWord+windowSize>=maxWord)?maxWord-1:currWord+windowSize;
        System.out.printf("Start Site:%d Word:%d End Site:%d Word:%d %n", startSite, startWord, endSite, endWord);
//        int[] bins=new int[101];
//        int[] countClose=new int[anchorAlignment.getSequenceCount()];
        for (int i = 0; i < anchorAlignment.getSequenceCount(); i++) {
            if (maskHighHets && highHet[i]) {
                continue;
            }
            long[] imj = anchorAlignment.getAllelePresenceForAllSites(i, 0).getBits();
            //long[] imj = anchorAlignment.getTaxaBitsNoClone(i, 0).getBits();
            long[] imn = anchorAlignment.getAllelePresenceForAllSites(i, 1).getBits();
            //long[] imn = anchorAlignment.getTaxaBitsNoClone(i, 1).getBits();
            for (int j = 0; j < i; j++) {
                if (maskHighHets && highHet[j]) {
                    continue;
                }
                long[] jmj = anchorAlignment.getAllelePresenceForAllSites(j, 0).getBits();
                //long[] jmj = anchorAlignment.getTaxaBitsNoClone(j, 0).getBits();
                long[] jmn = anchorAlignment.getAllelePresenceForAllSites(j, 1).getBits();
                //long[] jmn = anchorAlignment.getTaxaBitsNoClone(j, 1).getBits();
                int same = 0, diff = 0, hets = 0;
//                int minorSame=0;
                for (int w = startWord; w <= endWord; w++) {
                    long ihetMask = ~(imj[w] & imn[w]);
                    long jhetMask = ~(jmj[w] & jmn[w]);
                    long imjnh = imj[w] & ihetMask;
                    long imnnh = imn[w] & ihetMask;
                    long jmjnh = jmj[w] & jhetMask;
                    long jmnnh = jmn[w] & jhetMask;
                    same += Long.bitCount(imjnh & jmjnh) + Long.bitCount(imnnh & jmnnh);
                    diff += Long.bitCount(imjnh & jmnnh) + Long.bitCount(imnnh & jmjnh);
//                    same+=Long.bitCount(imj[w]&jmj[w])+Long.bitCount(imn[w]&jmn[w]);
//                    diff+=Long.bitCount(imj[w]&jmn[w])+Long.bitCount(imn[w]&jmj[w]);
//                    hets+=Long.bitCount(mj[w]&mn[w])+Long.bitCount(jmn[w]&jmj[w]);
//                    minorSame+=Long.bitCount(imn[w]&jmn[w]);
                }
//               matchInWin[j][i]=matchInWin[i][j]=same;
//               diffInWin[j][i]=diffInWin[i][j]=diff;
                matchInWin[j][i] = matchInWin[i][j] = (short) (same - hets);
                diffInWin[j][i] = diffInWin[i][j] = (short) (diff - hets);
//               int sum=same+diff-(2*hets);
//               double div=(double)diffInWin[j][i]/(double)sum;
//               if((div<0.02)&&(sum>100)) {
//                   countClose[i]++;
//                   countClose[j]++;
//               }
                //                if (div<0.03) System.out.printf("%s %s %d %d %d %g %n",anchorAlignment.getTaxaName(i),anchorAlignment.getTaxaName(j),matchInWin[j][i],diffInWin[j][i],sum, div);
//                   bins[(int)((diff*100)/sum)]++;
//               }
            }
        }
//        for (int i = 0; i < bins.length; i++) {
//            System.out.printf("%d %d %d%n",initialSite,i,bins[i]);
//
//        }
//        for (int i = 0; i < countClose.length; i++) {
//            System.out.printf("%d %s %d%n",initialSite, anchorAlignment.getTaxaName(i),countClose[i]);
//
//        }
    }

    public void writeAlignment(String outfile) {
        ExportUtils.writeToHapmap(impAlign, false, outfile, '\t', null);
    }

    public static void main(String[] args) {
        System.out.println("Running main method in FastImputation");

        String[] dargs = {
       //     "-hmp", "/Users/edbuckler/SolexaAnal/GBS/build110813/test/maize110812.cov10.fT1E1pLD.mgNoHet.c10.hmp.txt",
            "-hmp","/Volumes/LaCie/GEM-DH_Production_January2012_FINAL_chr10.hmp.txt",
     //       "-o", "/Users/edbuckler/SolexaAnal/GBS/build110813/test/maize110812.cov10.fT1E1pLD.mgNoHet.imp.c10.hmp.txt",
            "-o", "/Volumes/LaCie/GEM-DH_Production_January2012_FINAL_chr10.imp.hmp.txt",
            "sC", "10",
            "eC", "10"
        };
        if (args.length == 0) {
            args = dargs;
        }
        int sC = Integer.parseInt(args[5]);
        int eC = Integer.parseInt(args[7]);
        String outfile, anchorMapFile;
        for (int cNum = sC; cNum <= eC; cNum++) {
            outfile = args[3].replace("+", "" + cNum);
            anchorMapFile = args[1].replace("+", "" + cNum);
            System.out.println("Reading " + anchorMapFile);
            Alignment a = ImportUtils.readFromHapmap(anchorMapFile, null);
            System.out.printf("Read Alignment with %d taxa and %d sites %n", a.getSequenceCount(), a.getSiteCount());
            System.out.println("p1a:" + a.getBaseAsStringRow(1));
            //a = TBitAlignment.getInstance(a);
            a = ConvertSBitTBitPlugin.convertAlignment(a, ConvertSBitTBitPlugin.CONVERT_TYPE.tbit, null);
            System.out.println("TBA:" + a.getBaseAsStringRow(1));
            FastImputationBitFixedWindow fi = new FastImputationBitFixedWindow(a);
            System.out.println("Writing " + outfile);
            fi.writeAlignment(outfile);
        }
    }
}
