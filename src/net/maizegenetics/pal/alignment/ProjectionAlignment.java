/*
 * ProjectionAlignment
 */
package net.maizegenetics.pal.alignment;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;
import java.util.TreeMap;

import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.util.BitSet;
import net.maizegenetics.util.ExceptionUtils;
import net.maizegenetics.util.Utils;

/**
 * This class projects high Density markers on large group of taxa through a
 * look up table system. The lookup table generally needs be built through some
 * imputation approach.
 *
 * @author ed, terry, Gabriel Rodrigues Alves Margarido
 */
public class ProjectionAlignment extends AbstractAlignment implements MutableAlignment {

    private int[][] mySiteBreaks;  //temporary - saving not needed
    private int[][][] myHDTaxa;  //taxa ids should be saved
    private int[][] myPosBreaks;  //positions should be saved
    private final Alignment myBaseAlignment;  //high density marker alignment that is being projected.
    
    private int[][] cacheTaxonSiteBound, cacheTaxonDonors; 
    public long cacheUseCnt=0, lookupCnt=0;
    

    public ProjectionAlignment(Alignment hdAlign, IdGroup ldIDGroup) {
        super(ldIDGroup, hdAlign.getAlleleEncodings());
        myBaseAlignment = AlignmentUtils.optimizeForSites(hdAlign);
        mySiteBreaks = new int[getSequenceCount()][];
        myPosBreaks = new int[getSequenceCount()][];
        myHDTaxa = new int[getSequenceCount()][][];
        init();
    }
    
    public ProjectionAlignment(Alignment hdAlign, IdGroup ldIDGroup, int[][] myPosBreaks, 
            int[][][] myHDTaxa) {
        super(ldIDGroup, hdAlign.getAlleleEncodings());
        myBaseAlignment = AlignmentUtils.optimizeForSites(hdAlign);
        this.myPosBreaks = myPosBreaks;
        this.myHDTaxa = myHDTaxa;
        mySiteBreaks = new int[getSequenceCount()][];
        for (int taxon = 0; taxon < myHDTaxa.length; taxon++) {
            if(myHDTaxa[taxon]==null) continue;
            mySiteBreaks[taxon] = new int[myPosBreaks[taxon].length];
            for (int i = 0; i < myPosBreaks[taxon].length; i++) {
                int site = myBaseAlignment.getSiteOfPhysicalPosition(myPosBreaks[taxon][i], null);
                if (site < 0) {
                    site = -(site + 1);
                }
                this.mySiteBreaks[taxon][i] = site;
            }   
        }
        init();
    }
    
    public static ProjectionAlignment getInstance(String paFile, String baseHighDensityAlignmentFile) {
        return getInstance(paFile, ImportUtils.readFromHapmap(baseHighDensityAlignmentFile, null));
    } 
    
    
    public static ProjectionAlignment getInstance(String paFile, Alignment baseHighDensityAlignment) {
        BufferedReader br = null;
        String s;
        try {
            br = Utils.getBufferedReader(paFile);
            String[] sl=Utils.readLineSkipComments(br).split("\t");
            int baseTaxaCnt=Integer.parseInt(sl[0]);
            if(baseTaxaCnt!=baseHighDensityAlignment.getSequenceCount()) {
                System.err.println("Error in number of base taxa"); return null;
            }
            int taxaCnt=Integer.parseInt(sl[1]);
            IdGroup aIDG=new SimpleIdGroup(taxaCnt);
            for (int i = 0; i < baseTaxaCnt; i++) {
                sl=Utils.readLineSkipComments(br).split("\t");
                //change to hash map 
                int index=Integer.parseInt(sl[0]);
                if(!baseHighDensityAlignment.getFullTaxaName(index).equals(sl[1])) {
                    System.err.println("Names or order does not agree with base taxa"); return null;
                }
            }
            int[][] myPosBreaks = new int[taxaCnt][];
            int[][][] myHDTaxa = new int[taxaCnt][][];
            for (int i = 0; i < taxaCnt; i++) {
                sl=Utils.readLineSkipComments(br).split("\t");
                aIDG.setIdentifier(i, new Identifier(sl[0]));
                int breakTotal=sl.length-1;
                if(breakTotal==0) continue;  //no data
                myPosBreaks[i]=new int[sl.length-1];
                myHDTaxa[i]=new int[sl.length-1][2];
                for (int bp = 0; bp < myHDTaxa[i].length; bp++) {
                    String[] bptext=sl[bp+1].split(":");
                    myPosBreaks[i][bp]=Integer.parseInt(bptext[0]);
                    myHDTaxa[i][bp][0]=Integer.parseInt(bptext[1]);
                    myHDTaxa[i][bp][1]=Integer.parseInt(bptext[2]);
                }
            }
            return (new ProjectionAlignment(baseHighDensityAlignment, aIDG, myPosBreaks, myHDTaxa));
        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalArgumentException("Error reading Projection file: " + paFile + ": " + ExceptionUtils.getExceptionCauses(e));
        } finally {
            try {br.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        } 
        //return null;
    }
    
    private void init() {
        myNumSites=myBaseAlignment.getSiteCount();
        cacheTaxonSiteBound=new int[getSequenceCount()][2];
        cacheTaxonDonors=new int[getSequenceCount()][2];
        for (int i = 0; i < getSequenceCount(); i++) {
         //   translateTaxon(i,0);
            cacheNewTaxonSiteRange(i,0);
        }
//        System.out.println(Arrays.toString(mySiteBreaks[0]));
//        for (int i = 0; i < 16000; i++) {
//            System.out.println(i);
//            System.out.println("TOrig"+Arrays.toString(translateTaxon(0,i)));
//            System.out.println("TNew"+Arrays.toString(translateTaxonX(0,i)));
//            System.out.println("cacheB:"+Arrays.toString(cacheTaxonSiteBound[0]));
//        }
    }
    
    public void save(String outfile) {
        BufferedWriter bw = null;
        try {
            String fullFileName = Utils.addSuffixIfNeeded(outfile, ".pa.txt.gz", new String[]{".pa.txt", ".pa.txt.gz"});
            bw = Utils.getBufferedWriter(fullFileName);
            bw.write(myBaseAlignment.getSequenceCount()+"\t"+getSequenceCount()+"\n");
            bw.write("#Donor Haplotypes\n");
            for (int i = 0; i < myBaseAlignment.getSequenceCount(); i++) {
                bw.write(i+"\t"+myBaseAlignment.getFullTaxaName(i)+"\n");
            }
            bw.write("#Taxa Breakpoints\n");
            bw.write("#Block are defined position:donor1:donor2 (-1 means no hypothesis)\n");
            for (int i = 0; i < getSequenceCount(); i++) {
                bw.write(getFullTaxaName(i)+"\t");
                for (int p = 0; (myPosBreaks[i]!=null)&&(p < myPosBreaks[i].length); p++) {
                    bw.write(myPosBreaks[i][p]+":"+myHDTaxa[i][p][0]+":"+myHDTaxa[i][p][1]+"\t");
                }
                bw.write("\n");
            }   
        } catch (Exception e) {
            e.printStackTrace();
            throw new IllegalArgumentException("Error writing Projection file: " + outfile + ": " + ExceptionUtils.getExceptionCauses(e));
        } finally {
            try {bw.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        } 
    }
    


    public void setCompositionOfTaxon(int taxon, int[] posBreaks, int[][] hdTaxa) {
        if (posBreaks.length != hdTaxa.length) {
            throw new IllegalArgumentException("ProjectionAlignment: setCompositionOfTaxon: number of positions should equal number taxa.");
        }
        myPosBreaks[taxon] = posBreaks;
        myHDTaxa[taxon] = hdTaxa;
        mySiteBreaks[taxon] = new int[myPosBreaks[taxon].length];
        for (int i = 0; i < myPosBreaks[taxon].length; i++) {
            int site = myBaseAlignment.getSiteOfPhysicalPosition(posBreaks[i], null);
            if (site < 0) {
                site = -(site + 1);
            }
            this.mySiteBreaks[taxon][i] = site;
        }
    }
    
    public void setCompositionOfTaxon(int taxon, TreeMap<Integer,int[]> breakPoints) {
        breakPoints=removeBreakPointReduncdancy(breakPoints);
        int breaks=breakPoints.size();
        myPosBreaks[taxon]=new int[breaks];
        myHDTaxa[taxon] = new int[breaks][2];
        mySiteBreaks[taxon] = new int[breaks];
        int cnt=0;
        for (Map.Entry<Integer,int[]> bp : breakPoints.entrySet()) {
            myPosBreaks[taxon][cnt]=bp.getKey();
            myHDTaxa[taxon][cnt]=bp.getValue();
            int site = myBaseAlignment.getSiteOfPhysicalPosition( myPosBreaks[taxon][cnt], null);
            if (site < 0) {
                site = -(site + 1);
            }
            this.mySiteBreaks[taxon][cnt] = site;
            cnt++;
        }
    }
    
    private TreeMap<Integer,int[]> removeBreakPointReduncdancy(TreeMap<Integer,int[]> breakPoints) {
        int[] lastP=new int[]{-2,-2};
        ArrayList<Integer> posToRemove=new ArrayList<Integer>();
        for (Map.Entry<Integer,int[]> bp : breakPoints.entrySet()) {
            if(Arrays.equals(lastP,bp.getValue())) {
                  posToRemove.add(bp.getKey());
            } else {
                lastP=bp.getValue();
            }
        }
        for (Integer pos : posToRemove) {breakPoints.remove(pos);}
        return breakPoints;
    }
    
    public void setCompositionOfTaxon(String taxa, int[] posBreaks, int[][] hdTaxa) {
        int taxon = getIdGroup().whichIdNumber(taxa);
        System.out.printf("SetComp %s %d %n", taxa, taxon);
        setCompositionOfTaxon(taxon, posBreaks, hdTaxa);
    }

    public String getCompositionOfTaxon(int taxon) {
        StringBuilder sb = new StringBuilder(this.getIdGroup().getIdentifier(taxon).getFullName() + "\t");
        if (myPosBreaks[taxon] == null) {
            sb.append("NULL");
        } else {
            for (int i = 0; i < myPosBreaks[taxon].length; i++) {
                sb.append(myPosBreaks[taxon][i] + ":");
                sb.append(mySiteBreaks[taxon][i] + ":");
                sb.append(myHDTaxa[taxon][i] + "\t");
            }
        }
        return sb.toString();
    }

    public void reportPAComposition() {
        for (int i = 0; i < this.getSequenceCount(); i++) {
            System.out.println(getCompositionOfTaxon(i));
        }
    }

    private int[] translateTaxon(int taxon, int site) {
        if (mySiteBreaks[taxon] == null) {
            return null;
        }
        if((cacheTaxonSiteBound[taxon][0]<=site)&&(site<=cacheTaxonSiteBound[taxon][1])) {
            cacheUseCnt++;
            return cacheTaxonDonors[taxon];
        } else {
            lookupCnt++;
            return cacheNewTaxonSiteRange(taxon, site);
        }
//        int b = Arrays.binarySearch(mySiteBreaks[taxon], site);
//        if (b < 0) {
//            b = -(b + 2);  //this will not work if it does not start with zero.
//        }
//        return myHDTaxa[taxon][b];
    }
    
    private int[] translateTaxonX(int taxon, int site) {
        if (mySiteBreaks[taxon] == null) {
            return null;
        }
        int b = Arrays.binarySearch(mySiteBreaks[taxon], site);
        if (b < 0) {
            b = -(b + 2);  //this will not work if it does not start with zero.
        }
        return myHDTaxa[taxon][b];
    }
    
    private int[] cacheNewTaxonSiteRange(int taxon, int site){
        if (mySiteBreaks[taxon] == null) return null;
        int b = Arrays.binarySearch(mySiteBreaks[taxon], site);
        if (b < 0) {
            b = -(b + 2);  //this will not work if it does not start with zero.
        }
        cacheTaxonSiteBound[taxon][0]=mySiteBreaks[taxon][b];
        if((b+1)<mySiteBreaks[taxon].length) {
            cacheTaxonSiteBound[taxon][1]=mySiteBreaks[taxon][b+1];
        } else {
            cacheTaxonSiteBound[taxon][1]=myNumSites;
        }
        cacheTaxonDonors[taxon]=myHDTaxa[taxon][b];
        return myHDTaxa[taxon][b];
    }

    @Override
    public String getBaseAsString(int taxon, int site) { 
        return NucleotideAlignmentConstants.getNucleotideIUPAC(getBase(taxon, site));
    }

    @Override
    public String getDiploidAsString(int site, byte value) {
        return NucleotideAlignmentConstants.getNucleotideIUPAC(value);
    }

    @Override
    /**
     * This is the slow implementation of this. Most of these should be buffered
     * bit sets. Note the imputation of the taxon is likely to be the the same
     * for 64 or more sites in a row (potentially, 10,000s of sites in many
     * cases).
     */
    public byte getBase(int taxon, int site) {
        int[] t = translateTaxon(taxon, site);
        if((t==null)||(t[0]<0)) return Alignment.UNKNOWN_DIPLOID_ALLELE;
        byte b0, b1;
        b0=myBaseAlignment.getBase(t[0], site);
        if(t[0]==t[1]) {b1=b0;}
        else {b1=myBaseAlignment.getBase(t[1], site);}
        return AlignmentUtils.getUnphasedDiploidValueNoHets(b0, b1);
    }

    @Override
    public boolean isSBitFriendly() {
        return false;
    }

    @Override
    public boolean isTBitFriendly() {
        return false;
    }

//    @Override
//    public byte[] getBaseRow(int taxon) {
//
//        int numBreaks = mySiteBreaks[taxon].length;
//        byte[] result = new byte[myNumSites];
//        for (int i = 0; i < numBreaks - 1; i++) {
//            int[] hdTaxon = myHDTaxa[taxon][i];
//            for (int j = mySiteBreaks[taxon][i]; j < mySiteBreaks[taxon][i + 1]; j++) {
//               // result[j] = myBaseAlignment.getBase(hdTaxon, j);
//                result[j]=AlignmentUtils.getDiploidValue(myBaseAlignment.getBase(hdTaxon[0], j), myBaseAlignment.getBase(hdTaxon[1], j));
//            }
//        }
//
//        int[] hdTaxon = myHDTaxa[taxon][numBreaks - 1];
//        for (int j = mySiteBreaks[taxon][numBreaks - 1], n = getSiteCount(); j < n; j++) {
// //           result[j] = myBaseAlignment.getBase(hdTaxon, j);
//            result[j]=AlignmentUtils.getDiploidValue(myBaseAlignment.getBase(hdTaxon[0], j), myBaseAlignment.getBase(hdTaxon[1], j));
//        }
//
//        return result;
//    }

    @Override
    // TERRY - This Could be Optimized like getBaseRow()
    public byte[] getBaseRange(int taxon, int startSite, int endSite) {

        byte[] result = new byte[endSite - startSite];
        for (int i = startSite; i < endSite; i++) {
            result[i] = getBase(taxon, i);
        }
        return result;

    }

    @Override
    public byte getBase(int taxon, Locus locus, int physicalPosition) {
        return getBase(taxon, getSiteOfPhysicalPosition(physicalPosition, locus));
    }

    @Override
    public BitSet getAllelePresenceForAllSites(int taxon, int alleleNumber) {
        throw new UnsupportedOperationException();
    }

//    @Override
//    public BitSet getAllelePresenceForAllTaxa(int site, int alleleNumber) {
//        BitSet baseBitSet = myBaseAlignment.getAllelePresenceForAllTaxa(site, alleleNumber);
//        BitSet result = new OpenBitSet(getSequenceCount());
//        for (int i = 0, n = getSequenceCount(); i < n; i++) {
//            int index = translateTaxon(i, site);
//            if (baseBitSet.fastGet(index)) {
//                result.fastSet(i);
//            }
//        }
//        return result;
//    }

    @Override
    public long[] getAllelePresenceForSitesBlock(int taxon, int alleleNumber, int startBlock, int endBlock) {
        throw new UnsupportedOperationException();
    }

    @Override
    public int getIndelSize(int site) {
        return myBaseAlignment.getIndelSize(site);
    }

    @Override
    public boolean isIndel(int site) {
        return isIndel(site);
    }

    @Override
    public float getSiteScore(int seq, int site) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
        //return myBaseAlignment.getSiteScore(translateTaxon(seq, site), site);
    }

    @Override
    public boolean hasSiteScores() {
        return false;
    }

    @Override
    public SITE_SCORE_TYPE getSiteScoreType() {
        return Alignment.SITE_SCORE_TYPE.None;
    }

    @Override
    public boolean hasReference() {
        return myBaseAlignment.hasReference();
    }

    @Override
    public byte getReferenceAllele(int site) {
        return myBaseAlignment.getReferenceAllele(site);
    }

    @Override
    public byte[] getReference() {
        return myBaseAlignment.getReference();
    }

    @Override
    public byte[] getReference(int startSite, int endSite) {
        return myBaseAlignment.getReference(startSite, endSite);
    }

    @Override
    public boolean isPhased() {
        return false;
    }

    @Override
    public boolean isPositiveStrand(int site) {
        return myBaseAlignment.isPositiveStrand(site);
    }

    @Override
    public String getGenomeAssembly() {
        return myBaseAlignment.getGenomeAssembly();
    }

    @Override
    public GeneticMap getGeneticMap() {
        return myBaseAlignment.getGeneticMap();
    }

    @Override
    public int[] getPhysicalPositions() {
        return myBaseAlignment.getPhysicalPositions();
    }

    @Override
    public int getSiteCount() {
        return myBaseAlignment.getSiteCount();
    }

    @Override
    public int getPositionInLocus(int site) {
        return myBaseAlignment.getPositionInLocus(site);
    }

    @Override
    public Locus getLocus(int site) {
        return myBaseAlignment.getLocus(site);
    }

    @Override
    public Locus[] getLoci() {
        return myBaseAlignment.getLoci();
    }

    @Override
    public int getNumLoci() {
        return myBaseAlignment.getNumLoci();
    }

    @Override
    public int[] getLociOffsets() {
        return myBaseAlignment.getLociOffsets();
    }

    @Override
    public int getLocusSiteCount(Locus locus) {
        return myBaseAlignment.getLocusSiteCount(locus);
    }

    @Override
    public String[] getSNPIDs() {
        return myBaseAlignment.getSNPIDs();
    }

    @Override
    public String getSNPID(int site) {
        return myBaseAlignment.getSNPID(site);
    }

    @Override
    public int getSiteOfPhysicalPosition(int physicalPosition, Locus locus) {
        return myBaseAlignment.getSiteOfPhysicalPosition(physicalPosition, locus);
    }

    @Override
    public int getSiteOfPhysicalPosition(int physicalPosition, Locus locus, String snpID) {
        return myBaseAlignment.getSiteOfPhysicalPosition(physicalPosition, locus, snpID);
    }

    @Override
    public boolean retainsRareAlleles() {
        return myBaseAlignment.retainsRareAlleles();
    }

    @Override
    public String[][] getAlleleEncodings() {
        return myBaseAlignment.getAlleleEncodings();
    }

    @Override
    public String[] getAlleleEncodings(int site) {
        return myBaseAlignment.getAlleleEncodings(site);
    }

    @Override
    public String getBaseAsString(int site, byte value) {
        return myBaseAlignment.getBaseAsString(site, value);
    }

    @Override
    public int getMaxNumAlleles() {
        return myBaseAlignment.getMaxNumAlleles();
    }

//    @Override
//    public int getTotalGametesNotMissingForTaxon(int taxon) {
//
//        int numBreaks = mySiteBreaks[taxon].length;
//        int result = 0;
//        for (int i = 0; i < numBreaks - 1; i++) {
//            int hdTaxon = myHDTaxa[taxon][i];
//            for (int j = mySiteBreaks[taxon][i]; j < mySiteBreaks[taxon][i + 1]; j++) {
//                byte[] current = myBaseAlignment.getBaseArray(hdTaxon, j);
//                if (current[0] != Alignment.UNKNOWN_ALLELE) {
//                    result++;
//                }
//                if (current[1] != Alignment.UNKNOWN_ALLELE) {
//                    result++;
//                }
//            }
//        }
//
//        int hdTaxon = myHDTaxa[taxon][numBreaks - 1];
//        for (int j = mySiteBreaks[taxon][numBreaks - 1], n = getSiteCount(); j < n; j++) {
//            byte[] current = myBaseAlignment.getBaseArray(hdTaxon, j);
//            if (current[0] != Alignment.UNKNOWN_ALLELE) {
//                result++;
//            }
//            if (current[1] != Alignment.UNKNOWN_ALLELE) {
//                result++;
//            }
//        }
//
//        return result;
//
//    }
//
//    @Override
//    public int getTotalNotMissingForTaxon(int taxon) {
//
//        int numBreaks = mySiteBreaks[taxon].length;
//        int result = 0;
//        for (int i = 0; i < numBreaks - 1; i++) {
//            int hdTaxon = myHDTaxa[taxon][i];
//            for (int j = mySiteBreaks[taxon][i]; j < mySiteBreaks[taxon][i + 1]; j++) {
//                byte current = myBaseAlignment.getBase(hdTaxon, j);
//                if (current != Alignment.UNKNOWN_DIPLOID_ALLELE) {
//                    result++;
//                }
//            }
//        }
//
//        int hdTaxon = myHDTaxa[taxon][numBreaks - 1];
//        for (int j = mySiteBreaks[taxon][numBreaks - 1], n = getSiteCount(); j < n; j++) {
//            byte current = myBaseAlignment.getBase(hdTaxon, j);
//            if (current != Alignment.UNKNOWN_DIPLOID_ALLELE) {
//                result++;
//            }
//        }
//
//        return result;
//
//    }

    // TERRY - This Needs Work.
    public Object[][] getMajorMinorCounts() {

        if (myAlleleStates.length != 1) {
            return new Object[0][0];
        }

        long[][] counts = new long[16][16];

        if (myMaxNumAlleles >= 2) {
            for (int site = 0; site < myNumSites; site++) {
                byte indexI = myAlleles[site][0];
                byte indexJ = myAlleles[site][1];
                if (indexJ == UNKNOWN_ALLELE) {
                    indexJ = indexI;
                }
                counts[indexI][indexJ]++;
            }
        } else {
            for (int site = 0; site < myNumSites; site++) {
                byte indexI = myAlleles[site][0];
                counts[indexI][indexI]++;
            }
        }

        int numAlleles = 0;
        for (byte x = 0; x < 16; x++) {
            for (byte y = 0; y < 16; y++) {
                if (counts[x][y] != 0) {
                    numAlleles++;
                }
            }
        }

        Object[][] result = new Object[2][numAlleles];
        int nextResult = 0;
        for (byte x = 0; x < 16; x++) {
            for (byte y = 0; y < 16; y++) {
                if (counts[x][y] != 0) {
                    result[0][nextResult] = getBaseAsString(0, x) + ":" + getBaseAsString(0, y);
                    result[1][nextResult++] = counts[x][y];
                }
            }
        }

        boolean change = true;
        while (change) {

            change = false;

            for (int k = 0; k < numAlleles - 1; k++) {

                if ((Long) result[1][k] < (Long) result[1][k + 1]) {

                    Object temp = result[0][k];
                    result[0][k] = result[0][k + 1];
                    result[0][k + 1] = temp;

                    Object tempCount = result[1][k];
                    result[1][k] = result[1][k + 1];
                    result[1][k + 1] = tempCount;

                    change = true;
                }
            }

        }

        return result;
    }

//    @Override
//    public int[][] getAllelesSortedByFrequency(int site) {
//
//        int maxNumAlleles = myBaseAlignment.getMaxNumAlleles();
//        int numDataRows = myBaseAlignment.getTotalNumAlleles();
//        BitSet[] data = new BitSet[numDataRows];
//        for (int i = 0; i < numDataRows; i++) {
//            data[i] = getAllelePresenceForAllTaxa(site, i);
//        }
//        byte[] alleles = myBaseAlignment.getAlleles(site);
//
//        int[] counts = new int[16];
//        for (int i = 0; i < numDataRows; i++) {
//            byte indexI;
//            if ((retainsRareAlleles()) && (i == maxNumAlleles)) {
//                indexI = Alignment.RARE_ALLELE;
//            } else {
//                indexI = alleles[i];
//            }
//            counts[indexI] += (int) data[i].cardinality() * 2;
//            for (int j = i + 1; j < numDataRows; j++) {
//                byte indexJ;
//                if ((retainsRareAlleles()) && (j == maxNumAlleles)) {
//                    indexJ = Alignment.RARE_ALLELE;
//                } else {
//                    indexJ = alleles[j];
//                }
//                int ijHet = (int) OpenBitSet.intersectionCount(data[i], data[j]);
//                counts[indexI] -= ijHet;
//                counts[indexJ] -= ijHet;
//            }
//        }
//
//        int numAlleles = 0;
//        for (byte x = 0; x < Alignment.UNKNOWN_ALLELE; x++) {
//            if (counts[x] != 0) {
//                numAlleles++;
//            }
//        }
//
//        int current = 0;
//        int[][] result = new int[2][numAlleles];
//        for (byte x = 0; x < Alignment.UNKNOWN_ALLELE; x++) {
//            if (counts[x] != 0) {
//                result[0][current] = x;
//                result[1][current++] = counts[x];
//            }
//        }
//
//        boolean change = true;
//        while (change) {
//
//            change = false;
//
//            for (int k = 0; k < numAlleles - 1; k++) {
//
//                if (result[1][k] < result[1][k + 1]) {
//
//                    int temp = result[0][k];
//                    result[0][k] = result[0][k + 1];
//                    result[0][k + 1] = temp;
//
//                    int tempCount = result[1][k];
//                    result[1][k] = result[1][k + 1];
//                    result[1][k + 1] = tempCount;
//
//                    change = true;
//                }
//            }
//
//        }
//
//        return result;
//
//    }

    @Override
    public int getTotalNumAlleles() {
        return myBaseAlignment.getTotalNumAlleles();
    }

    @Override
    public void setBase(int taxon, int site, byte newBase) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void setBaseRange(int taxon, int startSite, byte[] newBases) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void addSite(int site) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void removeSite(int site) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void clearSiteForRemoval(int site) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void addTaxon(Identifier id) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void setTaxonName(int taxon, Identifier id) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void removeTaxon(int taxon) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void setPositionOfSite(int site, int position) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void setLocusOfSite(int site, Locus locus) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void setDepthForAlleles(int taxon, int site, short[] values) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void setCommonAlleles(int site, byte[] values) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void clean() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void setReferenceAllele(int site, byte diploidAllele) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    @Override
    public void setSNPID(int site, String name) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public boolean isDirty() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
}
