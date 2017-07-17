/*
 * TagsAtLocus
 */
package net.maizegenetics.gbs.maps;

import org.apache.commons.math.distribution.BinomialDistributionImpl;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import net.maizegenetics.gbs.tagdist.TagsByTaxa;
import net.maizegenetics.gbs.util.BaseEncoder;
import net.maizegenetics.pal.alignment.Alignment;
import java.util.List;
import org.biojava3.core.sequence.DNASequence;
import org.biojava3.core.sequence.compound.AmbiguityDNACompoundSet;
import org.biojava3.alignment.template.Profile;
import org.biojava3.alignment.Alignments;
import org.biojava3.alignment.SimpleGapPenalty;
import org.biojava3.alignment.SubstitutionMatrixHelper;
import org.biojava3.alignment.template.SequencePair;
import org.biojava3.alignment.template.SubstitutionMatrix;
import org.biojava3.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.pal.alignment.Locus;
import net.maizegenetics.pal.alignment.BitAlignment;
import net.maizegenetics.pal.alignment.AlignmentUtils;
import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;
import net.maizegenetics.prefs.TasselPrefs;
import net.maizegenetics.util.VCFUtil;

/**
 * Aligns and calls SNPs for all tags at a given locus.
 * 
 * @author jcg233, Gabriel Rodrigues Alves Margarido
 */
public class TagsAtLocus {

    ArrayList<SingleTagByTaxa> theTags = new ArrayList<SingleTagByTaxa>();
    private int minStartPosition;
    private int maxStartPosition;
    private int minTagLength;
    private int maxTagLength;
    private int chromosome;
    private byte strand;
    private int indexOfRef = Integer.MIN_VALUE;
    private int[] tagIndices = null;  // redirect from aligned tag indices to index in theTags
    private int[] positionsOfVariableSites;
    private byte[][] allelesAtVariableSitesByTag;
    private byte[] refCallsBySite = null;
    private int nTaxaCovered = Integer.MIN_VALUE;
    private int totalNReads = Integer.MIN_VALUE;
    private String status = "notSet";
    
    // For VCF output with depth (now used for custom SNP report in regular pipeline as well)
    private byte[][] myCommonAlleles = null;
    private short[][][] myAlleleDepthsInTaxa = null;
        
    private final static int maxSNPsPerLocus = 64;
    private final static int maxAlignmentSize = 10000;
    private final static int maxCountAtGeno = 500;
    private final static int maxNumAlleles = 3;
    private SubstitutionMatrix<NucleotideCompound> subMatrix = SubstitutionMatrixHelper.getNuc4_4();
    private SimpleGapPenalty gapPen = new SimpleGapPenalty((short) 5, (short) 2);
    private static int[] likelihoodRatioThreshAlleleCnt = null;  // index = sample size; value = min count of less tagged allele for likelihood ratio > 1
    // if less tagged allele has counts < likelihoodRatioThreshAlleleCnt[totalCount], call it a homozygote
    // where likelihood ratio = (binomial likelihood het) / (binomial likelihood all less tagged alleles are errors)

    static void setLikelihoodThresh(double errorRate) {   // initialize the likelihood ratio cutoffs for quantitative SNP calling
        likelihoodRatioThreshAlleleCnt = new int[maxCountAtGeno];
        System.out.println("\n\nInitializing the cutoffs for quantitative SNP calling likelihood ratio (pHet/pErr) >1\n");
        System.out.println("totalReadsForSNPInIndiv\tminLessTaggedAlleleCountForHet");
        for (int trials = 0; trials < 2; ++trials) {
            likelihoodRatioThreshAlleleCnt[trials] = 1;
        }
        int lastThresh = 1;
        for (int trials = 2; trials < likelihoodRatioThreshAlleleCnt.length; ++trials) {
            BinomialDistributionImpl binomHet = new BinomialDistributionImpl(trials, 0.5);
            BinomialDistributionImpl binomErr = new BinomialDistributionImpl(trials, errorRate);
            double LikeRatio;
            try {
                LikeRatio = binomHet.cumulativeProbability(lastThresh) / (1 - binomErr.cumulativeProbability(lastThresh) + binomErr.probability(lastThresh));
                while (LikeRatio <= 1.0) {
                    ++lastThresh;
                    LikeRatio = binomHet.cumulativeProbability(lastThresh) / (1 - binomErr.cumulativeProbability(lastThresh) + binomErr.probability(lastThresh));
                }
                likelihoodRatioThreshAlleleCnt[trials] = lastThresh;
                System.out.println(trials + "\t" + lastThresh);
            } catch (Exception e) {
                System.err.println("Error in the TagsAtLocus.BinomialDistributionImpl");
            }
        }
        System.out.println("\n");
    }

    public TagsAtLocus(int chromosome, byte strand, int startPosition, int tagLength, 
            boolean includeRefGenome, boolean fuzzyStartPositions, double errorRate) {
        this.chromosome = chromosome;
        this.strand = (includeRefGenome && fuzzyStartPositions) ? 1 : strand;
        this.minStartPosition = startPosition;
        this.maxStartPosition = startPosition;
        this.minTagLength = tagLength;
        this.maxTagLength = tagLength;
        positionsOfVariableSites = null;
        allelesAtVariableSitesByTag = null;
        if (likelihoodRatioThreshAlleleCnt == null) {
            setLikelihoodThresh(errorRate);
        }
    }

    public void addTag(int tagTOPMIndex, TagsOnPhysicalMapV3 theTOPM, TagsByTaxa theTBT, 
            boolean includeRefGenome, boolean fuzzyStartPositions) {
            SingleTagByTaxa singleTBT = new SingleTagByTaxa(tagTOPMIndex, theTOPM, theTBT, 
                includeRefGenome, fuzzyStartPositions);
        if (singleTBT.taxaWithTag > 0) {
            theTags.add(singleTBT);
            if (singleTBT.startPosition > minStartPosition) {
                maxStartPosition = singleTBT.startPosition;
            }
            if (singleTBT.tagLength < minTagLength) {
                minTagLength = singleTBT.tagLength;
            }
            if (singleTBT.tagLength > maxTagLength) {
                maxTagLength = singleTBT.tagLength;
            }
        }
    }
    
    public void addTag(int tagTOPMIndex, TagsOnPhysicalMap theTOPM, TagsByTaxa theTBT, 
            boolean includeRefGenome, boolean fuzzyStartPositions) {
        SingleTagByTaxa singleTBT = new SingleTagByTaxa(tagTOPMIndex, theTOPM, theTBT, 
                includeRefGenome, fuzzyStartPositions);
        if (singleTBT.taxaWithTag > 0) {
            theTags.add(singleTBT);
            if (singleTBT.startPosition > minStartPosition) {
                maxStartPosition = singleTBT.startPosition;
            }
            if (singleTBT.tagLength < minTagLength) {
                minTagLength = singleTBT.tagLength;
            }
            if (singleTBT.tagLength > maxTagLength) {
                maxTagLength = singleTBT.tagLength;
            }
        }
    }
    
    public void addRefTag(String refTag, int nLongsPerTag, String nullTag) {
        theTags.add(new SingleTagByTaxa(minStartPosition, strand, refTag, nLongsPerTag, nullTag));
        indexOfRef = theTags.size()-1;
    }

    public int getSize() {
        return theTags.size();
    }

    public int getChromosome() {
        return chromosome;
    }

    public byte getStrand() {
        return strand;
    }

    public int getMinStartPosition() {
        return minStartPosition;
    }

    public int getMaxStartPosition() {
        return maxStartPosition;
    }
    
    public int getMinTagLength() {
        return minTagLength;
    }
    
    public int getMaxTagLength() {
        return maxTagLength;
    }

    public void setMinStartPosition(int newMinStartPosition) {
        minStartPosition = newMinStartPosition;
    }

    public int getTOPMIndexOfTag(int tagIndex) {
        return theTags.get(tagIndex).tagTOPMIndex;
    }

    public int getTBTIndexOfTag(int tagIndex) {
        return theTags.get(tagIndex).tagTBTIndex;
    }

    public int getDivergenceOfTag(int tagIndex) {
        return theTags.get(tagIndex).divergence;
    }

    public byte getCallAtVariableSiteForTag(int site, int tagIndex) {
        return allelesAtVariableSitesByTag[site][tagIndex];
    }
    
    public byte getRefGeno(int site) {
        if (refCallsBySite == null || site > refCallsBySite.length-1) {
            return Alignment.UNKNOWN_DIPLOID_ALLELE;
        } else {
            return refCallsBySite[site];
        }
    }

    public int getNumberTaxaCovered() {
        if (theTags.size() < 1) {
            return 0;
        }
        if (nTaxaCovered == Integer.MIN_VALUE) {
            nTaxaCovered = 0;
            totalNReads = 0;
            boolean[] covered = new boolean[theTags.get(0).tagDist.length];  // initializes to false
            for (SingleTagByTaxa sTBT : theTags) {
                for (int tx = 0; tx < covered.length; ++tx) {
                    int reads = sTBT.tagDist[tx];
                    totalNReads += reads;
                    if (!covered[tx] && reads > 0) {
                        covered[tx] = true;
                    }
                }
            }
            for (int tx = 0; tx < covered.length; ++tx) {
                if (covered[tx]) {
                    ++nTaxaCovered;
                }
            }
            return nTaxaCovered; 
        } else {
            return nTaxaCovered;
        }
    }

    public int getTotalNReads() {
        if (theTags.size() < 1) {
            return 0;
        }
        if (totalNReads == Integer.MIN_VALUE) {
            getNumberTaxaCovered();
            return totalNReads; 
        } else {
            return totalNReads;
        }
    }

    private void assignRefTag() {
        int lengthOfRef = Integer.MIN_VALUE;
        int tagIndex = 0;
        for (SingleTagByTaxa sTBT : theTags) {
            if (sTBT.divergence == 0 && sTBT.tagLength > lengthOfRef) {
                indexOfRef = tagIndex;
                lengthOfRef = sTBT.tagLength;
            }
            ++tagIndex;
        }
    }
    
    public byte[][] getCommonAlleles() {
        return myCommonAlleles;
    }
    
    public short[][][] getAlleleDepthsInTaxa() {
        return myAlleleDepthsInTaxa;
    }

    // Qi's SNP caller for VCF output
    public byte[][] getSNPCallsVCF(boolean callBiallelicSNPsWithGap, boolean includeReferenceTag) {
        if (theTags.size() < 2) {
            status = "invariant";
            return null;
        }
        Alignment tagAlignment = getVariableSites();
        if (tagAlignment == null || tagAlignment.getSiteCount() < 1) {
            status = "invariant";
            return null;
        }
        int nSites = tagAlignment.getSiteCount();
        int nTaxa = theTags.get(0).tagDist.length;
        if (nTaxa < 1) {
            status = "noTaxa";
            return null;
        }
        status = "polymorphic";
        byte[][] callsBySite = new byte[nSites][nTaxa];
        if (includeReferenceTag) refCallsBySite = new byte[nSites];
        populateAllelesAtVariableSitesByTag(tagAlignment, nSites, includeReferenceTag, callBiallelicSNPsWithGap);
        positionsOfVariableSites = new int[nSites];
        myCommonAlleles = new byte[maxNumAlleles][nSites];
        myAlleleDepthsInTaxa = new short[maxNumAlleles][nSites][nTaxa];
        for (int s = 0; s < nSites; s++) {
            positionsOfVariableSites[s] = tagAlignment.getPositionInLocus(s);
            byte[] commonAlleles = getCommonAlleles(s, nTaxa, includeReferenceTag);
            int[][] alleleDepthsInTaxa = getAlleleDepthsInTaxa(commonAlleles, s, nTaxa, includeReferenceTag);
            setAlleleDepthsInTaxaForSite(s, alleleDepthsInTaxa, commonAlleles);
            for (int tx = 0; tx < nTaxa; tx++) {
                callsBySite[s][tx] = VCFUtil.resolveVCFGeno(commonAlleles, alleleDepthsInTaxa, tx);
            }
        }
        return callsBySite;
    }

    public byte[][] getSNPCallsQuant(boolean callBiallelicSNPsWithGap, boolean includeReferenceTag) {
        if (theTags.size() < 2) {
            status = "invariant";
            return null;
        }
        Alignment tagAlignment = this.getVariableSites();
        if (tagAlignment == null || tagAlignment.getSiteCount() < 1) {
            status = "invariant";
            return null;
        }
        int nSites = tagAlignment.getSiteCount();
        int nTaxa = theTags.get(0).tagDist.length;
        if (nTaxa < 1) {
            status = "noTaxa";  // this shouldn't happen but is here just as a check
            return null;
        }
        status = "polymorphic";
        byte[][] callsBySite = new byte[nSites][nTaxa];
        if (includeReferenceTag) refCallsBySite = new byte[nSites];
        populateAllelesAtVariableSitesByTag(tagAlignment, nSites, includeReferenceTag, callBiallelicSNPsWithGap);
        positionsOfVariableSites = new int[nSites];
        myCommonAlleles = new byte[maxNumAlleles][nSites];
        myAlleleDepthsInTaxa = new short[maxNumAlleles][nSites][nTaxa];
        for (int s = 0; s < nSites; s++) {
            positionsOfVariableSites[s] = tagAlignment.getPositionInLocus(s);
            byte[] commonAlleles = getCommonAlleles(s, nTaxa, includeReferenceTag); // NOTE: gap could be one of the common alleles (even if callBiallelicSNPsWithGap is false)
            int[][] alleleDepthsInTaxa = getAlleleDepthsInTaxa(commonAlleles, s, nTaxa, includeReferenceTag);
            setAlleleDepthsInTaxaForSite(s, alleleDepthsInTaxa, commonAlleles);
            for (int tx = 0; tx < nTaxa; tx++) {
                int count = 0;
                for (int a = 0; a < maxNumAlleles; a++) {
                    count += alleleDepthsInTaxa[a][tx];
                }
                if (count == 0) {
                    callsBySite[s][tx] = Alignment.UNKNOWN_DIPLOID_ALLELE;
                    continue;
                }
                // check for each possible homozygote
                boolean done = false;
                for (int a = 0; a < maxNumAlleles; a++) {
                    if ((count - alleleDepthsInTaxa[a][tx]) == 0) {
                        callsBySite[s][tx] = (byte) ((commonAlleles[a] << 4) | commonAlleles[a]);
                        done = true;
                        break;
                    }
                }
                if (done) {
                    continue;
                }
                callsBySite[s][tx] = resolveHetGeno(commonAlleles, alleleDepthsInTaxa, tx);
            }
        }
        return callsBySite;
    }

    public byte[][] getSNPCallsQuant(String refSeq, boolean callBiallelicSNPsWithGap) {
        // ToDo: UPDATE THIS FOR Tassel4 allele encoding
        if (theTags.size() < 2) {
            return null;
        }
        Alignment tagAlignment = this.getVariableSites(refSeq);
        if (tagAlignment == null || tagAlignment.getSiteCount() < 1) {
            return null;
        }
        int nSites = tagAlignment.getSiteCount();
        int nTaxa = theTags.get(0).tagDist.length;  // the number of taxa is the same for all tags
        if (nTaxa < 1) {
            return null;
        }
        byte[][] callsBySite = new byte[nSites][nTaxa];
        final int nAlignedTags = tagAlignment.getSequenceCount();
        tagIndices = new int[nAlignedTags];  // the reference sequence is not included
        allelesAtVariableSitesByTag = new byte[nSites][theTags.size()];
        for (int tg = 0; tg < nAlignedTags; tg++) {
            int indexInTheTags = Integer.parseInt(tagAlignment.getTaxaName(tg)); // taxaName in tagAlignment is set to indexInTheTags
            tagIndices[tg] = indexInTheTags;
            for (int s = 0; s < nSites; s++) {
                allelesAtVariableSitesByTag[s][tagIndices[tg]] = tagAlignment.getBase(tg, s);
            }
        }
        positionsOfVariableSites = new int[nSites];
        for (int s = 0; s < nSites; s++) {
            positionsOfVariableSites[s] = tagAlignment.getPositionInLocus(s);
            for (int tx = 0; tx < nTaxa; tx++) {
                int[] alleleCounts = new int[Byte.MAX_VALUE];
                for (int tg = 0; tg < nAlignedTags; tg++) {
                    int tagIndex = tagIndices[tg];
                    byte baseToAdd = allelesAtVariableSitesByTag[s][tagIndex];
                    if (baseToAdd == Alignment.UNKNOWN_DIPLOID_ALLELE && callBiallelicSNPsWithGap && maxStartPosition == minStartPosition) {
                        baseToAdd = NucleotideAlignmentConstants.GAP_DIPLOID_ALLELE;
                    }
                    alleleCounts[baseToAdd] += theTags.get(tagIndex).tagDist[tx];
                }
                callsBySite[s][tx] = resolveQuantGeno(alleleCounts);
            }
        }
        return callsBySite;
    }
    
    private void setAlleleDepthsInTaxaForSite(int site, int[][] alleleDepthsInTaxa, byte[] commonAlleles) {
        for (int a = 0; a < commonAlleles.length; a++) {
            myCommonAlleles[a][site] = commonAlleles[a];
        }
        for (int a = 0; a < alleleDepthsInTaxa.length; a++) {
            for (int tx = 0; tx < alleleDepthsInTaxa[a].length; tx++) {
                if (alleleDepthsInTaxa[a][tx] > 32767) { // max value is 32767
                    alleleDepthsInTaxa[a][tx] = 32767;
                }
                myAlleleDepthsInTaxa[a][site][tx] = (short)alleleDepthsInTaxa[a][tx];
            }
        }
    }

    private void populateAllelesAtVariableSitesByTag(Alignment tagAlignment, int nSites, boolean includeReferenceTag, boolean callBiallelicSNPsWithGap) {
        int nAlignedTags = tagAlignment.getSequenceCount();
        tagIndices = new int[nAlignedTags];
        allelesAtVariableSitesByTag = new byte[nSites][theTags.size()];
        for (int tg = 0; tg < nAlignedTags; tg++) {
            tagIndices[tg] = Integer.parseInt(tagAlignment.getTaxaName(tg).split("_")[0]);  // taxaName in tagAlignment is set to indexInTheTags_"refTag"|"no"
            for (int s = 0; s < nSites; s++) {
                if (includeReferenceTag && tagIndices[tg] == theTags.size()-1) {
                    refCallsBySite[s] = tagAlignment.getBase(tg, s); // diploid byte for the reference allele/geno
                } else {
                    byte allele = tagAlignment.getBaseArray(tg, s)[0]; // tags only have one base so the 1st allele (index [0]) sufffices
                    if (callBiallelicSNPsWithGap && allele == Alignment.UNKNOWN_ALLELE) {
                        allele = NucleotideAlignmentConstants.GAP_ALLELE;
                    }
                    allelesAtVariableSitesByTag[s][tagIndices[tg]] = allele;
                }
            }
        }
    }

    public int[] getPositionsOfVariableSites() {
        return positionsOfVariableSites;
    }
    
    public String getLocusReport(int minTaxaWithLocus, boolean[] varSiteKept) {
        int start, end, totalbp, refTag=Integer.MIN_VALUE;
        if (strand == -1) {
            end = minStartPosition;
            start = minStartPosition - maxTagLength + 1;
        } else {
            start = minStartPosition;
            end = minStartPosition + maxTagLength - 1;
        }
        totalbp = end - start + 1;
        int nVarSites=0, nVarSitesKept=0;
        String posVarSites = "", posVarsKept = "";
        if (status.equals("polymorphic")) {
            nVarSites = positionsOfVariableSites.length;
            for (int s = 0; s < nVarSites; s++) {
                posVarSites = s < nVarSites-1 ? posVarSites + positionsOfVariableSites[s] + ":"
                                              : posVarSites + positionsOfVariableSites[s];
                if (varSiteKept[s]) {
                    posVarsKept = posVarsKept + positionsOfVariableSites[s] + ":";
                    nVarSitesKept++;
                }
            }
            if (posVarsKept.length()>0) posVarsKept = posVarsKept.substring(0, posVarsKept.length()-1);
            else posVarsKept = "NA";
        } else {
            posVarSites = "NA";
            posVarsKept = "NA";
            if (status.equals("invariant") || status.contains("tooManyTags")) assignRefTag();
        }
        refTag = (indexOfRef==Integer.MIN_VALUE) ? 0 : 1;
        return
            chromosome +"\t"+
            start +"\t"+
            end +"\t"+
            strand +"\t"+
            totalbp +"\t"+
            theTags.size() +"\t"+
            this.getTotalNReads() +"\t"+
            this.getNumberTaxaCovered() +"\t"+
            minTaxaWithLocus+"\t"+
            status +"\t"+ 
            nVarSites +"\t"+ 
            posVarSites +"\t"+
            nVarSitesKept  +"\t"+
            posVarsKept +"\t"+
            refTag +"\t"+
            maxTagLength +"\t"+
            minTagLength +"\n"
        ;
    }

    private Alignment getVariableSites() {
        if (theTags.size() < 2) {
            status = "invariant";
            return null;
        }
        if (theTags.size() > maxAlignmentSize) {
            status = "tooManyTags(>"+maxAlignmentSize+")";
            return null;
        }
        if (indexOfRef == Integer.MIN_VALUE) this.assignRefTag();
        boolean checkReplicateProfiles = false;
        boolean printOutAlignments = true;
        List<DNASequence> lst = new ArrayList<DNASequence>();
        int tagIndex = 0;
        for (SingleTagByTaxa sTBT : theTags) {
            DNASequence ds = new DNASequence(sTBT.tagTrimmed);
            String refMark = (tagIndex == indexOfRef) ? "refTag" : "no";
            ds.setOriginalHeader(tagIndex + "_" + refMark);    // OriginalHeader set to indexInTheTags_'refTag'|'no'
            ds.setCompoundSet(AmbiguityDNACompoundSet.getDNACompoundSet());
            lst.add(ds);
            ++tagIndex;
        }
        Profile<DNASequence, NucleotideCompound> profile = Alignments.getMultipleSequenceAlignment(lst);
        if (checkReplicateProfiles) {
            System.out.printf("Clustal1:%d%n%s%n", minStartPosition, profile);
            Profile<DNASequence, NucleotideCompound> profile2 = Alignments.getMultipleSequenceAlignment(lst);
            System.out.printf("Clustal2:%d%n%s%n", minStartPosition, profile2);
        }
        String[] aseqs = new String[theTags.size()];
        String[] names = new String[theTags.size()];
        boolean refTagWithGaps = false;
        int[] positions = null;
        for (int i = 0; i < aseqs.length; i++) {
            aseqs[i] = profile.getAlignedSequence(i + 1).getSequenceAsString();
            names[i] = profile.getAlignedSequence(i + 1).getOriginalSequence().getOriginalHeader();
            if (names[i].split("_")[1].equals("refTag")) {  // names were set to indexInTheTags_"refTag"|"no"
                if (aseqs[i].contains("-")) {
                    refTagWithGaps = true;
                    positions = new int[aseqs[i].length()];
                    positions[0] = 0;
                    for (int site = 1; site < aseqs[i].length(); site++) {
                        positions[site] = (aseqs[i].charAt(site) == '-') ? (positions[site - 1]) : (positions[site - 1] + 1);
                    }
                }
            }
        }
        profile = null;
        Alignment aa = null;
        if (refTagWithGaps) {
            aa = BitAlignment.getNucleotideInstance(new SimpleIdGroup(names), aseqs, null, null, positions, 5, new Locus[]{Locus.UNKNOWN}, new int[]{0}, null, false, true);
        } else {
            aa = BitAlignment.getNucleotideInstance(new SimpleIdGroup(names), aseqs, null, null, null, 5, new Locus[]{Locus.UNKNOWN}, new int[]{0}, null, false, true);
        }
        Alignment faa = AlignmentUtils.removeSitesBasedOnFreqIgnoreMissing(aa, 0.000001, 1.0, 2);
//        if (printOutAlignments && refTagWithGaps) {
        if (printOutAlignments && (minStartPosition % 1000 == 0)) {
            String tagStr;
            System.out.println("\nHere is an example alignment for a TagLocus (1 out of every 1000 is displayed):");
            System.out.println("chr" + chromosome + "  pos:" + minStartPosition + "  strand:" + strand + "  All sites:");
            for (int tg = 0; tg < aseqs.length; tg++) {
                tagStr = aa.getBaseAsStringRow(tg);
                tagStr = tagStr.replaceAll(";", "");
                System.out.println(tagStr + " " + names[tg]);
            }
            System.out.println("chr" + chromosome + "  pos:" + minStartPosition + "  strand:" + strand + "  Polymorphic sites only:");
            for (int tg = 0; tg < aseqs.length; tg++) {
                tagStr = faa.getBaseAsStringRow(tg);
                tagStr = tagStr.replaceAll(";", "");
                System.out.println(tagStr + " " + names[tg]);
            }
            System.out.println();
        }
        if (faa.getSiteCount() > maxSNPsPerLocus) {
            status = "tooManyVariants(>"+maxSNPsPerLocus+")";
            return null;
        }
        if (faa.getSiteCount() < 1) {
            status = "noVarSitesInAlign";
            return null;
        }
        if (faa.getSequenceCount() < 2) {
            status = "onlyOneTagInAlign";
            return null;
        }
        return faa;
    }

    private Alignment getVariableSites(String refSeq) {
        if (theTags.size() < 2) {
            return null;
        }
        boolean printOutAlignments = true;
        int startRefGenIndex = 0, endRefGenIndex = 1, startTagIndex = 2, endTagIndex = 3; // relevant indices in alignStats[]  (alignedTagLen=4)
        DNASequence dsRefSeq = new DNASequence(refSeq);
        dsRefSeq.setCompoundSet(AmbiguityDNACompoundSet.getDNACompoundSet());
        int minRefGenIndex = Integer.MAX_VALUE, maxRefGenIndex = Integer.MIN_VALUE;
        ArrayList<SequencePair<DNASequence, NucleotideCompound>> pairwiseAligns = new ArrayList<SequencePair<DNASequence, NucleotideCompound>>();
        int tagIndex = 0;
        for (SingleTagByTaxa sTBT : theTags) {
            DNASequence ds = new DNASequence(sTBT.tagTrimmed);
            ds.setCompoundSet(AmbiguityDNACompoundSet.getDNACompoundSet());
            SequencePair<DNASequence, NucleotideCompound> psa = Alignments.getPairwiseAlignment(ds, dsRefSeq, PairwiseSequenceAlignerType.LOCAL, gapPen, subMatrix);
            int[] alignStats = getAlignStats(psa, printOutAlignments, tagIndex, sTBT.tagLength, sTBT.tagStrand);
            minRefGenIndex = adjustMinRefGenIndex(minRefGenIndex, alignStats[startRefGenIndex], alignStats[startTagIndex]);
            maxRefGenIndex = adjustMaxRefGenIndex(maxRefGenIndex, alignStats[endRefGenIndex], alignStats[endTagIndex], sTBT.tagLength, refSeq);
            pairwiseAligns.add(psa);
            ++tagIndex;
        }
        minStartPosition += minRefGenIndex - 1;
        if (printOutAlignments && minStartPosition > 10000000 && minStartPosition < 10100000) {
            System.out.println("minRefGenIndex:" + minRefGenIndex + "  maxRefGenIndex:" + maxRefGenIndex + "  ChrPositionAtMinRefGenIndex:" + minStartPosition + "\n");
        }
        String[] aseqs = new String[theTags.size()];  // omit the reference genome sequence
        String[] names = new String[theTags.size()];
        char[][] myAlign = getAlignment(pairwiseAligns, refSeq, minRefGenIndex, maxRefGenIndex, aseqs, names);
        if (printOutAlignments && minStartPosition > 10000000 && minStartPosition < 10100000) {
            writeAlignment(refSeq, myAlign, minRefGenIndex, maxRefGenIndex);
        }
        Alignment a = null;
        a = BitAlignment.getNucleotideInstance(new SimpleIdGroup(names), aseqs, null, null, null, TasselPrefs.getAlignmentMaxAllelesToRetain(), new Locus[]{Locus.UNKNOWN}, new int[]{0}, null, TasselPrefs.getAlignmentRetainRareAlleles(), true);
        Alignment fa = AlignmentUtils.removeSitesBasedOnFreqIgnoreMissing(a, 0.000001, 1.0, 2);
        if (printOutAlignments && minStartPosition > 10000000 && minStartPosition < 10100000) {
            System.out.println("chr" + chromosome + "  pos:" + minStartPosition + "  strand:" + strand + "  FA (alignment filtered for polymorphic sites):\n" + fa.toString());
        }
        if (fa.getSiteCount() > maxSNPsPerLocus * 5 || fa.getSiteCount() < 1 || fa.getSequenceCount() < 2) {
            return null;
        }
        return fa;
    }

    private int[] getAlignStats(SequencePair<DNASequence, NucleotideCompound> psa, boolean printOutAlignments, int tagIndex, int tagLength, byte tagStrand) {
        int[] alignStats = new int[5];
        int startRefGenIndex = 0, endRefGenIndex = 1, startTagIndex = 2, endTagIndex = 3, alignedTagLen = 4; // indices in alignStats[]
        alignStats[startRefGenIndex] = psa.getIndexInTargetAt(1);
        alignStats[endRefGenIndex] = psa.getIndexInTargetAt(psa.getLength());
        alignStats[startTagIndex] = psa.getIndexInQueryAt(1);
        alignStats[endTagIndex] = psa.getIndexInQueryAt(psa.getLength());
        alignStats[alignedTagLen] = alignStats[endTagIndex] - alignStats[startTagIndex] + 1;
        if (printOutAlignments && minStartPosition > 10000000 && minStartPosition < 10100000) {
            System.out.println("tagIndex:" + tagIndex
                    + "  startRefGenIndex:" + alignStats[startRefGenIndex]
                    + "  endRefGenIndex:" + alignStats[endRefGenIndex]
                    + "  tagLength:" + tagLength
                    + "  startTagIndex:" + alignStats[startTagIndex]
                    + "  endTagIndex:" + alignStats[endTagIndex]
                    + "  alignedTagLen:" + alignStats[alignedTagLen]
                    + "  originalStrand: " + tagStrand + "\n"
                    + psa);
        }
        return alignStats;
    }

    private int adjustMinRefGenIndex(int minRefGenIndex, int startRefGenIndex, int startTagIndex) {
        if (startRefGenIndex < minRefGenIndex) {
            minRefGenIndex = startRefGenIndex;
        }
        if (startTagIndex > 1 && startTagIndex < 4 && startRefGenIndex - startTagIndex + 1 < minRefGenIndex && startRefGenIndex - startTagIndex + 1 > 0) {
            // extend regional alignment if there was soft clipping of 1 or 2 bases and there is sufficient 5' refSeq available
            minRefGenIndex = startRefGenIndex - startTagIndex + 1;
        }
        return minRefGenIndex;
    }

    private int adjustMaxRefGenIndex(int maxRefGenIndex, int endRefGenIndex, int endTagIndex, int tagLength, String refSeq) {
        if (endRefGenIndex > maxRefGenIndex) {
            maxRefGenIndex = endRefGenIndex;
        }
        if (endTagIndex < tagLength
                && tagLength - endTagIndex < 3
                && endRefGenIndex + tagLength - endTagIndex > maxRefGenIndex
                && endRefGenIndex + tagLength - endTagIndex <= refSeq.length()) {
            // extend regional alignment if there was soft clipping of 1 or 2 bases and there is sufficient 3' refSeq available
            maxRefGenIndex = endRefGenIndex + tagLength - endTagIndex;
        }
        return maxRefGenIndex;
    }

    private char[][] getAlignment(ArrayList<SequencePair<DNASequence, NucleotideCompound>> pairwiseAligns,
            String refSeq, int minRefGenIndex, int maxRefGenIndex, String[] aseqs, String[] names) {
        int totAlignedLen = maxRefGenIndex - minRefGenIndex + 1;
        char[][] myAlign = new char[theTags.size()][totAlignedLen];  // omit the reference genome sequence
        for (int t = 0; t < myAlign.length; t++) {
            Arrays.fill(myAlign[t], 'N');
        }
        int tagIndex = 0;
        for (SingleTagByTaxa sTBT : theTags) {
            SequencePair<DNASequence, NucleotideCompound> psa = pairwiseAligns.get(tagIndex);
            int tagStart = psa.getIndexInQueryAt(1);
            int tagEnd = psa.getIndexInQueryAt(psa.getLength());
            int refSeqStart = psa.getIndexInTargetAt(1);
            int refSeqEnd = psa.getIndexInTargetAt(psa.getLength());
            if (tagStart > 1 && tagStart < 4 && refSeqStart - minRefGenIndex - tagStart + 1 > -1) {
                // extend tag start if there was soft clipping of 1 or 2 bases and there is sufficient 5' refSeq available
                for (int offset = tagStart - 1; offset > 0; offset--) {
                    myAlign[tagIndex][refSeqStart - minRefGenIndex - offset] = sTBT.tagTrimmed.charAt(tagStart - offset - 1);
                }
            }
            for (int i = 1; i <= psa.getLength(); i++) {
                char refBase = psa.getCompoundInTargetAt(i).getBase().charAt(0);
                if (refBase != '-') {
                    myAlign[tagIndex][psa.getIndexInTargetAt(i) - minRefGenIndex] = psa.getCompoundInQueryAt(i).getBase().charAt(0);
                }
            }
            int extension = sTBT.tagLength - tagEnd;
            if (extension > 0 && extension < 3 && refSeqEnd - minRefGenIndex + extension < refSeq.length()) {
                // extend tag end if there was soft clipping of 1 or 2 bases and there is sufficient 3' refSeq available
                for (int offset = sTBT.tagLength - tagEnd; offset > 0; offset--) {
                    myAlign[tagIndex][refSeqEnd - minRefGenIndex + offset] = sTBT.tagTrimmed.charAt(sTBT.tagLength - offset);
                }
            }
            aseqs[tagIndex] = new String(myAlign[tagIndex]);
            names[tagIndex] = tagIndex + "";
            ++tagIndex;
        }
        return myAlign;
    }

    private void writeAlignment(String refSeq, char[][] myAlign, int minRefGenIndex, int maxRefGenIndex) {
        System.out.println("All tags in the region aligned to the reference sequence (first line) (insertions relative to the reference excluded):");
        System.out.println(refSeq.substring(minRefGenIndex - 1, maxRefGenIndex));
        for (int tagIndex = 0; tagIndex < myAlign.length; tagIndex++) {
            for (int b = 0; b < myAlign[tagIndex].length; b++) {
                System.out.print(myAlign[tagIndex][b]);
            }
            System.out.print("\n");
        }
        System.out.print("\n");
    }

    private byte[] getCommonAlleles(int s, int nTaxa, boolean includeReferenceTag) {
        int[] alleleCounts = new int[16];
        int nTags = includeReferenceTag ? theTags.size()-1 : theTags.size();
        for (int tg = 0; tg < nTags; tg++) {
            byte baseToAdd = allelesAtVariableSitesByTag[s][tg];
            for (int tx = 0; tx < nTaxa; tx++) {
                alleleCounts[baseToAdd] += theTags.get(tg).tagDist[tx];
            }
        }
        byte[] commonAlleles = new byte[maxNumAlleles];
        int[][] sortedAlleleCounts = sortAllelesByCount(alleleCounts);
        // Ties between maxNumAlleles - 1 and maxNumAlleles are not handled
        for (int i = 0; i < maxNumAlleles; i++) {
            commonAlleles[i] = (byte) sortedAlleleCounts[0][i];
        }
        return commonAlleles;
    }

    private int[][] getAlleleDepthsInTaxa(byte[] commonAlleles, int s, int nTaxa, boolean includeReferenceTag) {
        int[][] allelesInTaxa = new int[maxNumAlleles][nTaxa];
        int nTags = includeReferenceTag ? theTags.size()-1 : theTags.size(); // skip the reference tag (=last tag, if present), as it has no tagDist[]
        for (int tg = 0; tg < nTags; tg++) {
            byte baseToAdd = allelesAtVariableSitesByTag[s][tg];
            for (int a = 0; a < maxNumAlleles; a++) {
                if (baseToAdd == commonAlleles[a]) {
                    for (int tx = 0; tx < nTaxa; tx++) {
                        allelesInTaxa[a][tx] += theTags.get(tg).tagDist[tx];
                    }
                }
            }
        }
        return allelesInTaxa;
    }
    
    
    private byte resolveHetGeno(byte[] alleles, int[][] allelesInTaxa, int tx) {
        int max = 0;
        byte maxAllele = Alignment.UNKNOWN_ALLELE;
        int nextMax = 0;
        byte nextMaxAllele = Alignment.UNKNOWN_ALLELE;
        for (int a = 0; a < maxNumAlleles; a++) {
            if (allelesInTaxa[a][tx] > max) {
                nextMax = max;
                nextMaxAllele = maxAllele;
                max = allelesInTaxa[a][tx];
                maxAllele = alleles[a];
            } else if (allelesInTaxa[a][tx] > nextMax) {
                nextMax = allelesInTaxa[a][tx];
                nextMaxAllele = alleles[a];
            }
        }
        int totCount = max + nextMax;
        if (totCount < maxCountAtGeno) {
            if (nextMax < likelihoodRatioThreshAlleleCnt[totCount]) {
                return (byte) ((maxAllele << 4) | maxAllele); // call it a homozygote
            } else {
                return (byte) ((maxAllele << 4) | nextMaxAllele); // call it a het
            }
        } else {
            if (nextMax / totCount < 0.1) {
                return (byte) ((maxAllele << 4) | maxAllele); // call it a homozygote
            } else {
                return (byte) ((maxAllele << 4) | nextMaxAllele); // call it a het
            }
        }
    }

    private byte resolveQuantGeno(int[] alleleCounts) {
        int[][] sortedAlleleCounts = sortAllelesByCount(alleleCounts);
        int a1Count = sortedAlleleCounts[1][0];
        if (a1Count == 0) {
            return Alignment.UNKNOWN_DIPLOID_ALLELE;
        }
        int a2Count = sortedAlleleCounts[1][1];  // What if a3Count = a2Count? -- this situation is not dealt with
        byte a1 = (byte) sortedAlleleCounts[0][0];
        if (a2Count == 0) {
            return a1;
        }
        byte a2 = (byte) sortedAlleleCounts[0][1];
        int totCount = a1Count + a2Count;
        if (totCount < maxCountAtGeno) {
            if (a2Count < likelihoodRatioThreshAlleleCnt[totCount]) {
                return a1;  // call it a homozygote
            } else {
                return AlignmentUtils.getDiploidValue(a1, a2);  // call it a het
                //return IUPACNucleotides.getDegerateSNPByteFromTwoSNPs(a1, a2); // call it a het
            }
        } else {
            if (a2Count / totCount < 0.1) {
                return a1;  // call it a homozygote
            } else {
                return AlignmentUtils.getDiploidValue(a1, a2);  // call it a het
                //return IUPACNucleotides.getDegerateSNPByteFromTwoSNPs(a1, a2); // call it a het
            }
        }
    }

    private int[][] sortAllelesByCount(int[] alleleCounts) {
        // note that 'N' is not included as an allele
        byte[] alleles = {NucleotideAlignmentConstants.A_ALLELE, NucleotideAlignmentConstants.C_ALLELE,
            NucleotideAlignmentConstants.G_ALLELE, NucleotideAlignmentConstants.T_ALLELE, NucleotideAlignmentConstants.GAP_ALLELE};
        int[][] result = new int[2][alleles.length]; // result[0][a]=allele; result[1][a]=count
        for (int i = 0; i < alleles.length; i++) {
            result[0][i] = alleles[i];
            result[1][i] = alleleCounts[alleles[i]];
        }
        boolean change = true;
        while (change) { // sort the alleles by descending frequency
            change = false;
            for (int k = 0; k < alleles.length - 1; k++) {
                if (result[1][k] < result[1][k + 1]) {
                    int temp = result[0][k];
                    result[0][k] = result[0][k + 1];
                    result[0][k + 1] = temp;
                    int tempCount = result[1][k];
                    result[1][k] = result[1][k + 1];
                    result[1][k + 1] = tempCount;
                    change = true;
                }
            }
        }
        return result;
    }

    private String padTagWithNs(SingleTagByTaxa tag, String refSeq) {
        StringBuilder sb = new StringBuilder();
        char[] nullBases;
        if (tag.tagStrand == -1) {
            nullBases = new char[tag.startPosition - minStartPosition - tag.tagTrimmed.length() + 1];
        } else {
            nullBases = new char[tag.startPosition - minStartPosition];
        }
        Arrays.fill(nullBases, 'N');
        sb.append(nullBases);
        sb.append(tag.tagTrimmed);
        if (tag.tagStrand == -1) {
            nullBases = new char[refSeq.length() - sb.length()];
        } else {
            nullBases = new char[refSeq.length() - sb.length()];
        }

        Arrays.fill(nullBases, 'N');
        sb.append(nullBases);
        return sb.toString();
    }
}

class SingleTagByTaxa {

    int tagTOPMIndex;
    int tagLength;
    int startPosition;
    byte tagStrand;
    int divergence;
    String tagTrimmed;
    int tagTBTIndex; //index in the TBT
    int taxaWithTag;
    short[] tagDist;  // observed count of the tag for each taxon

    SingleTagByTaxa(int tagTOPMIndex, TagsOnPhysicalMapV3 theTOPM, TagsByTaxa theTBT, boolean includeRefGenome, boolean fuzzyStartPositions) {
        tagStrand = Byte.MIN_VALUE;
        this.tagTOPMIndex = tagTOPMIndex;
        long[] tag = theTOPM.getTag(tagTOPMIndex);
        tagTBTIndex = theTBT.getTagIndex(tag);
        taxaWithTag = (tagTBTIndex > -1) ? theTBT.getNumberOfTaxaWithTag(tagTBTIndex) : 0;
        if (taxaWithTag > 0) {  // tags with 0 taxaWithTag will not be added to TagsAtLocus
            startPosition = theTOPM.getStartPosition(tagTOPMIndex);
            tagLength = theTOPM.getTagLength(tagTOPMIndex);
            divergence = theTOPM.getDivergence(tagTOPMIndex);
            tagTrimmed = BaseEncoder.getSequenceFromLong(tag).substring(0, tagLength);
            tagStrand = theTOPM.getStrand(tagTOPMIndex);
            if (includeRefGenome && fuzzyStartPositions) {
                if (tagStrand == -1) {
                    tagTrimmed = BaseEncoder.getReverseComplement(tagTrimmed);
                }
            } else if (tagLength < theTOPM.getTagSizeInLong() * BaseEncoder.chunkSize) {
                tagTrimmed = tagTrimmed
                        + theTOPM.getNullTag().substring(0, theTOPM.getTagSizeInLong() * BaseEncoder.chunkSize - tagLength).replace("A", "N");
            }
            tagDist = theTBT.getTaxaReadCountsForTag(tagTBTIndex);
        }
    }
    
    SingleTagByTaxa(int tagTOPMIndex, TagsOnPhysicalMap theTOPM, TagsByTaxa theTBT, boolean includeRefGenome, boolean fuzzyStartPositions) {
        tagStrand = Byte.MIN_VALUE;
        this.tagTOPMIndex = tagTOPMIndex;
        long[] tag = theTOPM.getTag(tagTOPMIndex);
        tagTBTIndex = theTBT.getTagIndex(tag);
        taxaWithTag = (tagTBTIndex > -1) ? theTBT.getNumberOfTaxaWithTag(tagTBTIndex) : 0;
        if (taxaWithTag > 0) {  // tags with 0 taxaWithTag will not be added to TagsAtLocus
            startPosition = theTOPM.getStartPosition(tagTOPMIndex);
            tagLength = theTOPM.getTagLength(tagTOPMIndex);
            divergence = theTOPM.getDivergence(tagTOPMIndex);
            tagTrimmed = BaseEncoder.getSequenceFromLong(tag).substring(0, tagLength);
            tagStrand = theTOPM.getStrand(tagTOPMIndex);
            if (includeRefGenome && fuzzyStartPositions) {
                if (tagStrand == -1) {
                    tagTrimmed = BaseEncoder.getReverseComplement(tagTrimmed);
                }
            } else if (tagLength < theTOPM.getTagSizeInLong() * BaseEncoder.chunkSize) {
                tagTrimmed = tagTrimmed
                        + theTOPM.getNullTag().substring(0, theTOPM.getTagSizeInLong() * BaseEncoder.chunkSize - tagLength).replace("A", "N");
            }
            tagDist = theTBT.getTaxaReadCountsForTag(tagTBTIndex);
        }
    }
    
    // this constructor can be used to add a refTag directly from the refererence genome
    SingleTagByTaxa(int startPosition, byte strand, String refTag, int nLongsPerTag, String nullTag) {
        tagTOPMIndex = Integer.MIN_VALUE;
        tagLength = (byte) refTag.length();
        this.startPosition = startPosition;
        tagStrand = strand;
        divergence = 0;
        tagTrimmed = refTag;
        if (tagLength < nLongsPerTag*BaseEncoder.chunkSize) {
            tagTrimmed = tagTrimmed + 
                    nullTag.substring(0,nLongsPerTag*BaseEncoder.chunkSize-tagLength).replace("A","N");
        }
        tagTBTIndex = Integer.MIN_VALUE;
        taxaWithTag = 0;
        tagDist = null;
    }
}
