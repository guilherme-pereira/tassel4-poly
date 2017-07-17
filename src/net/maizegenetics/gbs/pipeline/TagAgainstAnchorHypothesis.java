/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.pipeline;

import cern.jet.random.Binomial;
import edu.cornell.lassp.houle.RngPack.RandomJava;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import net.maizegenetics.gbs.maps.TagGeneticMappingInfo;
import net.maizegenetics.gbs.maps.TagMappingInfo;
import net.maizegenetics.gbs.maps.TagMappingInfoV3;
import net.maizegenetics.gbs.maps.TagMappingInfoV3.Aligner;
import net.maizegenetics.gbs.maps.TagsOnPhysicalMapV3;
import net.maizegenetics.gbs.tagdist.TagsByTaxaByte;
import net.maizegenetics.gbs.tagdist.TagsByTaxaByteHDF5TagGroups;
import net.maizegenetics.gbs.util.BaseEncoder;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.BitAlignment;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.util.OpenBitSet;

/**
 *
 * @author Fei Lu
 */
public class TagAgainstAnchorHypothesis {
    SimpleGenotypeSBit anchor;
    double[] anchorMaf;
    int[] chromosomeNumber;
    int[] chrStartIndex;
    int[] chrEndIndex;
    TagsByTaxaByteHDF5TagGroups tbt;
    int[] tbtRedirect;
    TagsOnPhysicalMapV3 topm;
    Aligner blockAligner;
    double pThresh = 0.000001;
    int minCount = 20;
    int minTagAlleleIntersection = 4;
    int testSiteNum = 100;
    
    Task[] jobs;
    static int threadNumPerCore = 32; //32 can be faster, but the admin will complain about too many threads
    int threadNum;

    
    /**
     * Constructor to run genetic mapping
     * @param hapMapHDF5 anchor map
     * @param tbtHDF5 TagsByTaxa of HDF5 format
     * @param blockFileS TBTTagBlockFile, used to block the position where the tag is aligned while scanning
     * @param outfileS output file of genetic mapping
     * @param pThresh P-value threshold, default: 1e-6
     * @param minCount minimum count when tag appear in taxa (In TBT file), default = 20, too low number lacks statistical power
     * @param coreNum default:-1, which means using all cores in a node. When the coreNum is set less than total core number, which means using coreNum cores, each core runs 1 thread
     * @param chunkSize number of tags in a chunk. This determines the time usage in a node
     */
    public TagAgainstAnchorHypothesis (String hapMapHDF5, String tbtHDF5, String topmFileS, String software, int testSiteNum, double pThresh, int minCount, int coreNum) {
        this.testSiteNum = testSiteNum;
        this.pThresh = pThresh;
        this.minCount = minCount;
        this.loadAnchorMap(hapMapHDF5);
        this.loadTBT(tbtHDF5);
        this.loadTOPM(topmFileS);
        this.loadBlockAligner(software);
        this.calculateThreadNum(coreNum);
        this.MTMapping();
    }

    /**
     * Pre-calculate number of chunks when qsub genetic mapping 
     * @param tbtHDF5
     * @param coreNum
     * @param tagNumPerThreadInChunk
     * @return 
     */
    public static int getChunkNum (String tbtHDF5, int chunkSize) {
        TagsByTaxaByteHDF5TagGroups tbt = new TagsByTaxaByteHDF5TagGroups (tbtHDF5);
        int tagNum = tbt.getTagCount();
        int chunkNum = 0;
        int left = tagNum % chunkSize;
        if (left == 0) chunkNum = tagNum / chunkSize;
        else chunkNum = tagNum / chunkSize + 1;
        System.out.println("TBT has " + tagNum + " tags");
        System.out.println("TBT will be devided into " + chunkNum + " chunks, " + chunkSize + " tags each");
        if (left != 0) System.out.println("The last chunk has " + left + " tags");
        System.out.println("The index of chunk are used submit parallel computation to different node");
        System.out.println("Tags in each chunk will be multi-threaded in one node");
        return chunkNum;
    }
    
    /**
     * MT genetic mapping
     * @param outfileS 
     */
    public void MTMapping () {
        int chunkNum = topm.getChunkNum();
        TagGeneticMappingInfo[][] gmChunk = new TagGeneticMappingInfo[topm.getMappingNum()][topm.getChunkSize()];
        String[] dataSetNames = topm.creatTagGeneticMappingInfoDatasets(0, topm.getMappingNum());
        for (int i = 0; i < chunkNum; i++) {
            int chunkStartTagIndex = i * topm.getChunkSize();
            int chunkEndTagIndex = chunkStartTagIndex + topm.getChunkSize();
            if (chunkEndTagIndex > topm.getTagCount()) chunkEndTagIndex = topm.getTagCount();
            int actualChunkSize = chunkEndTagIndex - chunkStartTagIndex;
            int[] chunkTopm2TbtIndex = new int[actualChunkSize];
            TagMappingInfoV3[][] tmiChunk = topm.getMappingInfoChunk(chunkStartTagIndex);
            int[][] blockChrPos = new int[actualChunkSize][];
            for (int j = 0; j < actualChunkSize; j++) {
                int topmTagIndex = chunkStartTagIndex+j;
                blockChrPos[j] = topm.getUniqueMappingOfAligner(topmTagIndex, blockAligner);
                long[] t = topm.getTag(topmTagIndex);
                int hit = tbt.getTagIndex(t);
                if (hit < 0) chunkTopm2TbtIndex[j] = Integer.MIN_VALUE;
                else chunkTopm2TbtIndex[j] = hit;   
            }
            int[] threadStartTagIndex;
            int[] threadEndTagIndex;
            int[] threadSize;
            int left = actualChunkSize % threadNum;
            int baseSize = actualChunkSize / threadNum;
            if (baseSize == 0) {
                threadStartTagIndex = new int[left];
                threadEndTagIndex = new int[left];
                for (int j = 0; j < threadStartTagIndex.length; j++) {
                    threadStartTagIndex[j] = chunkStartTagIndex + j;
                    threadEndTagIndex[j] = threadStartTagIndex[j]+1;
                }
            }
            else {
                threadSize = new int[threadNum];
                threadStartTagIndex = new int[threadNum];
                threadEndTagIndex = new int[threadNum];
                for (int j = 0; j < left; j++) {
                    threadSize[j] = baseSize + 1;
                }
                for (int j = left; j < threadNum; j++) {
                    threadSize[j] = baseSize;
                }
                threadStartTagIndex[0] = chunkStartTagIndex;
                threadEndTagIndex[0] = threadStartTagIndex[0] + threadSize[0];
                for (int j = 1; j < threadNum; j++) {
                    threadStartTagIndex[j] = threadEndTagIndex[j-1];
                    threadEndTagIndex[j] = threadStartTagIndex[j]+threadSize[j];
                }
            } 
            int actualThreadNum = threadStartTagIndex.length;                
            jobs = new Task[actualThreadNum];
            Thread[] mts = new Thread[actualThreadNum];
            long lastTimePoint = this.getCurrentTimeNano();
            for (int j = 0; j < actualThreadNum; j++) {
                jobs[j] = new Task (threadStartTagIndex[j], threadEndTagIndex[j], chunkStartTagIndex, chunkTopm2TbtIndex, tmiChunk, blockChrPos);
            }
            System.out.println("Loading this chunk to multiple threads took " + this.getTimeSpanSecond(lastTimePoint) + " seconds");
            System.out.println("Multiple threading mapping in progress...");
            lastTimePoint = this.getCurrentTimeNano();
            for (int j = 0; j < actualThreadNum; j++) {
                mts[j] = new Thread(jobs[j]);
                mts[j].start();
            }
            for (int j = 0; j < actualThreadNum; j++) {
                try {
                    mts[j].join();
                }
                catch (Exception e) {
                    System.out.println(e.toString());
                }
            }
            System.out.println((chunkStartTagIndex/topm.getChunkSize()+1) + " chunks are mapped. " + topm.getChunkNum() + " chunks in total");
            System.out.println("Each LD compirison took " + (double)this.getTimeSpanNano(lastTimePoint)/actualChunkSize/testSiteNum/topm.getMappingNum() + " nano seconds");
            System.out.println("Multiple threading mapping took " + this.getTimeSpanSecond(lastTimePoint) + " seconds");
            int cnt = 0;
            for (int j = 0; j < jobs.length; j++) {
                TagGeneticMappingInfo[][] sub = jobs[j].getResult();
                for (int k = 0; k < sub.length; k++) {
                    for (int u = 0; u < sub[0].length; u++) {
                        gmChunk[u][cnt+k] = sub[k][u];
                    }
                }
                cnt+=sub.length;
            }
            if (cnt < topm.getChunkSize()) {
                for (int j = cnt; j < topm.getChunkSize(); j++) {
                    for (int k = 0; k < topm.getMappingNum(); k++) {
                        gmChunk[k][j] = new TagGeneticMappingInfo();
                    }
                }
            }
            topm.writeTagGeneticMappingInfoDataSets(dataSetNames, gmChunk, i);
            topm.setIfHasGeneticMapping(true);
            System.out.println("Mapping result from chunk " + i + "(Index) was written\n");
            System.gc();
        }     
    }
    
    /**
     * Calculate the thread number used in this node
     * @param coreNum 
     */
    private void calculateThreadNum (int coreNum) {
        int numOfProcessors = Runtime.getRuntime().availableProcessors();
        if (coreNum < 0) {
            threadNum = numOfProcessors * threadNumPerCore;
        }
        else {
            if (coreNum == 0) {
                System.out.println("Core number = 0, This runs at least on 1 thread. Quit.");
                System.exit(0);
            }
            threadNumPerCore = 1;
            threadNum = coreNum * threadNumPerCore;
            System.out.println("TBT will be mapped by " + threadNum + " tasks");
            System.out.println("Each core runs 1 tasks, or 1 threads");
        }
        int acutualUseCoreNum = numOfProcessors;
        if (acutualUseCoreNum > threadNum) acutualUseCoreNum = threadNum;
        System.out.println("This node has " + numOfProcessors + " processors. Will use " + acutualUseCoreNum + " processors");
        System.out.println("TOPM will be mapped by " + threadNum + " threads");
        System.out.println("Each core runs " + threadNumPerCore + " threads");
        System.out.println("Each TOPM chunk has " + topm.getChunkSize() + " tags, which will be split and mapped by the " + threadNum + " threads");
        System.out.println("Each thread will map " + topm.getChunkSize()/threadNum + " tags");
        System.out.println("For each tags, " + topm.getMappingNum() + " hypothesis will be tested at " + testSiteNum+ " adjacent sites");
        int left = topm.getTagCount()%topm.getChunkSize();
        if (left != 0) {
            System.out.println("The last TOPM chunk has " + left + " tags");
        }
        System.out.println("");
	}
    
    /**
     * Class for conducting mapping on one thread
     */
    class Task implements Runnable {
        int tagStartIndex;
        int chunkStartTagIndex;
        TagsByTaxaByte subTBT =null;
        TagMappingInfoV3[][] tmiChunk = null;
        int[][] blockChrPos = null; 
        int[] chunkTopm2TbtIndex;
        int subTBTSize;
		TagGeneticMappingInfo[][] gmResult = null;
        HashMap<Integer, Integer> subtbt2ChunkTopmHash;
		Task (int tagStartIndex, int tagEndIndex, int chunkStartTagIndex, int[] subTopm2TbtIndex, TagMappingInfoV3[][] tmiChunk, int[][] blockChrPos) {
            this.chunkStartTagIndex = chunkStartTagIndex;
            this.chunkTopm2TbtIndex = subTopm2TbtIndex;
            this.tmiChunk = tmiChunk;
            this.blockChrPos = blockChrPos;
            this.buildHash(tagStartIndex, tagEndIndex);
			this.buildSubTBT(tagStartIndex, tagEndIndex);
		}
        
        private void buildHash (int tagStartIndex, int tagEndIndex) { 
            subtbt2ChunkTopmHash = new HashMap();
            int cnt = 0;
            int offset = tagStartIndex - chunkStartTagIndex;
            for (int i = 0; i < tagEndIndex-tagStartIndex; i++) {
                if (chunkTopm2TbtIndex[i+offset] < 0) continue;
                subtbt2ChunkTopmHash.put(cnt, i+offset);
                cnt++;
            }
            subTBTSize = cnt;
        }
        
        /**
         * Build TBTByte for each task/thread
         * @param tagStartIndex
         * @param tagEndIndex 
         */
        private void buildSubTBT (int tagStartIndex, int tagEndIndex) {
            this.tagStartIndex = tagStartIndex;
            int subTopmSize = tagEndIndex-tagStartIndex;
            this.populateResult(subTopmSize);
            if (subTBTSize == 0) return;
            
            long[][] tags = new long[tbt.getTagSizeInLong()][subTBTSize];
            byte[] tagLength = new byte[subTBTSize];
            byte[][] tagDist = new byte[tbt.getTaxaCount()][subTBTSize];
            String[] namesForTaxa = tbt.getTaxaNames();
            int tbtTagIndex;
            for (int i = 0; i < subTBTSize; i++) {
                tbtTagIndex = chunkTopm2TbtIndex[subtbt2ChunkTopmHash.get(i)];
                long[] t = tbt.getTag(tbtTagIndex);
                for (int j = 0; j < tags.length; j++) {
                    tags[j][i] = t[j];
                }
                tagLength[i] = (byte)tbt.getTagLength(tbtTagIndex);
                for (int j = 0; j < tbt.getTaxaCount(); j++) {
                    tagDist[j][i] = (byte)tbt.getReadCountForTagTaxon(tbtTagIndex, j);
                }
                             
            }
            subTBT = new TagsByTaxaByte (tags, tagLength, tagDist, namesForTaxa);
        }
        
        private void populateResult (int subTopmSize) {
            gmResult = new TagGeneticMappingInfo[subTopmSize][topm.getMappingNum()];
            for (int i = 0; i < gmResult.length; i++) {
                for (int j = 0; j < topm.getMappingNum(); j++) {
                    gmResult[i][j] = new TagGeneticMappingInfo();
                }
            }
        }
        
        /**
         * Return gmResult for output when all the tasks are done
         * @return 
         */
		public TagGeneticMappingInfo[][] getResult () {
			return gmResult;
		}

        @Override
		public void run() {
            if (subTBT == null) return;
            long lastTimePoint = getCurrentTimeNano();
			ArrayList<String> resultList = new ArrayList();
            int chunkTopmTagIndex;
            int[] chrPos;
            int blockChr = Integer.MIN_VALUE;
            int blockPos = Integer.MIN_VALUE;
            int blastPos = Integer.MIN_VALUE;
            TagMappingInfoV3 tmi;
            int bestAlignment = Integer.MIN_VALUE;
            int refDiv = Integer.MIN_VALUE;
            double[][] theResults=new double[1][];
            int[] siteIndexRange;
            ScanChromosome sc;
            int offset = tagStartIndex - chunkStartTagIndex;
            for (int i = 0; i < subTBT.getTagCount(); i++) {
                if(subTBT.getNumberOfTaxaWithTag(i) <minCount) continue;
                chunkTopmTagIndex = subtbt2ChunkTopmHash.get(i);
                chrPos = blockChrPos[chunkTopmTagIndex];
                if (chrPos == null) {
                    blockChr = Integer.MIN_VALUE;
                    blockPos = Integer.MIN_VALUE;
                }
                else {
                    blockChr = chrPos[0];
                    blockPos = chrPos[1];
                }        
                long[] testTag = subTBT.getTag(i);
                long[] testTagDist = getTagsInBits(subTBT, i, tbtRedirect, anchor.getTaxaNum());
                
                
                for (int j = 0; j < topm.getMappingNum(); j++) {
                    tmi = tmiChunk[j][chunkTopmTagIndex];
                    if (tmi.chromosome < 0) continue;
                    int siteIndex = anchor.getSiteIndex(tmi.chromosome, tmi.startPosition);
                    if (siteIndex == Integer.MIN_VALUE) continue;
                    siteIndexRange = anchor.getAdjacentSiteIndexRange(siteIndex, testSiteNum);
                    if (tmi.chromosome == blockChr) blastPos = blockPos; 
                    else blastPos = Integer.MIN_VALUE;
                    sc = new ScanChromosome(testTagDist, theResults, 0, siteIndexRange[0], siteIndexRange[1], pThresh, 1, blastPos);
                    sc.scan();
                    TagGeneticMappingInfo tgmi = new TagGeneticMappingInfo(theResults[0][3], (int)theResults[0][0], (int)theResults[0][2], (int)theResults[0][4], sc.maxSigPos-sc.minSigPos);
                    gmResult[this.subtbt2ChunkTopmHash.get(i)-offset][j] = tgmi;
                }
            }
		}
	}
    
    /**
     * Class used to scan SNPs in one chromosome
     */
    private class ScanChromosome {
        int chrIndex;
        int chrStartIndex;
        int chrEndIndex;
        OpenBitSet obsTdist;
        Binomial binomFunc=new Binomial(5, 0.5, new RandomJava());
        double[][] resultReport;
        double sigThreshold;
        int minSigPos=Integer.MIN_VALUE, maxSigPos=Integer.MIN_VALUE;
        int step=1;
        double bionomialThreshold=0.2;
        int blockPosition=Integer.MIN_VALUE;  //position to block from testing.
        int blockWindow=64;
        
        /**
         * 
         * @param tdist tag distribution if bits across all taxa
         * @param resultReport
         * @param chrIndex 
         * @param chrStartIndex SNP start index of this chromosome
         * @param chrEndIndex SNP end index of this chromosome
         * @param sigThreshold threadhold of P-value, default 1E-6
         * @param step step of scanning
         * @param blockPosition 
         */
        public ScanChromosome(long[] tdist, double[][] resultReport, int chrIndex, int chrStartIndex, int chrEndIndex, double sigThreshold, int step, int blockPosition) {
            obsTdist=new OpenBitSet(tdist,tdist.length);
            this.resultReport=resultReport;
            this.chrIndex = chrIndex;
            this.chrStartIndex = chrStartIndex;
            this.chrEndIndex = chrEndIndex;
            this.sigThreshold=sigThreshold;
            this.blockPosition=blockPosition;
            this.step=step;
        }

        public void scan() {
            long tests=0;
            int bestSite=-1, countSig=0;
            double bestP=2;
            for (int i = chrStartIndex; i < chrEndIndex; i+=step) {
                if(Math.abs(anchor.getPosition(i) -blockPosition)<blockWindow) continue;
                OpenBitSet obsMajor = anchor.obsMajor[i];
                OpenBitSet obsMinor = anchor.obsMinor[i];
                if(obsMinor.cardinality()>4) {
                    double p=fastTestSites(obsTdist, obsMajor, obsMinor, anchorMaf[i], binomFunc);
                    if(p<bestP) {bestP=p; bestSite=i;}
                    if(p<sigThreshold) {
                        countSig++;
                        if(minSigPos==Integer.MIN_VALUE) minSigPos=anchor.getPosition(i);
                        maxSigPos=anchor.getPosition(i);
                    }
                }
                tests++;
            }
            int chr=anchor.getChromosomeNumber(bestSite);
            double[] result={chr, bestSite, anchor.getPosition(bestSite),bestP, countSig};
            resultReport[chrIndex]=result;
        }
    }
    
    /**
     * fast test, using ratio to reduce the calculation of binomial and p (reduce 73%), but the time is only saved by 3%
     * @param obsTdist
     * @param obsMajor
     * @param obsMinor
     * @param maf
     * @param binomFunc
     * @return 
     */
    public double fastTestSites(OpenBitSet obsTdist, OpenBitSet obsMajor, OpenBitSet obsMinor, double maf, Binomial binomFunc) {
        double result=1;
        int tagMinorCount=0, tagMajorCount=0;
        tagMinorCount=(int)OpenBitSet.intersectionCount(obsTdist, obsMinor);
        tagMajorCount=(int)OpenBitSet.intersectionCount(obsTdist,obsMajor);

        int sumTagAllele = tagMinorCount + tagMajorCount;
        if(sumTagAllele<4) return result;
        
        double ratio;
        double minorProb;
        int cdfCount;
        ratio = (double)tagMinorCount/sumTagAllele;
        if (tagMinorCount < tagMajorCount) {
            minorProb = maf;
            cdfCount = tagMinorCount;
        }
        else {
            minorProb = 1-maf;
            cdfCount = tagMajorCount;
        }
        if ((ratio-minorProb) < -0.003) return result;
        binomFunc.setNandP(sumTagAllele,minorProb);
        try {
            result = binomFunc.cdf(cdfCount);
            //System.out.println(ratio+"\t"+minorProb+"\t"+gmResult);
  
        } catch (Exception e) {
            System.err.println("Error in the BinomialDistributionImpl");
        }
        return result;
    }
    
    /**
     * Test association using binomial test, collecting P-value
     * @param obsTdist
     * @param obsMajor
     * @param obsMinor
     * @param binomFunc
     * @return 
     */
    public double testSites(OpenBitSet obsTdist, OpenBitSet obsMajor, OpenBitSet obsMinor, double maf, Binomial binomFunc) {
        double result=1;
        int tagMinorCount=0, tagMajorCount=0;
        tagMinorCount=(int)OpenBitSet.intersectionCount(obsTdist, obsMinor);
        tagMajorCount=(int)OpenBitSet.intersectionCount(obsTdist,obsMajor);

        int sumTagAllele = tagMinorCount + tagMajorCount;
        if(sumTagAllele<4) return result;
    
        double minorProb = (tagMinorCount<tagMajorCount)? maf:1-maf;
        binomFunc.setNandP(sumTagAllele,minorProb);
        try {
            result = (tagMinorCount<tagMajorCount)?binomFunc.cdf(tagMinorCount):binomFunc.cdf(tagMajorCount);
  
        } catch (Exception e) {
            System.err.println("Error in the BinomialDistributionImpl");
        }
        return result;
    }
    
    /**
     * Test association using binomial test, collecting P-value
     * @param obsTdist
     * @param obsMajor
     * @param obsMinor
     * @param binomFunc
     * @return 
     */
    public static double testSites(OpenBitSet obsTdist, OpenBitSet obsMajor, OpenBitSet obsMinor, Binomial binomFunc) {
        double result=1;
        int minorCount=0, tagMinorCount=0, majorCount=0, tagMajorCount=0;
        tagMinorCount=(int)OpenBitSet.intersectionCount(obsTdist, obsMinor);
        tagMajorCount=(int)OpenBitSet.intersectionCount(obsTdist,obsMajor);
        minorCount=(int)obsMinor.cardinality();
        majorCount=(int)obsMajor.cardinality();

        int sumAllele = minorCount + majorCount;
        int sumTagAllele = tagMinorCount + tagMajorCount;
        if(sumTagAllele<4) return result;
    
        double minorProb = (tagMinorCount<tagMajorCount)?(double)minorCount/(double)sumAllele:(double)majorCount/(double)sumAllele;
        binomFunc.setNandP(sumTagAllele,minorProb);
        try {
            result = (tagMinorCount<tagMajorCount)?binomFunc.cdf(tagMinorCount):binomFunc.cdf(tagMajorCount);
  
        } catch (Exception e) {
            System.err.println("Error in the BinomialDistributionImpl");
        }
        return result;
    }
    
    /**
     * Return bit version of tag distribution across taxa
     * @param aTBT
     * @param tbtTagIndex
     * @param reDirect
     * @param anchorTaxa
     * @return 
     */
    private long[] getTagsInBits(TagsByTaxaByte aTBT, int tagIndex, int[] reDirect, int anchorTaxa) {
        int lgPerSite = (anchorTaxa / 64) + 1;
        long[] seq = new long[lgPerSite];
        for (int j = 0; j < aTBT.getTaxaCount(); j++) {
            if(reDirect[j]<0) continue;
            int index=reDirect[j]/64;
            int offset=reDirect[j]%64;
            if (aTBT.getReadCountForTagTaxon(tagIndex, j)>0) {  //reference alleles
                seq[index]=seq[index]|(1L<<offset);
            }
        }
        return seq;
	}
    
    private void redirect () {
        long lastTimePoint = this.getCurrentTimeNano();
        tbtRedirect = new int[tbt.getTaxaCount()];
        for (int i = 0; i < tbtRedirect.length; i++) {
            tbtRedirect[i] = anchor.getTaxonIndex(tbt.getTaxaName(i));
        }
        System.out.println("Taxa redirection took " + String.valueOf(this.getTimeSpanSecond(lastTimePoint)) + " seconds\n");
    }
    
    private void loadBlockAligner (String software) {
        blockAligner = Aligner.getAlignerFromName(software);
        if (blockAligner == null) {
            System.out.println("Please input correct aligner name, currently sopport Bowtie2, BWA and Blast");
            System.exit(1);
        }
    }
    
    private void loadTOPM (String topmFileS) {
        long lastTimePoint = this.getCurrentTimeNano();
        topm = new TagsOnPhysicalMapV3(topmFileS);
        System.out.println("Loading TOPM HDF5 took " + String.valueOf(this.getTimeSpanSecond(lastTimePoint)) + " seconds\n");
    }
    
    /**
     * load HDF5 TBT
     * @param tbtHDF5 
     */
    private void loadTBT (String tbtHDF5) {
        long lastTimePoint = this.getCurrentTimeNano();
        tbt = new TagsByTaxaByteHDF5TagGroups (tbtHDF5);
        System.out.println("Loading TBT HDF5 took " + String.valueOf(this.getTimeSpanSecond(lastTimePoint)) + " seconds");
        System.out.println("TBT has " + tbt.getTagCount() + " tags and " + tbt.getTaxaCount() + " taxa\n");
        this.redirect();
    }
    
    /**
     * load HDF5 genotype and convert to BitAlignment
     * @param hapMapHDF5 
     */
    private void loadAnchorMap (String hapMapHDF5) {
        System.out.println("Start loading anchor map");
        long lastTimePoint = this.getCurrentTimeNano();
        anchor = new SimpleGenotypeSBit(hapMapHDF5);
        System.out.println("Loading hapmap (SimpleGenotypeSBit) HDF5 took " + String.valueOf(this.getTimeSpanSecond(lastTimePoint)) + " seconds");
        System.out.println("The anchor map has " + anchor.getSiteNum() + " sites and " + anchor.getTaxaNum() + " taxa");
        chromosomeNumber = anchor.chromosomeNumber;
        chrStartIndex = anchor.chrStartIndex;
        chrEndIndex = anchor.chrEndIndex;
        anchorMaf = anchor.maf;
        System.gc();
        this.screenPrintGbMemoryCurrentUse();
        this.screenPrintGbMemoryAvailable();
        System.out.println();
    }
    
    private double getTimeSpanSecond (long lastTimePoint) {
        return (double)this.getTimeSpanNano(lastTimePoint)/1000000000;
    }
    
    private long getTimeSpanNano (long lastTimePoint) {
        return this.getCurrentTimeNano()- lastTimePoint;
    }
    
    private long getCurrentTimeNano () {
        return System.nanoTime();
    }
    
    private String getCurrentTimeHR() {
		SimpleDateFormat sd = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
        Date date = new Date();
        String str = sd.format(date);
        return str;
	}
    
    private double getGbMemoryAvailable () {
        return (double)(Runtime.getRuntime().maxMemory()- (Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory()))/1024/1024/1024;
    }
    
    private double getGbMemoryCurrentUse () {
        return (double)(Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory())/1024/1024/1024;
    }
    
    private void screenPrintTimeSpanSecond (long lastTimePoint) {
        System.out.println("Time span is " + String.valueOf(this.getTimeSpanSecond(lastTimePoint)) + " seconds");
    }
    
    private void screenPrintTimeSpanNano (long lastTimePoint) {
        System.out.println("Time span is " + String.valueOf(this.getTimeSpanNano(lastTimePoint)) + " ns");
    }
    
    private void screenPrintCurrentTimeHR () {
        System.out.println("Current time is " + this.getCurrentTimeHR());
    }
    
    private void screenPrintGbMemoryAvailable () {
        System.out.println("Available memory is " + String.valueOf(this.getGbMemoryAvailable()) + " GB");
    }
    
    private void screenPrintGbMemoryCurrentUse () {
        System.out.println("Current memory in use is " + String.valueOf(this.getGbMemoryCurrentUse()) + " GB");
    }
}
