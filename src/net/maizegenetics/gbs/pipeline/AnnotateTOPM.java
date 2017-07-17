/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.pipeline;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.TreeSet;
import net.maizegenetics.gbs.maps.PETagsOnPhysicalMapV3;
import net.maizegenetics.gbs.maps.TagGeneticMappingInfo;
import net.maizegenetics.gbs.maps.TagMappingInfoV3;
import net.maizegenetics.gbs.maps.TagMappingInfoV3.Aligner;
import net.maizegenetics.gbs.maps.TagsOnGeneticMap;
import net.maizegenetics.gbs.maps.TagsOnPhysicalMapV3;
import net.maizegenetics.gbs.tagdist.TagCounts;
import net.maizegenetics.gbs.tagdist.TagsByTaxa;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.gbs.util.BaseEncoder;
import net.maizegenetics.gbs.util.SAMUtils;
import net.maizegenetics.util.MultiMemberGZIPInputStream;

/**
 * Methods to annotate TOPM file, including adding mapping info from aligners, adding PE tag position and genetic position, model prediction for the best position
 * Code of mappingSource. 0: Bowtie2; 1: BWA; 2: BLAST; 3: BWAMEM; 4: PE one end; 5: PE the other end;
 * @author Fei Lu
 */
public class AnnotateTOPM {
    /**TOPM file that will be annotated*/
    TagsOnPhysicalMapV3 topm;
    /**Record Sam record. When tmiBuffers[0] is full, output to TOPM block. Substitute tmiBuffers[i] with tmiBuffers[i+1]. Size = numBuffers×maxMappingNum(-K option)×TOPM CHUNK_SIZE
     * Multiple buffers are used because the output in SAM is not exactly the order of input tag, roughly in the same order though
     */
    TagMappingInfoV3[][][] tmiBuffers = null;
    /**check if tmiBuffer[0] is full for output*/
    boolean[][] bufferLights = null;
    /**Number of buffers*/
    int bufferNum = 2;
    /**record how many mapping tags are imported in each buffer*/
    int[] lightCounts;
    /**Actual tag index of each tmiBuffers[i]*/
    int[] bufferStartTagIndex = null;
    /**Tag index range of the tmiBuffers*/
    int[] bufferTagIndexRange = null;
    /**When the second buffer has filled at least with this cutoff, update buffer*/
    int updateBufferCountCutoff;
    
    public static enum EvidenceType {
        GM((byte)1, 0),
        PE((byte)(1<<1), 1),
        SingleBestBowtie2((byte)(1<<2), 2),
        BestBWA((byte)(1<<3), 3),
        SingleBestBlast((byte)(1<<4), 4),
        SingleBestBWAMEM((byte)(1<<5), 5);
        private byte typeCode;
        private int index;
        EvidenceType (byte value, int indexInBestEvidence) {
            typeCode = value;
            index = indexInBestEvidence;
        }
        
        public byte getCode () {
            return typeCode;
        }
        
        public int getIndex () {
            return index;
        }
        
        public boolean getIfHasEvidence (EvidenceType type, boolean[] ifEvidence) {
            return ifEvidence[type.getIndex()];
        }
        
        public static byte code (boolean[] ifEvidence) {
            int evidenceByte = 0;
            for (int i = 0; i < ifEvidence.length; i++) {
                evidenceByte = evidenceByte << 1;
                if (ifEvidence[i] == true) evidenceByte = evidenceByte|1;              
            }
            return (byte)evidenceByte;
        }
        
        public static boolean[] deCode (byte evidenceByte) {
            boolean[] ifEvidence = new boolean[EvidenceType.getSize()];
            byte test = 1;
            for (int i = 0; i < ifEvidence.length; i++) {
                int value = (evidenceByte>>i)&test;
                if (value == 0) ifEvidence[ifEvidence.length-i-1] = false;
                else ifEvidence[ifEvidence.length-i-1] = true;
            }
            return ifEvidence;
        }
        
        public static int getSize () {
            return EvidenceType.values().length;
        }
        
    }
    
    /**
     * Constructor from a TOPM file
     * @param topm 
     */
    public AnnotateTOPM (TagsOnPhysicalMapV3 topm) {
        this.topm = topm;   
    }
    
    /**
     * Import prediction result to TOPM, contain hard coding on map index
     * @param indexDirS
     * @param predictDirS
     * @param priorityAligner 
     */
    public void annotateBestMappingImport (String indexDirS, String predictDirS, Aligner priorityAligner) {
        byte[] bestStrand = new byte[topm.getTagCount()];
        int[] bestChr = new int[topm.getTagCount()];
        int[] bestStartPos = new int[topm.getTagCount()];
        int[] bestEndPos = new int[topm.getTagCount()];
        byte[] bestDivergence = new byte[topm.getTagCount()];
        byte[] bestMapP = new byte[topm.getTagCount()];
        byte[] bestDcoP = new byte[topm.getTagCount()];
        byte[] multimaps = new byte[topm.getTagCount()];
        byte[] bestEvidence = new byte[topm.getTagCount()];
        byte[] bestMapIndex = new byte[topm.getTagCount()];
        int fileNum = new File(predictDirS).listFiles().length;
        for (int i = 0; i < fileNum; i++) {
            String predictFileS = String.valueOf(i)+".out.txt";
            predictFileS = new File(predictDirS, predictFileS).getAbsolutePath();
            String indexFileS = String.valueOf(i)+".index.txt";
            indexFileS = new File(indexDirS, indexFileS).getAbsolutePath();
            int currentTagIndex = -1;
            try {
                BufferedReader brp = new BufferedReader (new FileReader (predictFileS), 65536);
                BufferedReader bri = new BufferedReader (new FileReader (indexFileS), 65536);
                for (int j = 0; j < 5; j++) brp.readLine();
                bri.readLine();
                String temp;
                int tagIndex, mapIndex;
                String[] tem;
                ArrayList<Integer> mapIndexList = new ArrayList();
                Integer[] mapIndices;
                while ((temp = bri.readLine()) != null) {
                    tem = temp.split("\t");
                    tagIndex = Integer.parseInt(tem[0]);
                    mapIndex = Integer.parseInt(tem[1]);
                    if (tagIndex != currentTagIndex) {
                        if (currentTagIndex != -1) {
                            mapIndices = mapIndexList.toArray(new Integer[mapIndexList.size()]);
                            this.incorperateBestGeneticMapping(currentTagIndex, mapIndices, brp, bestStrand, bestChr, bestStartPos, bestEndPos, bestDivergence, bestMapP, bestDcoP, bestEvidence, bestMapIndex);
                        }
                        mapIndexList = new ArrayList();
                        currentTagIndex = tagIndex;
                    }
                    mapIndexList.add(mapIndex);
                }
                mapIndices = mapIndexList.toArray(new Integer[mapIndexList.size()]);
                this.incorperateBestGeneticMapping(currentTagIndex, mapIndices, brp, bestStrand, bestChr, bestStartPos, bestEndPos, bestDivergence, bestMapP, bestDcoP, bestEvidence, bestMapIndex);
            }
            catch (Exception e) {
                System.out.println(e.toString());
                System.out.println(currentTagIndex);
                System.exit(0);
            }
        }
        int[] priorityIndex = topm.getMappingIndicesOfAligner(priorityAligner);
        int[] alignerStartMapIndex = {0,5,10,15,20,25};
        int[][] alignerIndex = new int[alignerStartMapIndex.length][];
        alignerIndex[0] = topm.getMappingIndicesOfAligner(Aligner.Bowtie2);
        alignerIndex[1] = topm.getMappingIndicesOfAligner(Aligner.BWA);
        alignerIndex[2] = topm.getMappingIndicesOfAligner(Aligner.Blast);
        alignerIndex[3] = topm.getMappingIndicesOfAligner(Aligner.BWAMEM);
        alignerIndex[4] = topm.getMappingIndicesOfAligner(Aligner.PEEnd1);
        alignerIndex[5] = topm.getMappingIndicesOfAligner(Aligner.PEEnd2);
        for (int i = 0; i < topm.getTagCount(); i++) {
            if (bestStrand[i] == 0) {
                int numRank0 = this.getNumOfRank0(i, priorityIndex);
                TagMappingInfoV3 tmi = topm.getMappingInfo(i, priorityIndex[0]);
                if (numRank0 == 1) {
                    bestStrand[i] = tmi.strand;
                    bestChr[i] = tmi.chromosome;
                    bestStartPos[i] = tmi.startPosition;
                    bestEndPos[i] = tmi.endPosition;
                    bestDivergence[i] = tmi.divergence;
                    bestMapP[i] = tmi.mapP;
                    bestDcoP[i] = tmi.dcoP;
                    bestMapIndex[i] = (byte)priorityIndex[0];
                    multimaps[i] = 1;
                }
                else {
                    bestStrand[i] = Byte.MIN_VALUE;
                    bestChr[i] = Integer.MIN_VALUE;
                    bestStartPos[i] = Integer.MIN_VALUE;
                    bestEndPos[i] = Integer.MIN_VALUE;
                    bestDivergence[i] = Byte.MIN_VALUE;
                    bestMapP[i] = Byte.MIN_VALUE;
                    bestDcoP[i] = Byte.MIN_VALUE;
                    bestMapIndex[i] = Byte.MIN_VALUE;
                    if (tmi.strand == Byte.MIN_VALUE) multimaps[i] = 0;
                    else multimaps[i] = 99;
                }
            }
            boolean[] ifEvidence = EvidenceType.deCode(bestEvidence[i]);
            TagMappingInfoV3 tmi = topm.getMappingInfo(i, alignerIndex[4][0]);
            if (tmi.chromosome == bestChr[i] && tmi.startPosition == bestStartPos[i] && tmi.strand == bestStrand[i]) {
                ifEvidence[EvidenceType.PE.getIndex()] = true;
            }
            int numRank0 = this.getNumOfRank0(i, alignerIndex[0]);
            tmi = topm.getMappingInfo(i, alignerIndex[0][0]);
            if (numRank0 == 0 && tmi.chromosome == bestChr[i] && tmi.startPosition == bestStartPos[i] && tmi.strand == bestStrand[i]) {
                ifEvidence[EvidenceType.SingleBestBowtie2.getIndex()] = true;
            }
            tmi = topm.getMappingInfo(i, alignerIndex[1][0]);
            if (tmi.chromosome == bestChr[i] && tmi.startPosition == bestStartPos[i] && tmi.strand == bestStrand[i]) {
                ifEvidence[EvidenceType.BestBWA.getIndex()] = true;
            }
            tmi = topm.getMappingInfo(i, alignerIndex[2][0]);
            if (numRank0 == 0 && tmi.chromosome == bestChr[i] && tmi.startPosition == bestStartPos[i] && tmi.strand == bestStrand[i]) {
                ifEvidence[EvidenceType.SingleBestBlast.getIndex()] = true;
            }
            tmi = topm.getMappingInfo(i, alignerIndex[3][0]);
            if (numRank0 == 0 && tmi.chromosome == bestChr[i] && tmi.startPosition == bestStartPos[i] && tmi.strand == bestStrand[i]) {
                ifEvidence[EvidenceType.SingleBestBWAMEM.getIndex()] = true;
            }
            bestEvidence[i] = EvidenceType.code(ifEvidence);
        }
        topm.writeBestMappingDataSets(bestStrand, bestChr, bestStartPos, bestEndPos, bestDivergence, bestMapP, bestDcoP, multimaps, bestEvidence, bestMapIndex);
    }
    
    /**
     * Incorperate the alignment with best genetic mapping, contain hard coding on map index
     * @param tagIndex
     * @param mapIndices
     * @param br
     * @param bestStrand
     * @param bestChr
     * @param multimapsm
     * @param bestEvidence 
     */
    private void incorperateBestGeneticMapping (int tagIndex, Integer[] mapIndices, BufferedReader br, byte[] bestStrand, int[] bestChr, int[] bestStartPos, int[] bestEndPos, byte[] bestDivergence, byte[] bestMapP, byte[] bestDcoP, byte[] bestEvidence, byte[] bestMapIndices) {
        ArrayList<Integer> yMapIndexList = new ArrayList();
        try {
            for (int i = 0; i < mapIndices.length; i++) {
                String in = br.readLine();
                String[] temp;
                if (in.startsWith(" ")) {
                    temp = in.split("\\s+")[3].split(":");
                }
                else {
                    temp = in.split("\\s+")[2].split(":");
                }
                if (temp[1].equals("N")) continue;
                yMapIndexList.add(mapIndices[i]);
            }
            if (yMapIndexList.isEmpty()) return;
            double pBest = 1;
            int bestMapIndex = -1;
            for (int i = 0; i < yMapIndexList.size(); i++) {
                TagGeneticMappingInfo tgmi = topm.getGeneticMappingInfo(tagIndex, yMapIndexList.get(i));
                if (tgmi.p < 0) continue;
                if (tgmi.p < pBest) {
                    pBest = tgmi.p;
                    bestMapIndex = yMapIndexList.get(i);
                }
            }
            if (bestMapIndex == -1) return;
            if (bestMapIndex >= 25) bestMapIndex-=5;
            TagMappingInfoV3 tmi = topm.getMappingInfo(tagIndex, bestMapIndex);
            bestStrand[tagIndex] = tmi.strand;
            bestChr[tagIndex] = tmi.chromosome;
            bestStartPos[tagIndex] = tmi.startPosition;
            bestEndPos[tagIndex] = tmi.endPosition;
            bestDivergence[tagIndex] = tmi.divergence;
            bestMapP[tagIndex] = tmi.mapP;
            bestDcoP[tagIndex] = tmi.dcoP;
            bestEvidence[tagIndex] = (byte)(1<<EvidenceType.values().length-1);
            bestMapIndices[tagIndex] = (byte)bestMapIndex;
        }
        catch (Exception e) {
            System.out.println(e.toString());
            System.out.println(tagIndex);
            System.exit(1);
        }
    }
    
    /**
     * Command line version to predict if hypotheses are correct using weka RandomForest. Make sure weka is installed (in ClassPath) in local machine
     * Contain hard coding
     * @param modelFileS
     * @param tagCountFileS
     * @param inputDirS
     * @param indexDirS
     * @param outputDirS 
     */
    public void annotateBestMappingPredict (String modelFileS, String tagCountFileS, String inputDirS, String indexDirS, String outputDirS) {
        String[] attributeName = {"TagCount", "TagLength", "GC", "Chr", "Pos", "Source", "Rank", "Score", "SigSiteNum", "SigSiteRange", "P", "PRatio", "NumOfRank0", "NumOfAlign", "NumOfAlignAll", "Rorw"};
        TagCounts tc = new TagCounts(tagCountFileS, FilePacking.Bit);
        int fileSize = topm.getChunkSize();
        int fileNum;
        int left = topm.getTagCount()%fileSize;
        if (left == 0) {
            fileNum = topm.getTagCount()/fileSize;
        }
        else {
            fileNum = topm.getTagCount()/fileSize+1;
        }
        int[] alignerStartMapIndex = {0,5,10,15,20,25};
        int[][] alignerIndex = new int[alignerStartMapIndex.length][];
        alignerIndex[0] = topm.getMappingIndicesOfAligner(Aligner.Bowtie2);
        alignerIndex[1] = topm.getMappingIndicesOfAligner(Aligner.BWA);
        alignerIndex[2] = topm.getMappingIndicesOfAligner(Aligner.Blast);
        alignerIndex[3] = topm.getMappingIndicesOfAligner(Aligner.BWAMEM);
        alignerIndex[4] = topm.getMappingIndicesOfAligner(Aligner.PEEnd1);
        alignerIndex[5] = topm.getMappingIndicesOfAligner(Aligner.PEEnd2);
        String missing = "?";
        for (int i = 0; i < fileNum; i++) {
            String predictFileS = String.valueOf(i)+".pre.arff";
            predictFileS = new File(inputDirS, predictFileS).getAbsolutePath();
            String indexFileS = String.valueOf(i)+".index.txt";
            indexFileS = new File(indexDirS, indexFileS).getAbsolutePath();
            int startTagIndex = i * fileSize;
            int endTagIndex = startTagIndex + fileSize;
            if (endTagIndex > topm.getTagCount()) endTagIndex = topm.getTagCount();
            try {
                BufferedWriter bwp = new BufferedWriter (new FileWriter(predictFileS), 65536);
                BufferedWriter bwi = new BufferedWriter (new FileWriter(indexFileS), 65536);
                bwp.write("@relation trainingFinalSet\n\n");
                for (int j = 0; j < attributeName.length-1; j++) {
                    bwp.write("@attribute "+ attributeName[j]+ " numeric\n");
                }
                bwp.write("@attribute Rorw {Y,N}\n");
                bwp.write("\n@data\n");
                bwi.write("TagIndex\tMapIndex");
                bwi.newLine();
                for (int j = startTagIndex; j < endTagIndex; j++) {
                    int cnt = 0;
                    long[] tag = topm.getTag(j);
                    int index = tc.getTagIndex(tag);
                    if (index < 0) {
                        System.out.println("TagCount file and TOPM file don't match, program quits");
                        System.exit(0);
                    }
                    int readCount = tc.getReadCount(index);
                    String seq = BaseEncoder.getSequenceFromLong(tag);
                    for (int k = 0; k < topm.getTagLength(j); k++) {
                        if (seq.charAt(k) == 'G' || seq.charAt(k) == 'C') cnt++;
                    }
                    double gc = (double)cnt/topm.getTagLength(j);
                    int[] numOfRank0 = new int[alignerStartMapIndex.length];
                    int[] numOfAlign = new int[alignerStartMapIndex.length];
                    for (int k = 0; k < numOfRank0.length; k++) {
                        numOfRank0[k] = this.getNumOfRank0(j, alignerIndex[k]);
                        numOfAlign[k] = this.getNumOfAlign(j, alignerIndex[k]);
                    }
                    double pBest = this.getPBest(j);
                    if (pBest == 1) continue;
                    int numOfAlignAll = this.getNumOfAlignAll(j);
                    for (int k = 0; k < topm.getMappingNum(); k++) {
                        TagMappingInfoV3 tmi = topm.getMappingInfo(j, k);
                        TagGeneticMappingInfo tgmi = topm.getGeneticMappingInfo(j, k);
                        if (tmi.chromosome < 0) continue;
                        if (tgmi.p < 0) continue;
                        bwp.write(String.valueOf(this.boxcox(readCount, -0.181818))+","); 
                        bwp.write(String.valueOf(topm.getTagLength(j))+",");
                        bwp.write(String.valueOf(gc)+",");
                        bwp.write(String.valueOf(tmi.chromosome)+",");
                        bwp.write(String.valueOf(tmi.startPosition)+",");
                        bwp.write(String.valueOf(tmi.mappingSource)+",");
                        bwp.write(String.valueOf(tmi.mappingRank)+",");
                        if (tmi.mappingScore < 0) {
                            bwp.write(missing+",");
                        }
                        else {
                            bwp.write(String.valueOf(tmi.mappingScore)+",");
                        }
                        
                        if (tgmi.sigSiteNum < 0) {
                            bwp.write(missing+",");
                        }
                        else {
                            double boxValue = this.boxcox(tgmi.sigSiteNum, 0.101010); //262414 tags
                            bwp.write(String.valueOf(boxValue)+",");
                        }
                        if (tgmi.sigSiteRange < 0) {
                            bwp.write(missing+",");
                        }
                        else {
                            double boxValue = this.boxcox(tgmi.sigSiteRange, 0.424242); // 65536 tags and 262414 tags
                            bwp.write(String.valueOf(boxValue)+",");
                        }
                        if (tgmi.p < 0) {
                            bwp.write(missing+",");
                            bwp.write(missing+",");
                        }
                        else {
                            double boxValue = this.boxcox(minusLog10P(tgmi.p), -0.060606); // 262414 tags
                            bwp.write(String.valueOf(boxValue)+",");
                            if (tgmi.p == 0) {
                                bwp.write(String.valueOf(pBest/Double.MIN_VALUE)+",");
                            }
                            else {
                                bwp.write(String.valueOf(pBest/tgmi.p)+",");
                            }  
                        }
                        index = Arrays.binarySearch(alignerStartMapIndex, k);
                        if (index < 0) index = -index -2;
                        bwp.write(String.valueOf(numOfRank0[index])+",");
                        bwp.write(String.valueOf(numOfAlign[index])+",");
                        bwp.write(String.valueOf(numOfAlignAll)+",");
                        bwp.write("N,");
                        bwp.newLine();
                        bwi.write(String.valueOf(j)+"\t"+String.valueOf(k));
                        bwi.newLine();
                    }
                }
                bwp.flush();
                bwi.flush();
                bwp.close();
                bwi.close();
            }
            catch (Exception e) {
                System.out.println(e.toString());
                System.exit(1);
            }
        }
        for (int i = 0; i < fileNum; i++) {
            String predictFileS = String.valueOf(i)+".pre.arff";
            predictFileS = new File(inputDirS, predictFileS).getAbsolutePath();
            String outputFileS = String.valueOf(i)+".out.txt";
            outputFileS = new File(outputDirS, outputFileS).getAbsolutePath();
            try {
                Runtime run = Runtime.getRuntime();
                String cmd = "cmd /c java weka.classifiers.trees.RandomForest -p 0 -T " + predictFileS + " -l " + modelFileS + " > " +outputFileS;
                System.out.println(cmd);
                Process p = run.exec(cmd);
                p.waitFor();
                System.out.println("Prediction is made at " + outputFileS);
            }
            catch (Exception e) {
                System.out.println(e.toString());
                System.exit(1);
            }
        }
    }
    
    /**
     * Return the number of Rank0 in an aligner
     * @param tagIndex
     * @param mapIndex
     * @return 
     */
    public int getNumOfRank0 (int tagIndex, int[] mapIndex) {
        int cnt = 0;
        for (int i = 0; i < mapIndex.length; i++) {
            TagMappingInfoV3 tmi = topm.getMappingInfo(tagIndex, mapIndex[i]);
            if (tmi.chromosome < 0) continue;
            if (tmi.mappingRank != 0) continue;
            cnt++;
        }
        return cnt;
    }
    
    /**
     * Return number of alignment in an aligner
     * @param tagIndex
     * @param mapIndex
     * @return 
     */
    public int getNumOfAlign (int tagIndex, int[] mapIndex) {
        int cnt = 0;
        for (int i = 0; i < mapIndex.length; i++) {
            TagMappingInfoV3 tmi = topm.getMappingInfo(tagIndex, mapIndex[i]);
            if (tmi.chromosome < 0) continue;
            cnt++;
        }
        return cnt;
    }
    
    /**
     * Return the most significant p-value across all hypotheses
     * @param tagIndex
     * @return 
     */
    public double getPBest (int tagIndex) {
        double p = 1;
        for (int i = 0; i < topm.getMappingNum(); i++) {
            TagGeneticMappingInfo tgmi = topm.getGeneticMappingInfo(tagIndex, i);
            if (tgmi.p < 0) continue;
            if (tgmi.p < p) p = tgmi.p;
        }
        if (p == 0) p = Double.MIN_VALUE;
        return p;
    }
    
    /**
     * Return number of alignment hypotheses across all aligners, excluding PEEnd1 and PEEnd2 
     * @param tagIndex
     * @return 
     */
    public int getNumOfAlignAll (int tagIndex) {
        int cnt = 0;
        for (int i = 0; i < 20; i++) {
            TagMappingInfoV3 tmi = topm.getMappingInfo(tagIndex, i);
            if (tmi.chromosome < 0) continue;
            cnt++;
        }
        return cnt;
    }
    
    /**
     * Return -log10(p-value), if p-value doesn't exist(p < 0), reutrn Double.NEGATIVE_INFINITY
     * @param p
     * @return 
     */
    public double minusLog10P (double p) {
        if (p < 0) return Double.NEGATIVE_INFINITY;
        p = -Math.log10(p);
        if (p == Double.POSITIVE_INFINITY) p = -Math.log10(Double.MIN_VALUE);
        return p;
    }
    
    /**
     * Boxcox transform
     * @param y
     * @param lambda
     * @return 
     */
    private double boxcox (double y, double lambda) {
        if (lambda != 0) {
            return (Math.pow(y, lambda)-1)/lambda;
        }
        else {
            return Math.log(y);
        }
    }
    
    /**
     * Annotate the TOPM file using Bowtie2
     * @param samFileS
     * @param maxMappingNum is the max number of multiple alignment
     */
    public void annotateWithBowtie2 (String samFileS, int maxMappingNum) {
        String[] dataSetNames = topm.creatTagMappingInfoDatasets(topm.getMappingNum(), maxMappingNum);
        this.iniTMIBuffers(bufferNum, maxMappingNum);
        System.out.println("Reading SAM format tag alignment (Bowtie2) from: " + samFileS);
        System.out.println("Coverting SAM to TOPMHDF5...");
        byte mappingSource = Aligner.Bowtie2.getValue();
        int chunkCnt = 0;
        try {
            BufferedReader br;
            if (samFileS.endsWith(".gz")) {
                br = new BufferedReader(new InputStreamReader(new MultiMemberGZIPInputStream(new FileInputStream(new File(samFileS)))));
            } else {
                br = new BufferedReader(new FileReader(new File(samFileS)), 65536);
            }
            while(br.readLine().startsWith("@PG")==false) {};
            String inputStr = null;
            while((inputStr = br.readLine())!=null) {
                String[] temp =inputStr.split("\\s");
                int orientiation=Integer.parseInt(temp[1]);
                int chr = Integer.MIN_VALUE;
                byte strand = Byte.MIN_VALUE;
                int startPos = Integer.MIN_VALUE;
                int endPos = Integer.MIN_VALUE;
                short mappingScore = Short.MIN_VALUE;
                byte divergence = Byte.MIN_VALUE;
                String seqS = temp[9];
                if (orientiation == 4) {
                    
                }
                else if (orientiation == 16 || orientiation == 272) {
                    seqS = BaseEncoder.getReverseComplement(seqS);
                    chr = Integer.parseInt(temp[2]);
                    strand = -1;
                    int[] alignSpan = SAMUtils.adjustCoordinates(temp[5], Integer.parseInt(temp[3]));
                    startPos = alignSpan[1];
                    endPos = alignSpan[0];
                    mappingScore = Short.parseShort(temp[11].split(":")[2]);
                    if (temp[17].startsWith("NM")) {
                        divergence = Byte.parseByte(temp[17].split(":")[2]);
                    }
                    else {
                        divergence = Byte.parseByte(temp[16].split(":")[2]);
                    }
                }
                else {
                    chr = Integer.parseInt(temp[2]);
                    strand = 1;
                    int[] alignSpan = SAMUtils.adjustCoordinates(temp[5], Integer.parseInt(temp[3]));
                    startPos = alignSpan[0];
                    endPos = alignSpan[1];
                    mappingScore = Short.parseShort(temp[11].split(":")[2]);
                    if (temp[17].startsWith("NM")) {
                        divergence = Byte.parseByte(temp[17].split(":")[2]);
                    }
                    else {
                        divergence = Byte.parseByte(temp[16].split(":")[2]);
                    }
                }
                TagMappingInfoV3 theTMI = new TagMappingInfoV3(chr, strand, startPos, endPos, divergence, mappingSource, mappingScore);
                long[] seq = BaseEncoder.getLongArrayFromSeq(seqS,topm.getTagSizeInLong()*32);
                int tagIndex = topm.getTagIndex(seq);
                if (tagIndex < this.bufferTagIndexRange[0] || tagIndex >= this.bufferTagIndexRange[1]) {
                    System.out.println("The index of the tag from sam file is out of buffer range. Program quits.");
                    System.out.println("Please increase the buffer number");
                    System.exit(1);
                }
                int bufferIndex = Arrays.binarySearch(bufferStartTagIndex, tagIndex);
                if (bufferIndex < 0) bufferIndex = -bufferIndex-2;
                int bufferTagIndex = tagIndex % topm.getChunkSize();
                int mappingDatasetIndex = this.getMappingDatasetIndex(bufferIndex, bufferTagIndex);
                if (mappingDatasetIndex == Integer.MIN_VALUE) continue;
                tmiBuffers[bufferIndex][mappingDatasetIndex][bufferTagIndex] = theTMI;
                if (bufferLights[bufferIndex][bufferTagIndex] == false) {
                    lightCounts[bufferIndex]++;
                }
                bufferLights[bufferIndex][bufferTagIndex] = true;
                if (lightCounts[0] == topm.getChunkSize() && lightCounts[1] > this.updateBufferCountCutoff) {
                    this.saveTMIBufferToTOPM(tmiBuffers[0], dataSetNames, this.bufferStartTagIndex[0]/topm.getChunkSize(), mappingSource);
                    this.updateTMIBuffer();
                    System.out.println(++chunkCnt + " chunks are annotated. " + topm.getChunkNum() + " chunks in total");
                }   
            }
            this.saveTMIBufferToTOPM(tmiBuffers[0], dataSetNames, this.bufferStartTagIndex[0]/topm.getChunkSize(), mappingSource);
            topm.setMappingNum(topm.getMappingNum()+maxMappingNum);
            System.out.println(++chunkCnt + " chunks are annotated. " + topm.getChunkNum() + " chunks in total");
            br.close();
        } catch (Exception e) {
            
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    /**
     * Annotate the TOPM with BWA
     * @param samFileS
     * @param maxMappingNum 
     */
    public void annotateWithBWA (String samFileS, int maxMappingNum) {
        String[] dataSetNames = topm.creatTagMappingInfoDatasets(topm.getMappingNum(), maxMappingNum);
        this.iniTMIBuffers(bufferNum, maxMappingNum);
        System.out.println("Reading SAM format tag alignment (BWA) from: " + samFileS);
        System.out.println("Coverting SAM to TOPMHDF5...");
        byte mappingSource = Aligner.BWA.getValue();
        int chunkCnt = 0;
        try {
            BufferedReader br;
            if (samFileS.endsWith(".gz")) {
                br = new BufferedReader(new InputStreamReader(new MultiMemberGZIPInputStream(new FileInputStream(new File(samFileS)))));
            } else {
                br = new BufferedReader(new FileReader(new File(samFileS)), 65536);
            }
            String inputStr = null;
            while ((inputStr=br.readLine()).startsWith("@")) {};
            while(inputStr!=null) {
                String[] temp =inputStr.split("\\s");
                int orientiation=Integer.parseInt(temp[1]);
                int chr = Integer.MIN_VALUE;
                byte strand = Byte.MIN_VALUE;
                int startPos = Integer.MIN_VALUE;
                int endPos = Integer.MIN_VALUE;
                short mappingScore = Short.MIN_VALUE;
                byte divergence = Byte.MIN_VALUE;
                String seqS = temp[9];
                String XAString = null;
                if (temp[temp.length-1].startsWith("XA")) {
                    XAString = temp[temp.length-1].replaceFirst("XA:Z:", "");
                }
                if (orientiation == 4) {
                    
                }
                else if (orientiation == 16) {
                    seqS = BaseEncoder.getReverseComplement(seqS);
                    chr = Integer.parseInt(temp[2]);
                    strand = -1;
                    int[] alignSpan = SAMUtils.adjustCoordinates(temp[5], Integer.parseInt(temp[3]));
                    startPos = alignSpan[1];
                    endPos = alignSpan[0];
                    divergence = Byte.parseByte(temp[12].split(":")[2]);
                }
                else {
                    chr = Integer.parseInt(temp[2]);
                    strand = 1;
                    int[] alignSpan = SAMUtils.adjustCoordinates(temp[5], Integer.parseInt(temp[3]));
                    startPos = alignSpan[0];
                    endPos = alignSpan[1];
                    divergence = Byte.parseByte(temp[12].split(":")[2]);
                }
                TagMappingInfoV3 theTMI = new TagMappingInfoV3(chr, strand, startPos, endPos, divergence, mappingSource, mappingScore);
                long[] seq = BaseEncoder.getLongArrayFromSeq(seqS,topm.getTagSizeInLong()*32);
                int tagIndex = topm.getTagIndex(seq);
                if (tagIndex < this.bufferTagIndexRange[0] || tagIndex >= this.bufferTagIndexRange[1]) {
                    System.out.println("The index of the tag from sam file is out of buffer range. Program quits.");
                    System.out.println("Please increase the buffer number");
                    System.exit(1);
                    }  
                int bufferIndex = Arrays.binarySearch(bufferStartTagIndex, tagIndex);
                if (bufferIndex < 0) bufferIndex = -bufferIndex-2;
                int bufferTagIndex = tagIndex % topm.getChunkSize();
                int mappingDatasetIndex = this.getMappingDatasetIndex(bufferIndex, bufferTagIndex);
                if (mappingDatasetIndex == Integer.MIN_VALUE) continue;
                tmiBuffers[bufferIndex][mappingDatasetIndex][bufferTagIndex] = theTMI;
                if (bufferLights[bufferIndex][bufferTagIndex] == false) {
                    lightCounts[bufferIndex]++;
                }
                if (XAString != null) {
                    temp = XAString.split(";");
                    for (int i = 0; i < temp.length; i++) {
                        mappingDatasetIndex = this.getMappingDatasetIndex(bufferIndex, bufferTagIndex);
                        if (mappingDatasetIndex == Integer.MIN_VALUE) break;
                        String[] tem = temp[i].split(",");
                        chr = Integer.parseInt(tem[0]);
                        if (tem[1].startsWith("+")) {
                            strand = 1;
                            int[] alignSpan = SAMUtils.adjustCoordinates(tem[2], Integer.parseInt(tem[1].substring(1)));
                            startPos = alignSpan[1];
                            endPos = alignSpan[0];
                        }
                        else {
                            strand = -1;
                            int[] alignSpan = SAMUtils.adjustCoordinates(tem[2], Integer.parseInt(tem[1].substring(1)));
                            startPos = alignSpan[0];
                            endPos = alignSpan[1];
                        }
                        divergence = Byte.parseByte(tem[3]);
                        theTMI = new TagMappingInfoV3(chr, strand, startPos, endPos, divergence, mappingSource, mappingScore);
                        tmiBuffers[bufferIndex][mappingDatasetIndex][bufferTagIndex] = theTMI;
                    }
                }
                bufferLights[bufferIndex][bufferTagIndex] = true;
                if (lightCounts[0] == topm.getChunkSize() && lightCounts[1] > this.updateBufferCountCutoff) {
                    this.saveBWATMIBufferToTOPM(tmiBuffers[0], dataSetNames, this.bufferStartTagIndex[0]/topm.getChunkSize(), mappingSource);
                    this.updateTMIBuffer();
                    System.out.println(++chunkCnt + " chunks are annotated. " + topm.getChunkNum() + " chunks in total");
                }
                inputStr=br.readLine();
            }
            this.saveBWATMIBufferToTOPM(tmiBuffers[0], dataSetNames, this.bufferStartTagIndex[0]/topm.getChunkSize(), mappingSource);
            topm.setMappingNum(topm.getMappingNum()+maxMappingNum);
            System.out.println(++chunkCnt + " chunks are annotated. " + topm.getChunkNum() + " chunks in total");
            br.close();
        } catch (Exception e) {
            
            e.printStackTrace();
            System.exit(1);
        }       
    }
    
    /**
     * Annotate the TOPM file using bwa-mem
     * @param samFileS
     * @param maxMappingNum 
     */
    public void annotateWithBWAMEM (String samFileS, int maxMappingNum) {
        String[] dataSetNames = topm.creatTagMappingInfoDatasets(topm.getMappingNum(), maxMappingNum);
        this.iniTMIBuffers(bufferNum, maxMappingNum);
        System.out.println("Reading SAM format tag alignment (BWAMEM) from: " + samFileS);
        System.out.println("Coverting SAM to TOPMHDF5...");
        byte mappingSource = Aligner.BWAMEM.getValue();
        try {
            BufferedReader br;
            if (samFileS.endsWith(".gz")) {
                br = new BufferedReader(new InputStreamReader(new MultiMemberGZIPInputStream(new FileInputStream(new File(samFileS)))));
            } else {
                br = new BufferedReader(new FileReader(new File(samFileS)), 65536);
            }
            String inputStr = null;
            while ((inputStr=br.readLine()).startsWith("@")) {};
            int mapCnt = 0;
            int tagIndex = 0;
            int chunkCnt = 0;
            while(inputStr!=null) {
                String[] temp =inputStr.split("\\s");
                int orientiation=Integer.parseInt(temp[1]);
                int chr = Integer.MIN_VALUE;
                byte strand = Byte.MIN_VALUE;
                int startPos = Integer.MIN_VALUE;
                int endPos = Integer.MIN_VALUE;
                short mappingScore = Short.MIN_VALUE;
                byte divergence = Byte.MIN_VALUE;
                String seqS = temp[9];
                if (seqS.equals("*")) {
                    mapCnt++;
                }
                else {
                    mapCnt = 1;
                }
                if (mapCnt > maxMappingNum) {
                    inputStr=br.readLine();
                    continue;
                }
                if (orientiation == 4) {
                    
                }
                else if (orientiation == 16 || orientiation == 272) {
                    if (!seqS.equals("*")) seqS = BaseEncoder.getReverseComplement(seqS);
                    chr = Integer.parseInt(temp[2]);
                    strand = -1;
                    int[] alignSpan = SAMUtils.adjustCoordinates(temp[5], Integer.parseInt(temp[3]));
                    startPos = alignSpan[1];
                    endPos = alignSpan[0];
                    divergence = Byte.parseByte(temp[11].split(":")[2]);
                    mappingScore = Byte.parseByte(temp[12].split(":")[2]);
                }
                else if (orientiation == 0 || orientiation == 256) {
                    chr = Integer.parseInt(temp[2]);
                    strand = 1;
                    int[] alignSpan = SAMUtils.adjustCoordinates(temp[5], Integer.parseInt(temp[3]));
                    startPos = alignSpan[0];
                    endPos = alignSpan[1];
                    divergence = Byte.parseByte(temp[11].split(":")[2]);
                    mappingScore = Byte.parseByte(temp[12].split(":")[2]);
                }
                else {
                    //handle flag value of 2048 and 2064, use sequence in SAM can't find the tag in TOPM, since seqs are choped
                    //but seems they are all secondary alignment
                    //2048 and 2064 are very rare 0.04%, ignore for now.
                    //Todo, change tag identification from seq to index
                    inputStr=br.readLine();
                    continue;
                }
                TagMappingInfoV3 theTMI = new TagMappingInfoV3(chr, strand, startPos, endPos, divergence, mappingSource, mappingScore);
                if (!seqS.equals("*")) {
                    long[] seq = BaseEncoder.getLongArrayFromSeq(seqS,topm.getTagSizeInLong()*32);
                    tagIndex = topm.getTagIndex(seq);
                    if (tagIndex < this.bufferTagIndexRange[0] || tagIndex >= this.bufferTagIndexRange[1]) {
                        System.out.println("The index of the tag from sam file is out of buffer range. Program quits.");
                        System.out.println("Please increase the buffer number");
                        System.exit(1);
                    }
                }
                int bufferIndex = Arrays.binarySearch(bufferStartTagIndex, tagIndex);
                if (bufferIndex < 0) bufferIndex = -bufferIndex-2;
                int bufferTagIndex = tagIndex % topm.getChunkSize();
                int mappingDatasetIndex = this.getMappingDatasetIndex(bufferIndex, bufferTagIndex);
                if (mappingDatasetIndex == Integer.MIN_VALUE) {
                    inputStr=br.readLine();
                    continue;
                }
                tmiBuffers[bufferIndex][mappingDatasetIndex][bufferTagIndex] = theTMI;
                if (bufferLights[bufferIndex][bufferTagIndex] == false) {
                    lightCounts[bufferIndex]++;
                }
                bufferLights[bufferIndex][bufferTagIndex] = true;
                //since the MT output of bwamem are in the same order the fastq
                if (lightCounts[1] > this.updateBufferCountCutoff) {
                    this.saveBWATMIBufferToTOPM(tmiBuffers[0], dataSetNames, this.bufferStartTagIndex[0]/topm.getChunkSize(), mappingSource);
                    this.updateTMIBuffer();
                    System.out.println(++chunkCnt + " chunks are annotated. " + topm.getChunkNum() + " chunks in total");
                }
                inputStr=br.readLine();
            }
            this.saveTMIBufferToTOPM(tmiBuffers[0], dataSetNames, this.bufferStartTagIndex[0]/topm.getChunkSize(), mappingSource);
            topm.setMappingNum(topm.getMappingNum()+maxMappingNum);
            System.out.println(++chunkCnt + " chunks are annotated. " + topm.getChunkNum() + " chunks in total");
            br.close();
        } catch (Exception e) {
            
            e.printStackTrace();
            System.exit(1);
        }       
    }
    
    /**
     * Annotate the TOPM with BLAST from a directory, where slices of blast result are stored
     * @param blastDirS
     * @param maxMappingNum 
     */
    public void annotateWithBlastFromDir (String blastDirS, int maxMappingNum) {
        String[] dataSetNames = topm.creatTagMappingInfoDatasets(topm.getMappingNum(), maxMappingNum);
        this.iniTMIBuffers(2, maxMappingNum);
        System.out.println("Reading BLAST table format tag alignment (BLAST) from: " + blastDirS);
        System.out.println("Coverting BLAST to TOPMHDF5...");
        byte mappingSource = Aligner.Blast.getValue();
        File[] infiles = new File (blastDirS).listFiles();
        Arrays.sort(infiles);
        int chunkCnt = 0;
        try {
            BufferedReader br;
            for (int i = 0; i < infiles.length; i++) {
                System.out.println("Reading BLAST table format tag alignment (BLAST) from: " + infiles[i].getAbsolutePath());
                if (infiles[i].getName().endsWith("gz")) {
                    br = new BufferedReader(new InputStreamReader(new MultiMemberGZIPInputStream(new FileInputStream(infiles[i]))));
                }
                else {
                    br = new BufferedReader(new FileReader(infiles[i]), 65536);
                }
                String inputStr = null;
                while((inputStr = br.readLine())!=null) {
                    String[] temp =inputStr.split("\\s+");
                    int chr = Integer.parseInt(temp[1]);
                    byte strand = Byte.MIN_VALUE;
                    int startPos = Integer.parseInt(temp[8]);
                    int endPos = Integer.parseInt(temp[9]);
                    if (startPos < endPos) {
                        strand = 1;
                    }
                    else {
                        strand = -1;
                    }
                    short mappingScore = Short.parseShort(temp[11].replaceAll("\\..+", ""));
                    byte divergence = Byte.MIN_VALUE;
                    int tagIndex = Integer.parseInt(temp[0]);
                    if (tagIndex >= this.bufferStartTagIndex[1]) {
                        int n = (tagIndex - bufferStartTagIndex[1]) /topm.getChunkSize()+1;
                        for (int j = 0; j < n; j++) {
                            this.saveTMIBufferToTOPM(tmiBuffers[0], dataSetNames, this.bufferStartTagIndex[0]/topm.getChunkSize(), mappingSource);
                            this.updateTMIBuffer();
                        }  
                    }
                    int bufferIndex = Arrays.binarySearch(bufferStartTagIndex, tagIndex);
                    if (bufferIndex < 0) bufferIndex = -bufferIndex-2;
                    int bufferTagIndex = tagIndex % topm.getChunkSize();
                    int mappingDatasetIndex = this.getMappingDatasetIndex(bufferIndex, bufferTagIndex);
                    if (mappingDatasetIndex == Integer.MIN_VALUE) continue;
                    TagMappingInfoV3 theTMI = new TagMappingInfoV3(chr, strand, startPos, endPos, divergence, mappingSource, mappingScore);
                    tmiBuffers[bufferIndex][mappingDatasetIndex][bufferTagIndex] = theTMI;
                }
                br.close();
                System.gc();
            }
            this.saveTMIBufferToTOPM(tmiBuffers[0], dataSetNames, this.bufferStartTagIndex[0]/topm.getChunkSize(), mappingSource);
            topm.setMappingNum(topm.getMappingNum()+maxMappingNum);
            System.out.println(++chunkCnt + " chunks are annotated. " + topm.getChunkNum() + " chunks in total");
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    /**
     * Annotate the TOPM with BLAST
     * @param blastM8FileS
     * @param maxMappingNum 
     */
    public void annotateWithBLAST (String blastM8FileS, int maxMappingNum) {
        String[] dataSetNames = topm.creatTagMappingInfoDatasets(topm.getMappingNum(), maxMappingNum);
        this.iniTMIBuffers(2, maxMappingNum);
        System.out.println("Reading BLAST table format tag alignment (BLAST) from: " + blastM8FileS);
        System.out.println("Coverting BLAST to TOPMHDF5...");
        byte mappingSource = Aligner.Blast.getValue();
        int chunkCnt = 0;
        try {
            BufferedReader br;
            if (blastM8FileS.endsWith(".gz")) {
                br = new BufferedReader(new InputStreamReader(new MultiMemberGZIPInputStream(new FileInputStream(new File(blastM8FileS)))));
            } else {
                br = new BufferedReader(new FileReader(new File(blastM8FileS)), 65536);
            }
            String inputStr = null;
            while((inputStr = br.readLine())!=null) {
                String[] temp =inputStr.split("\\s+");
                int chr = Integer.parseInt(temp[1]);
                byte strand = Byte.MIN_VALUE;
                int startPos = Integer.parseInt(temp[8]);
                int endPos = Integer.parseInt(temp[9]);
                if (startPos < endPos) {
                    strand = 1;
                }
                else {
                    strand = -1;
                }
                short mappingScore = Short.parseShort(temp[11].replaceAll("\\..+", ""));
                byte divergence = Byte.MIN_VALUE;
                int tagIndex = Integer.parseInt(temp[0]);
                if (tagIndex >= this.bufferStartTagIndex[1]) {
                    int n = (tagIndex - bufferStartTagIndex[1]) /topm.getChunkSize()+1;
                    for (int j = 0; j < n; j++) {
                        this.saveTMIBufferToTOPM(tmiBuffers[0], dataSetNames, this.bufferStartTagIndex[0]/topm.getChunkSize(), mappingSource);
                        this.updateTMIBuffer();
                    }  
                }
                int bufferIndex = Arrays.binarySearch(bufferStartTagIndex, tagIndex);
                if (bufferIndex < 0) bufferIndex = -bufferIndex-2;
                int bufferTagIndex = tagIndex % topm.getChunkSize();
                int mappingDatasetIndex = this.getMappingDatasetIndex(bufferIndex, bufferTagIndex);
                if (mappingDatasetIndex == Integer.MIN_VALUE) continue;
                TagMappingInfoV3 theTMI = new TagMappingInfoV3(chr, strand, startPos, endPos, divergence, mappingSource, mappingScore);
                System.out.println(++chunkCnt + " chunks are annotated. " + topm.getChunkNum() + " chunks in total");
                tmiBuffers[bufferIndex][mappingDatasetIndex][bufferTagIndex] = theTMI;
            }
            this.saveTMIBufferToTOPM(tmiBuffers[0], dataSetNames, this.bufferStartTagIndex[0]/topm.getChunkSize(), mappingSource);
            topm.setMappingNum(topm.getMappingNum()+maxMappingNum);
            System.out.println(++chunkCnt + " chunks are annotated. " + topm.getChunkNum() + " chunks in total");
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    /**
     * Annotate the TOPM with PE, using bowtie2 sam file
     * @param PETOPMFileS
     * @param maxMappingNum 
     */
    public void annotateWithPE (String PETOPMFileS, int maxMappingNum) {
        byte forwardMappingSource = Aligner.PEEnd1.getValue(), backMappingSource = Aligner.PEEnd2.getValue();
        PETagsOnPhysicalMapV3 ptopm = new PETagsOnPhysicalMapV3(PETOPMFileS);
        String[] forwardDataSetNames = topm.creatTagMappingInfoDatasets(topm.getMappingNum(), maxMappingNum);
        topm.setMappingNum(topm.getMappingNum()+maxMappingNum);
        String[] backwardDataSetNames = topm.creatTagMappingInfoDatasets(topm.getMappingNum(), maxMappingNum);
        topm.setMappingNum(topm.getMappingNum()+maxMappingNum);
        TagMappingInfoV3[][] forwardBuffer;
        TagMappingInfoV3[][] backBuffer;
        for (int i = 0; i < topm.getChunkNum(); i++) {
            forwardBuffer = this.getPopulateTMIBuffer(maxMappingNum);
            backBuffer = this.getPopulateTMIBuffer(maxMappingNum);
            int startIndex = i*topm.getChunkSize();
            int endIndex = startIndex+topm.getChunkSize();
            if (endIndex > topm.getTagCount()) endIndex = topm.getTagCount();
            for (int j = startIndex; j < endIndex; j++) {
                long[] t = topm.getTag(j);
                int index = ptopm.getTagIndexWithLongestSeq(t);
                if (index == -1) continue;
                int max = ptopm.getMappingNum(index);
                if (max > maxMappingNum) max = maxMappingNum;
                for (int k = 0; k < max; k++) {
                    int chr = ptopm.getChr(index, k);
                    byte strand = ptopm.getStrand(index, k);
                    int startPos = ptopm.getStartPos(index, k);
                    short mappingScore = ptopm.getScore(index, k);
                    byte divergence = ptopm.getDivergence(index, k);
                    TagMappingInfoV3 theTMI = new TagMappingInfoV3(chr, strand, startPos, Integer.MIN_VALUE, divergence, forwardMappingSource, mappingScore);
                    forwardBuffer[k][j-startIndex] = theTMI;
                }
                max = ptopm.getMappingNum(ptopm.getPairIndex(index));
                if (max > maxMappingNum) max = maxMappingNum;
                for (int k = 0; k < max; k++) {
                    int chr = ptopm.getChr(ptopm.getPairIndex(index), k);
                    byte strand = ptopm.getStrand(ptopm.getPairIndex(index), k);
                    int startPos = ptopm.getStartPos(ptopm.getPairIndex(index), k);
                    short mappingScore = ptopm.getScore(ptopm.getPairIndex(index), k);
                    byte divergence = ptopm.getDivergence(ptopm.getPairIndex(index), k);
                    TagMappingInfoV3 theTMI = new TagMappingInfoV3(chr, strand, startPos, Integer.MIN_VALUE, divergence, backMappingSource, mappingScore);
                    backBuffer[k][j-startIndex] = theTMI;
                } 
            }
            this.saveTMIBufferToTOPM(forwardBuffer, forwardDataSetNames, i);
            this.saveTMIBufferToTOPM(backBuffer, backwardDataSetNames, i);
            System.out.println("Chunk " + i + "(index) with " + topm.getChunkSize() + " tags is annotated");
        }
    }
    
    /**
     * Annotate the TOPM with whole genome genetic mapping
     * @param TOGMFileS
     * @param maxMappingNum 
     */
    public void annotateWithGMGW (String TOGMFileS, int maxMappingNum) {
        TagsOnGeneticMap togm = new TagsOnGeneticMap(TOGMFileS, FilePacking.Text);
        TagGeneticMappingInfo[] gmChunk;
        String dataSetName = topm.creatTagGeneticMappingInfoGWDataset();
        for (int i = 0; i < topm.getChunkNum(); i++) {
            gmChunk = new TagGeneticMappingInfo[topm.getChunkSize()];
            for (int j = 0; j < topm.getChunkSize(); j++) {
                gmChunk[j] = new TagGeneticMappingInfo();
            }
            int startIndex = i*topm.getChunkSize();
            int endIndex = startIndex+topm.getChunkSize();
            if (endIndex > topm.getTagCount()) endIndex = topm.getTagCount();
            for (int j = startIndex; j < endIndex; j++) {
                long[] t = topm.getTag(j);
                int index = togm.getTagIndex(t);
                if (index < 0) continue;
                int chr = togm.getGChr(index);
                int position = togm.getGPos(index); //rough pos in GM
                TagGeneticMappingInfo tgmi = new TagGeneticMappingInfo(Double.NEGATIVE_INFINITY, chr, position, Integer.MIN_VALUE, Integer.MIN_VALUE);
                gmChunk[j-startIndex] = tgmi;     
            }
            topm.writeTagGeneticMappingInfoGWDataSet(dataSetName, gmChunk, i);
            if (i%100 == 0) System.out.println("Chunk " + i + "(index) with " + topm.getChunkSize() + " tags is annotated with genome wide genetic mapping");
        }
        topm.setIfHasGeneticMappingGW(true);
    }
    
    /**
     * Save mapping info (from BWA) in a buffer to TOPM
     * BWA doesn't provide mapping score, so the rank is just based on the order provided for multiple hits
     * @param tmiBuffer
     * @param dataSetNames
     * @param chunkIndex
     * @param mappingSource 
     */
    private void saveBWATMIBufferToTOPM (TagMappingInfoV3[][] tmiBuffer, String[] dataSetNames, int chunkIndex, byte mappingSource) {
        for (int i = 0; i < tmiBuffer.length; i++) {
            for (int j = 0; j < tmiBuffer[i].length; j++) {
                if (tmiBuffer[i][j] != null) continue;
                tmiBuffer[i][j] = new TagMappingInfoV3();
                tmiBuffer[i][j].setMappingSource(mappingSource);
            }
        }
        for (int i = 0; i < tmiBuffer.length; i++) {
            for (int j = 0; j < tmiBuffer[i].length; j++) {
                tmiBuffer[i][j].setMappingRank((byte)i);
            }
        }
        topm.writeTagMappingInfoDataSets(dataSetNames, tmiBuffer, chunkIndex);
    }
    
    
    /**
     * Save mapping info in a buffer to TOPM. This works for PE and Genetic mapping. When mappingSource == Byte.MIN_VALUE, it means the tag doesn't exist in PE list or Genetic mapping
     * @param tmiBuffer
     * @param dataSetNames
     * @param chunkIndex
     * @param mappingSource 
     */
    private void saveTMIBufferToTOPM (TagMappingInfoV3[][] tmiBuffer, String[] dataSetNames, int chunkIndex) {
        for (int i = 0; i < tmiBuffer[0].length; i++) {
            int sum = 0;
            for (int j = 0; j < tmiBuffer.length; j++) {
                if (tmiBuffer[j][i].mappingSource == Byte.MIN_VALUE) sum++;
            }
            if (sum == tmiBuffer.length) continue;
            TreeSet<Short> set = new TreeSet();
            for (int j = 0; j < tmiBuffer.length; j++) {
                set.add(tmiBuffer[j][i].mappingScore);
            }
            Short[] sA = set.toArray(new Short[set.size()]);
            byte[] rank = new byte[tmiBuffer.length];
            for (int j = 0; j < rank.length; j++) {
                rank[j] = (byte)(sA.length - Arrays.binarySearch(sA, tmiBuffer[j][i].mappingScore) -1);
                tmiBuffer[j][i].setMappingRank(rank[j]);
            }
        }
        topm.writeTagMappingInfoDataSets(dataSetNames, tmiBuffer, chunkIndex);
    }
    
    /**
     * Save mapping info in a buffer to TOPM. This works for bowtie2, blast and bwamem
     * @param tmiBuffer
     * @param dataSetNames
     * @param chunkIndex
     * @param mappingSource 
     */
    private void saveTMIBufferToTOPM (TagMappingInfoV3[][] tmiBuffer, String[] dataSetNames, int chunkIndex, byte mappingSource) {
        for (int i = 0; i < tmiBuffer.length; i++) {
            for (int j = 0; j < tmiBuffer[i].length; j++) {
                if (tmiBuffer[i][j] != null) continue;
                tmiBuffer[i][j] = new TagMappingInfoV3();
                tmiBuffer[i][j].setMappingSource(mappingSource);
            }
        }
        for (int i = 0; i < tmiBuffer[0].length; i++) {
            TreeSet<Short> set = new TreeSet();
            for (int j = 0; j < tmiBuffer.length; j++) {
                set.add(tmiBuffer[j][i].mappingScore);
            }
            Short[] sA = set.toArray(new Short[set.size()]);
            byte[] rank = new byte[tmiBuffer.length];
            for (int j = 0; j < rank.length; j++) {
                rank[j] = (byte)(sA.length - Arrays.binarySearch(sA, tmiBuffer[j][i].mappingScore) -1);
                tmiBuffer[j][i].setMappingRank(rank[j]);
            }
        }
        topm.writeTagMappingInfoDataSets(dataSetNames, tmiBuffer, chunkIndex);
    }
    
    /**
     * Update TMI buffers. After the buffers[0] is written to TOPM, buffers[i] moves to buffers[i-1]. This creats a new buffer at the end
     */
    private void updateTMIBuffer () {
        for (int i = 0; i < tmiBuffers.length-1; i++) {
            tmiBuffers[i] = tmiBuffers[i+1];
            bufferStartTagIndex[i] = bufferStartTagIndex[i+1];
            bufferLights[i] = bufferLights[i+1];
            lightCounts[i] = lightCounts[i+1];
        }
        tmiBuffers[tmiBuffers.length-1] = new TagMappingInfoV3[tmiBuffers[0].length][topm.getChunkSize()];
        bufferStartTagIndex[tmiBuffers.length-1]+=topm.getChunkSize();
        bufferLights[tmiBuffers.length-1] = new boolean[topm.getChunkSize()];
        lightCounts[tmiBuffers.length-1] = 0;
        this.calBufferTagIndexRange();
    }
    
    private TagMappingInfoV3[][] getPopulateTMIBuffer (int maxMappingNum) {
        TagMappingInfoV3[][] tmiBuffer = new TagMappingInfoV3[maxMappingNum][topm.getChunkSize()] ;
        for (int i = 0; i < tmiBuffer.length; i++) {
            for (int j = 0; j < tmiBuffer[i].length; j++) {
                tmiBuffer[i][j] = new TagMappingInfoV3();
            }
        }
        return tmiBuffer;
    }
    
    /**
     * Initialize TagMappingInfo buffers
     * @param bufferNum number of buffers
     * @param maxMappingNum maximum mapping information which will be held
     */
    private void iniTMIBuffers (int bufferNum, int maxMappingNum) {
        tmiBuffers = new TagMappingInfoV3[bufferNum][maxMappingNum][topm.getChunkSize()];
        bufferStartTagIndex = new int[bufferNum];
        for (int i = 0; i < bufferNum; i++) {
            bufferStartTagIndex[i] = i*topm.getChunkSize();
        }
        this.calBufferTagIndexRange();
        this.bufferLights = new boolean[bufferNum][topm.getChunkSize()];
        this.lightCounts = new int[bufferNum];
        updateBufferCountCutoff = (int)((double)topm.getChunkSize() * 0.2);
    }
    
    /**
     * Calculate the tag index range in these buffers
     */
    private void calBufferTagIndexRange () {
        this.bufferTagIndexRange = new int[2];
        bufferTagIndexRange[0] = this.bufferStartTagIndex[0];
        bufferTagIndexRange[1] = this.bufferStartTagIndex[tmiBuffers.length-1]+topm.getChunkSize();
    }
    
    
    private int getMappingDatasetIndex (int bufferIndex, int bufferTagIndex) {
        for (int i = 0; i < tmiBuffers[0].length; i++) {
            if (tmiBuffers[bufferIndex][i][bufferTagIndex] == null) {
                return i;
            }
        }
        return Integer.MIN_VALUE;
    }
}
