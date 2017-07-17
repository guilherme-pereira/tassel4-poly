/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package net.maizegenetics.gbs.pipeline;

import ch.systemsx.cisd.hdf5.*;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import net.maizegenetics.gbs.maps.PETagsOnPhysicalMap;
import net.maizegenetics.gbs.maps.PETagsOnPhysicalMapV3;
import net.maizegenetics.gbs.maps.TagGeneticMappingInfo;
import net.maizegenetics.gbs.maps.TagMappingInfoV3;
import net.maizegenetics.gbs.maps.TagMappingInfoV3.Aligner;
import net.maizegenetics.gbs.maps.TagsOnGeneticMap;
import net.maizegenetics.gbs.maps.TagsOnPhysMapHDF5;
import net.maizegenetics.gbs.maps.TagsOnPhysicalMap;
import net.maizegenetics.gbs.maps.TagsOnPhysicalMapV3;
import net.maizegenetics.gbs.pipeline.AnnotateTOPM.EvidenceType;
import net.maizegenetics.gbs.tagdist.PETagCounts;
import net.maizegenetics.gbs.tagdist.TagCounts;
import net.maizegenetics.gbs.tagdist.TagsByTaxa;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.gbs.tagdist.TagsByTaxaByteHDF5TagGroups;
import net.maizegenetics.gbs.util.BaseEncoder;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.BitAlignment;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.util.OpenBitSet;
import org.apache.log4j.Logger;

/**
 *
 * @author Fei Lu
 */
public class FeiPipelines {
    private Logger myLogger = Logger.getLogger(FeiPipelines.class);
    
    public FeiPipelines () {
        //this.pipelinePE();
        //this.testPipeline();
        this.GBSV3Pipeline();
        
    }
    
    public static void main (String[] args) {
        new FeiPipelines();
        
    }
    
    public void GBSV3Pipeline () {
        //this.addTripsicum();
        //this.mkV3Alignment(); //take some time to align using blast+
        //this.initializeV3TOPM();
        //this.annotateV3TOPM();
        //this.transformAnchorHDF5();
        //this.hypothesisGeneticMapping();
        //this.annotateTOPMHDF5WithGMGW();
        //this.predictHypothesis();
        //this.importPrediction();
        this.callSNPServer();
        //this.geneticMapping();
    
    }
    
    public void geneticMapping () {
        //String hapMapHDF5 = "/workdir/mingh/AllZeaGBSv27i3b.imp.hmp.h5";
        //String hapMapHDF5 = "/workdir/mingh/GBS27.small.imp.hmp.h5";
        //String hapMapHDF5 = "M:/GBSV3/genotype/GBS27.small.imp.hmp.h5";
        String hapMapHDF5 = "M:/GBSV3/genotype/GBS27.small.sBit.h5";
        //String tbtHDF5 = "/workdir/mingh/smallerTBTHDF5_mergedtaxa_pivot_20120921.h5";
        String tbtHDF5 = "M:/GBSV3/tbt/smallerTBTHDF5_mergedtaxa_pivot_20120921.h5";
        //String outfileS = "/workdir/mingh/outfile.txt";
        String outfileS = "M:/GBSV3/mappingResult/outfile.txt";
        //TagAgainstAnchorOriginal taa = new TagAgainstAnchorOriginal(hapMapHDF5, tbtHDF5, outfileS, 0.000001, 20, -1, 4096);
        //TagAgainstAnchor.getChunkNum(tbtHDF5, 1024);    
    }
    
    public void callSNPServer () {
        String topmFileS = "/workdir/mingh/v3.bowtie2.bwa.blast.bwamem.PE.GM1000Sites.GMGW.BestPos.topm.h5";
        String pedigreeFileS = "/workdir/mingh/AllZeaPedigree2012oct01C.txt";
        String tbtFileS = "/workdir/mingh/mergeTBTHDF5_mergedtaxa_pivot_20120921.h5";
        String referenceFileS = "/workdir/mingh/ZmB73_RefGen_v2.fa";
        String hapmapFileS = "/workdir/mingh/myGBSGenos.chr+.hmp.txt";
        
        String arguments = "-i " + tbtFileS + " -m " + topmFileS + " -o " + hapmapFileS + " -p " + pedigreeFileS + " -ref " + referenceFileS;
        arguments = arguments + " -mxSites 700000 -mnF 0.8 -mnMAF 0.001 -mnMAC 10 -mnLCov 0.1 -cF -sC 10 -eC 10";
        String[] args = arguments.split(" ");
        DiscoverySNPCallerV3Plugin d = new DiscoverySNPCallerV3Plugin ();
        d.setParameters(args);
		d.performFunction(null); 
    }
    
    public void importPrediction () {
        String topmFileS = "M:/GBSV3/topm/v3.bowtie2.bwa.blast.bwamem.PE.GM1000Sites.GMGW.BestPos.topm.h5";
        String predictDirS = "M:/GBSV3/wekaPrediction/output/";
        String indexDirS = "M:/GBSV3/wekaPrediction/index/";
        TagsOnPhysicalMapV3 topm = new TagsOnPhysicalMapV3(topmFileS);
        AnnotateTOPM anno = new AnnotateTOPM(topm);
        Aligner priorityAligner = Aligner.Bowtie2;
        anno.annotateBestMappingImport(indexDirS, predictDirS, priorityAligner);
    }
    
    public void predictHypothesis () {
        String topmFileS = "M:/GBSV3/topm/v3.bowtie2.bwa.blast.bwamem.PE.GM1000Sites.GMGW.BestPos.topm.h5";
        TagsOnPhysicalMapV3 topm = new TagsOnPhysicalMapV3(topmFileS);
        //String wekaFileS = "C:/Program Files/Weka-3-6/weka.jar";
        String modelFileS = "E:/Research/gbsv3/modelDevelop/262144Tags_RandomForest.model";
        String tagCountFileS = "M:/GBSV3/tagCount/gbsV3.cnt";
        String inputDirS = "M:/GBSV3/wekaPrediction/input/";
        String indexDirS = "M:/GBSV3/wekaPrediction/index/";
        String outputDirS = "M:/GBSV3/wekaPrediction/output/";
        AnnotateTOPM anno = new AnnotateTOPM(topm);
        anno.annotateBestMappingPredict(modelFileS, tagCountFileS, inputDirS, indexDirS, outputDirS);
    }
    
    public void annotateTOPMHDF5WithGMGW () {
        String inputFileS = "M:/GBSV3/topm/v3.topm.h5";
        TagsOnPhysicalMapV3 topm = new TagsOnPhysicalMapV3(inputFileS);
        AnnotateTOPM anno = new AnnotateTOPM (topm);
        String TOGMFileS = "M:/pav/PhyGenMapping/v1.togm.txt";
        anno.annotateWithGMGW(TOGMFileS, 1);
    }
    
    public void hypothesisGeneticMapping () {
        String topmFileS = "/workdir/mingh/v3.topm.h5";
        String hapMapHDF5 = "/workdir/mingh/GBS27.sBit.h5";
        String tbtHDF5 = "/workdir/mingh/mergeTBTHDF5_mergedtaxa_pivot_20120921.h5";
        TagAgainstAnchorHypothesis taah = new TagAgainstAnchorHypothesis (hapMapHDF5, tbtHDF5, topmFileS, "Bowtie2", 1000, 0.000001, 20, -1);
    }
    
    public void transformAnchorHDF5 () {
        String hapMapInputFileS = "/workdir/mingh/AllZeaGBSv27i3b.imp.hmp.h5";
        String sBitFileS = "/workdir/mingh/GBS27.sBit.h5";
        new SimpleGenotypeSBit(hapMapInputFileS, sBitFileS);
        //new SimpleGenotypeSBit(sBitFileS);
    }
    
    public void annotateV3TOPM () {
        String inputFileS = "M:/GBSV3/topm/v3.topm.h5";
        //String inputFileS = "/workdir/mingh/v3.topm.h5";
        TagsOnPhysicalMapV3 topm = new TagsOnPhysicalMapV3(inputFileS);
        AnnotateTOPM anno = new AnnotateTOPM (topm);
        String bowtie2SamFileS = "M:/GBSV3/alignment/v3.bowtie2.sam.gz";
        //anno.annotateWithBowtie2(bowtie2SamFileS, 5);
        String bwaSamFileS = "M:/GBSV3/alignment/v3.bwa.sam.gz";
        //anno.annotateWithBWA(bwaSamFileS, 5);
        //String blastDirS = "M:/GBSV3/alignment/blastResult/";
        String blastDirS = "/workdir/mingh/blastResult/";
        //anno.annotateWithBlastFromDir(blastDirS, 5);
        String bwamemSamFileS = "/workdir/mingh/v3.bwamem.sam.gz";
        //anno.annotateWithBWAMEM(bwamemSamFileS, 5);
        String PETOPMFileS = "M:/af/ptopm/PE.topm";
        anno.annotateWithPE(PETOPMFileS, 5);
        ////String TOGMFileS = "M:/pav/PhyGenMapping/v1.togm.txt";
        ////anno.annotateWithGMGW(TOGMFileS, 1);
    }
    
    public void initializeV3TOPM () {
        String inputFileS = "M:/GBSV3/tagCount/gbsV3.cnt";
        String outputFileS = "M:/GBSV3/topm/v3.topm.h5";
        TagCounts tc = new TagCounts(inputFileS, FilePacking.Bit);
        TagsOnPhysicalMapV3.createFile(tc, outputFileS);
    }
    
    /**
     * make alignment using bowtie2m, bwa and blast+, fasta file need to  be split to many small files to speed up blast (parallel)
     */
    public void mkV3Alignment () {
        String inputFileS = "M:/GBSV3/tagCount/gbsV3.cnt";
        String fastqFileS = "M:/GBSV3/alignment/v3.fastq";
        //new FeiUtils ().convertTagCount2Fastq(inputFileS, fastqFileS);
        
        String fastaFileS = "M:/GBSV3/alignment/v3.fasta";
        //TagCounts tc = new TagCounts(inputFileS, FilePacking.Bit);
        //tc.toFASTA(fastaFileS);
        
        String zipFastaFileS = "M:/GBSV3/alignment/v3.fasta.gz";
        String fastaDirS = "M:/GBSV3/alignment/splitFasta/";
        FeiUtils util = new FeiUtils();
        util.splitFastaFileS(zipFastaFileS, fastaDirS, 20000);  
    }
    
    public void addTripsicum () {
        String sourceDir = "M:/GBSV3/tagCount/source/";
        String mergeFileS = "M:/GBSV3/tagCount/gbsV3.cnt";
        String arguments = "-i " + sourceDir + " -o " + mergeFileS + " -c " + 10 ;
        String[] args = arguments.split(" ");
        MergeMultipleTagCountPlugin m = new MergeMultipleTagCountPlugin();
        m.setParameters(args);
		m.performFunction(null);  
    }
    
    public void testPipeline () {
        //this.mkSmallTagCount();
        //this.mkFastq();
        //this.mkFasta();
        //this.mkTOPM(); //old version
        //this.mkTOPMHDF5(); //old version
        //this.initializeTOPMHDF5();
        //this.annotateTOPMHDF5WithAligner();
        
        //this.mkSmallTOPM();
        //this.mkSmallTOPMWithBWAMEM(); //alternative method to get smallTOPM, fast before lab meeting
        //this.mkSmallAnchorHDF5Test();
        //this.transformAnchorHDF5Test();
        //this.mkSmallTBTHDF5Test();
        //this.mkTBTTagBlockTest();
        //this.geneticMappingOldTest();
        //this.geneticMappingTest(); //for evertying mapping
        //this.hypothesisGeneticMappingTest();
        //this.annotateTOPMHDF5WithGMGWTest();
        //this.outputTextMapTest();
        //this.outputTextMapTest2(); //output more tags for training model
        //this.predictHypothesisTest();
        //this.importPredictionTest();
        this.callSNPTest();
        //this.checkSNPCallv3();
        //this.checkSNPCallv2();
        //this.checkSNPCall();
        //this.trackTag();
        //this.callSNPServerTest();
    }
    
    public void callSNPServerTest () {
        String topmFileS = "/workdir/mingh/small65536.GM1000Site.GMGW.Best.V3.TOPM.h5";
        String pedigreeFileS = "/workdir/mingh/AllZeaPedigree2012oct01C.txt";
        String tbtFileS = "/workdir/mingh/small1024TBTHDF5_mergedtaxa_pivot_20120921.h5";
        String referenceFileS = "/workdir/mingh/ZmB73_RefGen_v2.fa";
        String hapmapFileS = "/workdir/mingh/myGBSGenos.chr+.hmp.txt";
        
        String arguments = "-i " + tbtFileS + " -m " + topmFileS + " -o " + hapmapFileS + " -p " + pedigreeFileS + " -ref " + referenceFileS;
        arguments = arguments + " -mxSites 700000 -mnF 0.8 -mnMAF 0.001 -mnMAC 10 -mnLCov 0.1 -cF -sC 10 -eC 10";
        String[] args = arguments.split(" ");
        DiscoverySNPCallerV3Plugin d = new DiscoverySNPCallerV3Plugin ();
        d.setParameters(args);
		d.performFunction(null); 
    }
    
    public void trackTag () {
        String tagSeq = "CAGCACGTCACCCATTCCGAAAATGGTTGTCGGGGTGCATACAAAGCACGAGTTTTTGCCACCG";
        String topmFileS = "M:/GBSV3/topm/v3.bowtie2.bwa.blast.bwamem.PE.GM1000Sites.GMGW.BestPos.topm.h5";
        TagsOnPhysicalMapV3 topm = new TagsOnPhysicalMapV3(topmFileS);
        long[] t = BaseEncoder.getLongArrayFromSeq(tagSeq);
        int index = topm.getTagIndex(t);
        int bestMapIndex = topm.getBestMapIndex(index);
        TagMappingInfoV3 tmi = topm.getMappingInfo(index, bestMapIndex);
        TagGeneticMappingInfo tgmi = topm.getGeneticMappingInfo(index, bestMapIndex);
        System.out.println("OK");
    }
    
    public void checkSNPCall () {
        String outfileV3 = "M:/outfileV3.txt";
        String outfileV2 = "M:/outfileV2.txt";
        String tagCountV3 = "M:/GBSV3/tagCount/gbsV3.cnt";
        String tagCountV2 = "N:/Zea/AllZeaBuild_2.X/02_TagCounts/02_MergedTagCounts/AllZeaMasterTags_c10_20120606.cnt";
        String updateV3 = "M:/updateV3.txt";
        String updateV2 = "M:/updateV2.txt";
        TagCounts tc3 = new TagCounts (tagCountV3, FilePacking.Byte);
        TagCounts tc2 = new TagCounts (tagCountV2, FilePacking.Byte);
        try {
            BufferedReader br = new BufferedReader (new FileReader(outfileV3), 65536);
            BufferedWriter bw = new BufferedWriter (new FileWriter(updateV3), 65536);
            String temp;
            while ((temp = br.readLine()) != null) {
                String[] tem = temp.split("\t");
                long[] t = new long[2];
                t[0] = Long.valueOf(tem[1]);
                t[1] = Long.valueOf(tem[2]);
                int index = tc3.getTagIndex(t);
                bw.write(temp+"\t"+tc3.getReadCount(index));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        try {
            BufferedReader br = new BufferedReader (new FileReader(outfileV2), 65536);
            BufferedWriter bw = new BufferedWriter (new FileWriter(updateV2), 65536);
            String temp;
            while ((temp = br.readLine()) != null) {
                String[] tem = temp.split("\t");
                long[] t = new long[2];
                t[0] = Long.valueOf(tem[1]);
                t[1] = Long.valueOf(tem[2]);
                int index = tc2.getTagIndex(t);
                bw.write(temp+"\t"+tc2.getReadCount(index));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void checkSNPCallv2 () {
        String topmFileS = "N:/Zea/AllZeaBuild_2.X/04_TOPM/AllZeaMasterTags_c10_20120703.topm";
        String outfile = "M:/outfileV2.txt";
        TagsOnPhysicalMap topm = new TagsOnPhysicalMap(topmFileS, true);
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(outfile), 65536);
            for (int i = 0; i < topm.getTagCount(); i++) {
                if (topm.getChromosome(i) != 10) continue;
                if (topm.getStrand(i) != 1) continue;
                if (topm.getStartPosition(i) != 70587) continue;
                long[] t = topm.getTag(i);
                bw.write(BaseEncoder.getSequenceFromLong(t)+"\t"+t[0]+"\t"+t[1]);
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void checkSNPCallv3 () {
        String topmFileS = "M:/GBSV3/topm/v3.bowtie2.bwa.blast.bwamem.PE.GM1000Sites.GMGW.BestPos.topm.h5";
        String outfile = "M:/outfileV3.txt";
        TagsOnPhysicalMapV3 topm = new TagsOnPhysicalMapV3(topmFileS);
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(outfile), 65536);
            for (int i = 0; i < topm.getTagCount(); i++) {
                if (topm.getChromosome(i) != 10) continue;
                if (topm.getStrand(i) != 1) continue;
                if (topm.getStartPosition(i) != 70587) continue;
                long[] t = topm.getTag(i);
                bw.write(BaseEncoder.getSequenceFromLong(t)+"\t"+t[0]+"\t"+t[1]);
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        } 
    }
    
    public void callSNPTest () {
        String topmFileS = "M:/GBStest/topm/small65536.GM1000Site.GMGW.Best.V3.TOPM.h5";
        String pedigreeFileS = "M:/GBStest/keyfile/AllZeaPedigree2012oct01C.txt";
        String tbtFileS = "M:/GBStest/tbt/small75000TBTHDF5_mergedtaxa_pivot_20120921.h5";
        String referenceFileS = "M:/Database/maizeReference/chr10.fasta";
        String hapmapFileS = "M:/GBStest/hapmap/myGBSGenos.chr+.hmp.txt";
        
        String arguments = "-i " + tbtFileS + " -m " + topmFileS + " -o " + hapmapFileS + " -p " + pedigreeFileS + " -ref " + referenceFileS;
        arguments = arguments + " -mxSites 10000 -mnF 0.8 -mnMAF 0.001 -mnMAC 10 -mnLCov 0.1 -cF -sC 10 -eC 10";
        String[] args = arguments.split(" ");
        DiscoverySNPCallerV3Plugin d = new DiscoverySNPCallerV3Plugin ();
        d.setParameters(args);
		d.performFunction(null); 
    }
    
    public void importPredictionTest () {
        String topmFileS = "M:/GBStest/topm/small65536.GM1000Site.GMGW.Best.V3.TOPM.h5";
        String predictDirS = "M:/GBStest/wekaPrediction/output/";
        String indexDirS = "M:/GBStest/wekaPrediction/index/";
        TagsOnPhysicalMapV3 topm = new TagsOnPhysicalMapV3(topmFileS);
        AnnotateTOPM anno = new AnnotateTOPM(topm);
        Aligner priorityAligner = Aligner.Bowtie2;
        anno.annotateBestMappingImport(indexDirS, predictDirS, priorityAligner);
        
    }
    
    public void predictHypothesisTest () {
        String topmFileS = "M:/GBStest/topm/small65536.GM1000Site.GMGW.V3.TOPM.h5";
        TagsOnPhysicalMapV3 topm = new TagsOnPhysicalMapV3(topmFileS);
        //String wekaFileS = "C:/Program Files/Weka-3-6/weka.jar";
        String modelFileS = "E:/Research/gbsv3/modelDevelop/262144Tags_RandomForest.model";
        String tagCountFileS = "M:/GBSV3/tagCount/gbsV3.cnt";
        String inputDirS = "M:/GBStest/wekaPrediction/input/";
        String indexDirS = "M:/GBStest/wekaPrediction/index/";
        String outputDirS = "M:/GBStest/wekaPrediction/output/";
        AnnotateTOPM anno = new AnnotateTOPM(topm);
        anno.annotateBestMappingPredict(modelFileS, tagCountFileS, inputDirS, indexDirS, outputDirS);
    }
    
    public void outputTextMapTest2 () {
        String topmFileS = "M:/GBSV3/topm/v3.bowtie2.bwa.blast.bwamem.PE.GM1000Sites.GMGW.topm.h5";
        String subTopmFileS = "M:/GBStest/topm/small100000.GM1000Site.GMGW.V3.TOPM.h5";
        int[] index = new int[100000];
        for (int i = 0; i < index.length; i++) index[i] = i;
        //TagsOnPhysicalMapV3 topm = new TagsOnPhysicalMapV3(topmFileS);
        //topm.writeSubTOPM(subTopmFileS, index);
        
        String textFileS = "E:/Research/gbsv3/hypothesisTest/small100000.GM.V3.TOPM.map.txt";
        String tagCountFileS = "M:/GBSV3/tagCount/gbsV3.cnt";
        TagsOnPhysicalMapV3 topm = new TagsOnPhysicalMapV3(subTopmFileS);
        topm.writeTextMap(tagCountFileS, textFileS);
    }
    
    public void outputTextMapTest () {
        String topmFileS = "M:/GBStest/topm/small65536.GM1000Site.GMGW.V3.TOPM.h5";
        String outfileS = "E:/Research/gbsv3/hypothesisTest/small65536.GM.V3.TOPM.map.txt";
        String tagCountFileS = "M:/GBSV3/tagCount/gbsV3.cnt";
        TagsOnPhysicalMapV3 topm = new TagsOnPhysicalMapV3(topmFileS);
        topm.writeTextMap(tagCountFileS, outfileS);
    }
    
    public void annotateTOPMHDF5WithGMGWTest () {
        String inputFileS = "M:/GBStest/topm/small65536.GM1000Site.GMGW.V3.TOPM.h5";
        TagsOnPhysicalMapV3 topm = new TagsOnPhysicalMapV3(inputFileS);
        AnnotateTOPM anno = new AnnotateTOPM (topm);
        String TOGMFileS = "M:/pav/PhyGenMapping/v1.togm.txt";
        anno.annotateWithGMGW(TOGMFileS, 1);
    }
    
    public void hypothesisGeneticMappingTest () {
        //on pc
        //String topmFileS = "M:/GBStest/topm/small65536.V3.TOPM - Copy.h5";
        //String hapMapHDF5 = "M:/GBStest/genotype/GBS27.small1024.sBit.h5"; //1024 sites
        //String tbtHDF5 = "M:/GBStest/tbt/small4096TBTHDF5_mergedtaxa_pivot_20120921.h5"; // 4096 tags
        //TagAgainstAnchorHypothesis taah = new TagAgainstAnchorHypothesis (hapMapHDF5, tbtHDF5, topmFileS, "Bowtie2", 0.000001, 20, 1);
        
        //on server
        String topmFileS = "/workdir/mingh/small65536.V3.TOPM.h5";
        String hapMapHDF5 = "/workdir/mingh/GBS27.sBit.h5"; //1024 sites
        String tbtHDF5 = "/workdir/mingh/mergeTBTHDF5_mergedtaxa_pivot_20120921.h5"; // 4096 tags
        TagAgainstAnchorHypothesis taah = new TagAgainstAnchorHypothesis (hapMapHDF5, tbtHDF5, topmFileS, "Bowtie2", 1000, 0.000001, 20, -1);
    }
    
    public void geneticMappingTest () {
        //on pc
        //String hapMapHDF5 = "M:/GBStest/genotype/GBS27.small1024.sBit.h5"; //1024 sites
        //String tbtHDF5 = "M:/GBStest/tbt/small4096TBTHDF5_mergedtaxa_pivot_20120921.h5"; // 4096 tags
        //String blockFileS = "M:/GBStest/tbtTagBlock/block.ttb";
        //String outfileS = "M:/GBStest/mappingResult/outfile.txt";
        
        //on server
        String hapMapHDF5 = "/workdir/mingh/GBS27.sBit.h5"; //1024 sites
        String tbtHDF5 = "/workdir/mingh/small192TBTHDF5_mergedtaxa_pivot_20120921.h5"; // 4096 tags
        String blockFileS = "/workdir/mingh/block.ttb";
        String outfileS = "/workdir/mingh/outfile.txt";
        TagAgainstAnchor taa = new TagAgainstAnchor(hapMapHDF5, tbtHDF5, blockFileS, outfileS, 0.000001, 20, -1, 4096);
    }
    
    public void geneticMappingOldTest () {
        String hapMapInputFileS = "M:/GBStest/genotype/GBS27.small1024.imp.hmp.h5";
        String tbtHDF5 = "M:/GBStest/tbt/small4096TBTHDF5_mergedtaxa_pivot_20120921.h5"; // 4096 tags
        String outfileS = "M:/GBStest/mappingResult/outfile.txt";
        new TagAgainstAnchorOld (hapMapInputFileS, tbtHDF5, outfileS, 0.000001, 20, -1, 256, false);
    }
    
    public void mkTBTTagBlockTest () {
        String tbtHDF5 = "M:/GBStest/tbt/small4096TBTHDF5_mergedtaxa_pivot_20120921.h5";
        String topmFileS = "M:/GBStest/topm/ini.topm.h5";
        String blockFileS = "M:/GBStest/tbtTagBlock/block.ttb";
        new TagMappingUtils().mkTBTTagBlockFile(tbtHDF5, topmFileS, blockFileS, "Bowtie2");
    }
    
    public void mkSmallTBTHDF5Test () {
        String inputTBTS = "M:/GBSV3/tbt/mergeTBTHDF5_mergedtaxa_pivot_20120921.h5";
        String outputTBTS = "M:/GBStest/tbt/small75000TBTHDF5_mergedtaxa_pivot_20120921.h5";
        TagsByTaxaByteHDF5TagGroups tbt = new TagsByTaxaByteHDF5TagGroups (inputTBTS);
        int[] tagIndex = new int[75000];
        for (int i = 0; i < tagIndex.length; i++) {
            //tagIndex[i] = (int)(tbt.getTagCount()*Math.random()); //random chosen
            tagIndex[i] = i; //first N chosen
        }
        Arrays.sort(tagIndex);
        tbt.writeDistFile(outputTBTS, tagIndex);
    }
    
    public void transformAnchorHDF5Test () {
        String hapMapInputFileS = "M:/GBStest/genotype/GBS27.small1024.imp.hmp.h5";
        String sBitFileS = "M:/GBStest/genotype/GBS27.small1024.sBit.h5";
        new SimpleGenotypeSBit(hapMapInputFileS, sBitFileS);
        //new SimpleGenotypeSBit(sBitFileS);
    }
    
    public void mkSmallAnchorHDF5Test () {
        String hapMapInputS = "M:/GBSV3/genotype/AllZeaGBSv27i3b.imp.hmp.h5";
        String hapMapOutputS = "M:/GBStest/genotype/GBS27.small1024.imp.hmp.h5";
        Alignment a = ImportUtils.readGuessFormat(hapMapInputS, true);
        
        int[] snpIndex = new int[1024];
        int startIndex = 0;
        for (int i = 0; i < snpIndex.length; i++) {
            snpIndex[i] = i+startIndex;
        }
        ExportUtils.writeToMutableHDF5(a, hapMapOutputS, snpIndex);
    }
    
    public void mkSmallTOPMWithBWAMEM () {
 //Step 1     
        String inputFileS = "M:/GBSV3/topm/v3.bowtie2.bwa.blast.topm.h5";
        String outputFileS = "M:/GBStest/topm/small65536.V3.TOPM.h5";
        TagsOnPhysicalMapV3 topm = new TagsOnPhysicalMapV3 (inputFileS);
        int[] index = new int[65536];
        for (int i = 0; i < index.length; i++) {
            index[i] = i;
        }
        topm.writeSubTOPM(outputFileS, index);
 
 /*//Step 2
        String inputFileS = "M:/GBStest/topm/small65536.V3.TOPM.h5";
        String fastqFileS = "M:/GBStest/alignment/addbwamem/65536.fq";
        FeiUtils util = new FeiUtils ();
        util.convertTOPM2Fastq(inputFileS, fastqFileS);
 */
 //Step 3 Alignment by BWA MEM
        String topmFileS = "M:/GBStest/topm/small65536.V3.TOPM.h5";
        String bwamemSamFileS = "M:/GBStest/alignment/addbwamem/65536.bwamem.sam";
        TagsOnPhysicalMapV3 topm2 = new TagsOnPhysicalMapV3(topmFileS);
        AnnotateTOPM anno = new AnnotateTOPM (topm2);
        anno.annotateWithBWAMEM(bwamemSamFileS, 5);
 
 //step 4
        topmFileS = "M:/GBStest/topm/small65536.V3.TOPM.h5";
        TagsOnPhysicalMapV3 topm3 = new TagsOnPhysicalMapV3(topmFileS);
        anno = new AnnotateTOPM (topm3);
        String PETOPMFileS = "M:/af/ptopm/PE.topm";
        anno.annotateWithPE(PETOPMFileS, 5);
 
    }
    
    public void mkSmallTOPM () {
        String inputFileS = "M:/GBSV3/topm/v3.bowtie2.bwa.blast.PE.topm.h5";
        String outputFileS = "M:/GBStest/topm/small65536.V3.TOPM.h5";
        TagsOnPhysicalMapV3 topm = new TagsOnPhysicalMapV3 (inputFileS);
        int[] index = new int[65536];
        for (int i = 0; i < index.length; i++) {
            index[i] = i;
        }
        topm.writeSubTOPM(outputFileS, index);
    }
    
    public void annotateTOPMHDF5WithAligner () {
        String inputFileS = "M:/GBStest/topm/ini.topm.h5";

        
        TagsOnPhysicalMapV3 topm = new TagsOnPhysicalMapV3(inputFileS);
        AnnotateTOPM anno = new AnnotateTOPM (topm);
        String bowtie2SamFileS = "M:/GBStest/alignment/bowtie2-K5.sam";
        anno.annotateWithBowtie2(bowtie2SamFileS, 5);
        String bwaSamFileS = "M:/GBStest/alignment/bwa-N5.sam";
        anno.annotateWithBWA(bwaSamFileS, 5);
        ////String blastFileS = "M:/GBStest/alignment/blast.m8.txt";
        ////anno.annotateWithBLAST(blastFileS, 5);
        String blastDirS = "M:/GBStest/alignment/blastOut/";
        anno.annotateWithBlastFromDir(blastDirS, 5);
        String bwaMemSamFileS = "M:/GBStest/alignment/bwa-mem-a.sam";
        anno.annotateWithBWAMEM(bwaMemSamFileS, 5);
        //String PETOPMFileS = "M:/af/ptopm/PE.topm";
        //anno.annotateWithPE(PETOPMFileS, 5);
    }
    
    public void initializeTOPMHDF5 () {
        String inputFileS = "M:/GBStest/tagCount/small.cnt";
        String outputFileS = "M:/GBStest/topm/ini.topm.h5";
        TagCounts tc = new TagCounts(inputFileS, FilePacking.Bit);
        TagsOnPhysicalMapV3.createFile(tc, outputFileS);
    }
    
    public void mkTOPMHDF5 () {
        String infileS = "M:/GBStest/topm/bowtie2.topm.bin";
        String outfileS = "M:/GBStest/topm/bowtie2.topm.h5";
        TagsOnPhysicalMap topm = new TagsOnPhysicalMap (infileS, false);
        TagsOnPhysMapHDF5.createFile(topm, outfileS, 4, 16);
    }
    
    public void mkTOPM () {
        String inputFileS = "M:/GBStest/alignment/bowtie2.sam";
        String outputFileS = "M:/GBStest/topm/bowtie2.topm.bin";
        new FeiUtils ().convertSAM2TOPM(inputFileS, outputFileS);
    }
    
    /**
     * Make FASTA file and do alignment using Blast, -m 8 -e 1e-10
     */
    public void mkFasta () {
        String inputFileS = "M:/GBStest/tagCount/small.cnt";
        String outputFileS = "M:/GBStest/alignment/small.fa";
        TagCounts tc = new TagCounts(inputFileS, FilePacking.Bit);
        tc.toFASTA(outputFileS);
    }
    
    /**
     * Make Fastq file and do alignment using Bowtie2 and BWA
     * Alignment files should be in the alignment folder
     */
    public void mkFastq () {
        String inputFileS = "M:/GBStest/tagCount/small.cnt";
        String outputFileS = "M:/GBStest/alignment/small.fq";
        new FeiUtils ().convertTagCount2Fastq(inputFileS, outputFileS);
    }
    
    public void mkSmallTagCount () {
        String inputFileS = "M:/GBStest/tagCount/434GFAAXX_s_4.cnt";
        String outputFileS = "M:/GBStest/tagCount/small.cnt";
        new FeiUtils ().mkSmallTagCountsFile(inputFileS, outputFileS, 10001, 500);
    }
    
    public void pipelinePE () {
        //parseQseq();
        //this.parseFastq();
        //this.checkPETagCounts(); //for checking, not in the pipeline.
        //this.mergePETagCounts();
        //this.contigPETagCounts();
        //this.mkPEstatistics();//for presentation, not included in pipeline.
        //this.mkFastaFileFromPE();
        //this.mkPEAlignment();
        
        //************************************
        //Deprecated
        //this.alignmentStep1();
        //this.alignmentStep2();
        //this.checkAlignmentOfPE(); //check alignment of longer sequence
        //************************************
    }
    
    public void checkAlignmentOfPE () {
        String TOGMFileS = "M:/pav/PhyGenMapping/v1.togm.txt";
        String topmFileS = "N:/Zea/AllZeaBuild_2.X/04_TOPM/2.6_production/02_MergedTOPM/AllZeaGBS_v2.6_MergedUnfiltProdTOPM_20130425.topm";
        String ptopmFileS = "M:/af/ptopm/merge.ptopm";
        String compareTableS = "E:/Research/af/alignmentImprovement/alignmentCompare.txt";
        FeiUtils fu = new FeiUtils ();
        fu.mkAlignmentCompareTable(TOGMFileS, topmFileS, ptopmFileS, compareTableS);
    }
    
    public void alignmentStep2 () {
        String fSamFileS = "M:/af/alignment/f.sam";
        String bSamFileS = "M:/af/alignment/b.sam";
        String contigSamFileS = "M:/af/alignment/c.sam";
        String infileS = "M:/af/mergePETagCounts/merge.con.pe.cnt";
        PETagCounts ptc = new PETagCounts (infileS, FilePacking.Bit);
        PETagsOnPhysicalMap topm = new PETagsOnPhysicalMap (ptc, fSamFileS, bSamFileS, contigSamFileS);
        String ptopmFileS = "M:/af/ptopm/merge.ptopm";
        topm.writeDistFile(ptopmFileS, FilePacking.Text);
    }
    
    public void alignmentStep1 () {
        String infileS = "M:/af/mergePETagCounts/merge.con.pe.cnt";
        PETagCounts ptc = new PETagCounts (infileS, FilePacking.Bit);
        String fastaFFileS = "M:/af/alignment/f.fasta.txt";
        String fastaBFileS = "M:/af/alignment/b.fasta.txt";
        String fastaCFileS = "M:/af/alignment/c.fasta.txt";
        ptc.mkFastaFile(fastaFFileS, fastaBFileS, fastaCFileS);
    }
    
    public void mkPEAlignment () {
        String fFastaFileS = "M:/af/alignment/f.fasta.txt";
        String bFastaFileS = "M:/af/alignment/b.fasta.txt";
        String fSamFileS = "M:/af/alignment/f.k5.sam.gz";
        String bSamFileS = "M:/af/alignment/b.k5.sam.gz";
        String a = null;
        PETagsOnPhysicalMapV3 ptopm = new PETagsOnPhysicalMapV3 (fFastaFileS, bFastaFileS, fSamFileS, bSamFileS);
        String PETOPM = "M:/af/ptopm/PE.topm";
        ptopm.writeBinaryFile(PETOPM);
        ptopm = new PETagsOnPhysicalMapV3(PETOPM);
    }
    
    public void mkFastaFileFromPE () {
        String infileS = "M:/af/mergePETagCounts/merge.con.pe.cnt";
        PETagCounts ptc = new PETagCounts (infileS, FilePacking.Bit);
        String fFastaFileS = "M:/af/alignment/f.fasta.txt";
        String bFastaFileS = "M:/af/alignment/b.fasta.txt";
        ptc.mkFastaFile(fFastaFileS, bFastaFileS);
    }
    
    public void mkPEstatistics () {
        String infileS = "M:/af/mergePETagCounts/merge.con.pe.cnt";
        PETagCounts ptc = new PETagCounts (infileS, FilePacking.Bit);
        System.out.println("Tag count is " + ptc.getTagCount());
        System.out.println("Read count is " + ptc.getTotalReadCount());
        System.out.println("Contig count is " + ptc.getContigCount());
        String staFileS = "E:/Research/af/PEstatistics/sta.txt";
        int maxCount = 1;
        for (int i = 0; i < ptc.getTagCount(); i++) {
            if (ptc.getReadCount(i) > maxCount) maxCount = ptc.getReadCount(i);
            //System.out.println(maxCount);            
        }
        System.out.println("Max read count is " + maxCount);
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(staFileS), 65536);
            bw.write("FLength\tBLength\tContigLength\tReadCount");
            bw.newLine();
            for (int i = 0; i < 10000; i++) {
                //if (ptc.getContigLengthInLong(i) == 0) continue; 
                bw.write(String.valueOf(ptc.getTagFLength(i))+"\t"+String.valueOf(ptc.getTagBLength(i))+"\t"+String.valueOf(ptc.getContigLength(i))+"\t"+String.valueOf(ptc.getReadCount(i)));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
     public void contigPETagCounts () {
        String infileS = "M:/af/mergePETagCounts/merge.pe.cnt";
        String outfileS = "M:/af/mergePETagCounts/merge.con.pe.cnt";
        String arguments = "-i " + infileS + " -o " + outfileS;
        String[] args = arguments.split(" ");
        ContigPETagCountPlugin m = new ContigPETagCountPlugin();
        m.setParameters(args);
		m.performFunction(null);
    }
    
    public void mergePETagCounts () {
        String inputDirS = "M:/af/PETagCounts/";
        String outfileS = "M:/af/mergePETagCounts/merge.pe.cnt";
        String arguments = "-i " + inputDirS + " -o " + outfileS;
        String[] args = arguments.split(" ");
        MergePETagCountPlugin m = new MergePETagCountPlugin();
        m.setParameters(args);
		m.performFunction(null); 
    }

    public void parseFastq () {
        String infile1 = "M:/af/Illumina/ImputationP15_1_1_fastq.txt";
        String infile2 = "M:/af/Illumina/ImputationP15_1_2_fastq.txt";
        
        String keyfile = "M:/af/key/ImputationP15_key.txt";
        String outputDirS = "M:/af/PETagCounts/";
        String arguments = "-iF " + infile1 + " -iB " + infile2 + " -k " + keyfile + " -e ApekI -l 8 -o " + outputDirS;
        String[] args = arguments.split(" ");
        FastqToPETagCountPlugin q = new FastqToPETagCountPlugin();
        q.setParameters(args);
		q.performFunction(null);     
    }
    
    public void checkPETagCounts () {
        String outputDirS = "M:/af/PETagCounts/";
        String txtDirS = "M:/af/text/";
        File[] files = new File (outputDirS).listFiles();
        for (int i = 0; i < files.length; i++) {
            PETagCounts p = new PETagCounts (files[i].getAbsolutePath(), FilePacking.Bit);
            String out = txtDirS + "/" + files[i].getName();
            p.writeDistFile(out, FilePacking.Text, 0);
        }
    }
     
    public void parseQseq () {
        String infile1 = "M:/af/Illumina/81546ABXX_8_1_qseq.txt";
        String infile2 = "M:/af/Illumina/81546ABXX_8_2_qseq.txt";
        String keyfile = "M:/af/key/81546ABXX_key.txt";
        String outputDirS = "M:/af/PETagCounts/";
        String arguments = "-iF " + infile1 + " -iB " + infile2 + " -k " + keyfile + " -e ApekI -o " + outputDirS;
        String[] args = arguments.split(" ");
        QseqToPETagCountPlugin q = new QseqToPETagCountPlugin();
        q.setParameters(args);
		q.performFunction(null);     
    }
    
}
