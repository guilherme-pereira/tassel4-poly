/*
 * TagsOnPhysMapHDF5
 */
package net.maizegenetics.gbs.maps;

import cern.colt.GenericSorting;
import ch.systemsx.cisd.hdf5.*;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import net.maizegenetics.gbs.maps.TagMappingInfoV3.Aligner;
import net.maizegenetics.gbs.tagdist.TagCounts;
import net.maizegenetics.gbs.tagdist.Tags;
import net.maizegenetics.gbs.tagdist.TagsByTaxa;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.gbs.util.BaseEncoder;
import net.maizegenetics.gbs.util.GBSHDF5Constants;
import org.apache.log4j.Logger;

/**
 * HDF5 version of TagsOnPhysical Map.  This is the preferred version of physical map as it uses less
 * memory, loads faster, and is more flexible with mapping positions.
 * <p>
 * Multiple mapping positions can be stored for each Tag.  For example, separate aligners could record their
 * positions in the {@link TagMappingInfoV3} objects.  Then the genetic mapping algorithm could be used to resolve,
 * which is the true mapping position.  MapPosition0 is used as the best mapping position, and used by the SNP caller.
 * The fields in {@link TagMappingInfoV3} may change for each aligner and genetic mapping, since part of the field might be missing different aligner and genetic mapping
 * <p>
 * TODO: createFile that includes a Locus filter, only exports positions on the same locus and position range
 * TODO: Resort map positions by quality (by model training)
 * 
 * 
 * @author Ed Buckler, Terry Casstevens, Fei Lu
 */
public class TagsOnPhysicalMapV3 extends AbstractTagsOnPhysicalMap implements TOPMInterface {

    private static final Logger myLogger = Logger.getLogger(TagsOnPhysicalMapV3.class);
    /**Shift 2^16*/
    private static final int BITS_TO_SHIFT_FOR_CHUNK = 16;
    /**65536 tags in a truck*/
    private static final int CHUNK_SIZE = 1 << BITS_TO_SHIFT_FOR_CHUNK;
    /**gzip format compression for TagMappingInfo*/
    private static HDF5GenericStorageFeatures genoFeatures = HDF5GenericStorageFeatures.createDeflation(5); //used by mapping object
    /**gzip format compression for other datasets*/
    private static HDF5IntStorageFeatures vectorFeatures = HDF5IntStorageFeatures.createDeflation(5); //used by vectors
    
    /**Number of physical positions from different aligner or aligner with different parameters*/
    protected int mappingNum = 0;
    /**Boolean value of if the file contains GM hypothesis test info*/
    protected boolean ifHasGM = false;
    /**Boolean value of if the file contains GM hypothesis test info*/
    protected boolean ifHasGMGW = false;
    /**Writer/Reader of TagsOnPhysicalMapV3*/
    private IHDF5Writer myHDF5 = null;
    /**Tag index (in whole tag list), for cachedTMI*/
    private int cachedTagIndex = -1;
    /**Tag mapping info index (in mappingNum), for cachedTMI*/
    private int cachedMapIndex = -1;
    /**Tag genetic mapping info index (in mappingNum), for cachedTGMI*/
    private int cachedGeneticMapIndex = -1;
    /**Current TMI*/
    private TagMappingInfoV3 cachedTMI = null;
    /**Current TGMI*/
    private TagGeneticMappingInfo cachedTGMI = null;
    /**Current TGMIGW*/
    private TagGeneticMappingInfo cachedTGMIGW = null;
    /**Chunk index where the tag belongs to. Max: TagCount>>BITS_TO_SHIFT_FOR_CHUNK+1*/
    private int cachedMappingChunkIndex = -1;
    /**Chunk index where the tag belongs to. Max: TagCount>>BITS_TO_SHIFT_FOR_CHUNK+1*/
    private int cachedGeneticMappingChunkIndex = -1;
    /**Chunk index where the tag belongs to. Max: TagCount>>BITS_TO_SHIFT_FOR_CHUNK+1*/
    private int cachedGeneticMappingGWChunkIndex = -1;
    /**All the pathes of mapping information, used to cache TMI block*/
    private String[] mapNames = null;
    /**All the pathes of genetic mapping information, used to cache TGMI block*/
    private String[] geneticMapNames = null;
    /**TMI chunk, 65536 * mappingNum */
    private TagMappingInfoV3[][] cachedTMIChunk = null;
    /**TGMI chunk, 65536 * mappingNum */
    private TagGeneticMappingInfo[][] cachedTGMIChunk = null;
     /**TGMI chunk, 65536 * mappingNum */
    private TagGeneticMappingInfo[] cachedTGMIGWChunk = null;
    /**Tag start index (in whole tag list) in current chunk*/
    private int mappingChunkStartTagIndex;
    /**Tag end index (in whole tag list) in current chunk, exclusive*/
    private int mappingChunkEndTagIndex;
    /**Tag start index (in whole tag list) in current chunk*/
    private int geneticMappingChunkStartTagIndex;
    /**Tag end index (in whole tag list) in current chunk, exclusive*/
    private int geneticMappingChunkEndTagIndex;
     /**Tag start index (in whole tag list) in current chunk*/
    private int geneticMappingGWChunkStartTagIndex;
    /**Tag end index (in whole tag list) in current chunk, exclusive*/
    private int geneticMappingGWChunkEndTagIndex;
    /**Current aligner being analyzed*/
    private Aligner currentAligner = null;
    /**mapping indices of currentAligner being analyzed*/
    private int[] mappingIndexOfAligner = null;
    private boolean cleanMap = true;
    private boolean cacheAllMappingBlocks = false;
    private HDF5CompoundType<TagMappingInfoV3> tmiType = null;
    private HDF5CompoundType<TagGeneticMappingInfo> tgmiType = null;
    private boolean hasDetailedMapping=false;
    
    int[] bestEndPosition;
    byte[] bestDivergence;
    byte[] bestMapP;
    byte[] bestDcoP;
    byte[] bestEvidence;
    byte[] bestMapIndex;
    /**
     * Initialize HDF5 TOPM from TagCounts file. The "MAXMAPPING"(set to 0), "TAGLENGTHINLONG", "tags" and "tagLength" are added.
     * @param inTags
     * @param newHDF5file 
     */
    public static void createFile (Tags inTags, String newHDF5file) {
        int tagLengthInLong = inTags.getTagSizeInLong();
        int tagCount = inTags.getTagCount();
        long[][] tags = new long[tagLengthInLong][tagCount];
        byte[] tagLength = new byte[tagCount];
        for (int i = 0; i < tagCount; i++) {
            long[] ct = inTags.getTag(i);
            for (int j = 0; j < tagLengthInLong; j++) {
                tags[j][i] = ct[j];
            }
            tagLength[i] = (byte) inTags.getTagLength(i);
        }
        IHDF5Writer h5 = null;
        try {
            myLogger.info("Creating HDF5 File: " + newHDF5file);
            IHDF5WriterConfigurator config = HDF5Factory.configure(new File(newHDF5file));
            config.overwrite();
            config.useUTF8CharacterEncoding();
            h5 = config.writer();
            h5.setIntAttribute(GBSHDF5Constants.ROOT, GBSHDF5Constants.MAXMAPPING, 0);
            h5.setBooleanAttribute(GBSHDF5Constants.ROOT, GBSHDF5Constants.IFHASGM, false);
            h5.setBooleanAttribute(GBSHDF5Constants.ROOT, GBSHDF5Constants.IFHASGMGW, false);
            h5.setIntAttribute(GBSHDF5Constants.ROOT, GBSHDF5Constants.TAGLENGTHINLONG, tagLengthInLong);
            h5.setIntAttribute(GBSHDF5Constants.ROOT, GBSHDF5Constants.TAGCOUNT, tagCount);
            h5.createLongMatrix(GBSHDF5Constants.TAGS, inTags.getTagSizeInLong(), tagCount, inTags.getTagSizeInLong(), tagCount, vectorFeatures);
            h5.writeLongMatrix(GBSHDF5Constants.TAGS, tags, vectorFeatures);
            tags=null;
            System.out.println("...Tags written");
            h5.createByteArray(GBSHDF5Constants.TAGLENGTH, tagCount, vectorFeatures);
            h5.writeByteArray(GBSHDF5Constants.TAGLENGTH, tagLength, vectorFeatures);
            tagLength=null;
            System.out.println("...Tags lengths written");
            System.gc();
            h5.flush();
            h5.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            h5.close();
            System.exit(1);
        }
    }
    
     /**
     * Constructor from a HDF5 TOPM file
     * @param theHDF5file 
     */
    public TagsOnPhysicalMapV3 (String theHDF5file) {
        System.out.println("Opening: " + theHDF5file);
        myHDF5 = HDF5Factory.open(theHDF5file);
        tmiType = myHDF5.compounds().getInferredType(TagMappingInfoV3.class);
        tgmiType = myHDF5.compounds().getInferredType(TagGeneticMappingInfo.class);
        this.tagLengthInLong = myHDF5.getIntAttribute(GBSHDF5Constants.ROOT, GBSHDF5Constants.TAGLENGTHINLONG);
        this.myNumTags = myHDF5.getIntAttribute(GBSHDF5Constants.ROOT, GBSHDF5Constants.TAGCOUNT);
        this.tags = myHDF5.readLongMatrix(GBSHDF5Constants.TAGS);
        this.tagLength = myHDF5.readByteArray(GBSHDF5Constants.TAGLENGTH);
        this.mappingNum = myHDF5.getIntAttribute(GBSHDF5Constants.ROOT, GBSHDF5Constants.MAXMAPPING);
        this.ifHasGM = myHDF5.getBooleanAttribute(GBSHDF5Constants.ROOT, GBSHDF5Constants.IFHASGM);
        this.ifHasGMGW = myHDF5.getBooleanAttribute(GBSHDF5Constants.ROOT, GBSHDF5Constants.IFHASGMGW);
        if (mappingNum != 0) {
            this.renameMapNames();
            this.getMappingInfo(0, 0);
        }
        if (ifHasGM) {
            this.renameGeneticMapNames();
            this.getGeneticMappingInfo(0, 0);
        }
        if (ifHasGMGW) {
            this.cacheGeneticMappingInfoGWChunk(0); 
            cachedTGMIGW = this.cachedTGMIGWChunk[0];
        }
       
        if(myHDF5.exists(GBSHDF5Constants.BEST_STRAND)) {
            bestStrand=myHDF5.readByteArray(GBSHDF5Constants.BEST_STRAND);
            bestChr= myHDF5.readIntArray(GBSHDF5Constants.BEST_CHR);
            bestStartPos= myHDF5.readIntArray(GBSHDF5Constants.BEST_STARTPOS);
            bestEndPosition= myHDF5.readIntArray(GBSHDF5Constants.BEST_ENDPOS);
            bestDivergence = myHDF5.readByteArray(GBSHDF5Constants.BEST_DIVERGENCE);
            bestMapP = myHDF5.readByteArray(GBSHDF5Constants.BEST_MAPP);
            bestDcoP = myHDF5.readByteArray(GBSHDF5Constants.BEST_DCOP);
            multimaps = myHDF5.readByteArray(GBSHDF5Constants.MULTIMAPS);
            bestEvidence = myHDF5.readByteArray(GBSHDF5Constants.BEST_EVIDENCE);
            bestMapIndex = myHDF5.readByteArray(GBSHDF5Constants.BEST_MAPINDEX);
            loadVariantsIntoMemory();
            populateChrAndVarPositions();
            initPhysicalSort();
        }    
        System.gc();
    }
    
    /**
     * Write text format map and genetic map for methods development
     * @param outputFileS 
     */
    public void writeTextMap (String tagCountFileS, String outputFileS) {
        TagCounts tc = new TagCounts(tagCountFileS, FilePacking.Byte);
        String[] title = {"AChr","AStartPos", "AEndPos", "Divergence", "Source", "Rank", "Score", "PValue", "GChr", "GPos", "SigSiteNum", "SigSiteRange"};
        try {
            BufferedWriter bw = new BufferedWriter (new FileWriter(outputFileS), 65536);
            bw.write("Tag\tTagLength\tTagCount\tGWGChr\tGWGPos\t");
            for (int i = 0; i < this.getMappingNum(); i++) {
                String mapIndexS = this.getThreeFigureString(i);
                for (int j = 0; j < title.length; j++) {
                    bw.write(title[j]+"-"+mapIndexS+"\t");
                }
            }
            bw.newLine();
            for (int i = 0; i < this.getTagCount(); i++) {
                long[] tag = this.getTag(i);
                bw.write(BaseEncoder.getSequenceFromLong(tag)+"\t"+String.valueOf(this.getTagLength(i))+"\t");
                int index = tc.getTagIndex(tag);
                bw.write(String.valueOf(tc.getReadCount(index))+"\t");
                TagMappingInfoV3 tmi;
                TagGeneticMappingInfo tgmi = this.getGeneticMappingInfoGW(i);
                bw.write(String.valueOf(tgmi.chromosome)+"\t"+String.valueOf(tgmi.position)+"\t");
                for (int j = 0; j < this.getMappingNum(); j++) {
                    tmi = this.getMappingInfo(i, j);
                    bw.write(String.valueOf(tmi.chromosome)+"\t");
                    bw.write(String.valueOf(tmi.startPosition)+"\t");
                    bw.write(String.valueOf(tmi.endPosition)+"\t");
                    bw.write(String.valueOf(tmi.divergence)+"\t");
                    bw.write(String.valueOf(tmi.mappingSource)+"\t");
                    bw.write(String.valueOf(tmi.mappingRank)+"\t");
                    bw.write(String.valueOf(tmi.mappingScore)+"\t");
                    tgmi = this.getGeneticMappingInfo(i, j);
                    bw.write(String.valueOf(tgmi.p)+"\t");
                    bw.write(String.valueOf(tgmi.chromosome)+"\t");
                    bw.write(String.valueOf(tgmi.position)+"\t");
                    bw.write(String.valueOf(tgmi.sigSiteNum)+"\t");
                    bw.write(String.valueOf(tgmi.sigSiteRange)+"\t");
                }
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            System.out.println(e.toString());
            System.exit(1);
        }
    }
    
    public void writeSubTOPM (String outputFileS, int[] tagIndex) {
        int tagLengthInLong = this.getTagSizeInLong();
        int tagCount = tagIndex.length;
        long[][] tags = new long[tagLengthInLong][tagCount];
        byte[] tagLength = new byte[tagCount];
        for (int i = 0; i < tagCount; i++) {
            long[] t = this.getTag(tagIndex[i]);
            for (int j = 0; j < tagLengthInLong; j++) {
                tags[j][i] = t[j];
            }
            tagLength[i] = (byte)this.getTagLength(tagIndex[i]);
        }
        IHDF5Writer h5 = null;
        try {
            myLogger.info("Creating HDF5 File: " + outputFileS);
            IHDF5WriterConfigurator config = HDF5Factory.configure(new File(outputFileS));
            config.overwrite();
            config.useUTF8CharacterEncoding();
            h5 = config.writer();
            h5.setIntAttribute(GBSHDF5Constants.ROOT, GBSHDF5Constants.MAXMAPPING, 0);
            h5.setBooleanAttribute(GBSHDF5Constants.ROOT, GBSHDF5Constants.IFHASGM, false);
            h5.setBooleanAttribute(GBSHDF5Constants.ROOT, GBSHDF5Constants.IFHASGMGW, false);
            h5.setIntAttribute(GBSHDF5Constants.ROOT, GBSHDF5Constants.TAGLENGTHINLONG, tagLengthInLong);
            h5.setIntAttribute(GBSHDF5Constants.ROOT, GBSHDF5Constants.TAGCOUNT, tagCount);
            h5.createLongMatrix(GBSHDF5Constants.TAGS, tagLengthInLong, tagCount, tagLengthInLong, tagCount, vectorFeatures);
            h5.writeLongMatrix(GBSHDF5Constants.TAGS, tags, vectorFeatures);
            tags=null;
            System.out.println("...Tags written");
            h5.createByteArray(GBSHDF5Constants.TAGLENGTH, tagCount, vectorFeatures);
            h5.writeByteArray(GBSHDF5Constants.TAGLENGTH, tagLength, vectorFeatures);
            tagLength=null;
            System.out.println("...Tags lengths written");
            System.gc();
            h5.flush();
            h5.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            h5.close();
            System.exit(1);
        }
        TagsOnPhysicalMapV3 topm = new TagsOnPhysicalMapV3(outputFileS);
        if (this.getMappingNum() != 0) {
            String[] mapNames = topm.creatTagMappingInfoDatasets(0, this.getMappingNum());
            TagMappingInfoV3[][] tmiBuffer = new TagMappingInfoV3[this.getMappingNum()][topm.getChunkSize()]; 
            int startTagIndex, endTagIndex, actualSize;
            for (int i = 0; i < topm.getChunkNum(); i++) {
                startTagIndex = i*topm.getChunkSize();
                endTagIndex = startTagIndex+topm.getChunkSize();
                actualSize = topm.getChunkSize();
                if (endTagIndex > topm.getTagCount()) {
                    actualSize = topm.getTagCount()-startTagIndex;
                    for (int j = actualSize; j < topm.getChunkSize(); j++) {
                        for (int k = 0; k < this.getMappingNum(); k++) {
                            tmiBuffer[k][j] = new TagMappingInfoV3();
                        }    
                    }
                }
                for (int j = 0; j < actualSize; j++) {
                    for (int k = 0; k < this.getMappingNum(); k++) {
                        tmiBuffer[k][j] = this.getMappingInfo(tagIndex[startTagIndex+j], k);
                    }
                }
                topm.writeTagMappingInfoDataSets(mapNames, tmiBuffer, i);
            }
            topm.setMappingNum(this.getMappingNum());
        }
        if (this.getIfHasGeneticMapping() == true) {
            String[] mapNames = topm.creatTagGeneticMappingInfoDatasets(0, this.getMappingNum());
            TagGeneticMappingInfo[][] tgmiBuffer = new TagGeneticMappingInfo[this.getMappingNum()][topm.getChunkSize()]; 
            int startTagIndex, endTagIndex, actualSize;
            for (int i = 0; i < topm.getChunkNum(); i++) {
                startTagIndex = i*topm.getChunkSize();
                endTagIndex = startTagIndex+topm.getChunkSize();
                actualSize = topm.getChunkSize();
                if (endTagIndex > topm.getTagCount()) {
                    actualSize = topm.getTagCount()-startTagIndex;
                    for (int j = actualSize; j < topm.getChunkSize(); j++) {
                        for (int k = 0; k < this.getMappingNum(); k++) {
                            tgmiBuffer[k][j] = new TagGeneticMappingInfo();
                        }    
                    }
                }
                for (int j = 0; j < actualSize; j++) {
                    for (int k = 0; k < this.getMappingNum(); k++) {
                        tgmiBuffer[k][j] = this.getGeneticMappingInfo(tagIndex[startTagIndex+j], k);
                    }
                }
                topm.writeTagGeneticMappingInfoDataSets(mapNames, tgmiBuffer, i);
            }
            topm.setIfHasGeneticMapping(this.getIfHasGeneticMapping());
        }
        if (this.getIfHasGeneticMappingGW() == true) {
            String dataSetName = topm.creatTagGeneticMappingInfoGWDataset();
            TagGeneticMappingInfo[] gmChunk;
            int startTagIndex, endTagIndex, actualSize;
            for (int i = 0; i < topm.getChunkNum(); i++) {
                gmChunk = new TagGeneticMappingInfo[topm.getChunkSize()];
                startTagIndex = i*topm.getChunkSize();
                endTagIndex = startTagIndex+topm.getChunkSize();
                actualSize = topm.getChunkSize();
                if (endTagIndex > topm.getTagCount()) {
                    actualSize = topm.getTagCount()-startTagIndex;
                    for (int j = actualSize; j < topm.getChunkSize(); j++) {
                        gmChunk[j] = new TagGeneticMappingInfo();  
                    }
                }
                for (int j = 0; j < actualSize; j++) {
                    gmChunk[j] = this.getGeneticMappingInfoGW(tagIndex[startTagIndex+j]);
                }
                topm.writeTagGeneticMappingInfoGWDataSet(dataSetName, gmChunk, i);
                if (i%100 == 0) System.out.println("Chunk " + i + "(index) with " + topm.getChunkSize() + " tags is annotated with genome wide genetic mapping");
            }
            topm.setIfHasGeneticMappingGW(true);
        }
    }
    
   
    
    /**
     * Creat datasets in HDF5 holding mapping information, which is used to annotate the TOPM with multiple alignment hypothesis
     * @param startTagIndex Start index of tag mapping information. This is essentially the current mappingNum
     * @param size the number of datasets which will be created
     * @return names of the datasets
     */
    public String[] creatTagMappingInfoDatasets (int startIndex, int size) {
        int chunkCount = this.getChunkNum(); 
        int chunkSize = this.getChunkSize();
        String[] dataSetNames = new String[size];
        for (int i = 0; i < size; i++) {
            dataSetNames[i] = GBSHDF5Constants.MAPBASE + this.getThreeFigureString(i+startIndex);
            myHDF5.compounds().createArray(dataSetNames[i], tmiType, chunkSize*chunkCount, chunkSize, genoFeatures);
            
        }
        System.out.println("Created new TagMappingInfo datasets. They are");
        for (int i = 0; i < dataSetNames.length; i++) {
            System.out.println(dataSetNames[i]);
        }
        return dataSetNames;
    }
    
    /**
     * Creat datasets in HDF5 holding genetic mapping information, which is used to test multiple alignment hypothesis
     * @param startTagIndex
     * @param size
     * @return 
     */
    public String[] creatTagGeneticMappingInfoDatasets (int startIndex, int size) {
        int chunkCount = this.getChunkNum(); 
        int chunkSize = this.getChunkSize();
        String[] dataSetNames = new String[size];
        for (int i = 0; i < size; i++) {
            dataSetNames[i] = GBSHDF5Constants.GENETICMAMMPING + this.getThreeFigureString(i+startIndex);
            myHDF5.compounds().createArray(dataSetNames[i], tgmiType, chunkSize*chunkCount, chunkSize, genoFeatures);
            
        }
        System.out.println("Created new TagGeneticMappingInfo datasets. They are");
        for (int i = 0; i < dataSetNames.length; i++) {
            System.out.println(dataSetNames[i]);
        }
        return dataSetNames;
    }
    
    /**
     * Creat dataset in HDF5 holding genome wide genetic mapping information, which is used to build training dataset to predict hypothesis genetic mapping
     * @return 
     */
    public String creatTagGeneticMappingInfoGWDataset () {
        int chunkCount = this.getChunkNum(); 
        int chunkSize = this.getChunkSize();
        myHDF5.compounds().createArray(GBSHDF5Constants.GENETICMAMMPINGGW, tgmiType, chunkSize*chunkCount, chunkSize, genoFeatures);
        System.out.println("Created new TagGeneticMappingInfoGW dataset: " + GBSHDF5Constants.GENETICMAMMPINGGW);
        return GBSHDF5Constants.GENETICMAMMPINGGW;
    }
    
    /**
     * Write best mapping positions from all hypotheses to HDF5 datasets, including strand, chr, pos and number of multimaps
     * @param bestStrand
     * @param bestChr
     * @param bestStartPos
     * @param multimaps 
     */
    public void writeBestMappingDataSets (byte[] bestStrand, int[] bestChr, int[] bestStartPos, int[] bestEndPos, byte[] bestDivergence, byte[] bestMapP, byte[] bestDcoP, byte[] multimaps, byte[] bestEvidence, byte[] bestMapIndex) {
        this.bestStrand = bestStrand;
        this.bestChr = bestChr;
        this.bestStartPos = bestStartPos;
        this.multimaps = multimaps;
        this.bestEvidence = bestEvidence;
        if (bestStrand.length != this.getTagCount() || bestChr.length != this.getTagCount() || bestStartPos.length != this.getTagCount() || multimaps.length != this.getTagCount()) {
            System.out.println("Size of best mapping arrays is not equal to tag count, program quits");
            System.exit(0);
        }
        myHDF5.createByteArray(GBSHDF5Constants.BEST_STRAND, this.getTagCount(), vectorFeatures);
        myHDF5.writeByteArray(GBSHDF5Constants.BEST_STRAND, bestStrand);
        myHDF5.createIntArray(GBSHDF5Constants.BEST_CHR, this.getTagCount(), vectorFeatures);
        myHDF5.writeIntArray(GBSHDF5Constants.BEST_CHR, bestChr);
        myHDF5.createIntArray(GBSHDF5Constants.BEST_STARTPOS, this.getTagCount(), vectorFeatures);
        myHDF5.writeIntArray(GBSHDF5Constants.BEST_STARTPOS, bestStartPos);
        myHDF5.createIntArray(GBSHDF5Constants.BEST_ENDPOS, this.getTagCount(), vectorFeatures);
        myHDF5.writeIntArray(GBSHDF5Constants.BEST_ENDPOS, bestEndPos);
        myHDF5.createByteArray(GBSHDF5Constants.BEST_DIVERGENCE, this.getTagCount(), vectorFeatures);
        myHDF5.writeByteArray(GBSHDF5Constants.BEST_DIVERGENCE, bestDivergence);
        myHDF5.createByteArray(GBSHDF5Constants.BEST_MAPP, this.getTagCount(), vectorFeatures);
        myHDF5.writeByteArray(GBSHDF5Constants.BEST_MAPP, bestMapP);
        myHDF5.createByteArray(GBSHDF5Constants.BEST_DCOP, this.getTagCount(), vectorFeatures);
        myHDF5.writeByteArray(GBSHDF5Constants.BEST_DCOP, bestDcoP);
        myHDF5.createByteArray(GBSHDF5Constants.MULTIMAPS, this.getTagCount(), vectorFeatures);
        myHDF5.writeByteArray(GBSHDF5Constants.MULTIMAPS, multimaps);
        myHDF5.createByteArray(GBSHDF5Constants.BEST_EVIDENCE, this.getTagCount(), vectorFeatures);
        myHDF5.writeByteArray(GBSHDF5Constants.BEST_EVIDENCE, bestEvidence);
        myHDF5.createByteArray(GBSHDF5Constants.BEST_MAPINDEX, this.getTagCount(), vectorFeatures);
        myHDF5.writeByteArray(GBSHDF5Constants.BEST_MAPINDEX, bestMapIndex);
        System.out.println("Best mapping positions from hypotheses are selected and saved to HDF5 TOPM");
    }
    
    /**
     * Write TMI buffer/chunk to HDF5 datasets, which is used to annotate the TOPM with multiple alignment hypothesis
     * @param dataSetNames
     * @param tmiChunk TMI chunk [dataSetNames.length]*[chunk_size]
     * @param chunkIndex index of this chunk
     */
    public void writeTagMappingInfoDataSets (String[] dataSetNames, TagMappingInfoV3[][] tmiChunk, int chunkIndex) {
        for (int i = 0; i < dataSetNames.length; i++) {   
            myHDF5.compounds().writeArrayBlock(dataSetNames[i], tmiType, tmiChunk[i], chunkIndex);
        }
    }
    
    /**
     * Write TGMI buffer/chunk to HDF5 datasets
     * @param dataSetNames
     * @param tmiChunk
     * @param chunkIndex 
     */
    public void writeTagGeneticMappingInfoDataSets (String[] dataSetNames, TagGeneticMappingInfo[][] tgmiChunk, int chunkIndex) {
        for (int i = 0; i < dataSetNames.length; i++) {   
            myHDF5.compounds().writeArrayBlock(dataSetNames[i], tgmiType, tgmiChunk[i], chunkIndex);
        }
    }
    
    /**
     * Write whole genome genetic mapping TGMI buffer/chunk to HDF5 datasets
     * @param dataSetName
     * @param tgmiChunk
     * @param chunkIndex 
     */
    public void writeTagGeneticMappingInfoGWDataSet (String dataSetName, TagGeneticMappingInfo[] tgmiChunk, int chunkIndex) {
        myHDF5.compounds().writeArrayBlock(dataSetName, tgmiType, tgmiChunk, chunkIndex);
    }
    
    /**
     * Set mappingNum attribute in HDF5
     * @param maxMapping 
     */
    public void setMappingNum (int mappingNum) {
        this.mappingNum = mappingNum;
        myHDF5.setIntAttribute(GBSHDF5Constants.ROOT, GBSHDF5Constants.MAXMAPPING, mappingNum);
        this.renameMapNames();
        System.out.println("TOPM maxMapping attibute was set to " + String.valueOf(mappingNum));
    }
    
    /**
     * Set if file has genetic mapping test result, true/false
     * @param value 
     */
    public void setIfHasGeneticMapping (boolean value) {
        this.ifHasGM = value;
        myHDF5.setBooleanAttribute(GBSHDF5Constants.ROOT, GBSHDF5Constants.IFHASGM, ifHasGM);
        this.renameGeneticMapNames();
        System.out.println("TOPM genetic mapping status was set to " + String.valueOf(value));
    }
    
    /**
     * Set if file has genome wide genetic mapping result, true/false
     * @param value 
     */
    public void setIfHasGeneticMappingGW (boolean value) {
        this.ifHasGMGW = value;
        myHDF5.setBooleanAttribute(GBSHDF5Constants.ROOT, GBSHDF5Constants.IFHASGMGW, ifHasGMGW);
        System.out.println("TOPM genome wide genetic mapping status was set to " + String.valueOf(value));
    }
    
    /**
     * Rename the mapNames based on the number of mapping
     */
    private void renameMapNames () {
        if (mappingNum == 0) {
            mapNames = null;
            return;
        }
        this.mapNames = new String[mappingNum];
        for (int i = 0; i < mappingNum; i++) {
            mapNames[i] = GBSHDF5Constants.MAPBASE + this.getThreeFigureString(i);
        }
    }
    
    private void renameGeneticMapNames () {
        if (this.getIfHasGeneticMapping() == false) {
            geneticMapNames = null;
            return;
        }
        geneticMapNames = new String[mappingNum];
        for (int i = 0; i < mappingNum; i++) {
            geneticMapNames[i] = GBSHDF5Constants.GENETICMAMMPING + this.getThreeFigureString(i);
        }
    }
    
    private boolean loadVariantsIntoMemory() {
        int howManyDef=0;
        int readBlock=this.getChunkSize();
        variantDefs=new byte[myNumTags][];
        variantOffsets=new byte[myNumTags][];
        if(!myHDF5.exists(GBSHDF5Constants.VARIANTDEF)) return false;
        for (int blockStep = 0; blockStep < myNumTags; blockStep+=readBlock) {
            int blockSize=(myNumTags-blockStep<readBlock)?myNumTags-blockStep:readBlock;
            byte[][] vd=myHDF5.readByteMatrixBlockWithOffset(GBSHDF5Constants.VARIANTDEF,blockSize,myMaxVariants,blockStep,0);
//            System.out.println(Arrays.toString(vd[0]));
            byte[][] vo=myHDF5.readByteMatrixBlockWithOffset(GBSHDF5Constants.VARIANTPOSOFF,blockSize,myMaxVariants,blockStep,0);
            for (int j = 0; j < blockSize; j++) {
                int cnt=0;
                for (byte bs : vd[j]) {if (bs!=TOPMInterface.BYTE_MISSING) cnt++;}
                if(cnt==0) continue;
                byte[] vdReDim=new byte[cnt];
                byte[] voReDim=new byte[cnt];
                for (int i = 0; i < cnt; i++) {
                    vdReDim[i]=vd[j][i];
                    voReDim[i]=vo[j][i];
                    howManyDef++;
                }
                variantDefs[blockStep+j]=vdReDim;
                variantOffsets[blockStep+j]=voReDim;
            }
            
            //byte[] vd=myHDF5.readByteArrayBlockWithOffset(null, i, i)
        }
        System.out.println("Real Variant Defs:"+howManyDef);
        return true;
    }
    
    private static boolean writeVariantsToHDF5(IHDF5Writer aHDF5, AbstractTagsOnPhysicalMap aTOPM) {
        int howManyDef=0;
        int readBlock=4096*16;
        
        int myNumTags=aTOPM.myNumTags;
        int myMaxVariants=aTOPM.myMaxVariants;
        aHDF5.createByteMatrix(GBSHDF5Constants.VARIANTDEF, myNumTags, myMaxVariants);
        aHDF5.createByteMatrix(GBSHDF5Constants.VARIANTPOSOFF, myNumTags, myMaxVariants);
//        variantDefs=new byte[myNumTags][];
//        variantOffsets=new byte[myNumTags][];
        if(!aHDF5.exists(GBSHDF5Constants.VARIANTDEF)) return false;
        byte[][] vd=new byte[readBlock][myMaxVariants];
        byte[][] vo=new byte[readBlock][myMaxVariants];
        for (int blockStep = 0; blockStep < myNumTags; blockStep+=readBlock) {
            int blockSize=(myNumTags-blockStep<readBlock)?myNumTags-blockStep:readBlock;
            vd=new byte[blockSize][myMaxVariants];
            vo=new byte[blockSize][myMaxVariants];
            for (int j = 0; j < blockSize; j++) {
                for (int v = 0; v < vo[0].length; v++) {
                    vd[j][v]=aTOPM.getVariantDef(blockStep+j, v);
                    vo[j][v]=aTOPM.getVariantPosOff(blockStep+j, v);
                }
            }
            aHDF5.writeByteMatrixBlockWithOffset(GBSHDF5Constants.VARIANTDEF,vd,blockStep,0);
            aHDF5.writeByteMatrixBlockWithOffset(GBSHDF5Constants.VARIANTPOSOFF,vo,blockStep,0);
            //byte[] vd=myHDF5.readByteArrayBlockWithOffset(null, i, i)
        }
        System.out.println("Real Variant Defs:"+howManyDef);
        return true;
    }
    
    /**
     * Load mapping information in a chunk to memory and reset tag start and end index
     * @param chunkIndex 
     */
    private synchronized void cacheMappingInfoChunk (int chunkIndex) {
        this.cachedTMIChunk = new TagMappingInfoV3[this.getMappingNum()][this.getChunkSize()];
        for (int i = 0; i < mappingNum; i++) {
            cachedTMIChunk[i] = myHDF5.compounds().readArrayBlock(mapNames[i], tmiType, this.getChunkSize(), chunkIndex);
        }
        this.cachedMappingChunkIndex = chunkIndex;
        this.mappingChunkStartTagIndex = chunkIndex*this.getChunkSize();
        this.mappingChunkEndTagIndex = mappingChunkStartTagIndex+this.getChunkSize();
        if (mappingChunkEndTagIndex > this.getTagCount()) mappingChunkEndTagIndex = this.getTagCount();
    }
    
    /**
     * Load genetic mapping information in a chunk to memory and reset tag start and end index
     * @param chunkIndex 
     */
    private synchronized void cacheGeneticMappingInfoChunk (int chunkIndex) {
        this.cachedTGMIChunk = new TagGeneticMappingInfo[this.getMappingNum()][this.getChunkSize()];
        for (int i = 0; i < mappingNum; i++) {
            cachedTGMIChunk[i] = myHDF5.compounds().readArrayBlock(geneticMapNames[i], tgmiType, this.getChunkSize(), chunkIndex);
        }
        this.cachedGeneticMappingChunkIndex = chunkIndex;
        this.geneticMappingChunkStartTagIndex = chunkIndex*this.getChunkSize();
        this.geneticMappingChunkEndTagIndex = geneticMappingChunkStartTagIndex+this.getChunkSize();
        if (geneticMappingChunkEndTagIndex > this.getTagCount()) geneticMappingChunkEndTagIndex = this.getTagCount();
    }
    /**
     * Load genome wide genetic mapping information in a chunk to memory and reset tag start and end index
     * @param chunkIndex 
     */
    private synchronized void cacheGeneticMappingInfoGWChunk (int chunkIndex) {
        this.cachedTGMIGWChunk = new TagGeneticMappingInfo[this.getChunkSize()];
        cachedTGMIGWChunk = myHDF5.compounds().readArrayBlock(GBSHDF5Constants.GENETICMAMMPINGGW, tgmiType, this.getChunkSize(), chunkIndex);
        this.cachedGeneticMappingGWChunkIndex = chunkIndex;
        this.geneticMappingGWChunkStartTagIndex = chunkIndex*this.getChunkSize();
        this.geneticMappingGWChunkEndTagIndex = geneticMappingGWChunkStartTagIndex+this.getChunkSize();
        if (geneticMappingGWChunkEndTagIndex > this.getTagCount()) geneticMappingGWChunkEndTagIndex = this.getTagCount();
    }
    
    /**
     * Update current cachedTMI
     * @param tagIndex
     * @param mapIndex 
     */
    private synchronized void cacheMappingInfo(int tagIndex, int mapIndex) {
        if (tagIndex == cachedTagIndex && mapIndex == cachedMapIndex) {
            return;
        }
        int chunkIndex = tagIndex >> BITS_TO_SHIFT_FOR_CHUNK;
        if (chunkIndex != this.cachedMappingChunkIndex) {
            this.cacheMappingInfoChunk(chunkIndex);
        }
        cachedTMI = this.cachedTMIChunk[mapIndex][tagIndex%this.getChunkSize()];
        cachedTagIndex = tagIndex;
        cachedMapIndex = mapIndex;
    }
    
    /**
     * Update current cachedTGMI
     * @param tagIndex
     * @param geneticMapIndex 
     */
    private synchronized void cacheGeneticMappingInfo(int tagIndex, int geneticMapIndex) {
        if (tagIndex == cachedTagIndex && geneticMapIndex == cachedGeneticMapIndex) {
            return;
        }
        int chunkIndex = tagIndex >> BITS_TO_SHIFT_FOR_CHUNK;
        if (chunkIndex != this.cachedGeneticMappingChunkIndex) {
            this.cacheGeneticMappingInfoChunk(chunkIndex); 
        }
        cachedTGMI = this.cachedTGMIChunk[geneticMapIndex][tagIndex%this.getChunkSize()];
        cachedTagIndex = tagIndex;
        cachedGeneticMapIndex = geneticMapIndex;
    }

    /**
     * Update current cachedTGMIGW (In memory)
     * @param tagIndex 
     */
    private synchronized void cacheGeneticMappingInfoGW(int tagIndex) {
        if (tagIndex == cachedTagIndex) {
            return;
        }
        int chunkIndex = tagIndex >> BITS_TO_SHIFT_FOR_CHUNK;
        if (chunkIndex != this.cachedGeneticMappingGWChunkIndex) {
            this.cacheGeneticMappingInfoGWChunk(chunkIndex); 
        }
        cachedTGMIGW = this.cachedTGMIGWChunk[tagIndex%this.getChunkSize()];
        cachedTagIndex = tagIndex;
    }

    /**
     * Return the total number of chunks
     * @return 
     */
    public int getChunkNum () {
        int num = this.getTagCount()/this.getChunkSize();
        if (this.getTagCount()% this.getChunkSize() == 0) return num;
        else return num+1;
    }
    
    /**
     * Return the chunk size (Number of tags in a chunk)
     * @return 
     */
    public int getChunkSize () {
        return this.CHUNK_SIZE;
    }
    
    public void getFileReadyForClosing() {
//        writeCachedVariantDefs();
//        writeCachedVariantOffsets();
//        saveCacheBackToFile();
    }

    /**
     * Return number of mapping result
     * @return 
     */
    public int getMappingNum () {
        return this.mappingNum; 
    }
    
    /**
     * Return if the file has genetic mapping test result
     * @return 
     */
    public boolean getIfHasGeneticMapping () {
        return this.ifHasGM;
    }
    
    public boolean getIfHasGeneticMappingGW () {
        return this.ifHasGMGW;
    }
    /**
     * Calculate TMI indices of an aligner
     * @param alignerName 
     */
    private synchronized void calMappingIndicesOfAligner (Aligner alignerName) {
        byte mappingSourceValue = alignerName.getValue();
        ArrayList<Integer> l = new ArrayList();
        for (int i = 0; i < cachedTMIChunk.length; i++) {
            if (cachedTMIChunk[i][0].mappingSource == mappingSourceValue) l.add(i);
        }
        mappingIndexOfAligner = new int[l.size()];
        for (int i = 0; i < mappingIndexOfAligner.length; i++) {
            mappingIndexOfAligner[i] = l.get(i);
        }
        this.currentAligner = alignerName;
    }
    
    /**
     * Return unique mapping (chr and startPosition) from an aligner, return null if it has multiple equally good position or doesn't align. This is used to block positions for genetic mapping
     * @param tagIndex
     * @param alignerName
     * @return 
     */
    public synchronized int[] getUniqueMappingOfAligner (int tagIndex, Aligner alignerName) {
        if (alignerName != this.currentAligner) this.calMappingIndicesOfAligner(alignerName);
        ArrayList<Integer> l = new ArrayList();
        int cnt = 0;
        for (int i = 0; i < mappingIndexOfAligner.length; i++) {
            TagMappingInfoV3 tempTMI = this.getMappingInfo(tagIndex, mappingIndexOfAligner[i]);
            if (tempTMI.mappingRank == 0) {
                if (tempTMI.chromosome < 0) return null;
                cnt++;
                if (cnt > 1) return null;
                l.add(tempTMI.chromosome);
                l.add(tempTMI.startPosition);
            }
        }
        if (l.isEmpty()) return null;
        int[] chrPos = new int[2];
        chrPos[0] = l.get(0);
        chrPos[1] = l.get(1);
        return chrPos;
    }
    
    @Override
    public int addVariant(int tagIndex, byte offset, byte base) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    public int getBestMapIndex (int tagIndex) {
        if(this.bestStrand == null) throw new IllegalStateException("Best mapping not present");
        return bestMapIndex[tagIndex];
    }
    
    /**
     * 
     * @param tagIndex
     * @return 
     */
    @Override
    public byte getDcoP(int tagIndex) {
        if(this.bestStrand == null) throw new IllegalStateException("Best mapping not present");
        return bestDcoP[tagIndex];
    }
    
    /**
     * @param tagIndex
     * @return 
     */
    @Override
    public byte getStrand (int tagIndex) {
        if(this.bestStrand == null) throw new IllegalStateException("Best mapping not present");
        return bestStrand[tagIndex];
    }
    
    /**
     * Blast doesn't have divergence, so it always return Byte.MIN_VALUE of Blast hits
     * @param tagIndex
     * @return 
     */
    @Override
    public byte getDivergence(int tagIndex) {
        if(this.bestStrand == null) throw new IllegalStateException("Best mapping not present");
        return bestDivergence[tagIndex];
    }

    /**
     * 
     * @param tagIndex
     * @return 
     */
    @Override
    public int getStartPosition (int tagIndex) {
        if(this.bestStrand == null) throw new IllegalStateException("Best mapping not present");
        return bestStartPos[tagIndex];
    }
    
    /**
     * EndPosition of PEEnd1 is probably not the EndPosition of the tag
     * @param index
     * @return 
     */
    @Override
    public int getEndPosition(int tagIndex) {
        if(this.bestStrand == null) throw new IllegalStateException("Best mapping not present");
        return bestEndPosition[tagIndex];
    }
    
    /**
     * Return the evidence
     * @param tagIndex
     * @return 
     */
    public byte getEvidence (int tagIndex) {
        if(this.bestStrand == null) throw new IllegalStateException("Best mapping not present");
        return bestEvidence[tagIndex];
    }
    
    @Override
    public byte getMapP(int tagIndex) {
        if(this.bestStrand == null) throw new IllegalStateException("Best mapping not present");
        return bestMapP[tagIndex];
    }
    
    /**
     * 
     * @param tagIndex
     * @return 
     */
    @Override
    public byte getMultiMaps (int tagIndex) {
        if(this.bestStrand == null) throw new IllegalStateException("Best mapping not present");
        return multimaps[tagIndex];
    }
    
    /**
     * Return the map indices of an aligner
     * @param alignerName
     * @return 
     */
    public int[] getMappingIndicesOfAligner (Aligner alignerName) {
        byte value = alignerName.getValue();
        ArrayList<Integer> list = new ArrayList();
        int tagIndex = 0;
        while (list.isEmpty()) {
            for (int i = 0; i < this.getMappingNum(); i++) {
                byte c = this.cachedTMIChunk[i][tagIndex].mappingSource;
                if (c == value) {
                    list.add(i);
                }
            }
            tagIndex++;
        }
        if (list.isEmpty()) return null;
        int[] indices = new int[list.size()];
        for (int i = 0; i < indices.length; i++) {
            indices[i] = list.get(i);
        }
        return indices;
    }
    
    @Override
    public int[] getPositionArray(int tagIndex) {
        int[] r = {this.bestChr[tagIndex], this.bestStrand[tagIndex], this.bestStartPos[tagIndex]};
        return r;
    }

    @Override
    public int getReadIndexForPositionIndex(int posIndex) {
        return indicesOfSortByPosition[posIndex];
    }
    
    /**
     * Return tag mapping information of a tag in one map
     * @param tagIndex
     * @param mapIndex
     * @return 
     */
    public TagMappingInfoV3 getMappingInfo (int tagIndex, int mapIndex) {
        this.cacheMappingInfo(tagIndex, mapIndex);
        return this.cachedTMI;
    }
    
    /**
     * Return mapping information of a while chunk, avoid issues while multi threads are trying to get TMI info, specifically for hypothesis genetic mapping
     * @param tagIndex
     * @return 
     */
    public TagMappingInfoV3[][] getMappingInfoChunk (int tagIndex) {
        if (this.getMappingNum() == 0) return null;
        this.getMappingInfo(tagIndex, 0);
        return this.cachedTMIChunk;
    }
    
    /**
     * Return tag genetic mapping information of a tag in one genetic map
     * @param tagIndex
     * @param geneticMapIndex
     * @return 
     */
    public TagGeneticMappingInfo getGeneticMappingInfo (int tagIndex, int geneticMapIndex) {
        this.cacheGeneticMappingInfo(tagIndex, geneticMapIndex);
        return this.cachedTGMI;
    }
    
    public TagGeneticMappingInfo getGeneticMappingInfoGW (int tagIndex) {
        this.cacheGeneticMappingInfoGW(tagIndex);
        return this.cachedTGMIGW;
    }
    
    /**
     * Example: Change 1 to 001, used in map names
     * @param number
     * @return 
     */
    private String getThreeFigureString (int number) {
        String s = String.valueOf(number);
        int length = s.length();
        for (int i = 0; i < 3-length; i++) {
            s = "0"+s;
        }
        return s;
    }
    
    @Override
    public int[] getUniquePositions(int chromosome) {
        if(myUniquePositions==null) populateChrAndVarPositions();
        return myUniquePositions[chromosome];
    }

    public void setMultimaps(int index, byte multimaps) {
        this.multimaps[index] = multimaps;
    }

    @Override
    public void setChromoPosition(int index, int chromosome, byte strand, int positionMin, int positionMax) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void setDivergence(int index, byte divergence) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void setMapP(int index, byte mapP) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void setMapP(int index, double mapP) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public synchronized void setVariantDef(int tagIndex, int variantIndex, byte def) {
        myHDF5.writeByteMatrixBlockWithOffset(GBSHDF5Constants.VARIANTDEF, new byte[][]{{def},}, tagIndex, variantIndex);
        variantDefs[tagIndex][variantIndex]=def;
    }

    @Override
    public synchronized void setVariantPosOff(int tagIndex, int variantIndex, byte offset) {
        myHDF5.writeByteMatrixBlockWithOffset(GBSHDF5Constants.VARIANTPOSOFF, new byte[][]{{offset},}, tagIndex, variantIndex);
        variantOffsets[tagIndex][variantIndex]=offset;
    }
    
    /**
     * Preferred method for setting variant information
     * @param tagIndex
     * @param defAndOffset Two dimension [0=definition, 1=offset][upto 16 bytes for each SNP]
     */
    public synchronized void setAllVariantInfo(int tagIndex, byte[][] defAndOffset) {
        myHDF5.writeByteMatrixBlockWithOffset(GBSHDF5Constants.VARIANTDEF, new byte[][]{defAndOffset[0]}, tagIndex, 0);
        myHDF5.writeByteMatrixBlockWithOffset(GBSHDF5Constants.VARIANTPOSOFF, new byte[][]{defAndOffset[1]}, tagIndex, 0);
        variantDefs[tagIndex]=defAndOffset[1];
        variantOffsets[tagIndex]=defAndOffset[1];
    }
    
    public long sortTable(boolean byHaplotype) {
        System.out.print("Starting Read Table Sort ...");
        if (byHaplotype == false) {
            //TODO change the signature at some time
            System.out.print("ERROR:  Position sorting has been eliminated ...");
            return -1;
        }
        long time = System.currentTimeMillis();
        GenericSorting.quickSort(0, tags[0].length, this, this);
        long totalTime = System.currentTimeMillis() - time;
        System.out.println("Done in " + totalTime + "ms");
        initPhysicalSort();
        return totalTime;
    }

    @Override
    public void clearVariants() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
}
