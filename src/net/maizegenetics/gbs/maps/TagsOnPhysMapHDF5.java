/*
 * TagsOnPhysMapHDF5
 */
package net.maizegenetics.gbs.maps;

import ch.systemsx.cisd.hdf5.*;
import java.io.File;
import java.util.Arrays;
import net.maizegenetics.gbs.util.GBSHDF5Constants;
import org.apache.log4j.Logger;

/**
 * HDF5 version of TagsOnPhysical Map.  This is the preferred version of physical map as it uses less
 * memory, loads faster, and is more flexible with mapping positions.
 * <p>
 * Multiple mapping positions can be stored for each Tag.  For example, separate aligners could record their
 * positions in the {@link TagMappingInfo} objects.  Then the genetic mapping algorithm could be used to resolve,
 * which is the true mapping position.  MapPosition0 is used as the best mapping position, and used by the SNP caller.
 *
 * <p>
 * TODO: createFile - needs to instantiate a TOPM just using a Tag Object
 * TODO: createFile that includes a Locus filter, only exports positions on the same locus and position range
 * TODO: Resort map positions by quality
 * 
 * 
 * @author Ed Buckler, Terry Casstevens
 */
public class TagsOnPhysMapHDF5 extends AbstractTagsOnPhysicalMap implements TOPMInterface {

    private static final Logger myLogger = Logger.getLogger(TagsOnPhysMapHDF5.class);
    private static final int NUM_UNITS_TO_CACHE_ON_GET = 64;
    private static final int BITS_TO_SHIFT_FOR_CHUNK = 16;
    private static final int CHUNK_SIZE = 1 << BITS_TO_SHIFT_FOR_CHUNK;
    private static HDF5GenericStorageFeatures genoFeatures = HDF5GenericStorageFeatures.createDeflation(5); //used by mapping object
    private static HDF5IntStorageFeatures vectorFeatures = HDF5IntStorageFeatures.createDeflation(5); //used by vectors
    
    private int maxMapping = 4;
    private IHDF5Writer myHDF5 = null;
    private int cachedMappingIndex = -1;
    private TagMappingInfo cachedTMI = null;
    private int cachedMappingBlock = -1;
    private TagMappingInfo[][] cachedTMIBlock = null;
    private boolean cleanMap = true;
    private boolean cacheAllMappingBlocks = false;
    private HDF5CompoundType<TagMappingInfo> tmiType = null;
    private boolean hasDetailedMapping=false;

    public static void createFile(AbstractTagsOnPhysicalMap inTags, String newHDF5file, int maxMapping, int maxVariants) {
        int tagLengthInLong = inTags.getTagSizeInLong();
        int tagCount = inTags.getTagCount();
     //   tagCount=tagCount;
        System.gc();
        long[][] tags;
        byte[] tagLength;
        if(inTags instanceof AbstractTagsOnPhysicalMap) {
            tags=((AbstractTagsOnPhysicalMap)inTags).getTagsArray();
            tagLength=((AbstractTagsOnPhysicalMap)inTags).getTagLengthArray();
        } else {
            tags = new long[tagLengthInLong][tagCount];
            tagLength = new byte[tagCount];
            for (int i = 0; i < tagCount; i++) {
                long[] ct = inTags.getTag(i);
                for (int j = 0; j < tagLengthInLong; j++) {
                    tags[j][i] = ct[j];
                }
                tagLength[i] = (byte) inTags.getTagLength(i);
            }
        }
        IHDF5Writer h5 = null;
        try {
            myLogger.info("Creating HDF5 File: " + newHDF5file);
            System.out.println("Creating HDF5 File: " + newHDF5file);
            IHDF5WriterConfigurator config = HDF5Factory.configure(new File(newHDF5file));
            config.overwrite();
            // config.dontUseExtendableDataTypes();
            config.useUTF8CharacterEncoding();
            h5 = config.writer();
            h5.setIntAttribute(GBSHDF5Constants.ROOT, GBSHDF5Constants.TAGCOUNT, tagCount);
            h5.setIntAttribute(GBSHDF5Constants.ROOT, GBSHDF5Constants.MAXVARIANTS, maxVariants);
            h5.setIntAttribute(GBSHDF5Constants.ROOT, GBSHDF5Constants.MAXMAPPING, maxMapping);
            h5.setIntAttribute(GBSHDF5Constants.ROOT, GBSHDF5Constants.TAGLENGTHINLONG, tagLengthInLong);
            //create tag matrix
            h5.createLongMatrix(GBSHDF5Constants.TAGS, inTags.getTagSizeInLong(), tagCount, inTags.getTagSizeInLong(), tagCount, vectorFeatures);
            h5.writeLongMatrix(GBSHDF5Constants.TAGS, tags, vectorFeatures);
            tags=null;
            System.gc();
            System.out.println("...Tags written");
            h5.createByteArray(GBSHDF5Constants.TAGLENGTH, tagCount, vectorFeatures);
            h5.writeByteArray(GBSHDF5Constants.TAGLENGTH, tagLength, vectorFeatures);
            tagLength=null;
            System.out.println("...Tags lengths written");
            
            //need to create branch between instance of Tags and TOPM
            
            
            //creating fast bestChr and bestPosition
            byte[] mmOut=new byte[tagCount];
            byte[] strandOut=new byte[tagCount];
            int[] chrOut=new int[tagCount];
            int[] posOut=new int[tagCount];
            
            for (int i = 0; i < tagCount; i++) {
                mmOut[i]=inTags.getMultiMaps(i);
                strandOut[i]=inTags.getStrand(i);
                chrOut[i]=inTags.getChromosome(i);
                posOut[i]=inTags.getStartPosition(i);
            }
            h5.createByteArray(GBSHDF5Constants.MULTIMAPS, tagCount);
            h5.writeByteArray(GBSHDF5Constants.MULTIMAPS, mmOut, vectorFeatures);
            h5.createByteArray(GBSHDF5Constants.BEST_STRAND, tagCount);
            h5.writeByteArray(GBSHDF5Constants.BEST_STRAND, strandOut, vectorFeatures);
            h5.createIntArray(GBSHDF5Constants.BEST_CHR, tagCount, vectorFeatures);
            h5.writeIntArray(GBSHDF5Constants.BEST_CHR, chrOut, vectorFeatures);
            h5.createIntArray(GBSHDF5Constants.BEST_STARTPOS, tagCount, vectorFeatures);
            h5.writeIntArray(GBSHDF5Constants.BEST_STARTPOS, posOut, vectorFeatures);
            mmOut=null;
            strandOut=null;
            chrOut=null;
            posOut=null;
            System.gc();
            System.out.println("...multimapping, strand, chr, position  written");
                        
            HDF5CompoundType<TagMappingInfo> tmiType = h5.compounds().getInferredType(TagMappingInfo.class);
            System.out.println("Chunk Size for Tags: " + CHUNK_SIZE);
            int numOfChunks = tagsToChunks(tagCount);
            System.out.println("Number of Chunks: " + numOfChunks);
            int numTagsPadded = numOfChunks * CHUNK_SIZE;
            for (int mi = 0; mi < 1; mi++) {
 //           for (int mi = 0; mi < maxMapping; mi++) {
 //               h5.compounds().createArray("map" + mi, tmiType, numTagsPadded, CHUNK_SIZE);
                h5.compounds().createArray(GBSHDF5Constants.MAPBASE + mi, tmiType, numTagsPadded, CHUNK_SIZE);
                TagMappingInfo[] thTMI = new TagMappingInfo[CHUNK_SIZE];
                int block = 0;
                for (int i = 0; i < numTagsPadded; i++) {
                    if ((mi == 0) && (i < tagCount)) {
                        thTMI[i % CHUNK_SIZE] = (new TagMappingInfo(inTags.getChromosome(i),
                                inTags.getStrand(i), inTags.getStartPosition(i),
                                inTags.getEndPosition(i), inTags.getDivergence(i)));
                    } else {
                        thTMI[i % CHUNK_SIZE] = new TagMappingInfo();
                    }
                    if ((i + 1) % CHUNK_SIZE == 0) {
//                        h5.compounds().writeArrayBlock("map" + mi, tmiType, thTMI, block);
                        h5.compounds().writeArrayBlock(GBSHDF5Constants.MAPBASE + mi, tmiType, thTMI, block);
                        thTMI = new TagMappingInfo[CHUNK_SIZE];
                        block++;
                        System.out.println("Tag locations written: " + (i + 1));
                    }
                }
                System.out.println("...map" + mi + " positions written");
            }
            
            h5.createByteMatrix(GBSHDF5Constants.VARIANTDEF, tagCount, maxVariants);
            int numVariants = inTags.getMaxNumVariants();
            if (numVariants > maxVariants) {
                throw new IllegalArgumentException("TagsOnPhysMapHDF5: createFile: max variants can't be less than original TOPM Variant Defs: " + numVariants);
            }
            writeVariantsToHDF5(h5, inTags, maxVariants);
            System.out.println("Variant offsets written");

        } finally {
            try {
                h5.close();
            } catch (Exception e) {
                // do nothing
            }
        }
    }

    private static int tagsToChunks(long tags) {
        return (int) (((tags - 1) >>> BITS_TO_SHIFT_FOR_CHUNK) + 1);
    }

    public TagsOnPhysMapHDF5(String filename) {
        this(filename, true);
    }
   
    public TagsOnPhysMapHDF5(String theHDF5file, boolean cacheAllMappingBlocks) {
        myMaxVariants = 16;
        this.cacheAllMappingBlocks = cacheAllMappingBlocks;
        System.out.println("Opening :" + theHDF5file);
        myHDF5 = HDF5Factory.open(theHDF5file);
        hasDetailedMapping=myHDF5.exists(GBSHDF5Constants.MAPBASE+"0");
        myNumTags = myHDF5.getIntAttribute(GBSHDF5Constants.ROOT, GBSHDF5Constants.TAGCOUNT);
        this.tagLengthInLong = myHDF5.getIntAttribute(GBSHDF5Constants.ROOT, GBSHDF5Constants.TAGLENGTHINLONG);
        this.tags = myHDF5.readLongMatrix(GBSHDF5Constants.TAGS);
        this.tagLength = myHDF5.readByteArray(GBSHDF5Constants.TAGLENGTH);
        this.multimaps = myHDF5.readByteArray(GBSHDF5Constants.MULTIMAPS);
        this.maxMapping = myHDF5.getIntAttribute(GBSHDF5Constants.ROOT, GBSHDF5Constants.MAXMAPPING);
        tmiType = myHDF5.compounds().getInferredType(TagMappingInfo.class);
        cachedTMIBlock = new TagMappingInfo[maxMapping][];
        if(hasDetailedMapping) cacheMappingInfo(0);
        if(!myHDF5.exists(GBSHDF5Constants.BEST_STRAND)) {populateBestMappings();}
        else {
            bestStrand=myHDF5.readByteArray("bestStrand");
            bestChr= myHDF5.readIntArray("bestChr");
            bestStartPos= myHDF5.readIntArray("bestStartPos");
        }
        loadVariantsIntoMemory();
        System.out.println(theHDF5file + " read with tags:" + myNumTags);
        populateChrAndVarPositions();
        initPhysicalSort();
        System.gc();
    }
    
    private boolean populateBestMappings() {
        bestStrand=new byte[myNumTags];
        bestChr=new int[myNumTags];
        bestStartPos=new int[myNumTags];
        if(myHDF5.exists(GBSHDF5Constants.MAPBASE+"0")==false) {
            Arrays.fill(bestStrand, TOPMInterface.BYTE_MISSING);
            Arrays.fill(bestChr, TOPMInterface.INT_MISSING);
            Arrays.fill(bestStartPos, TOPMInterface.INT_MISSING);
            return false;
        }
        for (int i = 0; i < myNumTags; i++) {
            int[] posArray=getPositionArray(i);
           // {cachedTMI.chromosome, cachedTMI.strand, cachedTMI.startPosition}
            bestStrand[i]=(byte)posArray[1];
            bestChr[i]=posArray[0];
            bestStartPos[i]=posArray[2];
        }
        myHDF5.createByteArray(GBSHDF5Constants.BEST_STRAND, myNumTags);
        myHDF5.writeByteArray(GBSHDF5Constants.BEST_STRAND, bestStrand, vectorFeatures);
        myHDF5.createIntArray(GBSHDF5Constants.BEST_CHR, myNumTags, vectorFeatures);
        myHDF5.writeIntArray(GBSHDF5Constants.BEST_CHR, bestChr, vectorFeatures);
        myHDF5.createIntArray(GBSHDF5Constants.BEST_STARTPOS, myNumTags, vectorFeatures);
        myHDF5.writeIntArray(GBSHDF5Constants.BEST_STARTPOS, bestStartPos, vectorFeatures);
        return true;
    }
    
    private boolean loadVariantsIntoMemory() {
        int howManyDef=0;
        int readBlock=4096*16;
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
    
    private static boolean writeVariantsToHDF5(IHDF5Writer aHDF5, AbstractTagsOnPhysicalMap aTOPM, int myMaxVariants) {
        int howManyDef=0;
        int readBlock=4096*16;
        
        int myNumTags=aTOPM.myNumTags;
//        int myMaxVariants=aTOPM.myMaxVariants;
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

    private void cacheMappingInfo(int index) {
        if (index == cachedMappingIndex) {
            return;
        }
        int block = index >> BITS_TO_SHIFT_FOR_CHUNK;
        if (cachedMappingBlock != block) {
            if (cleanMap == false) {
                saveCacheBackToFile();
            }
            for (int mi = 0; mi < maxMapping; mi++) {
                cachedTMIBlock[mi] = myHDF5.compounds().readArrayBlock(GBSHDF5Constants.MAPBASE+"0", tmiType, CHUNK_SIZE, block);
                cachedMappingBlock = block;
                if (cacheAllMappingBlocks == false) {
                    break;
                }
            }
        }
        this.cachedTMI = cachedTMIBlock[0][index % CHUNK_SIZE];
        this.cachedMappingIndex = index;
    }

    private void saveCacheBackToFile() {
        int block = cachedMappingIndex >> BITS_TO_SHIFT_FOR_CHUNK;
        if (cachedMappingBlock != block) {
            for (int mi = 0; mi < maxMapping; mi++) {
                if (cleanMap == false) {
                    myHDF5.compounds().writeArrayBlock(GBSHDF5Constants.MAPBASE+"0", tmiType, cachedTMIBlock[mi], block);
                }
                if (cacheAllMappingBlocks == false) {
                    break;
                }
            }
            if (cleanMap == false) {
                myHDF5.writeByteArray("multimaps", multimaps);
            }  //this could be made more efficient by just writing the block
            cleanMap = true;
        }
    }

    public void getFileReadyForClosing() {
//        writeCachedVariantDefs();
//        writeCachedVariantOffsets();
//        saveCacheBackToFile();
    }


    public TagMappingInfo getAlternateTagMappingInfo(int index, int mapIndex) {
        if(hasDetailedMapping) return null;
        if (cacheAllMappingBlocks == false) {
            cacheMappingInfo(index);
        }
        if (cachedMappingIndex != index) {
            cacheMappingInfo(index);
        }
        return cachedTMIBlock[mapIndex][index % CHUNK_SIZE];
    }

    public void setAlternateTagMappingInfo(int index, int mapIndex, TagMappingInfo theTMI) {
        if(!hasDetailedMapping) throw new IllegalStateException("Detailed mapping not present");
        if (cacheAllMappingBlocks == false) {
            cacheMappingInfo(index);
        }
        if (cachedMappingIndex != index) {
            cacheMappingInfo(index);
        }
        cachedTMIBlock[mapIndex][index % CHUNK_SIZE] = theTMI;
        if (multimaps[index] >= mapIndex) {
            multimaps[index] = (byte) (mapIndex + 1);
        }
        cleanMap = false;
    }

    public void swapTagMappingInfo(int index, int mapIndex, int mapIndex2) {
        if(!hasDetailedMapping) throw new IllegalStateException("Detailed mapping not present");
        if (cacheAllMappingBlocks == false) {
            cacheMappingInfo(index);
        }
        if (cachedMappingIndex != index) {
            cacheMappingInfo(index);
        }
        TagMappingInfo tempTMI = cachedTMIBlock[mapIndex][index % CHUNK_SIZE];
        cachedTMIBlock[mapIndex][index % CHUNK_SIZE] = cachedTMIBlock[mapIndex2][index % CHUNK_SIZE];
        cachedTMIBlock[mapIndex2][index % CHUNK_SIZE] = tempTMI;
        cleanMap = false;
    }

    @Override
    public int addVariant(int tagIndex, byte offset, byte base) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public byte getDcoP(int index) {
        if(!hasDetailedMapping) throw new IllegalStateException("Detailed mapping not present");
        if (cachedMappingIndex != index) {
            cacheMappingInfo(index);
        }
        return cachedTMI.dcoP;
    }

    @Override
    public byte getDivergence(int index) {
        if(!hasDetailedMapping) throw new IllegalStateException("Detailed mapping not present");
        if (cachedMappingIndex != index) {
            cacheMappingInfo(index);
        }
        return cachedTMI.divergence;
    }

    @Override
    public int getEndPosition(int index) {
        if(!hasDetailedMapping) throw new IllegalStateException("Detailed mapping not present");
        if (cachedMappingIndex != index) {
            cacheMappingInfo(index);
        }
        return cachedTMI.endPosition;
    }

    @Override
    public byte getMapP(int index) {
        if(!hasDetailedMapping) throw new IllegalStateException("Detailed mapping not present");
        if (cachedMappingIndex != index) {
            cacheMappingInfo(index);
        }
        return cachedTMI.mapP;
    }

    @Override
    public int[] getPositionArray(int index) {
        if (cachedMappingIndex != index) {
            cacheMappingInfo(index);
        }
        int[] r = {cachedTMI.chromosome, cachedTMI.strand, cachedTMI.startPosition};
        return r;
    }

    @Override
    public int getReadIndexForPositionIndex(int posIndex) {
        return indicesOfSortByPosition[posIndex];
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
    

    @Override
    public void clearVariants() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
}
