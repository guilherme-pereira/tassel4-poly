/*
 * TagsByTaxaByteHDF5TagGroups
 */
package net.maizegenetics.gbs.tagdist;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import ch.systemsx.cisd.hdf5.IHDF5WriterConfigurator;

import java.io.File;
import java.nio.ByteBuffer;

import java.util.Arrays;
import java.util.Random;
import java.util.TreeMap;

/**
 * Tags by Taxa file based on the HDF5 data structure.  This version is optimized
 * for rapid access of taxa within Tag (ie it buffers the taxa counts within one
 * tag).
 *
 * 
 * @author edbuckler
 */
public class TagsByTaxaByteHDF5TagGroups extends AbstractTagsByTaxa {

    static int chunkSize = 1 << 16;
    int tagCount = 0;
    IHDF5Writer h5 = null;
    TreeMap<String, Integer> taxaNameIndexTreeMap;
    private int bufferedTagIndex = Integer.MIN_VALUE;
    private byte[] bufferedTagDist = null;
    private boolean bufferChanged = false;

    public TagsByTaxaByteHDF5TagGroups(Tags inTags, String newHDF5file) {
        this.tagLengthInLong = inTags.getTagSizeInLong();
        this.tagCount = inTags.getTagCount();
        this.taxaNum = 0;
        this.tags = new long[tagLengthInLong][tagCount];
        this.tagLength = new byte[tagCount];
        for (int i = 0; i < tagCount; i++) {
            long[] ct = inTags.getTag(i);
            for (int j = 0; j < tagLengthInLong; j++) {
                tags[j][i] = ct[j];
            }
            tagLength[i] = (byte) inTags.getTagLength(i);
        }
        IHDF5WriterConfigurator config = HDF5Factory.configure(new File(newHDF5file));
        System.out.println("Creating HDF5 file: " + newHDF5file);
        config.overwrite();
        config.dontUseExtendableDataTypes();
        config.useUTF8CharacterEncoding();
        h5 = config.writer();
        h5.setIntAttribute("/", "tagCount", inTags.getTagCount());
        h5.setIntAttribute("/", "chunkSize", chunkSize);
        h5.setIntAttribute("/", "tagLengthInLong", tagLengthInLong);
        h5.setIntAttribute("/", "taxaNum", taxaNum);
        //create tag matrix
        h5.createLongMatrix("tags", inTags.getTagSizeInLong(), tagCount, inTags.getTagSizeInLong(), tagCount);
        h5.writeLongMatrix("tags", tags);
        h5.createByteArray("tagLength", tagCount);
        h5.writeByteArray("tagLength", tagLength);
        //create TBT matrix
        h5.createGroup("tbttg");
        int tagChunks = inTags.getTagCount() >> 16;
        if (inTags.getTagCount() % chunkSize > 0) {
            tagChunks++;
        }
        System.out.println(chunkSize);
        System.out.printf("tagChunks %d Div %g %n", tagChunks, (double) inTags.getTagCount() / (double) chunkSize);
        h5.setIntAttribute("tbttg/", "tagCount", inTags.getTagCount());
        h5.setIntAttribute("tbttg/", "tagChunks", tagChunks);

        for (int tc = 0; tc < tagChunks; tc++) {
            h5.createGroup("tbttg/c" + tc);
        }
        if (inTags instanceof TagsByTaxa) {
            TagsByTaxa inTBT = (TagsByTaxa) inTags;
            this.taxaNum = inTBT.getTaxaCount();
            h5.setIntAttribute("/", "taxaNum", taxaNum);
            h5.setIntAttribute("/", "tagCount", tagCount);
            h5.createStringVariableLengthArray("tbttg/taxaNames", taxaNum);
            h5.writeStringVariableLengthArray("tbttg/taxaNames", inTBT.getTaxaNames());

            if (inTBT instanceof TagsByTaxaByteHDF5TaxaGroups) {
                populateTBTMatrixTranspose((TagsByTaxaByteHDF5TaxaGroups) inTBT);
            } else {
                populateTBTMatrix(inTBT);
            }

//            for (int tg = 0; tg < tagCount; tg++) {
//                byte[] td=new byte[taxaNum];
//                for(int tx = 0; tx <taxaNum; tx++) td[tx]=(byte)inTBT.getReadCountForTagTaxon(tg, tx);
//                int chunk=tg>>16;
//                String d="tbttg/c"+chunk+"/"+tg;
//                byte[] deftc=encodeBySign(td);
//                h5.createByteArray(d, deftc.length);
//                h5.writeByteArray(d, deftc);    
//            }
            createFastTaxaMap();
        }
    }

    public void writeDistFile(String newHDF5file, int[] tagIndex) {
        int newTagCount = tagIndex.length;
        long[][] newTags = new long[tagLengthInLong][newTagCount];
        byte[] newTagLength = new byte[newTagCount];
        for (int i = 0; i < newTagCount; i++) {
            long[] ct = this.getTag(tagIndex[i]);
            for (int j = 0; j < tagLengthInLong; j++) {
                newTags[j][i] = ct[j];
            }
            newTagLength[i] = (byte) this.getTagLength(tagIndex[i]);
        }
        IHDF5WriterConfigurator config = HDF5Factory.configure(new File(newHDF5file));
        System.out.println("Creating HDF5 file: " + newHDF5file);
        config.overwrite();
        config.dontUseExtendableDataTypes();
        config.useUTF8CharacterEncoding();
        IHDF5Writer nH5 = config.writer();
        nH5.setIntAttribute("/", "tagCount", newTagCount);
        nH5.setIntAttribute("/", "chunkSize", chunkSize);
        nH5.setIntAttribute("/", "tagLengthInLong", tagLengthInLong);
        nH5.setIntAttribute("/", "taxaNum", taxaNum);
        //create tag matrix
        nH5.createLongMatrix("tags", this.getTagSizeInLong(), newTagCount, this.getTagSizeInLong(), newTagCount);
        nH5.writeLongMatrix("tags", newTags);
        nH5.createByteArray("tagLength", newTagCount);
        nH5.writeByteArray("tagLength", newTagLength);
        //create TBT matrix
        nH5.createGroup("tbttg");
        int tagChunks = newTagCount >> 16;
        if (newTagCount % chunkSize > 0) {
            tagChunks++;
        }
        System.out.println(newTagCount + " tags are written");
        System.out.println("chunk size is " + chunkSize);
        System.out.println("tag chunk number is " + tagChunks);
        nH5.setIntAttribute("tbttg/", "tagCount", newTagCount);
        nH5.setIntAttribute("tbttg/", "tagChunks", tagChunks);

        for (int tc = 0; tc < tagChunks; tc++) {
            nH5.createGroup("tbttg/c" + tc);
        }
        
        nH5.createStringVariableLengthArray("tbttg/taxaNames", taxaNum);
        nH5.writeStringVariableLengthArray("tbttg/taxaNames", this.getTaxaNames());
        
        for (int i = 0; i < newTagCount; i++) {
            byte[] td = new byte[taxaNum];
            for (int j = 0; j < taxaNum; j++) {
                td[j] = (byte) this.getReadCountForTagTaxon(tagIndex[i], j);
            }
            int chunk = i >> 16;
            String d = "tbttg/c" + chunk + "/" + i;
            byte[] deftc = encodeBySign(td);
            nH5.createByteArray(d, deftc.length);
            nH5.writeByteArray(d, deftc);
            if (i%1000 == 0) {
                System.out.println(i + " tag distributions are written");
            }
        }
        createFastTaxaMap();   
    }
    
    private void populateTBTMatrix(TagsByTaxa inTBT) {
        for (int tg = 0; tg < tagCount; tg++) {
            byte[] td = new byte[taxaNum];
            for (int tx = 0; tx < taxaNum; tx++) {
                td[tx] = (byte) inTBT.getReadCountForTagTaxon(tg, tx);
            }
            int chunk = tg >> 16;
            String d = "tbttg/c" + chunk + "/" + tg;
            byte[] deftc = encodeBySign(td);
            h5.createByteArray(d, deftc.length);
            h5.writeByteArray(d, deftc);
        }
    }

    private void populateTBTMatrixTranspose(TagsByTaxaByteHDF5TaxaGroups inTBT) {
        for (int tgC = 0; tgC < tagCount; tgC += chunkSize) {
            System.out.println("TagChunk:" + tgC);
            byte[][] td = new byte[chunkSize][taxaNum];
            int top = (tgC + chunkSize > tagCount) ? tagCount : (tgC + chunkSize);
            for (int tx = 0; tx < taxaNum; tx++) {
                for (int tg = tgC; tg < top; tg++) {
                    td[tg - tgC][tx] = (byte) inTBT.getReadCountForTagTaxon(tg, tx);
                }
            }
            for (int tg = tgC; tg < top; tg++) {
                int chunk = tgC >> 16;
                String d = "tbttg/c" + chunk + "/" + tg;
                byte[] deftc = encodeBySign(td[tg - tgC]);
                h5.createByteArray(d, deftc.length);
                h5.writeByteArray(d, deftc);
            }

        }
    }

    public TagsByTaxaByteHDF5TagGroups(String infile) {
        IHDF5WriterConfigurator config = HDF5Factory.configure(new File(infile));
        config.dontUseExtendableDataTypes();
        h5 = config.writer();
        tagCount = h5.getIntAttribute("/", "tagCount");
        chunkSize = h5.getIntAttribute("/", "chunkSize");
        tagLengthInLong = h5.getIntAttribute("/", "tagLengthInLong");
        taxaNum = h5.getIntAttribute("/", "taxaNum");
        this.tags = h5.readLongMatrix("tags");
        this.tagLength = h5.readByteArray("tagLength");
        createFastTaxaMap();
    }

    private void createFastTaxaMap() {
        taxaNames = h5.readStringArray("tbttg/taxaNames");
        taxaNameIndexTreeMap = new TreeMap<String, Integer>();
        for (int i = 0; i < taxaNames.length; i++) {
            taxaNameIndexTreeMap.put(taxaNames[i], i);
        }
    }

    public static byte[] createRandomDistribution(int arraySize, double proportionToFill, int maxValue) {
        byte[] result = new byte[arraySize];
        int fillNum = (int) (result.length * proportionToFill);
        Random r = new Random();
        for (int i = 0; i < fillNum; i++) {
            byte cnt = (byte) r.nextInt(maxValue);
            result[r.nextInt(arraySize)] = cnt;
        }
        return result;
    }

    public static byte[][] encodeBySign(byte[] source, int chunkSize) {
        int chunks = source.length / chunkSize;
        if (source.length % chunkSize != 0) {
            chunks++;
        }
        byte[][] result = new byte[chunks][];
        for (int i = 0; i < chunks; i++) {
            int s = i * chunkSize;
            int e = s + chunkSize;
            if (e > source.length) {
                e = source.length;
            }
            result[i] = encodeBySign(Arrays.copyOfRange(source, s, e));
        }
        return result;
    }

    public static byte[] encodeBySign(byte[] source) {
        ByteBuffer dest = ByteBuffer.allocate(source.length + 4);
        dest.putInt(source.length);
        byte runLength = 0;
        for (int i = 0; i < source.length; i++) {
            if (source[i] > 0) {
                if (runLength < 0) {
                    dest.put(runLength);
                }
                dest.put(source[i]);
                runLength = 0;
            } else {
                runLength--;
                if (runLength == Byte.MIN_VALUE) {
                    dest.put(runLength);
                    runLength = 0;
                }
            }
        }
        if (runLength < 0) {
            dest.put(runLength);
        }
        return Arrays.copyOf(dest.array(), dest.position());
    }

    public static byte[] decodeBySign(byte[][] srcCompChunk) {
        byte[][] resInChunk = new byte[srcCompChunk.length][];
        int totalLength = 0;
        for (int i = 0; i < srcCompChunk.length; i++) {
            resInChunk[i] = decodeBySign(srcCompChunk[i]);
            totalLength += resInChunk[i].length;
        }
        ByteBuffer result = ByteBuffer.allocate(totalLength);
        for (int i = 0; i < srcCompChunk.length; i++) {
            result.put(resInChunk[i]);
        }
        return result.array();
    }

    public static byte[] decodeBySign(byte[] source) {
        ByteBuffer srcB = ByteBuffer.wrap(source);
        int length = srcB.getInt();
        ByteBuffer dest = ByteBuffer.allocate(length);
        for (int i = 4; i < source.length; i++) {
            if (source[i] > 0) {
                dest.put(source[i]);
            } else {
                dest.position((-source[i]) + dest.position());
            }
        }
        return dest.array();
    }

    public TagsByTaxaByte convertToTBTByte() {
        byte[][] tagDist = new byte[this.taxaNum][this.tagCount];
        for (int tag = 0; tag < this.tagCount; tag++) {
            for (int taxon = 0; taxon < this.taxaNum; taxon++) {
                tagDist[taxon][tag] = (byte) getReadCountForTagTaxon(tag, taxon);
            }
        }
        return new TagsByTaxaByte(this.tags, this.tagLength, tagDist, this.getTaxaNames());
    }

    public static void main(String[] args) {
        String inTBTFile = "/Users/edbuckler/SolexaAnal/GBS/build20120110/tbt/434GFAAXX_s_3.tbt.byte";
        String path = "/Users/edbuckler/SolexaAnal/GBS/test/";
        //static String path="/Volumes/LaCie/";
        String file = "testToGZ.h5";
        TagsByTaxa inTBT = new TagsByTaxaByte(inTBTFile, FilePacking.Byte);
        long time = System.currentTimeMillis();
        TagsByTaxaByteHDF5TagGroups tHDF5 = new TagsByTaxaByteHDF5TagGroups(inTBT, path + file);
        tHDF5.getFileReadyForClosing();
        TagsByTaxaByteHDF5TagGroups rHDF5 = new TagsByTaxaByteHDF5TagGroups(path + file);
        int same = 0, diff = 0, count = 0;
        time = System.currentTimeMillis();
        int tags = 9;
        for (int i = 0; i < 10000; i += 11) {
            int taxon = i % inTBT.getTaxaCount();
            // taxon=15;
            if (i % 550 == 0) {
                tags = i % inTBT.getTagCount();
            }

            int newTaxonIndex = rHDF5.getIndexOfTaxaName(inTBT.getTaxaName(taxon));
            if (inTBT.getReadCountForTagTaxon(tags, taxon) == rHDF5.getReadCountForTagTaxon(tags, newTaxonIndex)) {
                same++;
            } else {
                diff++;
            }
            count++;
        }
        System.out.printf("Same %d Diff %d %n", same, diff);
        long duration = System.currentTimeMillis() - time;
        double rate = (double) duration / (double) count;
        System.out.printf("Rate %g %n", rate);
    }

    @Override
    public int getIndexOfTaxaName(String taxon) {
        Integer index;
        index = taxaNameIndexTreeMap.get(taxon);
        if (index == null) {
            return -1;
        }
        return index;
    }

    
    @Override
    public synchronized int getReadCountForTagTaxon(int tagIndex, int taxaIndex) {
        if (bufferedTagIndex != tagIndex) {
            bufferTagDist(tagIndex);
        }
        return bufferedTagDist[taxaIndex];
    }

    synchronized private void bufferTagDist(int tagIndex) {
        if (bufferChanged) {//save the old buffer 
            int chunk = bufferedTagIndex >> 16;
            String g = "tbttg/c" + chunk + "/" + bufferedTagIndex;
            h5.writeByteArray(g, bufferedTagDist);
        }
        int chunk = tagIndex >> 16;
        String g = "tbttg/c" + chunk + "/" + tagIndex;
        this.bufferedTagDist = decodeBySign(h5.readByteArray(g));
        this.bufferedTagIndex = tagIndex;
        bufferChanged = false;
    }

    @Override
    public void setMethodByRows(boolean rowSetMethod) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void setReadCountForTagTaxon(int tagIndex, int taxaIndex, int value) {
        if (bufferedTagIndex != taxaIndex) {
            bufferTagDist(taxaIndex);
        }
        bufferChanged = true;
        if (value > Byte.MAX_VALUE) {
            bufferedTagDist[taxaIndex] = Byte.MAX_VALUE;
        } else if (value < 0) {
            bufferedTagDist[taxaIndex] = 0;
        } else {
            bufferedTagDist[taxaIndex] = (byte) value;
        }
    }

    @Override
    public void initMatrices(int taxaNum, int tagNum) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void addTaxa(String[] addTaxaNames) {
        throw new UnsupportedOperationException("Is not supported, use other HDF5 class");
    }

    @Override
    public void getFileReadyForClosing() {
        bufferTagDist(0);
        h5.close();
    }
}
