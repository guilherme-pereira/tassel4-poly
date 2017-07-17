/*
 * TagsByTaxaByteHDF5TaxaGroups
 */
package net.maizegenetics.gbs.tagdist;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Writer;
import ch.systemsx.cisd.hdf5.IHDF5WriterConfigurator;

import java.io.File;

import java.nio.ByteBuffer;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.TreeMap;

/**
 * Tags by Taxa file based on the HDF5 data structure.  This version is optimized
 * for rapid access of tags within taxa (ie it buffers the tag counts within one
 * taxon).  It is good for adding, removing, and combining taxa
 * 
 * @author edbuckler
 */
public class TagsByTaxaByteHDF5TaxaGroups extends AbstractTagsByTaxa {

    static String path = "/Users/edbuckler/SolexaAnal/GBS/test/";
    //static String path="/Volumes/LaCie/";
    static String file = "testToGZ.h5";
    static int chunkSize = 1 << 16;
    int tagCount = 0;
    int tagChunks = 0;
    IHDF5Writer h5 = null;
    ArrayList<String> taxaDirList, taxaNameList;
    TreeMap<String, String> taxaNameDirTreeMap;
    private int bufferedTaxaIndex = Integer.MIN_VALUE;
    private int bufferedChunkIndex = Integer.MIN_VALUE;
    private byte[] bufferedTagDist = null;
    private boolean bufferChanged = false;

    public TagsByTaxaByteHDF5TaxaGroups(Tags inTags, String newHDF5file) {
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
        h5.createGroup("tbttx");
        tagChunks = inTags.getTagCount() >> 16;
        if (inTags.getTagCount() % chunkSize > 0) {
            tagChunks++;
        }
        System.out.println(chunkSize);
        System.out.printf("tagChunks %d Div %g %n", tagChunks, (double) inTags.getTagCount() / (double) chunkSize);
        h5.setIntAttribute("tbttx/", "tagCount", inTags.getTagCount());
        h5.setIntAttribute("tbttx/", "tagChunks", tagChunks);
        if (inTags instanceof TagsByTaxa) {
            TagsByTaxa inTBT = (TagsByTaxa) inTags;
            this.taxaNum = inTBT.getTaxaCount();
            h5.setIntAttribute("/", "taxaNum", taxaNum);
            for (int tx = 0; tx < inTBT.getTaxaCount(); tx++) {
                byte[] tc = new byte[inTBT.getTagCount()];
                for (int i = 0; i < tc.length; i++) {
                    tc[i] = (byte) inTBT.getReadCountForTagTaxon(i, tx);
                }
                this.addTaxon(inTBT.getTaxaName(tx), tc);
            }
        }
    }

    public TagsByTaxaByteHDF5TaxaGroups(String infile) {
        IHDF5WriterConfigurator config = HDF5Factory.configure(new File(infile));
        config.dontUseExtendableDataTypes();
        h5 = config.writer();
        tagCount = h5.getIntAttribute("/", "tagCount");
        chunkSize = h5.getIntAttribute("/", "chunkSize");
        tagLengthInLong = h5.getIntAttribute("/", "tagLengthInLong");
        taxaNum = h5.getIntAttribute("/", "taxaNum");
        tagChunks = h5.getIntAttribute("tbttx/", "tagChunks");
        this.tags = h5.readLongMatrix("tags");
        this.tagLength = h5.readByteArray("tagLength");
        createFastTaxaMap();
    }

    private void createFastTaxaMap() {
        ArrayList<String> tmpList = new ArrayList<String>(h5.getAllGroupMembers("tbttx/"));
        //      System.out.println(tmpList.toString());
        taxaNameDirTreeMap = new TreeMap<String, String>();
        for (String tx : tmpList) {
            String s = h5.getStringAttribute("tbttx/" + tx, "name");
            taxaNameDirTreeMap.put(s, "tbttx/" + tx);
        }
        this.taxaDirList = new ArrayList(taxaNameDirTreeMap.values());
        this.taxaNameList = new ArrayList(taxaNameDirTreeMap.keySet());
    }

    public boolean addTaxon(String taxonName, byte[] values) {
        synchronized (h5) {
            if (values.length != this.tagCount) {
                System.err.printf("Taxon (%s) does not have the right number of sites (%d)%n", taxonName, values.length);
                return false;
            }
            if (addTaxon(taxonName) == false) {
                return false;
            }
            String ldir = this.taxaNameDirTreeMap.get(taxonName);
            byte[][] defTag = encodeBySign(values, chunkSize);
            for (int c = 0; c < defTag.length; c++) {
                h5.writeByteArray(ldir + "/c" + c, defTag[c]);
            }
            return true;
        }
    }

    public boolean addTaxon(String taxonName) {
        synchronized (h5) {
            int thc = taxonName.hashCode();
            while (h5.isGroup("tbttx/" + thc)) {
                if (h5.getStringAttribute("tbttx/" + thc, "name").equals(taxonName)) {
                    System.err.printf("Taxon (%s) cannot be added already exist%n", taxonName);
                    return false;
                }
                thc++;
            }
            String lg = "tbttx/" + thc;
            h5.createGroup(lg);
            h5.setStringAttribute(lg, "name", taxonName);
            createFastTaxaMap();
            taxaNum = taxaNameList.size();
            h5.setIntAttribute("/", "taxaNum", taxaNum);
            return true;
        }
    }

    public boolean deleteTaxon(String taxonName) {
        synchronized (h5) {
            h5.delete(taxaNameDirTreeMap.get(taxonName));
            createFastTaxaMap();
            taxaNum = taxaNameList.size();
            h5.setIntAttribute("/", "taxaNum", taxaNum);
            return true;
        }
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
        ByteBuffer dest = ByteBuffer.allocate(source.length);
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

    public static void main(String[] args) {
        String inTBTFile = "/Users/edbuckler/SolexaAnal/GBS/build20120110/tbt/434GFAAXX_s_3.tbt.byte";
        TagsByTaxa inTBT = new TagsByTaxaByte(inTBTFile, FilePacking.Byte);
        long time = System.currentTimeMillis();
        TagsByTaxaByteHDF5TaxaGroups tHDF5 = new TagsByTaxaByteHDF5TaxaGroups(inTBT, path + file);
        TagsByTaxaByteHDF5TaxaGroups rHDF5 = new TagsByTaxaByteHDF5TaxaGroups(path + file);
        int same = 0, diff = 0, count = 0;
        time = System.currentTimeMillis();
        for (int i = 0; i < 1000000; i += 11) {
            int taxon = i % inTBT.getTaxaCount();
//            taxon=15;
            int tags = i % inTBT.getTagCount();
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
        int index = Collections.binarySearch(taxaNameList, taxon);
        return index;
    }

    @Override
    public int getReadCountForTagTaxon(int tagIndex, int taxaIndex) {
        int chunk = tagIndex >> 16;
        if ((bufferedChunkIndex != chunk) || (bufferedTaxaIndex != taxaIndex)) {
            bufferTagDist(chunk, taxaIndex);
        }
        int offset = tagIndex % chunkSize;
        return bufferedTagDist[offset];
    }

    public byte[] getReadCountDistributionForTaxon(int taxaIndex) {
        ByteBuffer bb = ByteBuffer.allocate(tagCount);
        for (int i = 0; i < tagChunks; i++) {
            bufferTagDist(i, taxaIndex);
            bb.put(bufferedTagDist);
        }
        return bb.array();
    }

    public int getNumberOfChunks() {
        return tagChunks;
    }

    public int getNumberOfTagsPerChunk() {
        return chunkSize;
    }

    synchronized private void bufferTagDist(int tagChunk, int taxaIndex) {
        if (bufferChanged) {//save the old buffer        
            String g = taxaDirList.get(taxaIndex) + "/c" + tagChunk;
            h5.writeByteArray(g, bufferedTagDist);
        }
        String g = taxaDirList.get(taxaIndex) + "/c" + tagChunk;
        this.bufferedTagDist = decodeBySign(h5.readByteArray(g));
        this.bufferedChunkIndex = tagChunk;
        this.bufferedTaxaIndex = taxaIndex;
        bufferChanged = false;
    }

    @Override
    public void setMethodByRows(boolean rowSetMethod) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void setReadCountForTagTaxon(int tagIndex, int taxaIndex, int value) {
        synchronized (h5) {
            int chunk = tagIndex >> 16;
            if ((bufferedChunkIndex != chunk) || (bufferedTaxaIndex != taxaIndex)) {
                bufferTagDist(chunk, taxaIndex);
            }
            int offset = tagIndex % chunkSize;
            bufferChanged = true;
            if (value > Byte.MAX_VALUE) {
                bufferedTagDist[offset] = Byte.MAX_VALUE;
            } else if (value < 0) {
                bufferedTagDist[offset] = 0;
            } else {
                bufferedTagDist[offset] = (byte) value;
            }
        }
    }

    @Override
    public void initMatrices(int taxaNum, int tagNum) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void addTaxa(String[] addTaxaNames) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void getFileReadyForClosing() {
        bufferTagDist(0, 0);
        h5.close();
    }

    @Override
    public int getTaxaCount() {
        return taxaNum;
    }

    @Override
    public String getTaxaName(int taxaIndex) {
        return taxaNameList.get(taxaIndex);
    }

    @Override
    public String[] getTaxaNames() {
        String[] array = taxaNameList.toArray(new String[taxaNameList.size()]);
        return array;
    }

    public TagsByTaxaByte convertToTBTByte() {
        byte[][] tagDist = new byte[this.taxaNum][this.tagCount];
        for (int taxon = 0; taxon < this.taxaNum; taxon++) {
            for (int tag = 0; tag < this.tagCount; tag++) {
                tagDist[taxon][tag] = (byte) getReadCountForTagTaxon(tag, taxon);
            }
        }
        return new TagsByTaxaByte(this.tags, this.tagLength, tagDist, this.getTaxaNames());
    }
}
