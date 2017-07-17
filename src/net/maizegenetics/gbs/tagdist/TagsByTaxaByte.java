/*
 * TagsByTaxaByte
 */
package net.maizegenetics.gbs.tagdist;

import java.io.File;

/**
 * A bit implementation of TagByTaxa, so only counts as high as 127 can be maintained.
 * @author ed
 */
public class TagsByTaxaByte extends AbstractTagsByTaxa {

    private byte[][] tagDist;

    public TagsByTaxaByte() {
    }

    public TagsByTaxaByte(String infile, FilePacking binary) {
        readDistFile(new File(infile), binary);
    }

    public TagsByTaxaByte(String[] taxaNames, Tags theDistinctReads) {
        this.taxaNames = taxaNames.clone();
        taxaNum = taxaNames.length;
        int tagNum = theDistinctReads.getTagCount();
        tagLengthInLong = theDistinctReads.getTagSizeInLong();
        initMatrices(taxaNum, tagNum);
        for (int i = 0; i < tagNum; i++) {
            long[] h = theDistinctReads.getTag(i);
            for (int j = 0; j < h.length; j++) {
                tags[j][i] = h[j];
            }
            tagLength[i] = (byte) theDistinctReads.getTagLength(i);
        }
        for (int i = 0; i < taxaNames.length; i++) {
            this.taxaNames[i] = taxaNames[i];
        }
    }

    public TagsByTaxaByte(long[][] tags, byte[] tagLength, byte[][] tagDist, String[] namesForTaxa) {
        this.tags = tags;
        this.tagLength = tagLength;
        this.tagDist = tagDist;
        this.taxaNames = namesForTaxa;
        this.taxaNum = namesForTaxa.length;
        this.tagLengthInLong = tags.length;
    }

    public TagsByTaxaByte(long[][] reads, byte[][] readDist, String[] namesForTaxa) {
        tags = reads;
        tagDist = readDist;
        taxaNames = namesForTaxa;
        taxaNum = namesForTaxa.length;
        //tagNum = reads[0].length;
    }

    @Override
    public void setReadCountForTagTaxon(int tagIndex, int taxaIndex, int value) {
        if (value > Byte.MAX_VALUE) {
            tagDist[taxaIndex][tagIndex] = Byte.MAX_VALUE;
        } else if (value < 0) {
            tagDist[taxaIndex][tagIndex] = 0;
        } else {
            tagDist[taxaIndex][tagIndex] = (byte) value;
        }
    }

    @Override
    public int getReadCountForTagTaxon(int tagIndex, int taxaIndex) {
        return tagDist[taxaIndex][tagIndex];
    }

    @Override
    public void initMatrices(int taxaNum, int tagNum) {
        taxaNames = new String[taxaNum];
        tags = new long[tagLengthInLong][tagNum];
        tagDist = new byte[taxaNum][tagNum];
        tagLength = new byte[tagNum];
    }

    @Override
    public void addTaxa(String[] addTaxaNames) {
        taxaNum = tagDist.length + addTaxaNames.length;
        byte[][] newTagDist = new byte[taxaNum][];
        String[] newTaxaNames = new String[taxaNum];
        for (int i = 0; i < newTagDist.length; i++) {
            if (i < tagDist.length) {
                newTagDist[i] = tagDist[i];
                newTaxaNames[i] = taxaNames[i];
            } else {
                newTagDist[i] = new byte[getTagCount()];
                newTaxaNames[i] = addTaxaNames[i - tagDist.length];
            }
        }
        tagDist = newTagDist;
        taxaNames = newTaxaNames;
    }

    public void setMethodByRows(boolean rowSetMethod) {
        throw (new NoSuchMethodError("This method is only implemented in classes that use a RandomAccessFile."));
    }

    public void getFileReadyForClosing() {
        throw (new NoSuchMethodError("This method is only implemented in classes that use a RandomAccessFile."));
    }
}
