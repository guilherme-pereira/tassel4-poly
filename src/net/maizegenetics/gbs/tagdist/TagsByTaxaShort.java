/*
 * TagsByTaxaShort
 */
package net.maizegenetics.gbs.tagdist;

import java.io.File;

/**
 * Tags by taxa that supports depths captured by a short (max 32767)
 *
 * @author jcg233
 */
public class TagsByTaxaShort extends AbstractTagsByTaxa {

    private short[][] tagDist;

    public TagsByTaxaShort() {
    }

    public TagsByTaxaShort(String infile, FilePacking binary) {
        readDistFile(new File(infile), binary);
    }

    public TagsByTaxaShort(String[] taxaNames, Tags theDistinctReads) {
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

    public TagsByTaxaShort(long[][] reads, short[][] readDist, String[] namesForTaxa) {
        tags = reads;
        tagDist = readDist;
        taxaNames = namesForTaxa;
        taxaNum = namesForTaxa.length;
    }

    @Override
    public void setReadCountForTagTaxon(int tagIndex, int taxaIndex, int value) {
        if (value > Short.MAX_VALUE) {
            tagDist[taxaIndex][tagIndex] = Short.MAX_VALUE;
        } else if (value < 0) {
            tagDist[taxaIndex][tagIndex] = 0;
        } else {
            tagDist[taxaIndex][tagIndex] = (short) value;
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
        tagDist = new short[taxaNum][tagNum];
        tagLength = new byte[tagNum];
    }

    @Override
    public void addTaxa(String[] addTaxaNames) {
        taxaNum = tagDist.length + addTaxaNames.length;
        short[][] newTagDist = new short[taxaNum][];
        String[] newTaxaNames = new String[taxaNum];
        for (int i = 0; i < newTagDist.length; i++) {
            if (i < tagDist.length) {
                newTagDist[i] = tagDist[i];
                newTaxaNames[i] = taxaNames[i];
            } else {
                newTagDist[i] = new short[getTagCount()];
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
