/*
 * TagsByTaxaBit
 */
package net.maizegenetics.gbs.tagdist;

import java.io.File;

import java.util.BitSet;

/**
 * A bit implementation of TagByTaxa, so only presence or absence of a tags is known.
 *
 * @deprecated While efficient for memory storage that lack of depth is a problem for certain algorithms.  HDF5 classes
 * are a better alternatives {@link TagsByTaxaByteHDF5TagGroups} and {@link TagsByTaxaByteHDF5TaxaGroups}
 * @author edbuckler
 */
@Deprecated
public class TagsByTaxaBit extends AbstractTagsByTaxa {

    private BitSet[] tagDist;

    public TagsByTaxaBit() {
    }

    public TagsByTaxaBit(String infile, FilePacking binary) {
        readDistFile(new File(infile), binary);
    }

    public TagsByTaxaBit(String[] taxaNames, Tags theDistinctReads) {
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

    public TagsByTaxaBit(long[][] reads, BitSet[] readDist, String[] taxaNames) {
        tags = reads;
        tagDist = readDist;
        this.taxaNames = taxaNames;
        taxaNum = taxaNames.length;
    }

    @Override
    public void setReadCountForTagTaxon(int tagIndex, int taxaIndex, int value) {
        tagDist[taxaIndex].set(tagIndex, (value != 0));
    }

    @Override
    public int getReadCountForTagTaxon(int tagIndex, int taxaIndex) {
        if (tagDist[taxaIndex].get(tagIndex)) {
            return 1;
        }
        return 0;
    }

    @Override
    public void initMatrices(int taxaNum, int tagNum) {
        taxaNames = new String[taxaNum];
        tags = new long[tagLengthInLong][tagNum];
        tagDist = new BitSet[taxaNum];
        for (int i = 0; i < taxaNum; i++) {
            tagDist[i] = new BitSet(tagNum);
        }
        tagLength = new byte[tagNum];
    }

    @Override
    public void addTaxa(String[] addTaxaNames) {
        taxaNum = tagDist.length + addTaxaNames.length;
        BitSet[] newTagDist = new BitSet[taxaNum];
        String[] newTaxaNames = new String[taxaNum];
        for (int i = 0; i < newTagDist.length; i++) {
            if (i < tagDist.length) {
                newTagDist[i] = tagDist[i];
                newTaxaNames[i] = taxaNames[i];
            } else {
                newTagDist[i] = new BitSet(getTagCount());
                newTaxaNames[i] = addTaxaNames[i - tagDist.length];
            }
        }
        tagDist = newTagDist;
        taxaNames = newTaxaNames;
    }

    public static void main(String[] args) {
        System.out.println("Running main method in TagsByTaxaBit");
        String infileFile = "/Users/edbuckler/SolexaAnal/GBS/dist/h10000HiSeq282_110214.dist.txt";
        String outfile = "/Users/edbuckler/SolexaAnal/GBS/test/h10000HiSeq282_110214.dist.bibin";
        String outfile2 = "/Users/edbuckler/SolexaAnal/GBS/test/dup_h10000HiSeq282_110214.dist.txt";
        TagsByTaxaBit tbtb = new TagsByTaxaBit(infileFile, FilePacking.Text);
        tbtb.writeDistFile(new File(outfile), FilePacking.Bit, -1);
        TagsByTaxaBit tbtb2 = new TagsByTaxaBit(outfile, FilePacking.Bit);
        tbtb2.writeDistFile(new File(outfile2), FilePacking.Text, -1);

    }

    public void setMethodByRows(boolean rowSetMethod) {
        throw (new NoSuchMethodError("This method is only implemented in classes that use a RandomAccessFile."));
    }

    public void getFileReadyForClosing() {
        throw (new NoSuchMethodError("This method is only implemented in classes that use a RandomAccessFile."));
    }
}
