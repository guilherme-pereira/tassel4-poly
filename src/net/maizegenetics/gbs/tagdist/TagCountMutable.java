package net.maizegenetics.gbs.tagdist;

/**
 * Resizeable tags count container.  Tags sequences are compressed in long, tags lengths are tracked,
 * and read counts are stored.  Has basic filtering methods.
 * Read counts can be modified.
 *
 *
 */
public class TagCountMutable extends TagCounts {

    int currentRows;

    public TagCountMutable(int tagLengthInLong, int maxSize) {
        currentRows = 0;
        this.tagLengthInLong = tagLengthInLong;
        initMatrices(maxSize);
    }

    public TagCountMutable(Tags origTagCount, int maxSize) {
        currentRows = 0;
        this.tagLengthInLong = origTagCount.getTagSizeInLong();
        initMatrices(maxSize);
        if (origTagCount instanceof TagCounts) {
            addReadCounts((TagCounts) origTagCount);
        } else {
            {
                addReadCounts(origTagCount, 0);
            }
        }
    }

    protected void initMatrices(int tagNum) {
        tags = new long[tagLengthInLong][tagNum];
        tagLength = new byte[tagNum];
        readCount = new int[tagNum];
    }

    public int getCurrentSize() {
        return currentRows;
    }

    public void addReadCount(long[] tag, int length, int count) {
        for (int j = 0; j < tagLengthInLong; j++) {
            this.tags[j][currentRows] = tag[j];
        }

        this.tagLength[currentRows] = (byte) length;
        this.readCount[currentRows] = count;
        currentRows++;
    }

    /**
     * This adds a series of TagCounts to the list, preserving the count that each read has.
     * Designed to make it easy to combine TagCounts from multiple barcodes, lanes or flowcells.
     * @param tagCountsToAdd
     */
    public void addReadCounts(TagCounts tagCountsToAdd) {
        for (int i = 0; i < tagCountsToAdd.getTagCount(); i++) {
            long[] ls = tagCountsToAdd.getTag(i);
            for (int j = 0; j < tagLengthInLong; j++) {
                this.tags[j][currentRows] = ls[j];
            }
            this.tagLength[currentRows] = (byte) tagCountsToAdd.getTagLength(i);
            this.readCount[currentRows] = tagCountsToAdd.getReadCount(i);
            currentRows++;
        }
    }

    public void addReadCounts(Tags tagsToAdd, int defaultCount) {
        for (int i = 0; i < tagsToAdd.getTagCount(); i++) {
            long[] ls = tagsToAdd.getTag(i);
            for (int j = 0; j < tagLengthInLong; j++) {
                this.tags[j][currentRows] = ls[j];
            }
            this.tagLength[currentRows] = (byte) tagsToAdd.getTagLength(i);
            this.readCount[currentRows] = defaultCount;
            currentRows++;
        }
    }

    public void collapseCounts() {
        this.sort();
        int collapsedRows = 0;
        int initSize = this.getCurrentSize();
        for (int i = 1; i < this.getSize(); i++) {
            if (tags[0][i] == 0) {
                readCount[i] = 0;
                continue;
            }
            if (compare(i, i - 1) == 0) {
                readCount[i] += readCount[i - 1];
                readCount[i - 1] = 0;
                for (int j = 0; j < tagLengthInLong; j++) {
                    tags[j][i - 1] = 0;
                }
                collapsedRows++;
                currentRows--;
            }
        }
        this.sort();
        System.out.println("Rows collapsed:" + collapsedRows);
        int retained = initSize - collapsedRows;
        System.out.println("Unique tags retained:" + retained);
    }

    public void removeRareTag(int minCutoff) {
        for (int i = 0; i < this.getCurrentSize(); i++) {
            if (readCount[i] < minCutoff) {
                tags[0][i] = 0;
                readCount[i] = 0;
            }
        }
        this.collapseCounts();
    }

    public int shrinkToCurrentRows() {
        for (int i = 0; i < tagLengthInLong; i++) {
            long[] t1 = new long[currentRows];
            System.arraycopy(tags[i], 0, t1, 0, currentRows);
            tags[i] = t1;
        }
        return currentRows;
    }

    @Override
    public int compare(int index1, int index2) {
        if ((tags[0][index1] == 0) || (tags[0][index2] == 0)) {
            if ((tags[0][index1] == tags[0][index2])) {
                return 0;
            }
            if (tags[0][index1] != 0) {
                return -1;
            }
            return 1;
        }
        for (int i = 0; i < tagLengthInLong; i++) {
            if (tags[i][index1] < tags[i][index2]) {
                return -1;
            }
            if (tags[i][index1] > tags[i][index2]) {
                return 1;
            }
        }
        return 0;
    }
}
