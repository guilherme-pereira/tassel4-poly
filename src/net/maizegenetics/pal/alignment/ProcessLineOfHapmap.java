/*
 * ProcessLineOfHapmap
 */
package net.maizegenetics.pal.alignment;

import net.maizegenetics.util.BitUtil;
import net.maizegenetics.util.OpenBitSet;

/**
 *
 * @author terry
 */
public class ProcessLineOfHapmap implements Runnable {

    private String[][] myTokens;
    private int mySite;
    private int myNumTaxa;
    private int myLineInFile;
    private boolean myComplete = false;
    private OpenBitSet[][] myData;
    private byte[][] myAlleles;
    private int[][] myAlleleMappings;
    private int myNumAlleles;
    private int myNumDataRows;
    private boolean myRetainRareAlleles;
    private boolean myIsSBit;
    private int myNumSitesToProcess;

    private ProcessLineOfHapmap(byte[][] alleles, OpenBitSet[][] data, boolean retainRareAlleles, String[][] tokens, int numSitesToProcess, int site, int numTaxa, int lineInFile, boolean isSBit) {
        setVariables(alleles, data, retainRareAlleles, tokens, numSitesToProcess, site, numTaxa, lineInFile, isSBit);
    }

    public static ProcessLineOfHapmap getInstance(byte[][] alleles, OpenBitSet[][] data, boolean retainRareAlleles, String[][] tokens, int numSitesToProcess, int site, int numTaxa, int lineInFile, boolean isSBit) {
        return new ProcessLineOfHapmap(alleles, data, retainRareAlleles, tokens, numSitesToProcess, site, numTaxa, lineInFile, isSBit);
    }

    private void setVariables(byte[][] alleles, OpenBitSet[][] data, boolean retainRareAlleles, String[][] tokens, int numSitesToProcess, int site, int numTaxa, int lineInFile, boolean isSBit) {

        myData = data;
        myTokens = tokens;
        myNumSitesToProcess = numSitesToProcess;
        if (myNumSitesToProcess > 64) {
            throw new IllegalStateException("ProcessLineOfHapmap: setVariables: Can't process more than 64 sites: " + myNumSitesToProcess);
        }
        mySite = site;
        myNumTaxa = numTaxa;
        myLineInFile = lineInFile;
        myComplete = false;
        myAlleles = alleles;
        myNumAlleles = myAlleles[0].length;
        myAlleleMappings = new int[myNumSitesToProcess][16];
        for (int s = 0; s < myNumSitesToProcess; s++) {
            for (int i = 0; i < 16; i++) {
                myAlleleMappings[s][i] = -1;
            }
        }
        myRetainRareAlleles = retainRareAlleles;
        myNumDataRows = myNumAlleles;
        if (myRetainRareAlleles) {
            myNumDataRows++;
        }

        if (myNumDataRows != myData.length) {
            throw new IllegalStateException("ProcessLineOfHapmap: setVariables: number of data rows: " + myNumDataRows + " should equal first dimension of data array: " + myData.length);
        }

        myIsSBit = isSBit;
    }

    @Override
    public void run() {
        try {
            if (myComplete) {
                throw new IllegalStateException("ImportUtils: ProcessLineOfHapmap: run: trying to run completed instance.");
            }
            myComplete = true;

            byte[][] data = new byte[myNumSitesToProcess][myNumTaxa];
            for (int s = 0; s < myNumSitesToProcess; s++) {
                for (int i = 0; i < myNumTaxa; i++) {
                    try {
                        data[s][i] = NucleotideAlignmentConstants.getNucleotideDiploidByte(myTokens[s][ImportUtils.NUM_HAPMAP_NON_TAXA_HEADERS + i]);
                    } catch (IndexOutOfBoundsException ex) {
                        throw new IllegalStateException("Number of Taxa: " + myNumTaxa + " does not match number of values at line in file: " + (myLineInFile + s) + " site: " + mySite + s);
                    } catch (Exception e) {
                        throw new IllegalStateException("Problem with line in file: " + (myLineInFile + s), e);
                    }
                }
            }

            setAlleles(data);
            if (myIsSBit) {
                setSBits(data);
            } else {
                setTBits(data);
            }

        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    private void setAlleles(byte[][] data) {
        for (int s = 0; s < myNumSitesToProcess; s++) {
            int[][] alleles = AlignmentUtils.getAllelesSortedByFrequency(data[s]);
            int resultSize = alleles[0].length;
            for (int i = 0; i < myNumAlleles; i++) {
                if (i < resultSize) {
                    myAlleles[mySite + s][i] = (byte) alleles[0][i];
                    myAlleleMappings[s][myAlleles[mySite + s][i]] = i;
                } else {
                    myAlleles[mySite + s][i] = Alignment.UNKNOWN_ALLELE;
                }
            }
        }
    }

    private void setSBits(byte[][] data) {

        int numLongs = BitUtil.bits2words(myNumTaxa);
        for (int s = 0; s < myNumSitesToProcess; s++) {

            long[][] bits = new long[myNumDataRows][numLongs];
            byte[] cb = new byte[2];
            for (int l = 0; l < numLongs - 1; l++) {
                long bitmask = 0x1L;
                for (int t = l * 64, n = l * 64 + 64; t < n; t++) {
                    cb[0] = (byte) ((data[s][t] >>> 4) & 0xf);
                    cb[1] = (byte) (data[s][t] & 0xf);
                    for (int i = 0; i < 2; i++) {
                        if (cb[i] != Alignment.UNKNOWN_ALLELE) {
                            int index = myAlleleMappings[s][cb[i]];
                            if (index == -1) {
                                if (myRetainRareAlleles) {
                                    bits[myNumAlleles][l] |= bitmask;
                                }
                            } else {
                                bits[index][l] |= bitmask;
                            }
                        }
                    }
                    bitmask = bitmask << 1;
                }
            }

            int lastLong = numLongs - 1;
            int numRemaining = myNumTaxa % 64;
            if (numRemaining == 0) {
                numRemaining = 64;
            }
            long bitmask = 0x1L;
            for (int t = lastLong * 64, n = lastLong * 64 + numRemaining; t < n; t++) {
                cb[0] = (byte) ((data[s][t] >>> 4) & 0xf);
                cb[1] = (byte) (data[s][t] & 0xf);
                for (int i = 0; i < 2; i++) {
                    if (cb[i] != Alignment.UNKNOWN_ALLELE) {
                        int index = myAlleleMappings[s][cb[i]];
                        if (index == -1) {
                            if (myRetainRareAlleles) {
                                bits[myNumAlleles][lastLong] |= bitmask;
                            }
                        } else {
                            bits[index][lastLong] |= bitmask;
                        }
                    }
                }
                bitmask = bitmask << 1;
            }

            for (int i = 0; i < myNumDataRows; i++) {
                myData[i][mySite + s] = new OpenBitSet(bits[i], numLongs);
            }

        }

    }

    private void setTBits(byte[][] data) {

        if (mySite % 64 != 0) {
            throw new IllegalStateException("ProcessLineOfHapmap: setTBits: starting site must begin a word: " + mySite);
        }
        int wordNum = mySite / 64;

        long[] bits = new long[myNumDataRows];

        for (int t = 0; t < myNumTaxa; t++) {
            for (int j = 0; j < myNumDataRows; j++) {
                bits[j] = 0;
            }

            byte[] cb = new byte[2];
            long bitmask = 0x1L;
            for (int s = 0; s < myNumSitesToProcess; s++) {
                cb[0] = (byte) ((data[s][t] >>> 4) & 0xf);
                cb[1] = (byte) (data[s][t] & 0xf);
                for (int i = 0; i < 2; i++) {
                    if (cb[i] != Alignment.UNKNOWN_ALLELE) {
                        int index = myAlleleMappings[s][cb[i]];
                        if (index == -1) {
                            if (myRetainRareAlleles) {
                                bits[myNumAlleles] |= bitmask;
                            }
                        } else {
                            bits[index] |= bitmask;
                        }
                    }
                }
                bitmask = bitmask << 1;
            }

            for (int i = 0; i < myNumDataRows; i++) {
                myData[i][t].setLong(wordNum, bits[i]);
            }

        }

    }
}
