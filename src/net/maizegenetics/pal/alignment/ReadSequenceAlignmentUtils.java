// ReadSequenceAlignmentUtils.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)
package net.maizegenetics.pal.alignment;

import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.pal.io.FormattedInput;
import net.maizegenetics.pal.io.InputSource;

import java.io.IOException;
import java.io.PushbackReader;
import net.maizegenetics.prefs.TasselPrefs;

/**
 * reads aligned sequence data from plain text files.<p>
 *
 * recognizes PHYLIP 3.4 INTERLEAVED,
 *              PHYLIP SEQUENTIAL,
 *              CLUSTAL and derived formats.<p>
 *
 * Other features:
 * - the dot as "copy character" is recognized,
 * - all base characters are capitalized,
 * - automatic data type estimation
 * - determination of corresponding base frequencies.
 *
 * @version $Id: ReadSequenceAlignmentUtils.java,v 1.1 2009/07/23 19:39:41 tcasstevens Exp $
 *
 * @author Korbinian Strimmer
 * @author Alexei Drummond
 */
public class ReadSequenceAlignmentUtils {

    private ReadSequenceAlignmentUtils() {
        // to prevent instantiation of this utility class.
    }

    /** read from stream, while letting the user set the maximum label (taxa) label length
     * as non-standard phylip formats now exist*/
    public static Alignment readBasicAlignments(PushbackReader input, int maxLabelLength)
            throws IOException {
        Alignment saa = null;
        saa = readFile(input, maxLabelLength);
        return saa;
    }

    /** read from file, while letting the user set the maximum label (taxa) label length
     * as non-standard phylip formats now exist*/
    public static Alignment readBasicAlignments(String file, int maxLabelLength)
            throws IOException {
        PushbackReader input = InputSource.openFile(file);
        Alignment saa = readBasicAlignments(input, maxLabelLength);
        input.close();
        return saa;
    }

    private static final boolean isType(PushbackReader in, String id) throws IOException {
        FormattedInput fi = FormattedInput.getInstance();
        for (int i = 0; i < id.length(); i++) {
            int c = fi.readNextChar(in);
            if (c != id.charAt(i)) {
                in.unread(c);
                return false;
            }
        }
        return true;
    }

    private static Alignment readFile(PushbackReader in, int maxLabelLength) throws IOException {
        return readPHYLIP(in, maxLabelLength);
    }

    // Read alignment (in PHYLIP 3.4 INTERLEAVED or PHYLIP SEQUENTIAL format)
    private static Alignment readPHYLIP(PushbackReader in, int maxLabelLength) {
        FormattedInput fi = FormattedInput.getInstance();
        IdGroup idGroup;
        int numSeqs = 0, numSites = 0, lineLength = 0;
        char[][] data = null;
        int c, pos = 0, seq = 0;

        try {
            // Parse PHYLIP header line
            numSeqs = fi.readInt(in);
            numSites = fi.readInt(in);

            String[] identifiers = new String[numSeqs];
            data = new char[numSeqs][numSites];


            // Determine whether sequences are in INTERLEAVED
            // or in sequential format
            String header = fi.readLine(in, false);

            boolean interleaved = true;

            if (header.length() > 0) {
                if (header.charAt(0) == 'S') {
                    interleaved = false;
                }
            }

            if (interleaved) // PHYLIP INTERLEAVED
            {

                // Reading data
                while (pos < numSites) {
                    // Go to next block
                    c = fi.readNextChar(in);
                    in.unread(c);

                    for (seq = 0; seq < numSeqs; seq++) {
                        lineLength = readSeqLineP(in, seq, pos, numSites, data, identifiers, fi, maxLabelLength, lineLength);
                    }
                    pos += lineLength;
                }
            } else // PHYLIP SEQUENTIAL
            {
                //System.out.println("PHYLIP SEQUENTIAL");

                for (seq = 0; seq < numSeqs; seq++) {
                    // Go to next block
                    c = fi.readNextChar(in);
                    in.unread(c);

                    // Read label
                    identifiers[seq] = fi.readLabel(in, maxLabelLength).toUpperCase();

                    // Read sequences
                    for (pos = 0; pos < numSites; pos++) {
                        data[seq][pos] = (char) fi.readNextChar(in);

                        if (data[0][pos] == '.') {
                            if (seq == 0) {
                                throw new IllegalArgumentException(
                                        "Copy character (.) in first sequence not allowed (pos. " + (pos + 1) + ")");
                            } else {
                                data[seq][pos] = data[0][pos];
                            }
                        }
                    }
                }
            }
            idGroup = new SimpleIdGroup(identifiers);
        } catch (IOException e) {
            throw new IllegalArgumentException("IO error after pos. " + (pos + 1) + ", seq. " + (seq + 1));
        }
        String[] s = new String[numSeqs];
        for (int i = 0; i < numSeqs; i++) {
            s[i] = (new String(data[i])).toUpperCase();
        }
        // SimpleAlignment saa = SimpleAlignment.getInstance(idGroup, s, new Nucleotides());
        String[] sites = new String[numSites];
        int[] positions = new int[numSites];
        for (int i = 0; i < numSites; i++) {
            positions[i] = i;
            sites[i] = Integer.toString(i);
        }
        return BitAlignment.getNucleotideInstance(idGroup, s, null, null, null, TasselPrefs.getAlignmentMaxAllelesToRetain(), new Locus[] {Locus.UNKNOWN}, new int[]{0}, null, TasselPrefs.getAlignmentRetainRareAlleles(), true);
    }

    private static int readSeqLineP(PushbackReader in, int s, int pos, int maxPos, char[][] data, String[] identifiers,
            FormattedInput fi, int maxLabelLength, int lineLength)
            throws IOException {
        if (pos == 0) {
            identifiers[s] = fi.readLabel(in, maxLabelLength).toUpperCase();
        }

        if (s == 0) {
            String thisLine = fi.readLine(in, false);

            if (thisLine.length() > maxPos - pos) {
                lineLength = maxPos - pos;
            } else {
                lineLength = thisLine.length();
            }

            for (int i = 0; i < lineLength; i++) {
                data[0][pos + i] = thisLine.charAt(i);
                if (data[0][pos + i] == '.') {
                    throw new IllegalArgumentException("Copy character (.) in first sequence not allowed (pos. " + (i + pos + 1) + ")");
                }
            }
        } else {
            for (int i = 0; i < lineLength; i++) {
                data[s][pos + i] = (char) fi.readNextChar(in);
                if (data[s][pos + i] == '.') {
                    data[s][pos + i] = data[0][pos + i];
                }
            }
            fi.nextLine(in);
        }
        return lineLength;
    }
}
