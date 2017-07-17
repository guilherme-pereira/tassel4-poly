// ReadDistanceMatrix.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)
package net.maizegenetics.pal.distance;

import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.pal.io.FormattedInput;
import net.maizegenetics.pal.io.InputSource;

import java.io.IOException;
import java.io.PushbackReader;

/**
 * reads pairwise distance matrices in PHYLIP format
 *  (full matrix)
 *
 * @version $Id: ReadDistanceMatrix.java,v 1.4 2009/07/02 20:26:13 pjbradbury Exp $
 *
 * @author Korbinian Strimmer
 * @author Alexei Drummond
 */
public class ReadDistanceMatrix {

    private ReadDistanceMatrix() {
    }

    /** read from stream */
    public static DistanceMatrix readDistanceMatrix(PushbackReader input)
            throws DistanceParseException {
        return readSquare(input);
    }

    /** read from file */
    public static DistanceMatrix readDistanceMatrix(String file)
            throws DistanceParseException, IOException {
        PushbackReader input = null;
        try {
            input = InputSource.openFile(file);
            return readSquare(input);
        } finally {
            input.close();
        }
    }

    // Read square matrix
    private static DistanceMatrix readSquare(PushbackReader in)
            throws DistanceParseException {
        FormattedInput fi = FormattedInput.getInstance();
        try {
            // Parse PHYLIP header line
            int numSeqs = fi.readInt(in);
            fi.nextLine(in);

            // Read distance and sequence names
            double[][] distance = new double[numSeqs][numSeqs];
            Identifier[] ids = new Identifier[numSeqs];
            for (int i = 0; i < numSeqs; i++) {
                ids[i] = new Identifier(fi.readWord(in).trim());
                for (int j = 0; j < numSeqs; j++) {
                    distance[i][j] = fi.readDouble(in);
                }
                fi.nextLine(in);
            }
            IdGroup idGroup = new SimpleIdGroup(ids);
            return new DistanceMatrix(distance, idGroup);
        } catch (IOException e) {
            throw new DistanceParseException("IO error");
        } catch (NumberFormatException e) {
            throw new DistanceParseException("Number format error");
        }
    }
}
