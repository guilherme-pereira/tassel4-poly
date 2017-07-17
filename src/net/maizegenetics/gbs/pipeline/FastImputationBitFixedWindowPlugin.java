/*
 * FastImputationBitFixedWindowPlugin
 */
package net.maizegenetics.gbs.pipeline;

import net.maizegenetics.util.ArgsEngine;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import java.awt.Frame;
import java.io.File;
import java.util.HashMap;
import javax.swing.ImageIcon;
import net.maizegenetics.baseplugins.ConvertSBitTBitPlugin;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.ImportUtils;

import org.apache.log4j.Logger;

/**
 * 
 * @deprecated Replaced by better methods {@link MinorWindowViterbiImputationPlugin}
 * @author Ed Buckler
 */
@Deprecated
public class FastImputationBitFixedWindowPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(FastImputationBitFixedWindowPlugin.class);
    private ArgsEngine myArgsEngine = null;
    private String myInputFile;
    private String myOutputFile;
    private boolean usePedigree = false;
    HashMap<String, Double> taxaFs = null;

    public FastImputationBitFixedWindowPlugin() {
        super(null, false);
    }

    public FastImputationBitFixedWindowPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    @Override
    public DataSet performFunction(DataSet input) {
        myLogger.info("Reading from input file: " + myInputFile);
        Alignment a = ImportUtils.readFromHapmap(myInputFile, null);
        myLogger.info("Performing imputation...");
        //a = TBitAlignment.getInstance(a);
        a = ConvertSBitTBitPlugin.convertAlignment(a, ConvertSBitTBitPlugin.CONVERT_TYPE.tbit, this);
        FastImputationBitFixedWindow fi;
        if (usePedigree) {
            boolean highHetTaxa[] = markHighHetTaxa(a);
            fi = new FastImputationBitFixedWindow(a, highHetTaxa);
        } else {
            fi = new FastImputationBitFixedWindow(a);
        }
        myLogger.info("Output file: " + myOutputFile);
        fi.writeAlignment(myOutputFile);
        return null;
    }

    private void printUsage() {
        myLogger.info(
                "\nUsage:\n"
                + "-hmp  Input Hapmap File\n"
                + "-o    Output File\n"
                + "-p    Pedigree file containing full sample names (or expected names after merging) & expected inbreeding\n"
                + "        coefficient (F) for each.  Only highly inbred taxa, with F >= 0.8 (e.g., S3 or more), will be used\n"
                + "        for imputation (default: calculate the heterozygosity for each taxon from the data)\n");
    }

    @Override
    public void setParameters(String[] args) {

        if (args.length == 0) {
            printUsage();
            throw new IllegalArgumentException("setParameters: no arguments given.");
        }

        if (myArgsEngine == null) {
            myArgsEngine = new ArgsEngine();
            myArgsEngine.add("-hmp", "--hmp-file", true);
            myArgsEngine.add("-o", "--output-file", true);
            myArgsEngine.add("-p", "--pedigree-file", true);
        }
        myArgsEngine.parse(args);

        myInputFile = myArgsEngine.getString("-hmp");
        if ((myInputFile == null) || myInputFile.length() == 0) {
            printUsage();
            throw new IllegalArgumentException("setParameters: no input file specified.");
        }

        myOutputFile = myArgsEngine.getString("-o");
        if ((myOutputFile == null) || myOutputFile.length() == 0) {
            printUsage();
            throw new IllegalArgumentException("setParameters: no output file specified.");
        }

        if (myArgsEngine.getBoolean("-p")) {
            String pedigreeFileStr = myArgsEngine.getString("-p");
            File pedigreeFile = new File(pedigreeFileStr);
            if (!pedigreeFile.exists() || !pedigreeFile.isFile()) {
                printUsage();
                throw new IllegalArgumentException("Can't find the pedigree input file (-p option: " + pedigreeFileStr + ").");
            }
            taxaFs = DiscoverySNPCallerPlugin.readTaxaFsFromFile(pedigreeFile);
            if (taxaFs == null) {
                throw new IllegalArgumentException("Problem reading the pedigree file. Progam aborted.");
            }
            usePedigree = true;
        }
    }

    private boolean[] markHighHetTaxa(Alignment a) {
        boolean[] highHet = new boolean[a.getSequenceCount()];
        int nHighHetTaxa = 0;
        for (int t = 0; t < a.getSequenceCount(); t++) {
            if (taxaFs.containsKey(a.getFullTaxaName(t))) {
                if (taxaFs.get(a.getFullTaxaName(t)) < 0.8) {
                    highHet[t] = true;
                    nHighHetTaxa++;
                }
            }
        }
        myLogger.info(nHighHetTaxa + " heterozygous taxa (with an Expected F < 0.8) were found in the input genotype file");
        return highHet;
    }

    @Override
    public ImageIcon getIcon() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getButtonName() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getToolTipText() {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
