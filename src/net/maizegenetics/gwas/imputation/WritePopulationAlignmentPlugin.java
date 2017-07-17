package net.maizegenetics.gwas.imputation;

import java.awt.Frame;
import java.util.Arrays;
import java.util.List;

import javax.swing.ImageIcon;

import org.apache.log4j.Logger;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.FilterAlignment;
import net.maizegenetics.pal.alignment.MutableAlignment;
import net.maizegenetics.pal.alignment.MutableNucleotideAlignment;
import net.maizegenetics.pal.alignment.MutableSingleEncodeAlignment;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.PluginEvent;

public class WritePopulationAlignmentPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(WritePopulationAlignmentPlugin.class);
    boolean mergeAlignments = true;
    boolean writeParentCalls = true;
    boolean writeNucleotides = true;
    boolean outputDiploid = false;
    double minSnpCoverage = 0.1;
    double maxMafForMono = 0.01;
    String baseFileName;

    public WritePopulationAlignmentPlugin(Frame parentFrame) {
        super(parentFrame, false);
    }

    @Override
    public DataSet performFunction(DataSet input) {
        List<Datum> theData = input.getDataOfType(PopulationData.class);

        if (theData.size() > 0) {
            if (writeParentCalls) {
                writeOutput(theData, false);
            }
            if (writeNucleotides) {
                writeOutput(theData, true);
            }
            fireDataSetReturned(new PluginEvent(input, input.getCreator().getClass()));
            return input;
        } else {
            return null;
        }

    }

    private void writeOutput(List<Datum> theData, boolean asNucleotides) {
        String filename;
        if (mergeAlignments) {
            if (asNucleotides) {
                filename = baseFileName + "nuc.hmp.txt";
            } else {
                filename = baseFileName + "parents.hmp.txt";
            }
            Alignment[] allOfTheAlignments = new Alignment[theData.size()];
            int count = 0;
            for (Datum datum : theData) {
                PopulationData family = (PopulationData) datum.getData();
                allOfTheAlignments[count++] = createOutputAlignment(family, asNucleotides);
            }
            Alignment alignment = MutableSingleEncodeAlignment.getInstance(allOfTheAlignments);
            ExportUtils.writeToHapmap(alignment, outputDiploid, filename, '\t', null);
        } else {
            for (Datum datum : theData) {
                PopulationData family = (PopulationData) datum.getData();
                String familyName = family.name.replace('/', '.');
                if (asNucleotides) {
                    filename = baseFileName + ".family." + familyName + "nuc.hmp.txt";
                } else {
                    filename = baseFileName + ".family." + familyName + "parents.hmp.txt";
                }
                ExportUtils.writeToHapmap(createOutputAlignment(family, asNucleotides), outputDiploid, filename, '\t', null);
            }
        }

    }

    private Alignment createOutputAlignment(PopulationData popdata, boolean asNucleotides) {
        Alignment out = null;

        if (!asNucleotides) {
            out = popdata.imputed;
        } else {
            //change the parent calls to original nucleotides
            Alignment outPoly = NucleotideImputationUtils.convertParentCallsToNucleotides(popdata);

            if (!Double.isNaN(minSnpCoverage) && !Double.isNaN(maxMafForMono)) {
                int nsnps = popdata.original.getSiteCount();
                double ngametes = 2 * popdata.original.getSequenceCount();
                int[] monomorphicSnps = new int[nsnps];
                int snpCount = 0;
                for (int s = 0; s < nsnps; s++) {
                    double coverage = popdata.original.getTotalGametesNotMissing(s) / ngametes;
                    if (!popdata.snpIndex.fastGet(s) && popdata.original.getMinorAlleleFrequency(s) <= maxMafForMono && coverage >= minSnpCoverage) {
                        monomorphicSnps[snpCount++] = s;
                    }
                }
                monomorphicSnps = Arrays.copyOf(monomorphicSnps, snpCount);
                Alignment fa = FilterAlignment.getInstance(popdata.original, monomorphicSnps);
                if (fa.getSiteCount() == 0) {	//If there are no monomorphic sites (e.g, have been pre-filtered), just return polymorphic ones
                    out = MutableSingleEncodeAlignment.getInstance(new Alignment[]{outPoly});
                } else { //Return both monomorphic and polymorphic sites
                    MutableAlignment outMono = MutableNucleotideAlignment.getInstance(fa);
                    // fill in all values with the major allele
                    nsnps = outMono.getSiteCount();
                    int ntaxa = outMono.getSequenceCount();
                    for (int s = 0; s < nsnps; s++) {
                        byte majorAllele = outMono.getMajorAllele(s);
                        byte major = (byte) ((majorAllele << 4) | majorAllele);
                        for (int t = 0; t < ntaxa; t++) {
                            outMono.setBase(t, s, major);
                        }
                    }
                    out = MutableSingleEncodeAlignment.getInstance(new Alignment[]{outPoly, outMono});
                }
            }
        }
        return out;
    }

    @Override
    public void setParameters(String[] args) {
        if (args == null || args.length == 0) {
            myLogger.info(getUsage());
            return;
        }

        int narg = args.length;
        for (int i = 0; i < narg; i++) {
            if (args[i].equals("-f") || args[i].equalsIgnoreCase("-file")) {
                baseFileName = args[++i];
            } else if (args[i].equals("-m") || args[i].equalsIgnoreCase("-merge")) {
                String val = args[++i];
                if (val.toUpperCase().startsWith("T")) {
                    mergeAlignments = true;
                } else {
                    mergeAlignments = false;
                }
            } else if (args[i].equals("-p") || args[i].equalsIgnoreCase("-parentCalls")) {
                String val = args[++i];
                if (val.toUpperCase().startsWith("T")) {
                    writeParentCalls = true;
                    writeNucleotides = false;
                } else {
                    writeParentCalls = false;
                    writeNucleotides = true;
                }
            } else if (args[i].equals("-o") || args[i].equalsIgnoreCase("-outputType")) {
                String val = args[++i];
                if (val.toUpperCase().startsWith("P")) {
                    writeParentCalls = true;
                    writeNucleotides = false;
                } else if (val.toUpperCase().startsWith("N")) {
                    writeParentCalls = false;
                    writeNucleotides = true;
                } else if (val.toUpperCase().startsWith("B")) {
                    writeParentCalls = true;
                    writeNucleotides = true;
                } else {
                    writeParentCalls = true;
                    writeNucleotides = false;
                }
            } else if (args[i].equals("-d") || args[i].equalsIgnoreCase("-diploid")) {
                String val = args[++i];
                if (val.toUpperCase().startsWith("T")) {
                    outputDiploid = true;
                } else {
                    outputDiploid = false;
                }
            } else if (args[i].equals("-c") || args[i].equalsIgnoreCase("-minCoverage")) {
                minSnpCoverage = Double.parseDouble(args[++i]);
            } else if (args[i].equals("-x") || args[i].equalsIgnoreCase("-maxMono")) {
                maxMafForMono = Double.parseDouble(args[++i]);
            } else if (args[i].equals("?")) {
                myLogger.info(getUsage());
            }
        }
    }

    public void setMergeAlignments(boolean mergeAlignments) {
        this.mergeAlignments = mergeAlignments;
    }

    public void setWriteParentCalls(boolean writeParentCalls) {
        this.writeParentCalls = writeParentCalls;
    }

    public void setOutputDiploid(boolean outputDiploid) {
        this.outputDiploid = outputDiploid;
    }

    public void setBaseFileName(String baseFileName) {
        this.baseFileName = baseFileName;
    }

    @Override
    public ImageIcon getIcon() {
        return null;
    }

    @Override
    public String getButtonName() {
        return "Write Populations";
    }

    @Override
    public String getToolTipText() {
        return null;
    }

    public String getUsage() {
        StringBuilder usage = new StringBuilder("The WritePopulationAlignmentPlugin requires the following parameter:\n");
        usage.append("-f or -file : The base file name for the ouput. .hmp.txt will be appended.\n");
        usage.append("The following parameters are optional:\n");
        usage.append("-m or -merge : if true families are merged into a single file, if false each family is output to a separate file (default = true)\n");
        usage.append("-o or -outputType : parents = output parent calls, nucleotides = output nucleotides, both = output both\n");
        usage.append("-d or -diploid : if true output is AA/CC/AC, if false output is A/C/M\n");
        usage.append("-c or -minCoverage : the minimum coverage for a monomorphic snp to be included in the nucleotide output (default = NaN, no monomorphic SNPs included)\n");
        usage.append("-x or -maxMono : the maximum minor allele frequency used to call monomorphic snps (default = NaN, none called)\n");
        usage.append("? : print the parameter list.\n");

        return usage.toString();
    }

	public void setWriteNucleotides(boolean writeNucleotides) {
		this.writeNucleotides = writeNucleotides;
	}

	public void setMinSnpCoverage(double minSnpCoverage) {
		this.minSnpCoverage = minSnpCoverage;
	}

	public void setMaxMafForMono(double maxMafForMono) {
		this.maxMafForMono = maxMafForMono;
	}
}
