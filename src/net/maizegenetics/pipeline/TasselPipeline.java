/*
 * TasselPipeline.java
 *
 * Created on June 11, 2009
 *
 */
package net.maizegenetics.pipeline;

import java.awt.Frame;

import java.io.BufferedReader;

import java.lang.reflect.Constructor;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Vector;
import java.util.regex.Pattern;

import net.maizegenetics.baseplugins.AbstractDisplayPlugin;
import net.maizegenetics.baseplugins.CombineDataSetsPlugin;
import net.maizegenetics.baseplugins.ConvertAlignmentCoordinatesPlugin;
import net.maizegenetics.baseplugins.ConvertSBitTBitPlugin;
import net.maizegenetics.baseplugins.CreateTreePlugin;
import net.maizegenetics.baseplugins.DistanceMatrixPlugin;
import net.maizegenetics.baseplugins.DistanceMatrixRangesPlugin;
import net.maizegenetics.baseplugins.ExportMultiplePlugin;
import net.maizegenetics.baseplugins.FileLoadPlugin;
import net.maizegenetics.baseplugins.FilterAlignmentPlugin;
import net.maizegenetics.baseplugins.FilterSiteNamePlugin;
import net.maizegenetics.baseplugins.FilterTaxaAlignmentPlugin;
import net.maizegenetics.baseplugins.FilterTraitsPlugin;
import net.maizegenetics.baseplugins.FixedEffectLMPlugin;
import net.maizegenetics.baseplugins.FlapjackLoadPlugin;
import net.maizegenetics.baseplugins.GenotypeImputationPlugin;
import net.maizegenetics.baseplugins.GenotypeSummaryPlugin;
import net.maizegenetics.baseplugins.IntersectionAlignmentPlugin;
import net.maizegenetics.baseplugins.KinshipPlugin;
import net.maizegenetics.baseplugins.LinkageDiseqDisplayPlugin;
import net.maizegenetics.baseplugins.LinkageDisequilibriumPlugin;
import net.maizegenetics.baseplugins.MLMPlugin;
import net.maizegenetics.baseplugins.MergeAlignmentsPlugin;
import net.maizegenetics.baseplugins.MergeAlignmentsSameSitesPlugin;
import net.maizegenetics.baseplugins.NumericalGenotypePlugin;
import net.maizegenetics.baseplugins.PlinkLoadPlugin;
import net.maizegenetics.baseplugins.SeparatePlugin;
import net.maizegenetics.baseplugins.SequenceDiversityPlugin;
import net.maizegenetics.baseplugins.SynonymizerPlugin;
import net.maizegenetics.baseplugins.TableDisplayPlugin;
import net.maizegenetics.baseplugins.UnionAlignmentPlugin;
import net.maizegenetics.baseplugins.genomicselection.RidgeRegressionEmmaPlugin;
import net.maizegenetics.gbs.maps.TagsOnPhysMapHDF5;
import net.maizegenetics.gbs.maps.TagsOnPhysicalMap;

import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.gui.LinkageDisequilibriumComponent;
import net.maizegenetics.pal.popgen.LinkageDisequilibrium.testDesign;

import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.pal.popgen.LinkageDisequilibrium.HetTreatment;

import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.plugindef.Plugin;
import net.maizegenetics.plugindef.PluginEvent;
import net.maizegenetics.plugindef.PluginListener;
import net.maizegenetics.plugindef.ThreadedPluginListener;

import net.maizegenetics.prefs.TasselPrefs;
import net.maizegenetics.progress.ProgressPanel;
import net.maizegenetics.tassel.DataTreePanel;
import net.maizegenetics.tassel.TASSELMainFrame;
import net.maizegenetics.util.ExceptionUtils;
import net.maizegenetics.util.Utils;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.apache.log4j.PropertyConfigurator;

/**
 *
 * @author terryc
 */
public class TasselPipeline implements PluginListener {

    private static final Logger myLogger = Logger.getLogger(TasselPipeline.class);
    private final TASSELMainFrame myMainFrame;
    private final Map<String, List> myForks = new LinkedHashMap<String, List>();
    private String myCurrentFork = null;
    private List<Plugin> myCurrentPipe = null;
    private Plugin myFirstPlugin = null;
    private final List<ThreadedPluginListener> myThreads = new ArrayList();
    private final Map<Plugin, Integer> myProgressValues = new HashMap<Plugin, Integer>();

    /**
     * Creates a new instance of TasselPipeline
     */
    public TasselPipeline(String args[], TASSELMainFrame frame) {

        myMainFrame = frame;

        java.util.Properties props = new java.util.Properties();
        props.setProperty("log4j.logger.net.maizegenetics", "INFO, stdout");
        props.setProperty("log4j.appender.stdout",
                "org.apache.log4j.ConsoleAppender");
        props.setProperty("log4j.appender.stdout.layout",
                "org.apache.log4j.TTCCLayout");
        PropertyConfigurator.configure(props);

        try {

            if ((args.length == 1) && (args[0].equalsIgnoreCase("-versionComment"))) {
                System.out.println("Version " + myMainFrame.version + " on " + myMainFrame.versionDate);
                return;
            }

            if ((args.length == 1) && (args[0].equalsIgnoreCase("-versionTag"))) {
                System.out.println("V" + myMainFrame.version);
                return;
            }

            myLogger.info("Tassel Version: " + myMainFrame.version + "  Date: " + myMainFrame.versionDate);
            myLogger.info("Max Available Memory Reported by JVM: " + Utils.getMaxHeapSizeMB() + " MB");

            parseArgs(args);

            if (myMainFrame != null) {
                ProgressPanel progressPanel = myMainFrame.getProgressPanel();
                if (progressPanel != null) {
                    Iterator itr = myForks.keySet().iterator();
                    while (itr.hasNext()) {
                        String key = (String) itr.next();
                        List<Plugin> current = myForks.get(key);
                        progressPanel.addPipelineSegment(current);
                    }
                }
            }

            for (int i = 0; i < myThreads.size(); i++) {
                ThreadedPluginListener current = (ThreadedPluginListener) myThreads.get(i);
                current.start();
            }

        } catch (Exception e) {
            System.exit(1);
        }

    }

    public static void main(String args[]) {
        if ((args.length >= 2) && (args[0].equalsIgnoreCase("-createXML"))) {
            String xmlFilename = args[1].trim();
            String[] temp = new String[args.length - 2];
            System.arraycopy(args, 2, temp, 0, temp.length);
            TasselPipelineXMLUtil.writeArgsAsXML(xmlFilename, temp);
        } else if ((args.length >= 2) && (args[0].equalsIgnoreCase("-translateXML"))) {
            String xmlFilename = args[1].trim();
            String[] result = TasselPipelineXMLUtil.readXMLAsArgs(xmlFilename);
            for (int i = 0; i < result.length; i++) {
                System.out.print(result[i]);
                System.out.print(" ");
            }
            System.out.println("");
        } else {
            TasselPrefs.setPersistPreferences(false);
            new TasselPipeline(args, null);
        }
    }

    public void parseArgs(String[] args) {

        if ((args.length >= 2) && (args[0].equalsIgnoreCase("-configFile"))) {
            String xmlFilename = args[1].trim();
            args = TasselPipelineXMLUtil.readXMLAsArgs(xmlFilename);
        }

        int index = 0;
        while (index < args.length) {

            try {

                String current = args[index++];
                current = current.replaceFirst("—", "-");

                if (!current.startsWith("-")) {
                    throw new IllegalArgumentException("TasselPipeline: parseArgs: expecting argument beginning with dash: " + current);
                }

                if (current.startsWith("-runfork")) {
                    String key = current.replaceFirst("-runfork", "-fork");
                    List specifiedPipe = (List) myForks.get(key);
                    if (specifiedPipe == null) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: unknown fork: " + current);
                    } else if (specifiedPipe.size() == 0) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: empty fork: " + current);
                    } else {
                        PluginEvent event = new PluginEvent(new DataSet((Datum) null, null));
                        ThreadedPluginListener thread = new ThreadedPluginListener((Plugin) specifiedPipe.get(0), event);
                        myThreads.add(thread);
                    }
                } else if (current.startsWith("-fork")) {
                    if ((myCurrentPipe != null) && (myCurrentPipe.size() != 0)) {
                        myCurrentPipe.get(myCurrentPipe.size() - 1).setThreaded(true);
                    }
                    myCurrentFork = current;
                    myCurrentPipe = new ArrayList<Plugin>();
                    myForks.put(myCurrentFork, myCurrentPipe);
                } else if (current.startsWith("-input")) {
                    String key = current.replaceFirst("-input", "-fork");
                    List specifiedPipe = (List) myForks.get(key);
                    if (specifiedPipe == null) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: unknown input: " + current);
                    } else {
                        Plugin lastCurrentPipe = null;
                        try {
                            lastCurrentPipe = (Plugin) myCurrentPipe.get(myCurrentPipe.size() - 1);
                        } catch (Exception e) {
                            throw new IllegalArgumentException("TasselPipeline: parseArgs: -input must come after plugin in current fork.");
                        }
                        Plugin endSpecifiedPipe = (Plugin) specifiedPipe.get(specifiedPipe.size() - 1);
                        lastCurrentPipe.receiveInput(endSpecifiedPipe);
                    }
                } else if (current.startsWith("-inputOnce")) {
                    String key = current.replaceFirst("-input", "-fork");
                    List specifiedPipe = (List) myForks.get(key);
                    if (specifiedPipe == null) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: unknown input: " + current);
                    } else {
                        CombineDataSetsPlugin combinePlugin = null;
                        try {
                            combinePlugin = (CombineDataSetsPlugin) myCurrentPipe.get(myCurrentPipe.size() - 1);
                        } catch (Exception e) {
                            throw new IllegalArgumentException("TasselPipeline: parseArgs: -inputOnce must follow -combine flag.");
                        }
                        Plugin endSpecifiedPipe = (Plugin) specifiedPipe.get(specifiedPipe.size() - 1);
                        combinePlugin.receiveDataSetOnceFrom(endSpecifiedPipe);
                    }
                } else if (current.startsWith("-combine")) {
                    current = current.replaceFirst("-combine", "-fork");
                    if ((myCurrentPipe != null) && (myCurrentPipe.size() != 0)) {
                        myCurrentPipe.get(myCurrentPipe.size() - 1).setThreaded(true);
                    }
                    myCurrentFork = current;
                    myCurrentPipe = new ArrayList<Plugin>();
                    myForks.put(myCurrentFork, myCurrentPipe);
                    integratePlugin(new CombineDataSetsPlugin(), false);
                } else if (current.equalsIgnoreCase("-t")) {
                    String traitFile = args[index++].trim();
                    loadFile(traitFile, FileLoadPlugin.TasselFileType.Numerical);
                } else if (current.equalsIgnoreCase("-s")) {
                    String inputFile = args[index++].trim();
                    loadFile(inputFile, FileLoadPlugin.TasselFileType.Sequence);
                } else if (current.equalsIgnoreCase("-p")) {
                    String inputFile = args[index++].trim();
                    loadFile(inputFile, FileLoadPlugin.TasselFileType.Polymorphism);
                } else if (current.equalsIgnoreCase("-a")) {
                    String inputFile = args[index++].trim();
                    loadFile(inputFile, FileLoadPlugin.TasselFileType.Annotated);
                } else if (current.equalsIgnoreCase("-k")) {
                    String kinshipFile = args[index++].trim();
                    loadFile(kinshipFile, FileLoadPlugin.TasselFileType.SqrMatrix);
                } else if (current.equalsIgnoreCase("-q")) {
                    String populationFile = args[index++].trim();
                    loadFile(populationFile, FileLoadPlugin.TasselFileType.Phenotype);
                } else if (current.equalsIgnoreCase("-h")) {
                    String hapFile = args[index++].trim();
                    loadFile(hapFile, FileLoadPlugin.TasselFileType.Hapmap);
                } else if (current.equalsIgnoreCase("-h5")) {
                    String hdf5File = args[index++].trim();
                    loadFile(hdf5File, FileLoadPlugin.TasselFileType.HDF5);
                } else if (current.equalsIgnoreCase("-r")) {
                    String phenotypeFile = args[index++].trim();
                    loadFile(phenotypeFile, FileLoadPlugin.TasselFileType.Phenotype);
                } else if (current.equalsIgnoreCase("-plink")) {
                    String pedFile = null;
                    String mapFile = null;
                    for (int i = 0; i < 2; i++) {
                        String fileType = args[index++].trim();
                        String filename = args[index++].trim();
                        if (fileType.equalsIgnoreCase("-ped")) {
                            pedFile = filename;
                        } else if (fileType.equalsIgnoreCase("-map")) {
                            mapFile = filename;
                        } else {
                            throw new IllegalArgumentException("TasselPipeline: parseArgs: -plink: unknown file type: " + fileType);
                        }
                    }
                    if ((pedFile == null) || (mapFile == null)) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: -plink must specify both ped and map files.");
                    }
                    PlinkLoadPlugin plugin = new PlinkLoadPlugin(myMainFrame, false);
                    plugin.setPedFile(pedFile);
                    plugin.setMapFile(mapFile);
                    integratePlugin(plugin, true);
                } else if (current.equalsIgnoreCase("-flapjack")) {
                    String genoFile = null;
                    String mapFile = null;
                    for (int i = 0; i < 2; i++) {
                        String fileType = args[index++].trim();
                        String filename = args[index++].trim();
                        if (fileType.equalsIgnoreCase("-geno")) {
                            genoFile = filename;
                        } else if (fileType.equalsIgnoreCase("-map")) {
                            mapFile = filename;
                        } else {
                            throw new IllegalArgumentException("TasselPipeline: parseArgs: -flapjack: unknown file type: " + fileType);
                        }
                    }
                    if ((genoFile == null) || (mapFile == null)) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: -flapjack must specify both ped and map files.");
                    }
                    FlapjackLoadPlugin plugin = new FlapjackLoadPlugin(myMainFrame, false);
                    plugin.setGenoFile(genoFile);
                    plugin.setMapFile(mapFile);
                    integratePlugin(plugin, true);
                } else if (current.equalsIgnoreCase("-fasta")) {
                    String fastaFile = args[index++].trim();
                    loadFile(fastaFile, FileLoadPlugin.TasselFileType.Fasta);
                } else if (current.equalsIgnoreCase("-geneticMap")) {
                    String geneticMapFile = args[index++].trim();
                    loadFile(geneticMapFile, FileLoadPlugin.TasselFileType.GeneticMap);
                } else if (current.equalsIgnoreCase("-table")) {
                    String tableFile = args[index++].trim();
                    loadFile(tableFile, FileLoadPlugin.TasselFileType.Table);
                } else if (current.equalsIgnoreCase("-vcf")) {
                    String vcfFile = args[index++].trim();
                    loadFile(vcfFile, FileLoadPlugin.TasselFileType.VCF);
                } else if (current.equalsIgnoreCase("-readSerialAlignment")) {
                    String file = args[index++].trim();
                    loadFile(file, FileLoadPlugin.TasselFileType.Serial);
                } else if (current.equalsIgnoreCase("-importGuess")) {
                    String file = args[index++].trim();
                    loadFile(file, FileLoadPlugin.TasselFileType.Unknown);
                } else if (current.equalsIgnoreCase("-convertTOPMtoHDF5")) {
                    String filename = args[index++].trim();
                    TagsOnPhysicalMap topm = null;
                    String h5Filename = null;
                    if (filename.endsWith(".topm.txt")) {
                        topm = new TagsOnPhysicalMap(filename, false);
                        h5Filename = filename.replaceAll(".topm.txt", ".topm.h5");
                    } else if (filename.endsWith(".topm.bin")) {
                        topm = new TagsOnPhysicalMap(filename, true);
                        h5Filename = filename.replaceAll(".topm.bin", ".topm.h5");
                    } else {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: -convertTOPMtoHDF5: Unknown file extension: " + filename);
                    }
                    TagsOnPhysMapHDF5.createFile(topm, h5Filename, 4, 8);
                } else if (current.equalsIgnoreCase("-taxaJoinStrict")) {
                    String temp = args[index++].trim();
                    if (temp.equalsIgnoreCase("false")) {
                        temp = TasselPrefs.TASSEL_IDENTIFIER_JOIN_TYPES.NonStrict.toString();
                    } else if (temp.equalsIgnoreCase("true")) {
                        temp = TasselPrefs.TASSEL_IDENTIFIER_JOIN_TYPES.Strict.toString();
                    }
                    TasselPrefs.TASSEL_IDENTIFIER_JOIN_TYPES type = null;
                    try {
                        type = TasselPrefs.TASSEL_IDENTIFIER_JOIN_TYPES.valueOf(temp);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: -taxaJoinStrict: Unknown type: " + temp + "  Should be: " + Arrays.toString(TasselPrefs.TASSEL_IDENTIFIER_JOIN_TYPES.values()));
                    }
                    TasselPrefs.putIDJoinStrict(type);
                } else if (current.equalsIgnoreCase("-taxaJoinNumLevels")) {
                    String temp = args[index++].trim();
                    int numLevels = 0;
                    try {
                        numLevels = Integer.parseInt(temp);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: Problem parsing Taxa Join Number Levels: " + temp);
                    }
                    TasselPrefs.putIDJoinNumLevels(numLevels);
                } else if (current.equalsIgnoreCase("-retainRareAlleles")) {
                    String temp = args[index++].trim();
                    boolean retain = true;
                    if (temp.equalsIgnoreCase("false")) {
                        retain = false;
                    } else if (temp.equalsIgnoreCase("true")) {
                        retain = true;
                    } else {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: -retainRareAlleles parameter must be true or false.");
                    }
                    TasselPrefs.putAlignmentRetainRareAlleles(retain);
                } else if (current.equalsIgnoreCase("-maxAllelesToRetain")) {
                    String temp = args[index++].trim();
                    int maxAlleles = 0;
                    try {
                        maxAlleles = Integer.parseInt(temp);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: Problem parsing max alleles to retain: " + temp);
                    }
                    TasselPrefs.putAlignmentMaxAllelesToRetain(maxAlleles);
                } else if (current.equalsIgnoreCase("-optimizeForTaxa")) {
                    FileLoadPlugin plugin = (FileLoadPlugin) findLastPluginFromCurrentPipe(new Class[]{FileLoadPlugin.class});
                    if (plugin == null) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No File Load step defined: " + current);
                    }
                    plugin.setIsFileCreatedSBit(false);
                } else if (current.equalsIgnoreCase("-convertToSiteOpt")) {
                    ConvertSBitTBitPlugin plugin = new ConvertSBitTBitPlugin(myMainFrame, false);
                    plugin.setType(ConvertSBitTBitPlugin.CONVERT_TYPE.sbit);
                    integratePlugin(plugin, true);
                } else if (current.equalsIgnoreCase("-convertToTaxaOpt")) {
                    ConvertSBitTBitPlugin plugin = new ConvertSBitTBitPlugin(myMainFrame, false);
                    plugin.setType(ConvertSBitTBitPlugin.CONVERT_TYPE.tbit);
                    integratePlugin(plugin, true);
                } else if (current.equalsIgnoreCase("-union")) {
                    UnionAlignmentPlugin plugin = new UnionAlignmentPlugin(myMainFrame, false);
                    integratePlugin(plugin, true);
                } else if (current.equalsIgnoreCase("-intersect")) {
                    IntersectionAlignmentPlugin plugin = new IntersectionAlignmentPlugin(myMainFrame, false);
                    integratePlugin(plugin, true);
                } else if (current.equalsIgnoreCase("-separate")) {
                    SeparatePlugin plugin = new SeparatePlugin(myMainFrame, false);
                    String temp = args[index].trim();
                    if (!temp.startsWith("-")) {
                        String[] chromosomes = temp.split(",");
                        plugin.setChromosomesToSeparate(chromosomes);
                        index++;
                    }
                    integratePlugin(plugin, true);
                } else if (current.equalsIgnoreCase("-synonymizer")) {
                    SynonymizerPlugin plugin = new SynonymizerPlugin(myMainFrame, false);
                    integratePlugin(plugin, true);
                } else if (current.equalsIgnoreCase("-mergeAlignments")) {
                    MergeAlignmentsPlugin plugin = new MergeAlignmentsPlugin(myMainFrame, false);
                    integratePlugin(plugin, true);
                } else if (current.equalsIgnoreCase("-mergeAlignmentsSameSites")) {
                    MergeAlignmentsSameSitesPlugin plugin = new MergeAlignmentsSameSitesPlugin(myMainFrame);
                    try {
                        for (int i = 0; i < 2; i++) {
                            String paraType = args[index++].trim();
                            String value = args[index++].trim();
                            if (paraType.equalsIgnoreCase("-input")) {
                                String[] files = value.split(",");
                                List<String> filenames = new ArrayList<String>();
                                for (int j = 0; j < files.length; j++) {
                                    filenames.add(files[j]);
                                }
                                plugin.setInputFiles(filenames);
                            } else if (paraType.equalsIgnoreCase("-output")) {
                                plugin.setOutputFile(value);
                            } else {
                                throw new IllegalArgumentException("TasselPipeline: parseArgs: -mergeAlignmentsSameSites: unknown descriptor: " + paraType);
                            }
                        }
                    } catch (IndexOutOfBoundsException e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: -mergeAlignmentsSameSites: not specified correctly.");
                    }
                    integratePlugin(plugin, true);
                } else if (current.equalsIgnoreCase("-excludeLastTrait")) {
                    FilterTraitsPlugin plugin = new FilterTraitsPlugin(myMainFrame, false);
                    ArrayList input = new ArrayList();
                    int[] excludeLast = new int[]{-1};
                    input.add(excludeLast);
                    plugin.setIncludeList(input);
                    integratePlugin(plugin, true);
                } else if (current.equalsIgnoreCase("-mlm")) {
                    MLMPlugin plugin = new MLMPlugin(myMainFrame, false);
                    integratePlugin(plugin, true);
                } else if (current.equalsIgnoreCase("-mlmVarCompEst")) {
                    MLMPlugin plugin = (MLMPlugin) findLastPluginFromCurrentPipe(new Class[]{MLMPlugin.class});
                    if (plugin == null) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No MLM step defined: " + current);
                    }
                    String method = args[index++].trim();
                    plugin.setVarCompEst(method);
                } else if (current.equalsIgnoreCase("-mlmCompressionLevel")) {
                    MLMPlugin plugin = (MLMPlugin) findLastPluginFromCurrentPipe(new Class[]{MLMPlugin.class});
                    if (plugin == null) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No MLM step defined: " + current);
                    }
                    String type = args[index++].trim();
                    if (type.equalsIgnoreCase("Optimum")) {
                        plugin.setCompressionType(MLMPlugin.CompressionType.Optimum);
                    } else if (type.equalsIgnoreCase("Custom")) {
                        plugin.setCompressionType(MLMPlugin.CompressionType.Custom);
                    } else if (type.equalsIgnoreCase("None")) {
                        plugin.setCompressionType(MLMPlugin.CompressionType.None);
                    } else {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: Unknown compression type: " + type);
                    }
                } else if (current.equalsIgnoreCase("-mlmCustomCompression")) {
                    MLMPlugin plugin = (MLMPlugin) findLastPluginFromCurrentPipe(new Class[]{MLMPlugin.class});
                    if (plugin == null) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No MLM step defined: " + current);
                    }
                    String temp = args[index++].trim();
                    double value = 0;
                    try {
                        value = Double.parseDouble(temp);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: Problem parsing custom compression: " + temp);
                    }
                    plugin.setCustomCompression(value);
                } else if (current.equalsIgnoreCase("-mlmOutputFile")) {
                    MLMPlugin plugin = (MLMPlugin) findLastPluginFromCurrentPipe(new Class[]{MLMPlugin.class});
                    if (plugin == null) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No MLM step defined: " + current);
                    }
                    String filename = args[index++].trim();
                    plugin.setOutputName(filename);
                } else if (current.equalsIgnoreCase("-mlmMaxP")) {
                    MLMPlugin plugin = (MLMPlugin) findLastPluginFromCurrentPipe(new Class[]{MLMPlugin.class});
                    if (plugin == null) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No MLM step defined: " + current);
                    }
                    String temp = args[index++].trim();
                    double maxP = 0;
                    try {
                        maxP = Double.parseDouble(temp);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: Problem parsing max P: " + temp);
                    }
                    plugin.setMaxp(maxP);
                } else if (current.equalsIgnoreCase("-glm")) {
                    FixedEffectLMPlugin plugin = new FixedEffectLMPlugin(myMainFrame, false);
                    integratePlugin(plugin, true);
                } else if (current.equalsIgnoreCase("-glmOutputFile")) {
                    FixedEffectLMPlugin plugin = (FixedEffectLMPlugin) findLastPluginFromCurrentPipe(new Class[]{FixedEffectLMPlugin.class});
                    if (plugin == null) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No GLM step defined: " + current);
                    }
                    String filename = args[index++].trim();
                    plugin.setOutputFile(filename);
                } else if (current.equalsIgnoreCase("-glmMaxP")) {
                    FixedEffectLMPlugin plugin = (FixedEffectLMPlugin) findLastPluginFromCurrentPipe(new Class[]{FixedEffectLMPlugin.class});
                    if (plugin == null) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No GLM step defined: " + current);
                    }
                    String temp = args[index++].trim();
                    double maxP = 0;
                    try {
                        maxP = Double.parseDouble(temp);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: Problem parsing max P: " + temp);
                    }
                    plugin.setMaxP(maxP);
                } else if (current.equalsIgnoreCase("-glmPermutations")) {
                    FixedEffectLMPlugin plugin = (FixedEffectLMPlugin) findLastPluginFromCurrentPipe(new Class[]{FixedEffectLMPlugin.class});
                    if (plugin == null) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No GLM step defined: " + current);
                    }
                    String temp = args[index++].trim();
                    int permutations = 0;
                    try {
                        permutations = Integer.parseInt(temp);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: Problem parsing number of permutations: " + temp);
                    }
                    plugin.setPermute(true);
                    plugin.setNumberOfPermutations(permutations);
                } else if (current.equalsIgnoreCase("-td_csv")) {
                    String csvFile = args[index++].trim();
                    getTableDisplayPlugin(csvFile, current);
                } else if (current.equalsIgnoreCase("-td_tab")) {
                    String tabFile = args[index++].trim();
                    getTableDisplayPlugin(tabFile, current);
                } else if (current.equalsIgnoreCase("-td_gui")) {
                    getTableDisplayPlugin(null, current);
                } else if (current.equalsIgnoreCase("-diversity")) {
                    SequenceDiversityPlugin plugin = new SequenceDiversityPlugin(myMainFrame, false);
                    integratePlugin(plugin, true);
                } else if (current.equalsIgnoreCase("-diversityStartBase")) {

                    SequenceDiversityPlugin plugin = null;
                    try {
                        plugin = (SequenceDiversityPlugin) myCurrentPipe.get(myCurrentPipe.size() - 1);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No SequenceDiversityPlugin step defined: " + current);
                    }

                    String str = args[index++].trim();
                    int start = -1;
                    try {
                        start = Integer.parseInt(str);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: Problem with Diversity Start Base number: " + str);
                    }
                    if (start < 0) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: Diversity Start Base can't be less than 0.");
                    }

                    plugin.setStartSite(start);

                } else if (current.equalsIgnoreCase("-diversityEndBase")) {

                    SequenceDiversityPlugin plugin = null;
                    try {
                        plugin = (SequenceDiversityPlugin) myCurrentPipe.get(myCurrentPipe.size() - 1);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No SequenceDiversityPlugin step defined: " + current);
                    }

                    String str = args[index++].trim();
                    int end = -1;
                    try {
                        end = Integer.parseInt(str);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: Problem with Diversity Start Base number: " + str);
                    }

                    plugin.setEndSite(end);

                } else if (current.equalsIgnoreCase("-diversitySlidingWin")) {
                    SequenceDiversityPlugin plugin = null;
                    plugin.setSlidingWindowAnalysis(true);
                } else if (current.equalsIgnoreCase("-diversitySlidingWinStep")) {

                    SequenceDiversityPlugin plugin = null;
                    try {
                        plugin = (SequenceDiversityPlugin) myCurrentPipe.get(myCurrentPipe.size() - 1);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No SequenceDiversityPlugin step defined: " + current);
                    }

                    String str = args[index++].trim();
                    int step = -1;
                    try {
                        step = Integer.parseInt(str);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: Problem with Diversity Sliding Win Step number: " + str);
                    }

                    plugin.setStepSize(step);
                    plugin.setSlidingWindowAnalysis(true);

                } else if (current.equalsIgnoreCase("-diversitySlidingWinSize")) {

                    SequenceDiversityPlugin plugin = null;
                    try {
                        plugin = (SequenceDiversityPlugin) myCurrentPipe.get(myCurrentPipe.size() - 1);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No SequenceDiversityPlugin step defined: " + current);
                    }

                    String str = args[index++].trim();
                    int size = -1;
                    try {
                        size = Integer.parseInt(str);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: Problem with Diversity Sliding Win Size number: " + str);
                    }

                    plugin.setWindowSize(size);
                    plugin.setSlidingWindowAnalysis(true);

                } else if (current.equalsIgnoreCase("-diversityTypeSites")) {

                    SequenceDiversityPlugin plugin = null;
                    try {
                        plugin = (SequenceDiversityPlugin) myCurrentPipe.get(myCurrentPipe.size() - 1);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No SequenceDiversityPlugin step defined: " + current);
                    }

                    Vector grp = new Vector();
                    String[] types = args[index++].trim().split(",");
                    for (int i = 0; i < types.length; i++) {
                        if (types[i].equalsIgnoreCase("ALL")) {
                            grp.add(new Integer(Alignment.POSITION_TYPE_ALL_GROUP));
                        }
                        // else if (types[i].equalsIgnoreCase("INDEL")) {
                        //     grp.add(new Integer(Alignment.POSITION_TYPE_INDEL_GROUP));
                        // }
                    }

                    plugin.setTypeOfSitesToAnalyze(grp);

                } else if (current.equalsIgnoreCase("-ld")) {
                    LinkageDisequilibriumPlugin plugin = new LinkageDisequilibriumPlugin(myMainFrame, false);
                    integratePlugin(plugin, true);
                } else if (current.equalsIgnoreCase("-ldPermNum")) {

                    LinkageDisequilibriumPlugin plugin = null;
                    try {
                        plugin = (LinkageDisequilibriumPlugin) myCurrentPipe.get(myCurrentPipe.size() - 1);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No LinkageDisequilibriumPlugin step defined: " + current);
                    }

                    String str = args[index++].trim();
                    int permNum = -1;
                    try {
                        permNum = Integer.parseInt(str);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: Problem with LD Permutation number: " + str);
                    }
                    if (permNum < 1) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: LD Permutation size can't be less than 1.");
                    }

                    plugin.setPermutationNumber(permNum);

                } else if (current.equalsIgnoreCase("-ldTestSite")) {

                    LinkageDisequilibriumPlugin plugin = null;
                    try {
                        plugin = (LinkageDisequilibriumPlugin) myCurrentPipe.get(myCurrentPipe.size() - 1);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No LinkageDisequilibriumPlugin step defined: " + current);
                    }

                    String str = args[index++].trim();
                    int testSite = -1;
                    try {
                        testSite = Integer.parseInt(str);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: Problem with LD Test Site number: " + str);
                    }
                    if (testSite < 0) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: LD Test Site can't be less than 0.");
                    }

                    plugin.setTestSite(testSite);

                } else if (current.equalsIgnoreCase("-ldWinSize")) {

                    LinkageDisequilibriumPlugin plugin = null;
                    try {
                        plugin = (LinkageDisequilibriumPlugin) myCurrentPipe.get(myCurrentPipe.size() - 1);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No LinkageDisequilibriumPlugin step defined: " + current);
                    }

                    String str = args[index++].trim();
                    int winSize = -1;
                    try {
                        winSize = Integer.parseInt(str);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: Problem with LD Window Size: " + str);
                    }
                    if (winSize < 1) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: LD Window Size can't be less than 1.");
                    }

                    plugin.setWinSize(winSize);

                } else if (current.equalsIgnoreCase("-ldRapidAnalysis")) {

                    LinkageDisequilibriumPlugin plugin = null;
                    try {
                        plugin = (LinkageDisequilibriumPlugin) myCurrentPipe.get(myCurrentPipe.size() - 1);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No LinkageDisequilibriumPlugin step defined: " + current);
                    }

                    String temp = args[index++].trim();
                    boolean rapid = true;
                    if (temp.equalsIgnoreCase("false")) {
                        rapid = false;
                    } else if (temp.equalsIgnoreCase("true")) {
                        rapid = true;
                    } else {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: LD Rapid Analysis parameter must be true or false.");
                    }

                    plugin.setRapidAnalysis(rapid);

                } else if (current.equalsIgnoreCase("-ldType")) {

                    LinkageDisequilibriumPlugin plugin = null;
                    try {
                        plugin = (LinkageDisequilibriumPlugin) myCurrentPipe.get(myCurrentPipe.size() - 1);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No LinkageDisequilibriumPlugin step defined: " + current);
                    }

                    String temp = args[index++].trim();
                    if (temp.equalsIgnoreCase("All")) {
                        plugin.setLDType(testDesign.All);
                    } else if (temp.equalsIgnoreCase("SlidingWindow")) {
                        plugin.setLDType(testDesign.SlidingWindow);
                    } else if (temp.equalsIgnoreCase("SiteByAll")) {
                        plugin.setLDType(testDesign.SiteByAll);
                    } else {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: LD Type parameter must be All, SlidingWindow, or SiteByAll.");
                    }

                } else if (current.equalsIgnoreCase("-ldHetTreatment")) {

                    LinkageDisequilibriumPlugin plugin = null;
                    try {
                        plugin = (LinkageDisequilibriumPlugin) myCurrentPipe.get(myCurrentPipe.size() - 1);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No LinkageDisequilibriumPlugin step defined: " + current);
                    }

                    String temp = args[index++].trim();
                    if (temp.equalsIgnoreCase("Haplotype")) {
                        plugin.setHetTreatment(HetTreatment.Haplotype);
                    } else if (temp.equalsIgnoreCase("Homozygous")) {
                        plugin.setHetTreatment(HetTreatment.Homozygous);
                    } else if (temp.equalsIgnoreCase("Genotype")) {
                        plugin.setHetTreatment(HetTreatment.Genotype);
                    } else {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: LD Het Treatment parameter must be Haplotype, Homozygous, or Genotype.");
                    }

                } else if (current.equalsIgnoreCase("-ldd")) {
                    String outputType = args[index++].trim();
                    getLinkageDiseqDisplayPlugin(outputType);
                } else if (current.equalsIgnoreCase("-ldplotsize")) {

                    LinkageDiseqDisplayPlugin plugin = null;
                    try {
                        plugin = (LinkageDiseqDisplayPlugin) myCurrentPipe.get(myCurrentPipe.size() - 1);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No LinkageDiseqDisplay step defined: " + current);
                    }

                    String str = args[index++].trim();
                    int plotSize = -1;
                    try {
                        plotSize = Integer.parseInt(str);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: Problem with LD Plot size number: " + str);
                    }
                    if (plotSize < 1) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: LD Plot size can't be less than 1.");
                    }

                    plugin.setImageSize(plotSize, plotSize);

                } else if (current.equalsIgnoreCase("-ldplotlabels")) {

                    LinkageDiseqDisplayPlugin plugin = null;
                    try {
                        plugin = (LinkageDiseqDisplayPlugin) myCurrentPipe.get(myCurrentPipe.size() - 1);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No LinkageDiseqDisplay step defined: " + current);
                    }

                    String temp = args[index++].trim();
                    boolean ldPlotLabels = true;
                    if (temp.equalsIgnoreCase("false")) {
                        ldPlotLabels = false;
                    } else if (temp.equalsIgnoreCase("true")) {
                        ldPlotLabels = true;
                    } else {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: LD Plot labels parameter must be true or false.");
                    }

                    plugin.setShowLabels(ldPlotLabels);

                } else if (current.equalsIgnoreCase("-o")) {

                    Plugin plugin = findLastPluginFromCurrentPipe(new Class[]{LinkageDiseqDisplayPlugin.class});

                    String temp = args[index++].trim();

                    if (plugin == null) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No LinkageDiseqDisplay step defined: " + current + " " + temp);
                    } else if (plugin instanceof LinkageDiseqDisplayPlugin) {
                        ((LinkageDiseqDisplayPlugin) plugin).setSaveFile(temp);
                    }

                } else if (current.equalsIgnoreCase("-ck")) {
                    KinshipPlugin plugin = new KinshipPlugin(myMainFrame, false);
                    integratePlugin(plugin, true);
                } else if (current.equalsIgnoreCase("-ckModelHets")) {
                    throw new IllegalArgumentException("TasselPipeline: parseArgs: -ckModelHets not needed in Tassel 4.0. It is designed to handle heterzygotes.");
                } else if (current.equalsIgnoreCase("-ckRescale")) {
                    throw new IllegalArgumentException("TasselPipeline: parseArgs: -ckRescale not needed in Tassel 4.0. It is designed to handle heterzygotes.");

                } else if (current.equalsIgnoreCase("-tree")) {

                    CreateTreePlugin plugin = new CreateTreePlugin(myMainFrame, false);
                    integratePlugin(plugin, true);

                    String temp = args[index++].trim();
                    if (temp.equalsIgnoreCase("Neighbor")) {
                        plugin.setNeighborJoining(true);
                    } else if (temp.equalsIgnoreCase("UPGMA")) {
                        plugin.setNeighborJoining(false);
                    } else {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: tree clustering method must be Neighbor or UPGMA: " + temp);
                    }

                } else if (current.equalsIgnoreCase("-treeSaveDistance")) {

                    CreateTreePlugin plugin = (CreateTreePlugin) findLastPluginFromCurrentPipe(new Class[]{CreateTreePlugin.class});
                    if (plugin == null) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No Create Tree step defined: " + current);
                    }

                    String temp = args[index++].trim();
                    boolean value = true;
                    if (temp.equalsIgnoreCase("false")) {
                        value = false;
                    } else if (temp.equalsIgnoreCase("true")) {
                        value = true;
                    } else {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: tree save distance matrix parameter must be true or false: " + value);
                    }

                    plugin.setReturnDistanceMatrix(value);

                } else if (current.equalsIgnoreCase("-gs")) {
                    RidgeRegressionEmmaPlugin plugin = new RidgeRegressionEmmaPlugin(myMainFrame, false);
                    integratePlugin(plugin, true);
                } else if (current.equalsIgnoreCase("-distanceMatrix")) {
                    DistanceMatrixPlugin plugin = new DistanceMatrixPlugin(myMainFrame, false);
                    integratePlugin(plugin, true);
                } else if (current.equalsIgnoreCase("-distMatrixRanges")) {
                    DistanceMatrixRangesPlugin plugin = new DistanceMatrixRangesPlugin(myMainFrame, false);
                    integratePlugin(plugin, true);
                } else if (current.equalsIgnoreCase("-distMatrixRangesLocus")) {
                    DistanceMatrixRangesPlugin plugin = null;
                    try {
                        plugin = (DistanceMatrixRangesPlugin) myCurrentPipe.get(myCurrentPipe.size() - 1);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No DistanceMatrixRangesPlugin step defined: " + current);
                    }

                    String str = args[index++].trim();
                    plugin.setLocus(str);
                } else if (current.equalsIgnoreCase("-distMatrixRangesTaxon")) {
                    DistanceMatrixRangesPlugin plugin = null;
                    try {
                        plugin = (DistanceMatrixRangesPlugin) myCurrentPipe.get(myCurrentPipe.size() - 1);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No DistanceMatrixRangesPlugin step defined: " + current);
                    }

                    String str = args[index++].trim();
                    plugin.setTaxon(str);
                } else if (current.equalsIgnoreCase("-distMatrixRangesPos")) {
                    DistanceMatrixRangesPlugin plugin = null;
                    try {
                        plugin = (DistanceMatrixRangesPlugin) myCurrentPipe.get(myCurrentPipe.size() - 1);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No DistanceMatrixRangesPlugin step defined: " + current);
                    }

                    String[] positions = args[index++].trim().split(",");
                    plugin.setPhysicalPositions(positions);
                } else if (current.equalsIgnoreCase("-distMatrixRangesPosFile")) {
                    DistanceMatrixRangesPlugin plugin = null;
                    try {
                        plugin = (DistanceMatrixRangesPlugin) myCurrentPipe.get(myCurrentPipe.size() - 1);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No DistanceMatrixRangesPlugin step defined: " + current);
                    }
                    String posFile = args[index++].trim();

                    List positions = new ArrayList();
                    BufferedReader br = null;
                    try {
                        br = Utils.getBufferedReader(posFile);
                        String inputline = br.readLine();
                        Pattern sep = Pattern.compile("\\s+");

                        while (inputline != null) {
                            inputline = inputline.trim();
                            String[] parsedline = sep.split(inputline);
                            for (int i = 0; i < parsedline.length; i++) {
                                if ((parsedline[i] != null) || (parsedline[i].length() != 0)) {
                                    positions.add(parsedline[i]);
                                }
                            }
                            inputline = br.readLine();
                        }
                    } finally {
                        br.close();
                    }

                    String[] positionArray = new String[positions.size()];
                    positionArray = (String[]) positions.toArray(positionArray);
                    plugin.setPhysicalPositions(positionArray);
                } else if (current.equalsIgnoreCase("-genotypeSummary")) {
                    GenotypeSummaryPlugin plugin = new GenotypeSummaryPlugin(myMainFrame, false);
                    String temp = args[index++].trim();
                    String[] types = temp.split(",");
                    plugin.setCaculateOverview(false);
                    plugin.setCalculateSiteSummary(false);
                    plugin.setCalculateTaxaSummary(false);
                    for (int i = 0; i < types.length; i++) {
                        if (types[i].equalsIgnoreCase("overall")) {
                            plugin.setCaculateOverview(true);
                        } else if (types[i].equalsIgnoreCase("site")) {
                            plugin.setCalculateSiteSummary(true);
                        } else if (types[i].equalsIgnoreCase("taxa")) {
                            plugin.setCalculateTaxaSummary(true);
                        } else if (types[i].equalsIgnoreCase("all")) {
                            plugin.setCaculateOverview(true);
                            plugin.setCalculateSiteSummary(true);
                            plugin.setCalculateTaxaSummary(true);
                        } else {
                            throw new IllegalArgumentException("TasselPipeline: parseArgs: -genotypeSummary illegal types: " + temp);
                        }
                    }
                    integratePlugin(plugin, true);
                } else if (current.equalsIgnoreCase("-export")) {
                    ExportMultiplePlugin plugin = new ExportMultiplePlugin(myMainFrame);
                    String temp = args[index].trim();
                    if (!temp.startsWith("-")) {
                        String[] filenames = temp.split(",");
                        plugin.setSaveFiles(filenames);
                        index++;
                    }
                    integratePlugin(plugin, false);
                } else if (current.equalsIgnoreCase("-exportType")) {

                    ExportMultiplePlugin plugin = (ExportMultiplePlugin) findLastPluginFromCurrentPipe(new Class[]{ExportMultiplePlugin.class});
                    if (plugin == null) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No Export step defined: " + current);
                    }

                    String type = args[index++].trim();
                    FileLoadPlugin.TasselFileType fileType = null;
                    try {
                        fileType = FileLoadPlugin.TasselFileType.valueOf(type);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: -exportType: Unknown type: " + type + "  Should be: " + Arrays.toString(FileLoadPlugin.TasselFileType.values()));
                    }
                    plugin.setAlignmentFileType(fileType);

                } else if (current.equalsIgnoreCase("-impute")) {
                    GenotypeImputationPlugin plugin = new GenotypeImputationPlugin(myMainFrame, false);
                    integratePlugin(plugin, true);
                } else if (current.equalsIgnoreCase("-imputeMethod")) {
                    GenotypeImputationPlugin plugin = (GenotypeImputationPlugin) findLastPluginFromCurrentPipe(new Class[]{GenotypeImputationPlugin.class});
                    if (plugin == null) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No Impute step defined: " + current);
                    }
                    String temp = args[index++].trim();
                    if (temp.equalsIgnoreCase(GenotypeImputationPlugin.ImpMethod.Length.toString())) {
                        plugin.setMethod(GenotypeImputationPlugin.ImpMethod.Length);
                    } else if (temp.equalsIgnoreCase(GenotypeImputationPlugin.ImpMethod.MajorAllele.toString())) {
                        plugin.setMethod(GenotypeImputationPlugin.ImpMethod.MajorAllele);
                    } else if (temp.equalsIgnoreCase(GenotypeImputationPlugin.ImpMethod.IBDProb.toString())) {
                        plugin.setMethod(GenotypeImputationPlugin.ImpMethod.IBDProb);
                    } else if (temp.equalsIgnoreCase(GenotypeImputationPlugin.ImpMethod.SimilarWindow.toString())) {
                        plugin.setMethod(GenotypeImputationPlugin.ImpMethod.SimilarWindow);
                    } else {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: Not defined impute method: " + temp);
                    }
                } else if (current.equalsIgnoreCase("-imputeMinLength")) {
                    GenotypeImputationPlugin plugin = (GenotypeImputationPlugin) findLastPluginFromCurrentPipe(new Class[]{GenotypeImputationPlugin.class});
                    if (plugin == null) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No Impute step defined: " + current);
                    }
                    String temp = args[index++].trim();
                    int minLength = 0;
                    try {
                        minLength = Integer.parseInt(temp);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: Problem parsing impute min length: " + temp);
                    }
                    plugin.setMinLength(minLength);
                } else if (current.equalsIgnoreCase("-imputeMaxMismatch")) {
                    GenotypeImputationPlugin plugin = (GenotypeImputationPlugin) findLastPluginFromCurrentPipe(new Class[]{GenotypeImputationPlugin.class});
                    if (plugin == null) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No Impute step defined: " + current);
                    }
                    String temp = args[index++].trim();
                    int maxMismatch = 0;
                    try {
                        maxMismatch = Integer.parseInt(temp);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: Problem parsing impute max mismatch: " + temp);
                    }
                    plugin.setMaxMisMatch(maxMismatch);
                } else if (current.equalsIgnoreCase("-imputeMinProb")) {
                    GenotypeImputationPlugin plugin = (GenotypeImputationPlugin) findLastPluginFromCurrentPipe(new Class[]{GenotypeImputationPlugin.class});
                    if (plugin == null) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No Impute step defined: " + current);
                    }
                    String temp = args[index++].trim();
                    double minProb = 0;
                    try {
                        minProb = Double.parseDouble(temp);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: Problem parsing impute min probability: " + temp);
                    }
                    plugin.setMinProb(minProb);
                } else if (current.equalsIgnoreCase("-filterAlign")) {
                    FilterAlignmentPlugin plugin = new FilterAlignmentPlugin(myMainFrame, false);
                    integratePlugin(plugin, true);
                } else if (current.equalsIgnoreCase("-filterAlignMinCount")) {
                    FilterAlignmentPlugin plugin = (FilterAlignmentPlugin) findLastPluginFromCurrentPipe(new Class[]{FilterAlignmentPlugin.class});
                    if (plugin == null) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No Filter Alignment step defined: " + current);
                    }
                    String temp = args[index++].trim();
                    int minCount = 0;
                    try {
                        minCount = Integer.parseInt(temp);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: Problem parsing filter alignment min count: " + temp);
                    }
                    plugin.setMinCount(minCount);
                } else if (current.equalsIgnoreCase("-filterAlignMinFreq")) {
                    FilterAlignmentPlugin plugin = (FilterAlignmentPlugin) findLastPluginFromCurrentPipe(new Class[]{FilterAlignmentPlugin.class});
                    if (plugin == null) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No Filter Alignment step defined: " + current);
                    }
                    String temp = args[index++].trim();
                    double minFreq = 0;
                    try {
                        minFreq = Double.parseDouble(temp);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: Problem parsing filter alignment min frequency: " + temp);
                    }
                    plugin.setMinFreq(minFreq);
                } else if (current.equalsIgnoreCase("-filterAlignMaxFreq")) {
                    FilterAlignmentPlugin plugin = (FilterAlignmentPlugin) findLastPluginFromCurrentPipe(new Class[]{FilterAlignmentPlugin.class});
                    if (plugin == null) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No Filter Alignment step defined: " + current);
                    }
                    String temp = args[index++].trim();
                    double maxFreq = 0;
                    try {
                        maxFreq = Double.parseDouble(temp);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: Problem parsing filter alignment max frequency: " + temp);
                    }
                    plugin.setMaxFreq(maxFreq);
                } else if (current.equalsIgnoreCase("-filterAlignStart")) {
                    FilterAlignmentPlugin plugin = (FilterAlignmentPlugin) findLastPluginFromCurrentPipe(new Class[]{FilterAlignmentPlugin.class});
                    if (plugin == null) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No Filter Alignment step defined: " + current);
                    }
                    String temp = args[index++].trim();
                    int start = 0;
                    try {
                        start = Integer.parseInt(temp);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: Problem parsing filter alignment start: " + temp);
                    }
                    plugin.setStart(start);
                } else if (current.equalsIgnoreCase("-filterAlignEnd")) {
                    FilterAlignmentPlugin plugin = (FilterAlignmentPlugin) findLastPluginFromCurrentPipe(new Class[]{FilterAlignmentPlugin.class});
                    if (plugin == null) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No Filter Alignment step defined: " + current);
                    }
                    String temp = args[index++].trim();
                    int end = 0;
                    try {
                        end = Integer.parseInt(temp);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: Problem parsing filter alignment end: " + temp);
                    }
                    plugin.setEnd(end);
                } else if (current.equalsIgnoreCase("-filterAlignStartPos")) {
                    FilterAlignmentPlugin plugin = (FilterAlignmentPlugin) findLastPluginFromCurrentPipe(new Class[]{FilterAlignmentPlugin.class});
                    if (plugin == null) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No Filter Alignment step defined: " + current);
                    }
                    String temp = args[index++].trim();
                    int startPos = 0;
                    try {
                        startPos = Integer.parseInt(temp);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: Problem parsing filter alignment start physical position: " + temp);
                    }
                    plugin.setStartPos(startPos);
                } else if (current.equalsIgnoreCase("-filterAlignEndPos")) {
                    FilterAlignmentPlugin plugin = (FilterAlignmentPlugin) findLastPluginFromCurrentPipe(new Class[]{FilterAlignmentPlugin.class});
                    if (plugin == null) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No Filter Alignment step defined: " + current);
                    }
                    String temp = args[index++].trim();
                    int endPos = 0;
                    try {
                        endPos = Integer.parseInt(temp);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: Problem parsing filter alignment end physical position: " + temp);
                    }
                    plugin.setEndPos(endPos);
                } else if (current.equalsIgnoreCase("-filterAlignLocus")) {
                    FilterAlignmentPlugin plugin = (FilterAlignmentPlugin) findLastPluginFromCurrentPipe(new Class[]{FilterAlignmentPlugin.class});
                    if (plugin == null) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No Filter Alignment step defined: " + current);
                    }
                    String temp = args[index++].trim();
                    plugin.setLocusStr(temp);
                } else if (current.equalsIgnoreCase("-filterAlignExtInd")) {
                    FilterAlignmentPlugin plugin = (FilterAlignmentPlugin) findLastPluginFromCurrentPipe(new Class[]{FilterAlignmentPlugin.class});
                    if (plugin == null) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No Filter Alignment step defined: " + current);
                    }
                    plugin.setExtractIndels(true);
                } else if (current.equalsIgnoreCase("-filterAlignRemMinor")) {
                    FilterAlignmentPlugin plugin = (FilterAlignmentPlugin) findLastPluginFromCurrentPipe(new Class[]{FilterAlignmentPlugin.class});
                    if (plugin == null) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No Filter Alignment step defined: " + current);
                    }
                    plugin.setFilterMinorSNPs(true);
                } else if (current.equalsIgnoreCase("-filterAlignSliding")) {
                    FilterAlignmentPlugin plugin = (FilterAlignmentPlugin) findLastPluginFromCurrentPipe(new Class[]{FilterAlignmentPlugin.class});
                    if (plugin == null) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No Filter Alignment step defined: " + current);
                    }
                    plugin.setDoSlidingHaps(true);
                } else if (current.equalsIgnoreCase("-filterAlignHapLen")) {
                    FilterAlignmentPlugin plugin = (FilterAlignmentPlugin) findLastPluginFromCurrentPipe(new Class[]{FilterAlignmentPlugin.class});
                    if (plugin == null) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No Filter Alignment step defined: " + current);
                    }
                    String temp = args[index++].trim();
                    int hapLen = 0;
                    try {
                        hapLen = Integer.parseInt(temp);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: Problem parsing filter alignment haplotype length: " + temp);
                    }
                    plugin.setWinSize(hapLen);
                } else if (current.equalsIgnoreCase("-filterAlignStepLen")) {
                    FilterAlignmentPlugin plugin = (FilterAlignmentPlugin) findLastPluginFromCurrentPipe(new Class[]{FilterAlignmentPlugin.class});
                    if (plugin == null) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: No Filter Alignment step defined: " + current);
                    }
                    String temp = args[index++].trim();
                    int stepLen = 0;
                    try {
                        stepLen = Integer.parseInt(temp);
                    } catch (Exception e) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: Problem parsing filter alignment step length: " + temp);
                    }
                    plugin.setStepSize(stepLen);
                } else if (current.equalsIgnoreCase("-numericalGenoTransform")) {
                    NumericalGenotypePlugin plugin = new NumericalGenotypePlugin();

                    String temp = args[index++].trim();
                    if (temp.equalsIgnoreCase(NumericalGenotypePlugin.TRANSFORM_TYPE.collapse.toString())) {
                        plugin.setTransformType(NumericalGenotypePlugin.TRANSFORM_TYPE.collapse);
                    } else if (temp.equalsIgnoreCase(NumericalGenotypePlugin.TRANSFORM_TYPE.separated.toString())) {
                        plugin.setTransformType(NumericalGenotypePlugin.TRANSFORM_TYPE.separated);
                    } else {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: Not defined genotype transform type: " + temp);
                    }

                    integratePlugin(plugin, true);
                } else if (current.equalsIgnoreCase("-includeTaxa")) {
                    FilterTaxaAlignmentPlugin plugin = new FilterTaxaAlignmentPlugin(myMainFrame, false);
                    String[] taxa = args[index++].trim().split(",");
                    Identifier[] ids = new Identifier[taxa.length];
                    for (int i = 0; i < taxa.length; i++) {
                        ids[i] = new Identifier(taxa[i]);
                    }
                    plugin.setIdsToKeep(new SimpleIdGroup(ids));
                    integratePlugin(plugin, true);
                } else if (current.equalsIgnoreCase("-includeTaxaInFile")) {
                    FilterTaxaAlignmentPlugin plugin = new FilterTaxaAlignmentPlugin(myMainFrame, false);
                    String taxaListFile = args[index++].trim();

                    List taxa = new ArrayList();
                    BufferedReader br = null;
                    try {
                        br = Utils.getBufferedReader(taxaListFile);
                        String inputline = br.readLine();
                        Pattern sep = Pattern.compile("\\s+");

                        while (inputline != null) {
                            inputline = inputline.trim();
                            String[] parsedline = sep.split(inputline);
                            for (int i = 0; i < parsedline.length; i++) {
                                if ((parsedline[i] != null) || (parsedline[i].length() != 0)) {
                                    taxa.add(parsedline[i]);
                                }
                            }
                            inputline = br.readLine();
                        }
                    } finally {
                        br.close();
                    }

                    Identifier[] ids = new Identifier[taxa.size()];
                    for (int i = 0; i < taxa.size(); i++) {
                        ids[i] = new Identifier((String) taxa.get(i));
                    }
                    plugin.setIdsToKeep(new SimpleIdGroup(ids));
                    integratePlugin(plugin, true);
                } else if (current.equalsIgnoreCase("-excludeTaxa")) {
                    FilterTaxaAlignmentPlugin plugin = new FilterTaxaAlignmentPlugin(myMainFrame, false);
                    String[] taxa = args[index++].trim().split(",");
                    Identifier[] ids = new Identifier[taxa.length];
                    for (int i = 0; i < taxa.length; i++) {
                        ids[i] = new Identifier(taxa[i]);
                    }
                    plugin.setIdsToRemove(new SimpleIdGroup(ids));
                    integratePlugin(plugin, true);
                } else if (current.equalsIgnoreCase("-excludeTaxaInFile")) {
                    FilterTaxaAlignmentPlugin plugin = new FilterTaxaAlignmentPlugin(myMainFrame, false);
                    String taxaListFile = args[index++].trim();

                    List taxa = new ArrayList();
                    BufferedReader br = null;
                    try {
                        br = Utils.getBufferedReader(taxaListFile);
                        String inputline = br.readLine();
                        Pattern sep = Pattern.compile("\\s+");

                        while (inputline != null) {
                            inputline = inputline.trim();
                            String[] parsedline = sep.split(inputline);
                            for (int i = 0; i < parsedline.length; i++) {
                                if ((parsedline[i] != null) || (parsedline[i].length() != 0)) {
                                    taxa.add(parsedline[i]);
                                }
                            }
                            inputline = br.readLine();
                        }
                    } finally {
                        br.close();
                    }

                    Identifier[] ids = new Identifier[taxa.size()];
                    for (int i = 0; i < taxa.size(); i++) {
                        ids[i] = new Identifier((String) taxa.get(i));
                    }
                    plugin.setIdsToRemove(new SimpleIdGroup(ids));
                    integratePlugin(plugin, true);
                } else if (current.equalsIgnoreCase("-includeSiteNames")) {
                    FilterSiteNamePlugin plugin = new FilterSiteNamePlugin(myMainFrame, false);
                    String[] names = args[index++].trim().split(",");
                    plugin.setSiteNamesToKeep(names);
                    integratePlugin(plugin, true);
                } else if (current.equalsIgnoreCase("-includeSiteNamesInFile")) {
                    FilterSiteNamePlugin plugin = new FilterSiteNamePlugin(myMainFrame, false);
                    String siteNameListFile = args[index++].trim();

                    List siteNames = new ArrayList();
                    BufferedReader br = null;
                    try {
                        br = Utils.getBufferedReader(siteNameListFile);
                        String inputline = br.readLine();
                        Pattern sep = Pattern.compile("\\s+");

                        while (inputline != null) {
                            inputline = inputline.trim();
                            String[] parsedline = sep.split(inputline);
                            for (int i = 0; i < parsedline.length; i++) {
                                if ((parsedline[i] != null) || (parsedline[i].length() != 0)) {
                                    siteNames.add(parsedline[i]);
                                }
                            }
                            inputline = br.readLine();
                        }
                    } finally {
                        br.close();
                    }

                    String[] siteNameArray = new String[siteNames.size()];
                    siteNameArray = (String[]) siteNames.toArray(siteNameArray);
                    plugin.setSiteNamesToKeep(siteNameArray);
                    integratePlugin(plugin, true);
                } else if (current.equalsIgnoreCase("-excludeSiteNames")) {
                    FilterSiteNamePlugin plugin = new FilterSiteNamePlugin(myMainFrame, false);
                    String[] sites = args[index++].trim().split(",");
                    plugin.setSiteNamesToRemove(sites);
                    integratePlugin(plugin, true);
                } else if (current.equalsIgnoreCase("-excludeSiteNamesInFile")) {
                    FilterSiteNamePlugin plugin = new FilterSiteNamePlugin(myMainFrame, false);
                    String siteNameListFile = args[index++].trim();

                    List siteNames = new ArrayList();
                    BufferedReader br = null;
                    try {
                        br = Utils.getBufferedReader(siteNameListFile);
                        String inputline = br.readLine();
                        Pattern sep = Pattern.compile("\\s+");

                        while (inputline != null) {
                            inputline = inputline.trim();
                            String[] parsedline = sep.split(inputline);
                            for (int i = 0; i < parsedline.length; i++) {
                                if ((parsedline[i] != null) || (parsedline[i].length() != 0)) {
                                    siteNames.add(parsedline[i]);
                                }
                            }
                            inputline = br.readLine();
                        }
                    } finally {
                        br.close();
                    }

                    String[] names = new String[siteNames.size()];
                    for (int i = 0; i < siteNames.size(); i++) {
                        names[i] = (String) siteNames.get(i);
                    }
                    plugin.setSiteNamesToRemove(names);
                    integratePlugin(plugin, true);
                } else if (current.equalsIgnoreCase("-newCoordinates")) {
                    ConvertAlignmentCoordinatesPlugin plugin = new ConvertAlignmentCoordinatesPlugin(myMainFrame);
                    String mapFile = args[index++].trim();
                    plugin.setMapFilename(mapFile);
                    integratePlugin(plugin, true);
                } else {

                    try {
                        Plugin plugin = null;
                        String possibleClassName = current.substring(1);
                        List<String> matches = Utils.getFullyQualifiedClassNames(possibleClassName);
                        for (String match : matches) {
                            try {
                                Class currentMatch = Class.forName(match);
                                Constructor constructor = currentMatch.getConstructor(Frame.class);
                                plugin = (Plugin) constructor.newInstance(myMainFrame);
                                break;
                            } catch (NoSuchMethodException nsme) {
                                myLogger.warn("Self-describing Plugins should implement this constructor: " + current);
                                myLogger.warn("public Plugin(Frame parentFrame) {");
                                myLogger.warn("   super(parentFrame, false);");
                                myLogger.warn("}");
                            } catch (Exception e) {
                                // do nothing
                            }
                        }

                        if (plugin == null) {
                            try {
                                Class possibleClass = Class.forName(possibleClassName);
                                Constructor constructor = possibleClass.getConstructor(Frame.class);
                                plugin = (Plugin) constructor.newInstance(myMainFrame);
                            } catch (NoSuchMethodException nsme) {
                                myLogger.warn("Self-describing Plugins should implement this constructor: " + current);
                                myLogger.warn("public Plugin(Frame parentFrame) {");
                                myLogger.warn("   super(parentFrame, false);");
                                myLogger.warn("}");
                            } catch (Exception e) {
                                // do nothing
                            }
                        }

                        if (plugin != null) {
                            List pluginArgs = new ArrayList();
                            if (index == args.length) {
                                throw new IllegalArgumentException("TasselPipeline: parseArgs: No -endPlugin flag specified.");
                            }
                            String temp = args[index++].trim();
                            while (!temp.equalsIgnoreCase("-endPlugin")) {
                                pluginArgs.add(temp);
                                if (index == args.length) {
                                    throw new IllegalArgumentException("TasselPipeline: parseArgs: No -endPlugin flag specified.");
                                }
                                temp = args[index++].trim();
                            }

                            String[] result = new String[pluginArgs.size()];
                            result = (String[]) pluginArgs.toArray(result);

                            plugin.setParameters(result);
                            integratePlugin(plugin, true);
                        } else {
                            throw new IllegalArgumentException("TasselPipeline: parseArgs: Unknown parameter: " + current);
                        }

                    } catch (UnsupportedOperationException usoe) {
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: this plugin is not self-described: " + current);
                    } catch (Exception e) {
                        ExceptionUtils.logExceptionCauses(e, myLogger, Level.ERROR);
                        throw new IllegalArgumentException("TasselPipeline: parseArgs: Unknown parameter: " + current);
                    }

                }
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }

        }

        if (myFirstPlugin != null) {
            //((AbstractPlugin) myFirstPlugin).trace(0);
            tracePipeline();
        } else {
            myLogger.warn("parseArgs: no arguments specified.");
        }

    }

    private void tracePipeline() {

        for (int i = 0; i < myThreads.size(); i++) {
            Plugin current = (Plugin) ((ThreadedPluginListener) myThreads.get(i)).getPluginListener();
            ((AbstractPlugin) current).trace(0);
        }

    }

    public FileLoadPlugin loadFile(String filename, FileLoadPlugin.TasselFileType fileType) {

        myLogger.info("loadFile: " + filename);

        FileLoadPlugin plugin = new FileLoadPlugin(myMainFrame, false);
        if (fileType == null) {
            plugin.setTheFileType(FileLoadPlugin.TasselFileType.Unknown);
        } else {
            plugin.setTheFileType(fileType);
        }

        plugin.setOpenFiles(new String[]{filename});

        integratePlugin(plugin, true);

        return plugin;

    }

    public TableDisplayPlugin getTableDisplayPlugin(String filename, String flag) {

        TableDisplayPlugin plugin = null;

        if (flag.equalsIgnoreCase("-td_gui")) {
            plugin = new TableDisplayPlugin(myMainFrame, true);
            integratePlugin(plugin, false);
        } else if (flag.equalsIgnoreCase("-td_tab")) {
            filename = Utils.addSuffixIfNeeded(filename, ".txt");
            myLogger.info("getTableDisplayPlugin: " + filename);
            plugin = new TableDisplayPlugin(myMainFrame, false);
            plugin.setDelimiter("\t");
            plugin.setSaveFile(filename);
            integratePlugin(plugin, false);
        } else if (flag.equalsIgnoreCase("-td_csv")) {
            filename = Utils.addSuffixIfNeeded(filename, ".csv");
            myLogger.info("getTableDisplayPlugin: " + filename);
            plugin = new TableDisplayPlugin(myMainFrame, false);
            plugin.setDelimiter(",");
            plugin.setSaveFile(filename);
            integratePlugin(plugin, false);
        }

        return plugin;

    }

    public LinkageDiseqDisplayPlugin getLinkageDiseqDisplayPlugin(String type) {

        LinkageDiseqDisplayPlugin plugin = new LinkageDiseqDisplayPlugin(myMainFrame, true);

        if (type.equalsIgnoreCase("gui")) {
            plugin = new LinkageDiseqDisplayPlugin(null, true);
            plugin.setBlockSchematic(false);
            plugin.setLowerCorner(LinkageDisequilibriumComponent.P_VALUE);
            plugin.setUpperCorner(LinkageDisequilibriumComponent.RSQUARE);
        } else {
            plugin = new LinkageDiseqDisplayPlugin(null, false);
            plugin.setBlockSchematic(false);
            plugin.setLowerCorner(LinkageDisequilibriumComponent.P_VALUE);
            plugin.setUpperCorner(LinkageDisequilibriumComponent.RSQUARE);

            if (type.equalsIgnoreCase("png")) {
                plugin.setOutformat(AbstractDisplayPlugin.Outformat.png);
            } else if (type.equalsIgnoreCase("gif")) {
                plugin.setOutformat(AbstractDisplayPlugin.Outformat.gif);
            } else if (type.equalsIgnoreCase("bmp")) {
                plugin.setOutformat(AbstractDisplayPlugin.Outformat.bmp);
            } else if (type.equalsIgnoreCase("jpg")) {
                plugin.setOutformat(AbstractDisplayPlugin.Outformat.jpg);
            } else if (type.equalsIgnoreCase("svg")) {
                plugin.setOutformat(AbstractDisplayPlugin.Outformat.svg);
            } else {
                throw new IllegalArgumentException("TasselPipeline: getLinkageDiseqDisplayPlugin: unknown output type: " + type);
            }
        }

        integratePlugin(plugin, false);

        return plugin;

    }

    private void integratePlugin(Plugin plugin, boolean displayDataTree) {

        if (myFirstPlugin == null) {
            myFirstPlugin = plugin;
        }

        if (displayDataTree) {
            plugin.addListener(this);
        }

        if (myCurrentPipe == null) {
            myCurrentPipe = new ArrayList();

        }

        if (myCurrentPipe.size() == 0) {
            myCurrentPipe.add(plugin);
        } else {
            plugin.receiveInput(((Plugin) myCurrentPipe.get(myCurrentPipe.size() - 1)));
            myCurrentPipe.add(plugin);
        }

    }

    private Plugin findLastPluginFromAll(Class[] types) {

        if ((myCurrentPipe != null) && (myCurrentPipe.size() != 0)) {
            for (int i = myCurrentPipe.size() - 1; i >= 0; i--) {
                Plugin current = (Plugin) myCurrentPipe.get(i);
                if (matchType(types, current)) {
                    return current;
                }

            }
        }

        List keys = new ArrayList(myForks.keySet());
        for (int i = keys.size() - 1; i >= 0; i--) {
            List currentPipe = (List) myForks.get(keys.get(i));
            for (int j = currentPipe.size() - 1; j >= 0; j--) {
                Plugin current = (Plugin) currentPipe.get(j);
                if (matchType(types, current)) {
                    return current;
                }

            }
        }

        return null;

    }

    private Plugin findLastPluginFromCurrentPipe(Class[] types) {

        if ((myCurrentPipe != null) && (myCurrentPipe.size() != 0)) {
            for (int i = myCurrentPipe.size() - 1; i >= 0; i--) {
                Plugin current = (Plugin) myCurrentPipe.get(i);
                if (matchType(types, current)) {
                    return current;
                }

            }
        }

        return null;

    }

    private boolean matchType(Class[] types, Object test) {

        for (int i = 0; i < types.length; i++) {
            if (types[i].isInstance(test)) {
                return true;
            }

        }

        return false;

    }

    /**
     * Returns Tassel data set after complete.
     *
     * @param event event
     */
    public void dataSetReturned(PluginEvent event) {
        DataSet tds = (DataSet) event.getSource();
        if ((tds != null) && (tds.getSize() != 0) && (myMainFrame != null)) {
            myMainFrame.getDataTreePanel().addDataSet(tds, DataTreePanel.NODE_TYPE_DEFAULT);
        }

    }

    /**
     * Returns progress of execution.
     *
     * @param event event
     */
    public void progress(PluginEvent event) {

        if (myMainFrame == null) {

            DataSet ds = (DataSet) event.getSource();
            if (ds != null) {
                List percentage = ds.getDataOfType(Integer.class);
                Plugin plugin = ds.getCreator();
                Integer lastValue = (Integer) myProgressValues.get(plugin);
                if (lastValue == null) {
                    lastValue = new Integer(0);
                }

                if (percentage.size() > 0) {
                    Datum datum = (Datum) percentage.get(0);
                    Integer percent = (Integer) datum.getData();
                    if (percent >= lastValue) {
                        myLogger.info(ds.getCreator().getClass().getName() + ": progress: " + percent.intValue() + "%");
                        lastValue = lastValue + 10;
                        myProgressValues.put(plugin, lastValue);
                    }
                }
            }

        }

    }
}
