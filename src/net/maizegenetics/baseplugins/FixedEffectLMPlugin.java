package net.maizegenetics.baseplugins;

import java.awt.Container;
import java.awt.Dialog;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;
import java.util.TreeSet;
import java.util.regex.Pattern;

import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JTextField;

import org.apache.log4j.Logger;

import net.maizegenetics.gui.ReportDestinationDialog;
import net.maizegenetics.jGLiM.LinearModelUtils;
import net.maizegenetics.jGLiM.dm.CovariateModelEffect;
import net.maizegenetics.jGLiM.dm.FactorModelEffect;
import net.maizegenetics.jGLiM.dm.ModelEffect;
import net.maizegenetics.jGLiM.dm.ModelEffectUtils;
import net.maizegenetics.jGLiM.dm.SweepFastLinearModel;
import net.maizegenetics.jGLiM.dm.SweepFastNestedModel;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrixFactory;
import net.maizegenetics.pal.alignment.GeneticMap;
import net.maizegenetics.pal.alignment.MarkerPhenotype;
import net.maizegenetics.pal.alignment.MarkerPhenotypeAdapter;
import net.maizegenetics.pal.alignment.MarkerPhenotypeAdapterUtils;
import net.maizegenetics.pal.alignment.Phenotype;
import net.maizegenetics.pal.alignment.SimplePhenotype;
import net.maizegenetics.pal.alignment.Trait;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.pal.report.SimpleTableReport;
import net.maizegenetics.pal.report.TableReport;
import net.maizegenetics.plugindef.AbstractPlugin;
import net.maizegenetics.plugindef.DataSet;
import net.maizegenetics.plugindef.Datum;

public class FixedEffectLMPlugin extends AbstractPlugin {

    private static final Logger myLogger = Logger.getLogger(FixedEffectLMPlugin.class);
    private GeneticMap myMap;
    private boolean hasMap;
    private boolean writeOutputToFile = false;
    private String outputName = null;
    private boolean filterOutput = false;
    private double maxp = 1;
    private boolean reportBLUEs = true;
    private boolean permute = false;
    private int numberOfPermutations = 1000;
    private double[][] minP;
    private Random randomizer = new Random();

    public FixedEffectLMPlugin(Frame parentFrame, boolean isInteractive) {
        super(parentFrame, isInteractive);
    }

    public DataSet performFunction(DataSet input) {

        try {
            List<Datum> maps = input.getDataOfType(GeneticMap.class);
            if (maps.size() > 0) {
                myMap = (GeneticMap) maps.get(0).getData();
                hasMap = true;
            } else {
                myMap = null;
                hasMap = false;
            }

            List<Datum> datasets = input.getDataOfType(new Class[]{Phenotype.class, MarkerPhenotype.class});
            if (datasets.size() < 1) {
                String msg = "Error in performFunction: No appropriate dataset selected.";
                myLogger.error(msg);
                if (isInteractive()) {
                    JOptionPane.showMessageDialog(getParentFrame(), msg, "Error in GLM", JOptionPane.ERROR_MESSAGE);
                }

                return null;
            }

            if (isInteractive()) {

                //glm options: permutations
                GLMOptionsDialog options = new GLMOptionsDialog(getParentFrame());
                permute = options.permute();
                numberOfPermutations = options.getPermutationNumber();
                options.dispose();

                //choose the report destination: GUI or text file
                ReportDestinationDialog rdd = new ReportDestinationDialog();
                rdd.setLocationRelativeTo(getParentFrame());
                rdd.setVisible(true);
                if (!rdd.isOkayChecked()) {
                    return null;
                }
                writeOutputToFile = rdd.wasUseFileChecked();
                if (writeOutputToFile) {
                    outputName = rdd.getOutputFileName();
                }
                filterOutput = rdd.wasRestrictOutputChecked();
                if (filterOutput) {
                    maxp = rdd.getMaxP();
                }
                rdd.dispose();
            }

            LinkedList<Datum> results = new LinkedList<Datum>();
            for (Datum dataset : datasets) {
                LinkedList<Datum> result = null;
                try {
                    result = processDatum(dataset);
                } catch (Exception e) {
                    myLogger.error(e.toString());
                    e.printStackTrace();
                    if (isInteractive()) {
                        JOptionPane.showMessageDialog(getParentFrame(), e.toString(), "Error in GLM", JOptionPane.ERROR_MESSAGE);
                    }
                }
                if (result != null) {
                    results.addAll(result);
                }
                fireDataSetReturned(new DataSet(result, this));
            }
            return new DataSet(results, this);
        } finally {
            fireProgress(100);
        }

    }

    private LinkedList<Datum> processDatum(Datum dataset) {
        ArrayList<Object[]> markerTestResults = new ArrayList<Object[]>();
        ArrayList<Object[]> alleleEstimateResults = new ArrayList<Object[]>();
        boolean hasAlleleNames = false;
        String blank = new String("");

        //build the base model
        MarkerPhenotypeAdapter theAdapter;
        if (dataset.getDataType().equals(MarkerPhenotype.class)) {
            theAdapter = new MarkerPhenotypeAdapter((MarkerPhenotype) dataset.getData());
        } else {
            theAdapter = new MarkerPhenotypeAdapter((Phenotype) dataset.getData());
        }

        //numberof markers
        int numberOfMarkers = theAdapter.getNumberOfMarkers();

        //are there markers? If not calculate BLUEs for taxa
        if (numberOfMarkers == 0) {
            return calculateBLUEsFromPhenotypes(theAdapter, dataset.getName());
        }

        //numbers of various model components
        int numberOfCovariates = theAdapter.getNumberOfCovariates();
        int numberOfFactors = theAdapter.getNumberOfFactors();
        int numberOfPhenotypes = theAdapter.getNumberOfPhenotypes();

        //calculate total iterations
        int expectedIterations = numberOfPhenotypes * numberOfMarkers;
        int iterationsSofar = 0;
        int percentFinished = 0;

        //open files for output if that was chosen
        File tempFile = null;
        File ftestFile = null;
        File blueFile = null;
        BufferedWriter ftestWriter = null;
        BufferedWriter BLUEWriter = null;
        String ftestHeader = "Trait\tMarker\tLocus\tLocus_pos\tChr\tChr_pos\tmarker_F\tmarker_p\tperm_p\tmarkerR2\tmarkerDF\tmarkerMS\terrorDF\terrorMS\tmodelDF\tmodelMS";
        String BLUEHeader = "Trait\tMarker\tObs\tLocus\tLocus_pos\tChr\tChr_pos\tAllele\tEstimate";
        if (writeOutputToFile) {
            //create the file objects; avoid overwriting permanent files
            String outputbase = outputName;
            if (outputbase.endsWith(".txt")) {
                int index = outputbase.lastIndexOf(".");
                outputbase = outputbase.substring(0, index);
            }
            String datasetNameNoSpace = dataset.getName().trim().replaceAll("\\ ", "_");

            ftestFile = new File(outputbase + "_" + datasetNameNoSpace + "_ftest.txt");
            int count = 0;
            while (ftestFile.exists()) {  //avoids overwriting existing files by appending a number
                count++;
                ftestFile = new File(outputbase + "_" + datasetNameNoSpace + "_ftest" + count + ".txt");
            }

            blueFile = new File(outputbase + "_" + datasetNameNoSpace + "_BLUEs.txt");
            count = 0;
            while (blueFile.exists()) {
                count++;
                blueFile = new File(outputbase + "_" + datasetNameNoSpace + "_BLUEs" + count + ".txt");
            }

            tempFile = new File(outputbase + "_" + datasetNameNoSpace + "_ftest.tmp");

            try {
                if (permute) {
                    ftestWriter = new BufferedWriter(new FileWriter(tempFile));
                    ftestWriter.write(ftestHeader);
                    ftestWriter.newLine();
                } else {
                    ftestWriter = new BufferedWriter(new FileWriter(ftestFile));
                    ftestWriter.write(ftestHeader);
                    ftestWriter.newLine();
                }

                if (reportBLUEs) {
                    BLUEWriter = new BufferedWriter(new FileWriter(blueFile));
                    BLUEWriter.write(BLUEHeader);
                    BLUEWriter.newLine();
                }
            } catch (IOException e) {
                myLogger.error("Failed to open file for output");
                myLogger.error(e);
                return null;
            }
        }

        //initialize minP
        if (permute) {
            minP = new double[numberOfPhenotypes][numberOfPermutations];
            for (int i = 0; i < numberOfPermutations; i++) {
                for (int j = 0; j < numberOfPhenotypes; j++) {
                    minP[j][i] = 1;
                }
            }
        }

        //cycle through the phenotypes
        //notation: 
        //X is the design matrix without the markers, rows of X will be deleted if marker data is missing
        //Xm is the design matrix with markers
        //y is the data

        for (int ph = 0; ph < numberOfPhenotypes; ph++) {

            //get phenotype data
            double[] phenotypeData = theAdapter.getPhenotypeValues(ph);

            //keep track of missing rows
            boolean[] missing = theAdapter.getMissingPhenotypes(ph);

            //get factors
            ArrayList<String[]> factorList = MarkerPhenotypeAdapterUtils.getFactorList(theAdapter, ph, missing);

            //get covariates
            ArrayList<double[]> covariateList = MarkerPhenotypeAdapterUtils.getCovariateList(theAdapter, ph, missing);

            //if permutation analysis has been requested, generate a set of permuted data
            double[][] permutedData = null;
            if (permute) {
                permutedData = permuteData(phenotypeData, missing, factorList, covariateList, theAdapter);
            }

            //cycle through the markers
            for (int m = 0; m < numberOfMarkers; m++) {
                Object[] markerData = theAdapter.getMarkerValue(ph, m);
                boolean[] finalMissing = new boolean[missing.length];
                System.arraycopy(missing, 0, finalMissing, 0, missing.length);
                MarkerPhenotypeAdapterUtils.updateMissing(finalMissing, theAdapter.getMissingMarkers(ph, m));

                //identify the rows with no missing data
                int[] nonmissingRows = MarkerPhenotypeAdapterUtils.getNonMissingIndex(finalMissing);

                //number of observations in this data set
                int numberOfObs = nonmissingRows.length;

                //the phenotype data
                double[] y = new double[numberOfObs];
                for (int i = 0; i < numberOfObs; i++) {
                    y[i] = phenotypeData[nonmissingRows[i]];
                }

                //keep track of X matrix columns, so that the position of marker effect estimates in beta can be determined
                int firstMarkerAlleleEstimate = 1; //for the mean

                //first, add the mean to the model
                ArrayList<ModelEffect> modelEffects = new ArrayList<ModelEffect>();
                FactorModelEffect meanEffect = new FactorModelEffect(new int[numberOfObs], false);
                meanEffect.setID("mean");
                modelEffects.add(meanEffect);

                //add factors
                if (numberOfFactors > 0) {
                    for (int f = 0; f < numberOfFactors; f++) {
                        String[] afactor = factorList.get(f);
                        String[] factorLabels = new String[numberOfObs];
                        for (int i = 0; i < numberOfObs; i++) {
                            factorLabels[i] = afactor[nonmissingRows[i]];
                        }
                        FactorModelEffect fme = new FactorModelEffect(ModelEffectUtils.getIntegerLevels(factorLabels), true, theAdapter.getFactorName(f));
                        modelEffects.add(fme);
                        firstMarkerAlleleEstimate += fme.getNumberOfLevels() - 1;
                    }
                }

                //add covariates
                if (numberOfCovariates > 0) {
                    for (int c = 0; c < numberOfCovariates; c++) {
                        double[] covar = new double[numberOfObs];
                        double[] covariateData = covariateList.get(c);
                        for (int i = 0; i < numberOfObs; i++) {
                            covar[i] = covariateData[nonmissingRows[i]];
                        }
                        modelEffects.add(new CovariateModelEffect(covar, theAdapter.getCovariateName(c)));
                        firstMarkerAlleleEstimate++;
                    }
                }

                //add marker
                ModelEffect markerEffect;
                boolean markerIsDiscrete = theAdapter.isMarkerDiscrete(m);
                ArrayList<Object> alleleNames = new ArrayList<Object>();

                if (markerIsDiscrete) {
                    Object[] markers = new Object[numberOfObs];
                    for (int i = 0; i < numberOfObs; i++) {
                        markers[i] = markerData[nonmissingRows[i]];
                    }

                    int[] markerLevels = ModelEffectUtils.getIntegerLevels(markers, alleleNames);
                    markerEffect = new FactorModelEffect(markerLevels, true, theAdapter.getMarkerName(m));
                    hasAlleleNames = true;
                } else {
                    double[] markerdbl = new double[numberOfObs];
                    for (int i = 0; i < numberOfObs; i++) {
                        markerdbl[i] = ((Double) markerData[nonmissingRows[i]]).doubleValue();
                    }
                    markerEffect = new CovariateModelEffect(markerdbl, theAdapter.getMarkerName(m));
                }

                int[] alleleCounts = markerEffect.getLevelCounts();

                modelEffects.add(markerEffect);
                int markerEffectNumber = modelEffects.size() - 1;

                //if taxa are replicated, add taxa:marker to use as error term
                //if the full list has duplicates, check to see if it still has
                Identifier[] taxaSublist = new Identifier[numberOfObs];
                Identifier[] taxa = theAdapter.getTaxa(ph);
                for (int i = 0; i < numberOfObs; i++) {
                    taxaSublist[i] = taxa[nonmissingRows[i]];
                }
                boolean areTaxaReplicated = containsDuplicates(taxaSublist);

                double[] markerSSdf = null, errorSSdf = null, modelSSdf = null;
                double F, p;
                double[] beta = null;
                if (areTaxaReplicated && markerIsDiscrete) {
                    ModelEffect taxaEffect = new FactorModelEffect(ModelEffectUtils.getIntegerLevels(taxaSublist), true);
                    modelEffects.add(taxaEffect);
                    SweepFastNestedModel sfnm = new SweepFastNestedModel(modelEffects, y);
                    double[] taxaSSdf = sfnm.getTaxaInMarkerSSdf();
                    double[] residualSSdf = sfnm.getErrorSSdf();
                    markerSSdf = sfnm.getMarkerSSdf();
                    errorSSdf = sfnm.getErrorSSdf();
                    modelSSdf = sfnm.getModelcfmSSdf();
                    F = markerSSdf[0] / markerSSdf[1] / taxaSSdf[0] * taxaSSdf[1];
                    try {
                        p = LinearModelUtils.Ftest(F, markerSSdf[1], taxaSSdf[1]);
                    } catch (Exception e) {
                        p = Double.NaN;
                    }
                    beta = sfnm.getBeta();
                    int markerdf = (int) markerSSdf[1];
                    if (permute && markerdf > 0) {
                        updatePermutationPValues(ph, permutedData, nonMissingIndex(missing, finalMissing), getXfromModelEffects(modelEffects), sfnm.getInverseOfXtX(), markerdf);
                    }
                } else {
                    SweepFastLinearModel sflm = new SweepFastLinearModel(modelEffects, y);
                    modelSSdf = sflm.getModelcfmSSdf();
                    markerSSdf = sflm.getMarginalSSdf(markerEffectNumber);
                    errorSSdf = sflm.getResidualSSdf();
                    F = markerSSdf[0] / markerSSdf[1] / errorSSdf[0] * errorSSdf[1];
                    try {
                        p = LinearModelUtils.Ftest(F, markerSSdf[1], errorSSdf[1]);
                    } catch (Exception e) {
                        p = Double.NaN;
                    }
                    beta = sflm.getBeta();

                    int markerdf = (int) markerSSdf[1];
                    if (permute && markerdf > 0) {
                        updatePermutationPValues(ph, permutedData, nonMissingIndex(missing, finalMissing), getXfromModelEffects(modelEffects), sflm.getInverseOfXtX(), markerdf);
                    }
                }

                if (!filterOutput || p < maxp) {
                    String traitname = theAdapter.getPhenotypeName(ph);
                    if (traitname == null) {
                        traitname = blank;
                    }
                    String marker = theAdapter.getMarkerName(m);
                    if (marker == null) {
                        marker = blank;
                    }
                    String locus = theAdapter.getLocusName(m);
                    Integer site = new Integer(theAdapter.getLocusPosition(m));
                    String chrname = "";
                    Double chrpos = Double.NaN;
                    if (hasMap) {
                        int ndx = -1;
                        ndx = myMap.getMarkerIndex(marker);
                        if (ndx > -1) {
                            chrname = myMap.getChromosome(ndx);
                            chrpos = myMap.getGeneticPosition(ndx);
                        }
                    }

                    Object[] result = new Object[16];
                    int col = 0;
                    result[col++] = traitname;
                    result[col++] = marker;
                    result[col++] = locus;
                    result[col++] = site;
                    result[col++] = chrname;
                    result[col++] = chrpos;
                    result[col++] = new Double(F);
                    result[col++] = new Double(p);
                    result[col++] = Double.NaN; //placeholder for permutation p value
                    result[col++] = new Double(markerSSdf[0] / (modelSSdf[0] + errorSSdf[0]));
                    result[col++] = new Double(markerSSdf[1]);
                    result[col++] = new Double(markerSSdf[0] / markerSSdf[1]);
                    result[col++] = new Double(errorSSdf[1]);
                    result[col++] = new Double(errorSSdf[0] / errorSSdf[1]);
                    result[col++] = new Double(modelSSdf[1]);
                    result[col++] = new Double(modelSSdf[0] / modelSSdf[1]);

                    if (writeOutputToFile) {
                        StringBuilder sb = new StringBuilder();
                        sb.append(result[0]);
                        for (int i = 1; i < 16; i++) {
                            sb.append("\t").append(result[i]);
                        }
                        try {
                            ftestWriter.write(sb.toString());
                            ftestWriter.newLine();
                        } catch (IOException e) {
                            myLogger.error("Failed to write output to ftest file. Ending prematurely");
                            try {
                                ftestWriter.flush();
                                BLUEWriter.flush();
                            } catch (Exception e1) {
                            }
                            myLogger.error(e);
                            return null;
                        }
                    } else {
                        markerTestResults.add(result);
                    }

                    //output allele estimates
                    //{"Trait", "Marker", "Locus", "Locus.Postion", "Chr", "Chr.pos", "Allele", "Estimate"}
                    int numberOfMarkerAlleles = alleleNames.size();
                    if (numberOfMarkerAlleles == 0) {
                        numberOfMarkerAlleles++;
                    }

                    for (int i = 0; i < numberOfMarkerAlleles; i++) {
                        result = new Object[9];
                        result[0] = traitname;
                        result[1] = marker;
                        result[2] = new Integer(alleleCounts[i]);
                        result[3] = locus;
                        result[4] = site;
                        result[5] = chrname;
                        result[6] = chrpos;

                        if (numberOfMarkerAlleles == 1) {
                            result[7] = "";
                        } else {
                            result[7] = alleleNames.get(i);
                        }

                        if (i == numberOfMarkerAlleles - 1) {
                            result[8] = 0.0;
                        } else {
                            result[8] = beta[firstMarkerAlleleEstimate + i];
                        }

                        if (writeOutputToFile) {
                            StringBuilder sb = new StringBuilder();
                            sb.append(result[0]);
                            for (int j = 1; j < 9; j++) {
                                sb.append("\t").append(result[j]);
                            }
                            try {
                                BLUEWriter.write(sb.toString());
                                BLUEWriter.newLine();
                            } catch (IOException e) {
                                myLogger.error("Failed to write output to ftest file. Ending prematurely");
                                try {
                                    ftestWriter.flush();
                                    BLUEWriter.flush();
                                } catch (Exception e1) {
                                }
                                myLogger.error(e);
                                return null;
                            }
                        } else {
                            alleleEstimateResults.add(result);
                        }

                    }
                }

                int tmpPercent = ++iterationsSofar * 100 / expectedIterations;
                if (tmpPercent > percentFinished) {
                    percentFinished = tmpPercent;
                    fireProgress(percentFinished);
                }
            }
        }

        fireProgress(0);

        //close the output files, if open
        if (writeOutputToFile) {
            try {
                ftestWriter.close();
                BLUEWriter.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        //if permuting, set up the null distributions, fill in perm_p
        HashMap<String, Integer> traitnameMap = new HashMap<String, Integer>();
        if (permute) {
            for (int ph = 0; ph < numberOfPhenotypes; ph++) {
                Arrays.sort(minP[ph]);
                traitnameMap.put(theAdapter.getPhenotypeName(ph), ph);
            }

            if (writeOutputToFile) {
                try {
                    BufferedReader tempReader = new BufferedReader(new FileReader(tempFile));
                    ftestWriter = new BufferedWriter(new FileWriter(ftestFile));

                    //write the header
                    ftestWriter.write(tempReader.readLine());
                    ftestWriter.newLine();

                    String input;
                    String[] data;
                    Pattern tab = Pattern.compile("\t");
                    while ((input = tempReader.readLine()) != null) {
                        data = tab.split(input);
                        String trait = data[0];
                        double pval = Double.parseDouble(data[7]);
                        int ph = traitnameMap.get(trait);
                        int ndx = Arrays.binarySearch(minP[ph], pval);
                        if (ndx < 0) {
                            ndx = -ndx - 1;
                        }
                        if (ndx == 0) {
                            ndx = 1;
                        }
                        data[8] = Double.toString((double) ndx / (double) numberOfPermutations);
                        ftestWriter.write(data[0]);
                        for (int i = 1; i < data.length; i++) {
                            ftestWriter.write("\t");
                            ftestWriter.write(data[i]);
                        }
                        ftestWriter.newLine();
                    }
                    tempReader.close();
                    ftestWriter.close();
                    tempFile.delete();
                } catch (IOException e) {
                    myLogger.error(e);
                }
            } else {
                for (Object[] result : markerTestResults) {
                    String trait = result[0].toString();
                    double pval = (Double) result[7];
                    int ph = traitnameMap.get(trait);
                    int ndx = Arrays.binarySearch(minP[ph], pval);
                    if (ndx < 0) {
                        ndx = -ndx - 1;
                    }
                    if (ndx == 0) {
                        ndx = 1;
                    }
                    result[8] = new Double((double) ndx / (double) numberOfPermutations);
                }
            }
        }


        //columns without meaningful information should be deleted from the data
        //check to see if the adapter has marker names
        //only report chr and chr.pos if there is a genetic map

        //set colnames

        String[] columnLabels = new String[]{"Trait", "Marker", "Locus", "Locus_pos", "Chr", "Chr_pos", "marker_F", "marker_p", "perm_p",
            "markerR2", "markerDF", "markerMS", "errorDF", "errorMS", "modelDF", "modelMS"};
        boolean hasMarkerNames = theAdapter.hasMarkerNames();

        LinkedList<Integer> outputList = new LinkedList<Integer>();
        outputList.add(0);
        if (hasMarkerNames) {
            outputList.add(1);
        }
        outputList.add(2);
        outputList.add(3);
        if (hasMap) {
            outputList.add(4);
            outputList.add(5);
        }
        outputList.add(6);
        outputList.add(7);
        if (permute) {
            outputList.add(8);
        }
        for (int i = 9; i < 16; i++) {
            outputList.add(i);
        }

        LinkedList<Datum> resultset = new LinkedList<Datum>();
        int nrows = markerTestResults.size();
        Object[][] table = new Object[nrows][];

        int numberOfColumns = outputList.size();
        String[] colnames = new String[numberOfColumns];
        int count = 0;
        for (Integer ndx : outputList) {
            colnames[count++] = columnLabels[ndx];
        }

        for (int i = 0; i < nrows; i++) {
            table[i] = new Object[numberOfColumns];
            Object[] result = markerTestResults.get(i);
            count = 0;
            for (Integer ndx : outputList) {
                table[i][count++] = result[ndx];
            }
        }

        //the table name and comments
        StringBuilder tableName = new StringBuilder("GLM_marker_test_");
        tableName.append(dataset.getName());
        StringBuilder comments = new StringBuilder("Tests of Marker-Phenotype Association");
        comments.append("GLM: fixed effect linear model\n");
        comments.append("Data set: ").append(dataset.getName());
        comments.append("\nmodel: trait = mean");
        for (int i = 0; i < theAdapter.getNumberOfFactors(); i++) {
            comments.append(" + ");
            comments.append(theAdapter.getFactorName(i));
        }
        for (int i = 0; i < theAdapter.getNumberOfCovariates(); i++) {
            comments.append(" + ");
            comments.append(theAdapter.getCovariateName(i));
        }
        comments.append(" + marker");

        if (writeOutputToFile) {
            comments.append("\nOutput written to " + ftestFile.getPath());
        }

        TableReport markerTestReport = new SimpleTableReport("Marker Test", colnames, table);
        resultset.add(new Datum(tableName.toString(), markerTestReport, comments.toString()));

        //add BLUE table
        int[] outputIndex;
        columnLabels = new String[]{"Trait", "Marker", "Obs", "Locus", "Locus_pos", "Chr", "Chr_pos", "Allele", "Estimate"};
        if (hasAlleleNames) {
            if (hasMarkerNames && hasMap) {
                outputIndex = new int[]{0, 1, 2, 3, 4, 5, 6, 7, 8};
            } else if (hasMarkerNames) {
                outputIndex = new int[]{0, 1, 2, 3, 4, 7, 8};
            } else if (hasMap) {
                outputIndex = new int[]{0, 2, 3, 4, 5, 6, 7, 8};
            } else {
                outputIndex = new int[]{0, 2, 3, 4, 7, 8};
            }
        } else {
            if (hasMarkerNames && hasMap) {
                outputIndex = new int[]{0, 1, 2, 3, 4, 5, 6, 8};
            } else if (hasMarkerNames) {
                outputIndex = new int[]{0, 1, 2, 3, 4, 8};
            } else if (hasMap) {
                outputIndex = new int[]{0, 2, 3, 4, 5, 6, 8};
            } else {
                outputIndex = new int[]{0, 2, 3, 4, 8};
            }
        }

        nrows = alleleEstimateResults.size();
        table = new Object[nrows][];
        numberOfColumns = outputIndex.length;
        colnames = new String[numberOfColumns];
        for (int j = 0; j < numberOfColumns; j++) {
            colnames[j] = columnLabels[outputIndex[j]];
        }
        for (int i = 0; i < nrows; i++) {
            table[i] = new Object[numberOfColumns];
            Object[] result = alleleEstimateResults.get(i);
            for (int j = 0; j < numberOfColumns; j++) {
                table[i][j] = result[outputIndex[j]];
            }
        }

        //the table name and comments
        tableName = new StringBuilder("GLM allele estimates for ");
        tableName.append(dataset.getName());
        comments = new StringBuilder("Marker allele effect estimates\n");
        comments.append("GLM: fixed effect linear model\n");
        comments.append("Data set: ").append(dataset.getName());
        comments.append("\nmodel: trait = mean");
        for (int i = 0; i < theAdapter.getNumberOfFactors(); i++) {
            comments.append(" + ");
            comments.append(theAdapter.getFactorName(i));
        }
        for (int i = 0; i < theAdapter.getNumberOfCovariates(); i++) {
            comments.append(" + ");
            comments.append(theAdapter.getCovariateName(i));
        }
        comments.append(" + marker");

        if (writeOutputToFile) {
            comments.append("\nOutput written to " + blueFile.getPath());
        }

        TableReport alleleEstimateReport = new SimpleTableReport("Allele Estimates", colnames, table);
        resultset.add(new Datum(tableName.toString(), alleleEstimateReport, comments.toString()));
        fireProgress(0);
        return resultset;
    }

    private void updatePermutationPValues(int phenotype, double[][] permutedData, int[] nonMissingIndex, DoubleMatrix X, DoubleMatrix G, int markerdf) {
        int nobs = nonMissingIndex.length;
        int rX = X.numberOfColumns();
        double dferror = nobs - rX;

        //determine L, the contrast
        DoubleMatrix L = DoubleMatrixFactory.DEFAULT.make(rX, markerdf, 0);
        for (int i = 0; i < markerdf; i++) {
            L.set(rX - markerdf + i, i, 1);
        }

        for (int p = 0; p < numberOfPermutations; p++) {
            DoubleMatrix y = DoubleMatrixFactory.DEFAULT.make(nobs, 1);
            for (int i = 0; i < nobs; i++) {
                y.set(i, 0, permutedData[p][nonMissingIndex[i]]);
            }
            DoubleMatrix yty = y.crossproduct(y);
            DoubleMatrix Xty = X.crossproduct(y);
            DoubleMatrix beta = G.mult(Xty);
            double errorSS = yty.minus(Xty.crossproduct(beta)).get(0, 0);
            DoubleMatrix Ltbeta = L.crossproduct(beta);
            DoubleMatrix invLtGL = L.crossproduct(G.mult(L));
            invLtGL.invert();
            double dfL = L.numberOfColumns();
            double F = (Ltbeta.crossproduct(invLtGL.mult(Ltbeta)).get(0, 0)) / errorSS / dfL * dferror;
            double pval;
            try {
                pval = LinearModelUtils.Ftest(F, dfL, dferror);
            } catch (Exception e) {
                pval = 1;
            }
            minP[phenotype][p] = Math.min(minP[phenotype][p], pval);
        }
    }

    private double[][] permuteData(double[] phenotypeData, boolean[] missing, ArrayList<String[]> factorList,
            ArrayList<double[]> covariateList, MarkerPhenotypeAdapter theAdapter) {
        int[] nonMissingIndex = MarkerPhenotypeAdapterUtils.getNonMissingIndex(missing);
        int nobs = nonMissingIndex.length;
        double[][] permutedData = new double[numberOfPermutations][nobs];

        //fit a model without markers
        double[] y = new double[nobs];
        for (int i = 0; i < nobs; i++) {
            y[i] = phenotypeData[nonMissingIndex[i]];
        }

        //first, add the mean to the model
        ArrayList<ModelEffect> modelEffects = new ArrayList<ModelEffect>();
        FactorModelEffect meanEffect = new FactorModelEffect(new int[nobs], false);
        meanEffect.setID("mean");
        modelEffects.add(meanEffect);

        //add factors
        int numberOfFactors = 0;
        if (factorList != null) {
            numberOfFactors = factorList.size();
        }
        if (numberOfFactors > 0) {
            for (int f = 0; f < numberOfFactors; f++) {
                String[] afactor = factorList.get(f);
                String[] factorLabels = new String[nobs];
                for (int i = 0; i < nobs; i++) {
                    factorLabels[i] = afactor[nonMissingIndex[i]];
                }
                FactorModelEffect fme = new FactorModelEffect(ModelEffectUtils.getIntegerLevels(factorLabels), true, theAdapter.getFactorName(f));
                modelEffects.add(fme);
            }
        }

        //add covariates
        int numberOfCovariates = 0;
        if (covariateList != null) {
            numberOfCovariates = covariateList.size();
        }
        if (numberOfCovariates > 0) {
            for (int c = 0; c < numberOfCovariates; c++) {
                double[] covar = new double[nobs];
                double[] covariateData = covariateList.get(c);
                for (int i = 0; i < nobs; i++) {
                    covar[i] = covariateData[nonMissingIndex[i]];
                }
                modelEffects.add(new CovariateModelEffect(covar, theAdapter.getCovariateName(c)));
            }
        }

        //create and solve the model
        SweepFastLinearModel theModel = new SweepFastLinearModel(modelEffects, y);
        DoubleMatrix yhat = theModel.getPredictedValues();
        DoubleMatrix r = theModel.getResiduals();
        int nrows = r.numberOfRows();
        double[] res = new double[nrows];
        for (int i = 0; i < nrows; i++) {
            res[i] = r.get(i, 0);
        }

        //permute residuals and add to yhat
        permutedData = new double[numberOfPermutations][nobs];
        for (int p = 0; p < numberOfPermutations; p++) {
            LinearModelUtils.shuffle(res, randomizer);
            for (int i = 0; i < nobs; i++) {
                permutedData[p][i] = yhat.get(i, 0) + res[i];
            }
        }

        return permutedData;
    }

    private DoubleMatrix getXfromModelEffects(ArrayList<ModelEffect> someEffects) {
        int n = someEffects.size();
        DoubleMatrix[][] mats = new DoubleMatrix[1][n];
        for (int i = 0; i < n; i++) {
            mats[0][i] = someEffects.get(i).getX();
        }
        return DoubleMatrixFactory.DEFAULT.compose(mats);
    }

    private int[] nonMissingIndex(boolean[] originalMissing, boolean[] finalMissing) {
        //count final not missing
        int finalCount = 0;
        for (boolean miss : finalMissing) {
            if (!miss) {
                finalCount++;
            }
        }
        int[] index = new int[finalCount];
        int n = originalMissing.length;
        int count = 0;
        finalCount = 0;
        for (int i = 0; i < n; i++) {
            if (!finalMissing[i]) {
                index[finalCount++] = count;
            }
            if (!originalMissing[i]) {
                count++;
            }
        }
        return index;
    }

    private LinkedList<Datum> calculateBLUEsFromPhenotypes(MarkerPhenotypeAdapter mpa, String datasetName) {
        if (isInteractive()) {
            String msg = "The data set you have selected does not contain any marker data. Do you want to calculate BLUEs (best linear unbiased estimates) of the phenotypes?";
            String title = "Calculate BLUEs";

            int action = JOptionPane.showConfirmDialog(getParentFrame(), msg, title, JOptionPane.YES_NO_OPTION);
            if (action != JOptionPane.YES_OPTION) {
                return null;
            }
        }


        LinkedList<Datum> theResults = new LinkedList<Datum>();
        LinkedList<Object[]> anovaResults = new LinkedList<Object[]>();
        LinkedList<double[]> blueList = new LinkedList<double[]>();
        ArrayList<ArrayList<Object>> taxaListList = new ArrayList<ArrayList<Object>>();

        //numbers of various model components
        int numberOfCovariates = mpa.getNumberOfCovariates();
        int numberOfFactors = mpa.getNumberOfFactors();
        int numberOfPhenotypes = mpa.getNumberOfPhenotypes();

        //for each phenotype do the following:
        for (int ph = 0; ph < numberOfPhenotypes; ph++) {

            //get phenotype data
            double[] phenotypeData = mpa.getPhenotypeValues(ph);

            //keep track of missing rows
            boolean[] missing = mpa.getMissingPhenotypes(ph);

            //get factors
            ArrayList<String[]> factorList = MarkerPhenotypeAdapterUtils.getFactorList(mpa, ph, missing);

            //get covariates
            ArrayList<double[]> covariateList = MarkerPhenotypeAdapterUtils.getCovariateList(mpa, ph, missing);

            //identify the rows with no missing data
            int[] nonmissingRows = MarkerPhenotypeAdapterUtils.getNonMissingIndex(missing);

            //number of observations in this data set
            int numberOfObs = nonmissingRows.length;

            //the phenotype data
            double[] y = new double[numberOfObs];
            for (int i = 0; i < numberOfObs; i++) {
                y[i] = phenotypeData[nonmissingRows[i]];
            }

            //first, add the mean to the model
            ArrayList<ModelEffect> modelEffects = new ArrayList<ModelEffect>();
            FactorModelEffect meanEffect = new FactorModelEffect(new int[numberOfObs], false);
            meanEffect.setID("mean");
            modelEffects.add(meanEffect);

            //now add the taxa
            Identifier[] alltaxa = mpa.getTaxa(ph);
            Identifier[] taxa = new Identifier[numberOfObs];
            ArrayList<Object> taxaIds = new ArrayList<Object>();
            for (int i = 0; i < numberOfObs; i++) {
                taxa[i] = alltaxa[nonmissingRows[i]];
            }
            int[] taxaLevels = ModelEffectUtils.getIntegerLevels(taxa, taxaIds);
            taxaListList.add(taxaIds);
            FactorModelEffect taxaEffect = new FactorModelEffect(taxaLevels, true, "Taxa");
            modelEffects.add(taxaEffect);

            //add factors
            if (numberOfFactors > 0) {
                for (int f = 0; f < numberOfFactors; f++) {
                    String[] afactor = factorList.get(f);
                    String[] factorLabels = new String[numberOfObs];
                    for (int i = 0; i < numberOfObs; i++) {
                        factorLabels[i] = afactor[nonmissingRows[i]];
                    }
                    FactorModelEffect fme = new FactorModelEffect(ModelEffectUtils.getIntegerLevels(factorLabels), true, mpa.getFactorName(f));
                    modelEffects.add(fme);
                }
            }

            //add covariates
            if (numberOfCovariates > 0) {
                for (int c = 0; c < numberOfCovariates; c++) {
                    double[] covar = new double[numberOfObs];
                    double[] covariateData = covariateList.get(c);
                    for (int i = 0; i < numberOfObs; i++) {
                        covar[i] = covariateData[nonmissingRows[i]];
                    }
                    modelEffects.add(new CovariateModelEffect(covar, mpa.getCovariateName(c)));
                }
            }

            //solve the model
            SweepFastLinearModel sflm = new SweepFastLinearModel(modelEffects, y);
            double[] taxaSSdf = sflm.getMarginalSSdf(1);
            double[] modelSSdf = sflm.getFullModelSSdf();
            double[] errorSSdf = sflm.getResidualSSdf();
            double[] beta = sflm.getBeta();
            double F, p;
            F = taxaSSdf[0] / taxaSSdf[1] / errorSSdf[0] * errorSSdf[1];
            try {
                p = LinearModelUtils.Ftest(F, taxaSSdf[1], errorSSdf[1]);
            } catch (Exception e) {
                p = Double.NaN;
            }

            //record {"Trait", "F", "p", "taxaDF", "taxaMS", "errorDF", "errorMS", "modelDF", "modelMS"}
            Object[] result = new Object[9];
            result[0] = mpa.getPhenotypeName(ph);
            result[1] = new Double(F);
            result[2] = new Double(p);
            result[3] = new Double(taxaSSdf[1]);
            result[4] = new Double(taxaSSdf[0] / taxaSSdf[1]);
            result[5] = new Double(errorSSdf[1]);
            result[6] = new Double(errorSSdf[0] / errorSSdf[1]);
            result[7] = new Double(modelSSdf[1]);
            result[8] = new Double(modelSSdf[0] / modelSSdf[1]);
            anovaResults.add(result);

            // calculate the BLUEs for taxa
            // first calculate the average of the level estimates for all other factors (including the mean)
            double overallMean = beta[0]; //the mean
            int nEffects = modelEffects.size();
            int start = 0;
            for (int i = 1; i < nEffects; i++) {
                ModelEffect me = modelEffects.get(i);
                if (me instanceof FactorModelEffect && !me.getID().equals("Taxa")) {
                    FactorModelEffect fme = (FactorModelEffect) me;
                    int nLevels = fme.getNumberOfLevels();
                    int nEstimates;
                    if (fme.getRestricted()) {
                        nEstimates = nLevels - 1;
                    } else {
                        nEstimates = nLevels;
                    }
                    double factorMean = 0;
                    for (int j = 0; j < nEstimates; j++) {
                        factorMean += beta[j + start];
                    }
                    factorMean /= nLevels;
                    overallMean += factorMean;
                    start += nEstimates;
                } else {
                    start += me.getNumberOfLevels();
                }
            }

            int n = taxaIds.size();
            double[] blues = new double[n];
            for (int i = 0; i < n - 1; i++) {
                blues[i] = beta[i + 1] + overallMean;
            }
            blues[n - 1] = overallMean;
            blueList.add(blues);
            taxaListList.add(taxaIds);

        }

        //add the anova table to the results
        String[] anovaColumnLabels = new String[]{"Trait", "F", "p", "taxaDF", "taxaMS", "errorDF", "errorMS", "modelDF", "modelMS"};
        Object[][] table = new Object[anovaResults.size()][];
        anovaResults.toArray(table);
        String datumName = "Phenotype ANOVA from " + datasetName;
        StringBuilder datumComments = new StringBuilder("ANOVA for Phenotypes using GLM\n");
        datumComments.append("Data set: ").append(datasetName);
        datumComments.append("\nmodel: trait = mean + taxa");
        for (int i = 0; i < mpa.getNumberOfFactors(); i++) {
            datumComments.append(" + ");
            datumComments.append(mpa.getFactorName(i));
        }
        for (int i = 0; i < mpa.getNumberOfCovariates(); i++) {
            datumComments.append(" + ");
            datumComments.append(mpa.getCovariateName(i));
        }

        SimpleTableReport str = new SimpleTableReport(datumName, anovaColumnLabels, table);
        Datum theAnova = new Datum(datumName, str, datumComments.toString());
        theResults.add(theAnova);

        //add the BLUEs to the results
        //the individual phenotypes must first be merged to create a unified table with taxa as rows and phenotypes as columns

        TreeSet<Identifier> taxaSet = new TreeSet<Identifier>();
        for (ArrayList<Object> list : taxaListList) {
            for (Object taxon : list) {
                taxaSet.add((Identifier) taxon);
            }
        }

        HashMap<Identifier, Integer> taxaMap = new HashMap<Identifier, Integer>();
        int count = 0;
        for (Identifier taxon : taxaSet) {
            taxaMap.put(taxon, count++);
        }
        String[] blueColumnLabels = new String[numberOfPhenotypes + 1];
        blueColumnLabels[0] = "Taxa";
        for (int i = 0; i < numberOfPhenotypes; i++) {
            blueColumnLabels[i + 1] = mpa.getPhenotypeName(i);
        }

        int nrows = taxaSet.size();
        double[][] blues = new double[nrows][numberOfPhenotypes];
        for (int r = 0; r < nrows; r++) {
            for (int c = 0; c < numberOfPhenotypes; c++) {
                blues[r][c] = Double.NaN;
            }
        }

        LinkedList<Trait> traitList = new LinkedList<Trait>();
        for (int c = 0; c < numberOfPhenotypes; c++) {
            traitList.add(new Trait(mpa.getPhenotypeName(c), false, Trait.TYPE_DATA));
            double[] pheno = blueList.get(c);
            int n = pheno.length;
            ArrayList<Object> taxaList = taxaListList.get(c);
            for (int i = 0; i < n; i++) {
                int ndx = taxaMap.get(taxaList.get(i));
                if (ndx > -1) {
                    blues[ndx][c] = pheno[i];
                }
            }
        }

        Identifier[] taxaIds = new Identifier[taxaSet.size()];
        taxaSet.toArray(taxaIds);
        Phenotype thePhenotype = new SimplePhenotype(new SimpleIdGroup(taxaIds), traitList, blues);

        theResults.add(new Datum("BLUEs_" + datasetName, thePhenotype, "BLUEs calculated from " + datasetName));
        return theResults;
    }

    @Override
    public String getButtonName() {
        return "GLM";
    }

    @Override
    public ImageIcon getIcon() {
        URL imageURL = FixedEffectLMPlugin.class.getResource("images/LinearAssociation.gif");
        if (imageURL == null) {
            return null;
        } else {
            return new ImageIcon(imageURL);
        }
    }

    @Override
    public String getToolTipText() {
        return "Use fixed effect model to test associations";
    }

    private boolean containsDuplicates(Identifier[] ids) {
        HashSet<Identifier> idSet = new HashSet<Identifier>();
        for (Identifier id : ids) {
            idSet.add(id);
        }
        if (idSet.size() < ids.length) {
            return true;
        }
        return false;
    }

    public void setOutputFile(String filename) {
        writeOutputToFile = true;
        outputName = filename;
    }

    public void setRestrictOutput(boolean restrict) {
        filterOutput = restrict;
    }

    public void setMaxP(double value) {
        filterOutput = true;
        maxp = value;
    }

    public void setPermute(boolean permute) {
        this.permute = permute;
    }

    public void setNumberOfPermutations(int numberOfPermutations) {
        this.numberOfPermutations = numberOfPermutations;
    }

    /**
     * Sets the randomizer seed so that permutation results are reproducible for
     * testing. The same seed will reproduce the same sequence of pseudo-random
     * numbers.
     *
     * @param seed
     */
    public void setRandomizer(long seed) {
        randomizer.setSeed(seed);
    }

    class GLMOptionsDialog extends JDialog {

        JCheckBox chkPermute = new JCheckBox("Use permutation test for markers", false);
        JTextField txtPermutationNumber = new JTextField(8);

        GLMOptionsDialog(Frame parentFrame) {
            super(parentFrame);
            setLocationRelativeTo(parentFrame);
            Container content = getContentPane();
            content.setLayout(new GridBagLayout());
            setModalityType(Dialog.ModalityType.APPLICATION_MODAL);
            setTitle("GLM Options");
            txtPermutationNumber.setText("1000");
            txtPermutationNumber.setEnabled(false);
            JLabel lblPermNumber = new JLabel("Number of Permutations: ");

            JButton btnOK = new JButton("OK");
            btnOK.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent arg0) {
                    setVisible(false);
                }
            });

            chkPermute.addActionListener(new ActionListener() {
                @Override
                public void actionPerformed(ActionEvent evt) {
                    if (chkPermute.isSelected()) {
                        txtPermutationNumber.setEnabled(true);
                    } else {
                        txtPermutationNumber.setEnabled(false);
                    }
                }
            });

            GridBagConstraints gbc = new GridBagConstraints();
            gbc.gridx = 0;
            gbc.gridy = 0;
            gbc.gridwidth = 2;
            gbc.anchor = GridBagConstraints.WEST;
            gbc.insets = new Insets(10, 8, 5, 8);
            content.add(chkPermute, gbc);

            gbc.gridy++;
            gbc.gridwidth = 1;
            content.add(lblPermNumber, gbc);

            gbc.gridx++;
            content.add(txtPermutationNumber, gbc);

            gbc.gridx = 0;
            gbc.gridy++;
            gbc.gridwidth = 2;
            gbc.anchor = GridBagConstraints.CENTER;
            content.add(btnOK, gbc);

            pack();
            setVisible(true);

        }

        int getPermutationNumber() {
            return Integer.parseInt(txtPermutationNumber.getText());
        }

        boolean permute() {
            return chkPermute.isSelected();
        }
    }
}
