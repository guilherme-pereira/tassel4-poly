package net.maizegenetics.stats.MLM;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.TreeSet;

import javax.swing.JOptionPane;

import org.apache.log4j.Logger;


import net.maizegenetics.baseplugins.MLMPlugin;
import net.maizegenetics.jGLiM.LinearModelUtils;
import net.maizegenetics.jGLiM.dm.FactorModelEffect;
import net.maizegenetics.jGLiM.dm.ModelEffectUtils;
import net.maizegenetics.jGLiM.dm.SweepFast;
import net.maizegenetics.jGLiM.SymmetricMatrixInverterDM;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrix;
import net.maizegenetics.matrixalgebra.Matrix.DoubleMatrixFactory;
import net.maizegenetics.pal.alignment.GeneticMap;
import net.maizegenetics.pal.alignment.MarkerPhenotype;
import net.maizegenetics.pal.alignment.MarkerPhenotypeAdapter;
import net.maizegenetics.pal.alignment.MarkerPhenotypeAdapterUtils;
import net.maizegenetics.pal.alignment.Phenotype;
import net.maizegenetics.pal.distance.DistanceMatrix;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.pal.report.SimpleTableReport;
import net.maizegenetics.pal.tree.UPGMATree;
import net.maizegenetics.plugindef.Datum;
import net.maizegenetics.stats.EMMA.EMMAforDoubleMatrix;

public class CompressedMLMusingDoubleMatrix {

    private static final Logger myLogger = Logger.getLogger(CompressedMLMusingDoubleMatrix.class);
    private final boolean useCompression;
    private final boolean useP3D;
    private final double compression;
//    private DoubleMatrix ZKZ;
    private final MarkerPhenotypeAdapter theAdapter;
    private final MLMPlugin parentPlugin;
    private final DistanceMatrix kinshipMatrix;
    private double resvar, genvar, lnlk;
    private boolean testMarkers = true;
    private final ArrayList<Object[]> compressionResults = new ArrayList<Object[]>();
    private final ArrayList<Object[]> resultsMain = new ArrayList<Object[]>();
    private final ArrayList<Object[]> resultsAlleles = new ArrayList<Object[]>();
    private final ArrayList<Object[]> resultsBlups = new ArrayList<Object[]>();
    private boolean[] prevMissing = null;
    private SymmetricMatrixInverterDM Vminus = null;
    private String datasetName;
    private final String[] headerMain;
    private final String[] headerAlleles;
    private final String[] headerBlups = new String[]{};
    private final String[] headerCompression = new String[]{"Trait", "# groups", "Compression", "-2LnLk", "Var_genetic", "Var_error"};
    private GeneticMap myGeneticMap;

    public CompressedMLMusingDoubleMatrix(MLMPlugin parentPlugin, Datum dataset, DistanceMatrix kinshipMatrix, boolean useCompression, boolean useP3D, double compression, GeneticMap map) {
        this.parentPlugin = parentPlugin;
        if (dataset.getData().getClass().equals(MarkerPhenotype.class)) {
            theAdapter = new MarkerPhenotypeAdapter((MarkerPhenotype) dataset.getData());
        } else if (dataset.getData() instanceof Phenotype) {
            theAdapter = new MarkerPhenotypeAdapter((Phenotype) dataset.getData());
//			testMarkers = false;
        } else {
            theAdapter = null;
        }
        this.kinshipMatrix = kinshipMatrix;
        this.useCompression = useCompression;
        this.useP3D = useP3D;
        this.compression = compression;
        datasetName = dataset.getName();
        myGeneticMap = map;
        if (myGeneticMap == null) {
            headerMain = new String[]{"Trait", "Marker", "Locus", "Site", "df", "F", "p", "errordf", "markerR2", "Genetic Var", "Residual Var", "-2LnLikelihood"};
            headerAlleles = new String[]{"Trait", "Marker", "Locus", "Site", "Allele", "Effect", "Obs"};
        } else {
            headerMain = new String[]{"Trait", "Marker", "Chr", "Pos", "Locus", "Site", "df", "F", "p", "errordf", "markerR2", "Genetic Var", "Residual Var", "-2LnLikelihood"};
            headerAlleles = new String[]{"Trait", "Marker", "Chr", "Pos", "Locus", "Site", "Allele", "Effect", "Obs"};
        }
    }

    public List<Datum> solve() {

        int numberOfMarkers = theAdapter.getNumberOfMarkers();
        int numberOfPhenotypes = theAdapter.getNumberOfPhenotypes();

        //calculate total iterations
        int expectedIterations = numberOfPhenotypes * numberOfMarkers;
        int iterationsSofar = 0;

        //open output files and write headers, if that is where data should go
        BufferedWriter mainOutputWriter = null;
        BufferedWriter alleleOutputWriter = null;
        BufferedWriter blupOutputWriter = null;
        if (parentPlugin.isWriteOutputToFile()) {
            String outputbase = parentPlugin.getOutputName();
            if (outputbase.endsWith(".txt")) {
                int index = outputbase.lastIndexOf(".");
                outputbase = outputbase.substring(0, index);
            }
            String datasetNameNoSpace = datasetName.trim().replaceAll("\\ ", "_");
            File mainFile = new File(outputbase + "_" + datasetNameNoSpace + "_stats.txt");
            File alleleFile = new File(outputbase + "_" + datasetNameNoSpace + "_effects.txt");
            File blupFile = new File(outputbase + "_" + datasetNameNoSpace + "_blups.txt");
            try {
                mainOutputWriter = new BufferedWriter(new FileWriter(mainFile));
                mainOutputWriter.write(getTabbedStringFromArray(headerMain));
                mainOutputWriter.newLine();

                alleleOutputWriter = new BufferedWriter(new FileWriter(alleleFile));
                alleleOutputWriter.write(getTabbedStringFromArray(headerAlleles));
                alleleOutputWriter.newLine();

                //blups are not implemented

            } catch (IOException e) {
                myLogger.error("Failed to open file for output.");
                myLogger.error(e);
                return null;
            }
        }

        //cycle through the phenotypes
        for (int ph = 0; ph < numberOfPhenotypes; ph++) {
            //get phenotype data
            double[] phenotypeData = theAdapter.getPhenotypeValues(ph);

            //get the taxa
            Identifier[] theTaxa = theAdapter.getTaxa(ph);

            //keep track of missing rows
            boolean[] missing = theAdapter.getMissingPhenotypes(ph);

            //get factors
            ArrayList<String[]> factorList = MarkerPhenotypeAdapterUtils.getFactorList(theAdapter, ph, missing);

            //get covariates
            ArrayList<double[]> covariateList = MarkerPhenotypeAdapterUtils.getCovariateList(theAdapter, ph, missing);

            //update missing for taxa not in the kinship matrix or the distance matrix.
            //Create kinship and distance matrices with taxa in phenotype
            IdGroup nonmissingIds = updateMissingWithKinship(missing, theTaxa);
            DistanceMatrix kin = new DistanceMatrix(kinshipMatrix, nonmissingIds);

            //calculate the number of nonmissing observations
            int totalObs = missing.length;
            int nonMissingObs = 0;
            for (boolean m : missing) {
                if (!m) {
                    nonMissingObs++;
                }
            }

            //create phenotype matrix
            int count = 0;
            DoubleMatrix y = DoubleMatrixFactory.DEFAULT.make(nonMissingObs, 1);
            Identifier[] nonmissingTaxa = new Identifier[nonMissingObs];
            for (int i = 0; i < totalObs; i++) {
                if (!missing[i]) {
                    y.set(count, 0, phenotypeData[i]);
                    nonmissingTaxa[count] = theTaxa[i];
                    count++;
                }
            }

            //create the Z matrix
            count = 0;
            DoubleMatrix Z = DoubleMatrixFactory.DEFAULT.make(nonMissingObs, kin.getIdCount());
            for (int i = 0; i < totalObs; i++) {
                if (!missing[i]) {
                    Z.set(count++, kin.whichIdNumber(theTaxa[i]), 1);
                }
            }

            //fixed effects matrix
            DoubleMatrix fixed = LinearModelUtils.createFixedEffectsArray(factorList, covariateList, missing);

            //fit data without markers
            DoubleMatrix[] zk = computeZKZ(y, fixed, Z, kin, theAdapter.getPhenotypeName(ph));
//            EMMAforDoubleMatrix emlm = new EMMAforDoubleMatrix(y, fixed, zk[0], zk[1], 0, Double.NaN);
            EMMAforDoubleMatrix emlm = new EMMAforDoubleMatrix(y, fixed, zk[1], zk[0], 0, Double.NaN);
            emlm.solve();
            genvar = emlm.getVarRan();
            resvar = emlm.getVarRes();
            lnlk = emlm.getLnLikelihood();
            int baseModeldf = emlm.getDfModel();

            //record the results
            Object[] tableRow;
            if (myGeneticMap != null) {
                //{"Trait","Marker","Chr","Pos","Locus","Site","df","F","p","errordf","MarkerR2","Genetic Var","Residual Var", "-2LnLikelihood"}
                tableRow = new Object[]{theAdapter.getPhenotypeName(ph),
                            "None",
                            "",
                            "",
                            "",
                            "",
                            new Integer(0),
                            new Double(Double.NaN),
                            new Double(Double.NaN),
                            new Double(nonMissingObs - baseModeldf),
                            new Double(Double.NaN),
                            new Double(genvar),
                            new Double(resvar),
                            new Double(-2 * lnlk)};
            } else {
                //{"Trait","Marker","Locus","Site","df","F","p","errordf","MarkerR2","Genetic Var","Residual Var", "-2LnLikelihood"}
                tableRow = new Object[]{theAdapter.getPhenotypeName(ph),
                            "None",
                            "",
                            "",
                            new Integer(0),
                            new Double(Double.NaN),
                            new Double(Double.NaN),
                            new Double(nonMissingObs - baseModeldf),
                            new Double(Double.NaN),
                            new Double(genvar),
                            new Double(resvar),
                            new Double(-2 * lnlk)};
            }

            if (parentPlugin.isWriteOutputToFile()) {
                if (!writeObjectToFile(mainOutputWriter, tableRow)) {
                    return null;
                }
            } else {
                resultsMain.add(tableRow);
            }

            //the BLUPs
            //not implemented

            if (useP3D) {
            	DoubleMatrix ZKZ = zk[0].mult(zk[1]).tcrossproduct(zk[0]);
                Vminus = new SymmetricMatrixInverterDM(calculateV(ZKZ, genvar, resvar));
            }

            //iterate markers
            if (testMarkers) {
                for (int m = 0; m < numberOfMarkers; m++) {
                    Object[] markerData = theAdapter.getMarkerValue(ph, m);
                    boolean[] missingMarkers = theAdapter.getMissingMarkers(ph, m);
                    boolean[] finalMissing = new boolean[missing.length];
                    System.arraycopy(missing, 0, finalMissing, 0, missing.length);
                    MarkerPhenotypeAdapterUtils.updateMissing(finalMissing, missingMarkers);
                    boolean[] additionalMissing = new boolean[nonMissingObs];

                    //only data for which missing=false are in the Z matrix
                    //the block below finds the rows of Z that have no marker data.
                    //Those rows/columns will need to be removed from ZKZ or from V, depending on the analysis method.
                    int missingCount = 0;
                    for (int i = 0; i < totalObs; i++) {
                        if (!missing[i]) {
                            if (finalMissing[i]) {
                                additionalMissing[missingCount] = true;
                            } else {
                                additionalMissing[missingCount] = false;
                            }
                            missingCount++;
                        }
                    }

                    //adjust y for missing data
                    int numberOfRowsKept = 0;
                    for (boolean miss : additionalMissing) {
                        if (!miss) {
                            numberOfRowsKept++;
                        }
                    }
                    int[] keepIndex = new int[numberOfRowsKept];
                    count = 0;
                    for (int i = 0; i < nonMissingObs; i++) {
                        if (!additionalMissing[i]) {
                            keepIndex[count++] = i;
                        }
                    }

                    DoubleMatrix ymarker = y.getSelection(keepIndex, null);
                    count = 0;
                    for (int i = 0; i < totalObs; i++) {
                        if (!finalMissing[i]) {
                            ymarker.set(count++, 0, phenotypeData[i]);
                        }
                    }

                    //adjust the fixed effects
                    DoubleMatrix fixed2 = fixed.getSelection(keepIndex, null);

                    //add marker data to fixed effects
                    ArrayList<Object> markerIds = new ArrayList<Object>();
                    boolean markerIsDiscrete = theAdapter.isMarkerDiscrete(m);
                    int nAlleles;
                    int markerdf;
                    DoubleMatrix X;
                    int[] alleleCounts = null;

                    if (markerIsDiscrete) {
                        Object[] finalMarkerData = new Object[numberOfRowsKept];
                        count = 0;
                        int n = markerData.length;
                        for (int i = 0; i < n; i++) {
                            if (!finalMissing[i]) {
                                finalMarkerData[count++] = markerData[i];
                            }
                        }
                        FactorModelEffect markerEffect = new FactorModelEffect(ModelEffectUtils.getIntegerLevels(finalMarkerData, markerIds), true);
                        X = fixed2.concatenate(markerEffect.getX(), false);
                        nAlleles = markerEffect.getNumberOfLevels();
                        alleleCounts = markerEffect.getLevelCounts();
                        markerdf = nAlleles - 1;
                    } else {
                        double[] finalMarkerValues = new double[numberOfRowsKept];
                        count = 0;
                        int n = markerData.length;
                        for (int i = 0; i < n; i++) {
                            if (!finalMissing[i]) {
                                finalMarkerValues[count++] = (Double) markerData[i];
                            }
                        }
                        X = fixed2.concatenate(DoubleMatrixFactory.DEFAULT.make(numberOfRowsKept, 1, finalMarkerValues), false);
                        nAlleles = 1;
                        alleleCounts = new int[]{numberOfRowsKept};
                        markerdf = 1;
                    }

                    CompressedMLMResult result = new CompressedMLMResult();
                    //need to add marker information to result once Alignment is stable

                    if (useP3D) {
                        testMarkerUsingP3D(result, ymarker, X, Vminus.getInverse(additionalMissing), markerdf);
                    } else {
                    	DoubleMatrix Zsel = zk[0].getSelection(keepIndex, null);
                        testMarkerUsingEMMA(result, ymarker, X, zk[1], Zsel, nAlleles);
                        markerdf = result.modeldf - baseModeldf;
                    }

                    //if the results are to be filtered on pmax check for that condition
                    boolean recordTheseResults = true;
                    if (parentPlugin.isFilterOutput() && result.p > parentPlugin.getMaxp()) {
                        recordTheseResults = false;
                    }

                    if (recordTheseResults) {
                        //add result to main
                        //{"Trait","Marker","Chr","Pos","Locus","Site","df","F","p","errordf","MarkerR2","Genetic Var","Residual Var", "-2LnLikelihood"};
                        String markername = theAdapter.getMarkerName(m);
                        String chr = "";
                        String pos = "";
                        String locus = theAdapter.getLocusName(m);
                        String site = Integer.toString(theAdapter.getLocusPosition(m));
                        if (myGeneticMap != null) {
                            int ndx = myGeneticMap.getMarkerIndex(markername);
                            if (ndx >= 0) {
                                chr = myGeneticMap.getChromosome(ndx);
                                pos = Double.toString(myGeneticMap.getGeneticPosition(ndx));
                            }
                        }
                        double errordf = (double) (ymarker.numberOfRows() - result.modeldf);

                        if (myGeneticMap != null) {
                            tableRow = new Object[]{theAdapter.getPhenotypeName(ph),
                                        markername,
                                        chr,
                                        pos,
                                        locus,
                                        site,
                                        new Integer(markerdf),
                                        new Double(result.F),
                                        new Double(result.p),
                                        new Double(errordf),
                                        new Double(result.r2),
                                        new Double(genvar),
                                        new Double(resvar),
                                        new Double(-2 * lnlk)};
                        } else {
                            tableRow = new Object[]{theAdapter.getPhenotypeName(ph),
                                        markername,
                                        locus,
                                        site,
                                        new Integer(markerdf),
                                        new Double(result.F),
                                        new Double(result.p),
                                        new Double(errordf),
                                        new Double(result.r2),
                                        new Double(genvar),
                                        new Double(resvar),
                                        new Double(-2 * lnlk)};
                        }
                        if (parentPlugin.isWriteOutputToFile()) {
                            if (!writeObjectToFile(mainOutputWriter, tableRow)) {
                                return null;
                            }
                        } else {
                            resultsMain.add(tableRow);
                        }

                        //add result to alleles
                        //"Trait","Marker","Chr","Pos","Locus","Site","Allele","Effect"

                        if (!markerIsDiscrete) {
                            if (myGeneticMap != null) {
                                tableRow = new Object[]{theAdapter.getPhenotypeName(ph),
                                            markername,
                                            chr,
                                            pos,
                                            locus,
                                            site,
                                            "",
                                            result.beta.get(result.beta.numberOfRows() - 1, 0),
                                            numberOfRowsKept
                                        };
                            } else {
                                tableRow = new Object[]{theAdapter.getPhenotypeName(ph),
                                            markername,
                                            locus,
                                            site,
                                            "",
                                            result.beta.get(result.beta.numberOfRows() - 1, 0),
                                            numberOfRowsKept
                                        };
                            }
                            //record the results
                            if (parentPlugin.isWriteOutputToFile()) {
                                if (!writeObjectToFile(alleleOutputWriter, tableRow)) {
                                    return null;
                                }
                            } else {
                                resultsAlleles.add(tableRow);
                            }

                        } else if (nAlleles > 1) {
                            for (int a = 0; a < nAlleles; a++) {
                                Double estimate;
                                if (a < nAlleles - 1) {
                                    estimate = result.beta.get(result.beta.numberOfRows() - nAlleles + 1 + a, 0);
                                } else {
                                    estimate = 0.0;
                                }
                                if (myGeneticMap != null) {
                                    tableRow = new Object[]{theAdapter.getPhenotypeName(ph),
                                                markername,
                                                chr,
                                                pos,
                                                locus,
                                                site,
                                                markerIds.get(a),
                                                estimate,
                                                alleleCounts[a]
                                            };
                                } else {
                                    tableRow = new Object[]{theAdapter.getPhenotypeName(ph),
                                                markername,
                                                locus,
                                                site,
                                                markerIds.get(a),
                                                estimate,
                                                alleleCounts[a]
                                            };
                                }
                                //record the results
                                if (parentPlugin.isWriteOutputToFile()) {
                                    if (!writeObjectToFile(alleleOutputWriter, tableRow)) {
                                        return null;
                                    }
                                } else {
                                    resultsAlleles.add(tableRow);
                                }
                            }
                        }

                    }
                    iterationsSofar++;
                    int progress = (int) ((double) iterationsSofar / (double) expectedIterations * 100);
                    parentPlugin.updateProgress(progress);
                }
            }

        }

        parentPlugin.updateProgress(0);
        if (parentPlugin.isWriteOutputToFile()) {
            try {
                mainOutputWriter.close();
                alleleOutputWriter.close();
            } catch (IOException e) {
            }
        }

        if (parentPlugin.isWriteOutputToFile()) {
            if (useCompression) {
                writeCompressionResultsToFile();
            }
            return null;
        }

        return formatResults();
    }

    private String getTabbedStringFromArray(Object[] array) {
        StringBuffer sb = new StringBuffer();
        sb.append(array[0]);
        int n = array.length;
        for (int i = 1; i < n; i++) {
            sb.append("\t").append(array[i]);
        }
        return sb.toString();
    }

    private boolean writeObjectToFile(BufferedWriter bw, Object[] row) {
        try {
            bw.write(getTabbedStringFromArray(row));
            bw.newLine();
        } catch (IOException e) {
            String msg = "Unable to write main output to file in MLM. Application ending.";
            myLogger.error(msg);
            e.printStackTrace();
            if (parentPlugin.isInteractive()) {
                JOptionPane.showInternalMessageDialog(parentPlugin.getParentFrame(), msg, "IO Error", JOptionPane.ERROR_MESSAGE);
            }
            return false;
        }
        return true;
    }

    public List<Datum> formatResults() {
        LinkedList<Datum> output = new LinkedList<Datum>();

        //generate comments
        StringBuilder options = new StringBuilder();
        options.append("Compression = ").append(useCompression);
        options.append("Compression = ").append(useCompression);
        if (useCompression) {
            options.append(", compression level = ").append(compression).append("\n");
        }
        if (useP3D) {
            options.append("P3D = ").append(useP3D).append(". Variance components were estimated only for the model without any markers.\n");
        } else {
            options.append("P3D = ").append(useP3D).append(". Variance components were estimated for each marker.\n");
        }

        StringBuilder model = new StringBuilder();
        model.append("Model: trait = mean + ");
        int nFactors = theAdapter.getNumberOfFactors();
        for (int f = 0; f < nFactors; f++) {
            model.append(theAdapter.getFactorName(f)).append(" + ");
        }
        int nCovar = theAdapter.getNumberOfCovariates();
        for (int c = 0; c < nCovar; c++) {
            model.append(theAdapter.getCovariateName(c)).append(" + ");
        }
        model.append("marker\n");

        if (resultsMain.size() > 0) {
            Object[][] theTable = new Object[resultsMain.size()][];
            resultsMain.toArray(theTable);
            String reportName = "MLM_statistics_for_" + datasetName;
            StringBuilder comment = new StringBuilder();
            comment.append("MLM statistics for compressed MLM\n");
            comment.append("Dataset: ").append(datasetName).append("\n");
            comment.append(options).append(model);
            SimpleTableReport str = new SimpleTableReport(reportName, headerMain, theTable);
            output.add(new Datum(reportName, str, comment.toString()));
        }

        if (resultsAlleles.size() > 0) {
            Object[][] theTable = new Object[resultsAlleles.size()][];
            resultsAlleles.toArray(theTable);
            String reportName = "MLM_effects_for_" + datasetName;
            StringBuilder comment = new StringBuilder();
            comment.append("MLM SNP effect estimates\n");
            comment.append("Dataset: ").append(datasetName).append("\n");
            comment.append(options).append(model);
            SimpleTableReport str = new SimpleTableReport(reportName, headerAlleles, theTable);
            output.add(new Datum(reportName, str, comment.toString()));
        }

        if (resultsBlups.size() > 0) {
            Object[][] theTable = new Object[resultsBlups.size()][];
            resultsBlups.toArray(theTable);
            String reportName = "MLM_BLUPs_for_" + datasetName;
            StringBuilder comment = new StringBuilder();
            comment.append("MLM BLUPs\n");
            comment.append("Dataset: ").append(datasetName).append("\n");
            comment.append(options).append(model);
            SimpleTableReport str = new SimpleTableReport(reportName, headerBlups, theTable);
            output.add(new Datum(reportName, str, comment.toString()));
        }

        if (compressionResults.size() > 0) {
            Object[][] theTable = new Object[compressionResults.size()][];
            compressionResults.toArray(theTable);
            String reportName = "MLM_compression_for_" + datasetName;
            StringBuilder comment = new StringBuilder();
            comment.append("MLM compression report\n");
            comment.append("Dataset: ").append(datasetName).append("\n");
            comment.append(options).append(model);
            SimpleTableReport str = new SimpleTableReport(reportName, headerCompression, theTable);
            output.add(new Datum(reportName, str, comment.toString()));
        }

        return output;
    }

    public boolean writeCompressionResultsToFile() {
        try {
            String outputbase = parentPlugin.getOutputName();
            if (outputbase.endsWith(".txt")) {
                int index = outputbase.lastIndexOf(".");
                outputbase = outputbase.substring(0, index);
            }
            String datasetNameNoSpace = datasetName.trim().replaceAll("\\ ", "_");
            BufferedWriter bw = new BufferedWriter(new FileWriter(outputbase + "_" + datasetNameNoSpace + "_compression.txt"));
            bw.write(getTabbedStringFromArray(headerCompression));
            bw.newLine();
            for (Object[] row : compressionResults) {
                bw.write(getTabbedStringFromArray(row));
                bw.newLine();
            }
            bw.close();
        } catch (IOException e) {
            String msg = "Unable to write compression results to file in MLM.";
            myLogger.error(msg);
            e.printStackTrace();
            if (parentPlugin.isInteractive()) {
                JOptionPane.showInternalMessageDialog(parentPlugin.getParentFrame(), msg, "IO Error", JOptionPane.ERROR_MESSAGE);
            }
            return false;
        }
        return true;
    }

    /**
     * Computes ZKZ. If compression is specified then the compressed ZKZ is calculated along with compressed versions of Z and K.
     * @param data	the phenotype data. Needed for optimizing compression.
     * @param X	the incidence matrix specifying all fixed effects other than markers
     * @param Z the kinship incidence matrix
     * @param kin	the genetic similarity matrix
     * @param traitname	the name of the phenotype passed in data
     * @return	an array containing the Z matrix as its first element and the K matrix as its second element. If compression is specified, then both are the compressed versions.
     */
    public DoubleMatrix[] computeZKZ(DoubleMatrix data, DoubleMatrix X, DoubleMatrix Z, DistanceMatrix kin, String traitname) {
    	DoubleMatrix[] zkMatrices = new DoubleMatrix[2]; 
    	CompressedDoubleMatrix.kinshipMethod kinmethod = CompressedDoubleMatrix.kinshipMethod.avg;
    	
        //Kmatrix
        int nkin = kin.getSize();
        int nrow = nkin;
        int ncol = nrow;

        DoubleMatrix K = DoubleMatrixFactory.DEFAULT.make(nrow, ncol);
        for (int r = 0; r < nrow; r++) {
            for (int c = 0; c < ncol; c++) {
                K.set(r, c, kin.getDistance(r, c));
            }
        }

        if (!useCompression) {
        	zkMatrices[0] = Z;
        	zkMatrices[1] = K;
        } else if (Double.isNaN(compression)) {
            //are taxa replicated?
            //sum columns of Z. If any sum > 1, then yes
            int n = Z.numberOfColumns();
            int count = 0;
            boolean taxaReplicated = false;
            while (count < n && !taxaReplicated) {
                if (Z.columnSum(count++) > 1.5) {
                    taxaReplicated = true;
                }
            }

            DistanceMatrix distance = calculateDistanceFromKin(kin);
            CompressedDoubleMatrix cm = new CompressedDoubleMatrix(kin, new UPGMATree(distance));
            EMMAforDoubleMatrix emlm = new EMMAforDoubleMatrix(data, X, K, Z, 0, Double.NaN);
            
            emlm.solve();
            double bestlnlk = emlm.getLnLikelihood();
            int bestCompression = nkin;

            double exponent = 1;
            double base = 0.98;
            double maxexponent = Math.log(1 / ((double) nkin)) / Math.log(base);
            parentPlugin.updateProgress((int) (exponent * 100 / maxexponent));
            //int g = (int) (nkin * Math.pow(base, exponent));
            int g = (int) (nkin);
            while (g > 1 || (g == 1 && taxaReplicated)) {
                cm.setNumberOfGroups(g);

                DoubleMatrix compressedZ = cm.getCompressedZ(Z);
                DoubleMatrix compressedK = cm.getCompressedMatrix(kinmethod);
                try {
                    emlm = new EMMAforDoubleMatrix(data, X, compressedK, compressedZ, 0, Double.NaN);
                    emlm.solve();

                    //output number of groups, compression level (= number of taxa / number of groups), -2L, genvar, resvar
                    compressionResults.add(new Object[]{traitname, g,
                                ((double) nkin) / ((double) g),
                                -2 * emlm.getLnLikelihood(),
                                emlm.getVarRan(),
                                emlm.getVarRes()});

                    if (Double.isNaN(bestlnlk) || emlm.getLnLikelihood() > bestlnlk) {
                        bestlnlk = emlm.getLnLikelihood();
                        bestCompression = g;
                        resvar = emlm.getVarRes();
                        genvar = emlm.getVarRan();
                    }
                } catch (Exception e) {
                    System.out.println("Compression failed for g = " + g);
                }

                int prev = g;
                while (g == prev) {
                    exponent++;
                    int prog = (int) (exponent * 100 / maxexponent);
                    prog = Math.min(prog, 100);
                    parentPlugin.updateProgress(prog);
                    g = (int) (nkin * Math.pow(base, exponent));
                }
            }

            //for g = 1 use GLM to estimate beta and errvar
            if (!taxaReplicated) {
                SweepFast sweep = new SweepFast(X, data);
                sweep.XTXSweepSetDmin();
                n = X.numberOfColumns();
                double ssres = sweep.getResidualSS();
                double errordf = (double) (data.numberOfRows() - n);
                double errvar = ssres / errordf;
                double lnlk = (errordf * Math.log(2 * Math.PI * errvar) + errordf);

                compressionResults.add(new Object[]{traitname, g,
                            ((double) nkin) / ((double) g),
                            lnlk,
                            new Double(0.0),
                            errvar});

                if (Double.isNaN(bestlnlk) || emlm.getLnLikelihood() > bestlnlk) {
                    bestlnlk = emlm.getLnLikelihood();
                    bestCompression = g;
                    resvar = emlm.getVarRes();
                    genvar = 0;
                }

            }

            cm.setNumberOfGroups(bestCompression);
            zkMatrices[0] = cm.getCompressedZ(Z);
            zkMatrices[1] = cm.getCompressedMatrix(kinmethod);
            parentPlugin.updateProgress(0);

        } else {
            DistanceMatrix distance = calculateDistanceFromKin(kin);
            CompressedDoubleMatrix cm = new CompressedDoubleMatrix(kin, new UPGMATree(distance));
            int g = (int) Math.round(nkin / compression);
            cm.setNumberOfGroups(g);
            zkMatrices[0] = cm.getCompressedZ(Z);
            zkMatrices[1] = cm.getCompressedMatrix(kinmethod);
        }
        
        return zkMatrices;
    }

    public void testMarkerUsingEMMA(CompressedMLMResult result, DoubleMatrix y, DoubleMatrix X, DoubleMatrix K, DoubleMatrix Z, int nAlleles) {
        EMMAforDoubleMatrix emlm = new EMMAforDoubleMatrix(y, X, K, Z, nAlleles, Double.NaN);
        emlm.solve();
        result.beta = emlm.getBeta();
        double[] Fp = emlm.getMarkerFp();
        result.F = Fp[0];
        result.p = Fp[1];
        result.modeldf = emlm.getDfModel();
        genvar = emlm.getVarRan();
        resvar = emlm.getVarRes();
        lnlk = emlm.getLnLikelihood();

        calculateRsquare(X, y, emlm.getInvH(), result, nAlleles - 1);
    }

    public void testMarkerUsingP3D(CompressedMLMResult result, DoubleMatrix y, DoubleMatrix X, DoubleMatrix invV, int markerdf) {
        //calculate beta
        DoubleMatrix invXVX = X.crossproduct(invV).mult(X);
        invXVX.invert();
        result.beta = invXVX.mult(X.crossproduct(invV.mult(y)));

        //test for markerdf = 0
        if (markerdf == 0) {
            result.F = Double.NaN;
            result.p = Double.NaN;
            result.r2 = 0.0;
        } else {
            //calculate F test, p-value of F test
            int nparm = result.beta.numberOfRows();
            DoubleMatrix M = DoubleMatrixFactory.DEFAULT.make(markerdf, nparm, 0);
            for (int i = 0; i < markerdf; i++) {
                M.set(i, nparm - markerdf + i, 1);
            }
            DoubleMatrix Mb = M.mult(result.beta);
            DoubleMatrix invMiM = M.mult(invXVX.tcrossproduct(M));
            try {
                invMiM.invert();
                result.F = Mb.crossproduct(invMiM.mult(Mb)).get(0, 0) / markerdf;
            } catch (Exception ex) {
                result.F = Double.NaN;
            }
            try {
                result.p = LinearModelUtils.Ftest(result.F, markerdf, y.numberOfRows() - nparm);
            } catch (Exception e) {
                result.p = Double.NaN;
            }

            calculateRsquare(X, y, invV, result, markerdf);
        }

    }

    private void calculateRsquare(DoubleMatrix X, DoubleMatrix y, DoubleMatrix invV, CompressedMLMResult result, int markerdf) {
        //calculate R2
        //from Buse(1973) Am. Stat. 27:106-108.
        //R^2 = ymarker'*inverseV*ymarker / (y-mean)'*inverseV*(y-mean)
        //where ymarker = yhat(full model) - yhat(model without marker)
        //as Xm*betam where Xm is the columns of X due to the marker and betam is the portion of beta due to the markers adjusted for the marker mean

        int dimX = X.numberOfColumns();
        int dimXreduced = dimX - markerdf;
        int[] colsToKeep = new int[dimXreduced];
        for (int i = 0; i < dimXreduced; i++) {
            colsToKeep[i] = i;
        }
        DoubleMatrix Xreduced = X.getSelection(null, colsToKeep);

        //calculate reduced beta
        DoubleMatrix invXVX = Xreduced.crossproduct(invV).mult(Xreduced);
        invXVX.invert();
        DoubleMatrix betaReduced = invXVX.mult(Xreduced.crossproduct(invV.mult(y)));

        //calculate yhat = yhatFull - yhatReduced
        DoubleMatrix yhat = X.mult(result.beta);
        DoubleMatrix yhatReduced = Xreduced.mult(betaReduced);
        yhat.minusEquals(yhatReduced);

        //calculate ydev = y - mean
        double sum = 0;
        int n = y.numberOfRows();
        for (int i = 0; i < n; i++) {
            sum += y.get(i, 0);
        }
        double mean = sum / n;

        DoubleMatrix ydev = y.scalarAdd(-mean);
        double numerator = yhat.crossproduct(invV).mult(yhat).get(0, 0);
        double denominator = ydev.crossproduct(invV).mult(ydev).get(0, 0);
        result.r2 = numerator / denominator;
    }

    public DoubleMatrix calculateV(DoubleMatrix ZKZ, double genvar, double resvar) {
        DoubleMatrix V = ZKZ.scalarMult(genvar);
        int n = V.numberOfRows();
        for (int i = 0; i < n; i++) {
            V.set(i, i, V.get(i, i) + resvar);
        }
        return V;
    }

    public DistanceMatrix calculateDistanceFromKin(DistanceMatrix kin) {
        int n = kin.getSize();
        double max = kin.getDistance(0, 0);
        for (int i = 0; i < n; i++) {
            max = Math.max(max, kin.getDistance(i, i));
        }

        double constant;
        if (max > 2) {
            constant = max;
        } else if (max > 1) {
            constant = 2;
        } else {
            constant = 1;
        }

        DistanceMatrix distanceMatrix = new DistanceMatrix(kin);
        for (int r = 0; r < n; r++) {
            distanceMatrix.setDistance(r, r, constant - kin.getDistance(r, r));
            for (int c = r + 1; c < n; c++) {
                double newval = constant - kin.getDistance(r, c);
                distanceMatrix.setDistance(r, c, newval);
                distanceMatrix.setDistance(c, r, newval);
            }
        }
        return distanceMatrix;
    }

    /**
     * @param missing a boolean[] equal true when a value is missing in that row
     * @param phenotypeTaxa the taxa
     * @return an IdGroup with the taxa that are in both the kinship matrix and the phenotype
     */
    public IdGroup updateMissingWithKinship(boolean[] missing, Identifier[] phenotypeTaxa) {
        int n = missing.length;
        TreeSet<Identifier> kinSet = new TreeSet<Identifier>();
        for (int i = 0; i < n; i++) {
            int col = kinshipMatrix.whichIdNumber(phenotypeTaxa[i]);
            if (!missing[i]) {
                if (col == -1) {
                    missing[i] = true;
                } else {
                    kinSet.add(phenotypeTaxa[i]);
                }
            }
        }

        Identifier[] taxa = new Identifier[kinSet.size()];
        kinSet.toArray(taxa);
        return new SimpleIdGroup(taxa);
    }

    class CompressedMLMResult {

        DoubleMatrix beta = null;
        double F = Double.NaN;
        double p = Double.NaN;
        double r2 = Double.NaN;
        int modeldf;
        int ngroups;
    }

    public void setTestMarkers(boolean testMarkers) {
        this.testMarkers = testMarkers;
    }

    public void setMyGeneticMap(GeneticMap myGeneticMap) {
        this.myGeneticMap = myGeneticMap;
    }
}
