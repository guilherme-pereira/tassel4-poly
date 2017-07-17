package net.maizegenetics.jGLiM;
/*
 * jGLiM: Java for General Linear Models
 * for more information: http://www.maizegenetics.net
 *
 * Copyright (C) 2005 Peter Bradbury
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 */

//package net.maizegenetics.jGLiM;

import cern.colt.matrix.DoubleMatrix1D;

/**
 * Created using IntelliJ IDEA.
 * Author: Peter Bradbury
 * Date: Jan 27, 2005
 * Time: 4:41:15 PM
 */
public interface LinearModel {
    double getSSModel();
    double getSSError();
    int getdfModel();
    int getdfError();
    double getMSModel();
    double getMSError();
    double getRsq();
    double getF();
    double getp();
    int getN();
    DoubleMatrix1D getData();
    DoubleMatrix1D getBeta();
    int getAnalysisType();
    void setAnalysisType(int analysisType);
    Statistic getStatistic(int statisticType);
    int getNumberOfPermutations();
    void setNumberOfPermutations(int numberOfPermutations);
    double[] getNullFDistribution();
    double getPermutationp();
    void setModelForError(LinearModel lm);
    LinearModel getModelForError();
    ProgressIndicator getProgressIndicator();
    void setProgressIndicator(ProgressIndicator progress);
    
    void calculate();
    void recalculate(DoubleMatrix1D permutedData);
    void permute();
    void permute(int numberOfPermutations);
    void permute(int[] shuffleIndex);
     
    /**
     * Use when the dfError calculated by the model is not desired error for the F net.maizegenetics.jGLiM.test.
     *
     * @param dfError
     */
    void setdfError(int dfError);
    void setSSError(double SSError);
}
