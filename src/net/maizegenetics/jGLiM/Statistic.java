package net.maizegenetics.jGLiM;


public interface Statistic {
    int STATISTIC_SS_MODEL = 1;
    int STATISTIC_SS_ERROR = 2;
    int STATISTIC_DF_MODEL = 3;
    int STATISTIC_DF_ERROR = 4;
    int STATISTIC_MS_MODEL = 5;
    int STATISTIC_MS_ERROR = 6;
    int STATISTIC_F = 7;
    int STATISTIC_P = 8;
    int STATISTIC_PERMUTE_P = 9;
    int STATISTIC_NULL_FDIST = 10;
    int STATISTIC_LSESTIMATES = 11;
    int STATISTIC_RSQ = 12;
    int STATISTIC_PERMUTATION_NUMBER = 13;
    int STATISTIC_ANALYSIS_TYPE = 14;
    int STATISTIC_VARIANCE_ADD_EST = 15;
    int STATISTIC_VARIANCE_RES_EST = 16;
    int STATISTIC_VARIANCE_TOT_EST = 17;
    int STATISTIC_VARIANCE_ADD_STD = 18;
    int STATISTIC_VARIANCE_RES_STD = 19;
    int STATISTIC_VARIANCE_TOT_STD = 20;
    int STATISTIC_VARIANCE_ADD_Z = 21;
    int STATISTIC_VARIANCE_RES_Z = 22;
    int STATISTIC_VARIANCE_TOT_Z = 23;
    int STATISTIC_VARIANCE_ADD_P = 24;
    int STATISTIC_VARIANCE_RES_P = 25;
    int STATISTIC_VARIANCE_TOT_P = 26;
    int STATISTIC_H2_EST = 27;
    int STATISTIC_H2_STD = 28;
    int STATISTIC_H2_Z = 29;
    int STATISTIC_H2_P = 30;

    String getName();

    int getType();

    Object getValue();
}