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

/**
 * Created using IntelliJ IDEA.
 * Author: Peter Bradbury
 * Date: Jan 28, 2005
 * Time: 9:48:35 AM
 */
public class BasicStatistic implements Statistic {
    //fields
    private String name;
    private int type;
    private Object value;

    //constructors
    public BasicStatistic(Object value, int type, String name) {
        this.type = type;
        this.value = value;
        if (name == null) name = getStatisticNameFromType(type);
        else this.name = name;
    }

    public BasicStatistic(Object value, int type) {
        this.type = type;
        this.value = value;
        name = getStatisticNameFromType(type);
    }

    //methods
    public static String getStatisticNameFromType(int type) {

        switch(type) {
            case STATISTIC_SS_MODEL:
                return "SSModel";
            case STATISTIC_SS_ERROR:
                return "SSError";
            case STATISTIC_DF_MODEL:
                return "dfModel";
            case STATISTIC_DF_ERROR:
                return "dfError";
            case STATISTIC_MS_MODEL:
                return "MSModel";
            case STATISTIC_MS_ERROR:
                return "MSError";
            case STATISTIC_F:
                return "F";
            case STATISTIC_P:
                return "p";
            case STATISTIC_PERMUTE_P:
                return "p_perm";
            case STATISTIC_NULL_FDIST:
                return "Fdist";
            case STATISTIC_LSESTIMATES:
                return "LSMean";
            case STATISTIC_RSQ:
                return "Rsq";
            case STATISTIC_PERMUTATION_NUMBER:
                return "N_permutations";
            case STATISTIC_ANALYSIS_TYPE:
                return "Analysis Type";
            default:
                return "";
        }
    }

    public String getName() {
        return name;
    }

    public int getType() {
        return type;
    }

    public Object getValue() {
        return value;
    }

    public boolean equals(Object obj) {
        if (obj instanceof Statistic) {
            Statistic otherStatistic = (Statistic) obj;
            if (otherStatistic.getName() != name) return false;
            if (otherStatistic.getType() != type) return false;
            if (otherStatistic.getValue() != value) return false;
            return true;
        }
        return false;
    }

    public String toString() {
        return name + " = " + value.toString();
    }
}
