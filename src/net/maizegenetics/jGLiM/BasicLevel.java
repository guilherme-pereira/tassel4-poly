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

import java.util.*;


/**
 * Created using IntelliJ IDEA.
 * Author: Peter Bradbury
 * Date: Dec 29, 2004
 * Time: 9:41:55 AM
 */
public class BasicLevel implements Level{
    //fields
    private Comparable[] levelValues;
    private int index;
    private int random;

    //constructors
    public BasicLevel(Comparable[] levelValues) {
        this.levelValues = levelValues;
    }

    public BasicLevel(Level[] levels) {
        int nlevels = 0;
        for (int i = 0; i < levels.length; i++) {
            nlevels = +getNumberOfSublevels();
        }
        levelValues = new Comparable[nlevels];
        int cnt = 0;
        for (int i = 0; i < levels.length; i++) {
            for (int j = 0; j < levels[i].getNumberOfSublevels(); j++) {
                levelValues[cnt++] = levels[i].getSublevel(j);
            }
        }
    }

    public BasicLevel(Collection c) {
        LinkedList<Comparable> levelList = new LinkedList<Comparable>();
        Iterator it = c.iterator();
        while (it.hasNext()) {
            Level level = (Level) it.next();
            for (int i = 0; i < level.getNumberOfSublevels(); i++) {
                levelList.add(level.getSublevel(i));
            }
        }
        levelValues = new Comparable[levelList.size()];
        levelList.toArray(levelValues);
    }

    //methods
    public int compareTo(Level otherLevel) {
        int result;
        for (int i = 0; i < levelValues.length; i++) {
            result = levelValues[i].compareTo(otherLevel.getSublevel(i));
            if (result != 0) return result;
        }
        return 0;
    }

    public boolean equals(Object o) {
        if (!(o instanceof Level)) return false;
        Level otherLevel = (Level) o;
        if (levelValues.length != otherLevel.getNumberOfSublevels()) return false;
        for (int i = 0; i < levelValues.length; i++) {
            if (!(levelValues[i].equals(otherLevel.getSublevel(i)))) return false;
        }
        return true;
    }

    public int hashCode() {
        int hashval = 0;
        for (int i = 0; i < levelValues.length; i++) {
            hashval += levelValues[i].hashCode();
        }
        return hashval;
    }

    public String toString() {
        StringBuffer out = new StringBuffer(levelValues[0].toString());

        for (int i = 1; i < levelValues.length; i++) {
            out.append(":");
            out.append(levelValues[i].toString());
        }
        return out.toString();
    }

    public int getNumberOfSublevels() {
        return levelValues.length;
    }

    public Comparable[] getSublevels() {
        return levelValues;
    }

    public Comparable getSublevel(int sublevel) {
        return levelValues[sublevel];
    }

    public boolean contains(Comparable sublevel) {
        for (int i = 0; i < levelValues.length; i++) {
            if (levelValues[i].equals(sublevel)) return true;
        }
        return false;
    }

    public int getIndex() {
        return index;
    }

    public void setIndex(int index) {
        this.index = index;
    }

    public int getRandom() {
        return random;
    }

    public void setRandom(int random) {
        this.random = random;
    }

    public static Comparator indexComparator() {
        return new Comparator() {
            public int compare(Object o1, Object o2) {
                BasicLevel basicLevel1 = (BasicLevel) o1;
                BasicLevel basicLevel2 = (BasicLevel) o2;
                int result = basicLevel1.compareTo(basicLevel2);
                if (result != 0) return result;
                return basicLevel1.getIndex() - basicLevel2.getIndex();
            }
        };
    }

    public static Comparator randomComparator() {
        return new Comparator() {
            public int compare(Object o1, Object o2) {
                BasicLevel basicLevel1 = (BasicLevel) o1;
                BasicLevel basicLevel2 = (BasicLevel) o2;
                int result = basicLevel1.compareTo(basicLevel2);
                if (result != 0) return result;
                if (basicLevel1.getRandom() > basicLevel2.getRandom()) return 1;
                if (basicLevel1.getRandom() < basicLevel2.getRandom()) return -1;
                return 0;
            }
        };
    }
}
