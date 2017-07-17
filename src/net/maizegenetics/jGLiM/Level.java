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
 * Date: Dec 29, 2004
 * Time: 9:33:56 AM
 */
public interface Level extends Comparable<Level> {
    int getNumberOfSublevels();
    Comparable getSublevel(int sublevel);
    Comparable[] getSublevels();
    boolean contains(Comparable sublevel);

}
