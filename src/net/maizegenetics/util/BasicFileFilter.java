package net.maizegenetics.util;

import javax.swing.filechooser.FileFilter;
import java.io.File;


/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2004</p>
 * <p>Company: </p>
 * @author not attributable
 * @version 1.0
 */


public class BasicFileFilter extends FileFilter
{
  private String extension;

  public BasicFileFilter( String extension ) { this.extension = extension; }

  public String getDescription() { return "." + extension + " files"; }

  public boolean accept(File f) {
    if (f.isDirectory())
      return true;

    String name = f.getName();
    int pos = name.lastIndexOf( '.' );
    if (pos == -1) return false;

    String the_extension = name.substring( pos+1 );

    if (the_extension.equals( extension )) return true;

    return false;
  }

  public String getExtension() {
    return extension;
  }
}

