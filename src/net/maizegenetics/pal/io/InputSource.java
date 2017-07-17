// InputSource.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)


package net.maizegenetics.pal.io;

import java.io.*;

import java.net.URL;


/**
 * convenience class to open input streams
 * linked to files, stdin, and strings
 *
 * @version $Id: InputSource.java,v 1.2 2009/05/03 06:00:10 tcasstevens Exp $
 *
 * @author Korbinian Strimmer
 */
public class InputSource extends PushbackReader
{
	//
	// Public stuff
	//

	/**
	 * open file for reading
	 *
	 * @param name file name
	 *
	 * @return input stream
	 */
	public static InputSource openFile(String name)
		throws FileNotFoundException
	{

        try {

        if (name.startsWith("http")) {
            URL url = new URL(name);
            return new InputSource(
			new BufferedReader(
			new InputStreamReader(url.openStream())));
        } else {
            return new InputSource(
			new BufferedReader(
			new FileReader(name)));
        }

        } catch (Exception e) {
            e.printStackTrace();
        }

        return null;
		
	}

	/**
	 * open standard input
	 *
	 * @return input stream
	 */			
	public static InputSource openStdIn()
	{
		return
			new InputSource(
			new BufferedReader(
			new InputStreamReader(System.in)));
	}

	/**
	 * "open" string for reading
	 *
	 * @param input string serving as source
	 *
	 * @return input stream
	 */
	public static InputSource openString(String input)
	{
		return new InputSource(new StringReader(input));
	}
	
	
	//
	// Private stuff
	//
	
	
	// Private constructor

	private InputSource(Reader in)
	{
		super(in);
	}
}
