// ReadTree.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)


package net.maizegenetics.pal.tree;

import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.io.FormattedInput;
import net.maizegenetics.pal.io.InputSource;

import java.io.IOException;
import java.io.PushbackReader;


/**
 * constructs a tree reading in a New Hampshire treefile, taking care
 * for internal labels and branch lengths and binary/nonbinary and
 * rooted/unrooted trees
 *
 * @version $Id: ReadTree.java,v 1.1 2007/01/12 03:26:17 tcasstevens Exp $
 *
 * @author Korbinian Strimmer
 * @author Alexei Drummond
 */
public class ReadTree extends SimpleTree
{
	//
	// Public stuff
	//

	/**
	 * read tree from input stream
	 *
	 * @param input input stream
	 */
	public ReadTree(PushbackReader input) throws TreeParseException
	{
		super();

		readNH(input);

		// node heights should be populated as well - AD
		NodeUtils.lengths2Heights(getRoot());

		createNodeList();
	}

	/**
	 * read tree from file
	 *
	 * @param file name of file
	 */
	public ReadTree(String file) throws TreeParseException, IOException
	{
		super();

		PushbackReader input = InputSource.openFile(file);
		readNH(input);
		input.close();

		// node heights should be populated as well - AD
		NodeUtils.lengths2Heights(getRoot());
		createNodeList();
	}


	//
	// Private stuff
	//

	private FormattedInput fi = FormattedInput.getInstance();

	// Construct tree by reading a New Hampshire tree
	private void readNH(PushbackReader input, Node currentNode)
		throws TreeParseException
	{
		try
		{
			int c = fi.readNextChar(input);
			if (c == '(')
			{
				int count = 0;
				do
				{
					Node newNode = NodeFactory.createNode();
					currentNode.addChild(newNode);
					readNH(input, newNode);
					count++;

					c = fi.readNextChar(input);
				}
				while (c == ',');

				if (c != ')')	{
					throw new TreeParseException("Missing closing bracket");
				}

				if (count < 2){
					throw new TreeParseException("Node with single child enountered");
				}

			}
			else
			{
				input.unread(c);
			}

			// Read label (any length)
			currentNode.setIdentifier(new Identifier(fi.readLabel(input, -1)));

			// Read distance
			c = fi.readNextChar(input);

			if (c == ':')
			{
				currentNode.setBranchLength(fi.readDouble(input,true));
			}
			else
			{
				input.unread(c);
			}
		}

		catch (IOException e) {
			throw new TreeParseException("IO error");
		}
		catch (NumberFormatException e) {
			throw new TreeParseException("Error while parsing number");
		}
	}

	// Construct tree by reading a New Hampshire tree
	private void readNH(PushbackReader input) throws TreeParseException
	{
		try
		{
			readNH(input, getRoot());

			// Drop terminating semicolon
			int c = fi.readNextChar(input);
			if (c != ';')
			{
				throw new TreeParseException("Missing terminating semicolon");
			}
		}

			catch (IOException e)
			{
				throw new TreeParseException();
			}
	}
}
