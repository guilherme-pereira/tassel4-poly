// CladeSystem.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)

package net.maizegenetics.pal.tree;

import net.maizegenetics.pal.ids.IdGroup;

import java.io.PrintWriter;
import java.io.StringWriter;

/**
 * data structure for a set of splits 
 *
 * @version $Id: CladeSystem.java,v 1.1 2007/01/12 03:26:16 tcasstevens Exp $
 *
 * @author Alexei Drummond
 */
public class CladeSystem {
	//
	// Public stuff
	//

	/**
	 * @param idGroup  sequence labels
	 * @param size     number of clades
	 */
	public CladeSystem(IdGroup idGroup, int size) {
		this.idGroup = idGroup;
		clades = new boolean[size][idGroup.getIdCount()];
	}

	/** get number of clades */
	public int getCladeCount() {		
		return clades.length;
	}

	/** get number of labels */
	public int getLabelCount() {		
		return clades[0].length;
	}

	/** get clade array */
	public boolean[][] getCladeArray() {		
		return clades;
	}

	/** get clade */
	public boolean[] getClade(int i) {		
		return clades[i];
	}

	/** get idGroup */
	public IdGroup getIdGroup() {		
		return idGroup;
	}

	/**
	  + test whether a clade is contained in this clade system
	  * (assuming the same leaf order)
	  *
	  * @param clade clade
	  */
	public boolean hasClade(boolean[] clade) {
		for (int i = 0; i < clades.length; i++)
		{
			if (SplitUtils.isSame(clade, clades[i])) return true;
		}
			
		return false;
	}


	/** print clade system */
	public String toString() {
		StringWriter sw = new StringWriter();
		PrintWriter pw = new PrintWriter(sw);
		
		for (int i = 0; i < getLabelCount(); i++)
		{
			pw.println(idGroup.getIdentifier(i));
		}
		pw.println();
		
		for (int i = 0; i < getCladeCount(); i++)
		{
			for (int j = 0; j < getLabelCount(); j++)
			{
				if (clades[i][j] == true)
					pw.print('*');
				else
					pw.print('.');
			}
			
			pw.println();
		}

		return sw.toString();
	}
	
	// ********************************************************************
	// STATIC METHODS
	// ********************************************************************
	
	/**
	 * @return all clade systems for a group of trees.
	 */
	public static CladeSystem[] getCladeSystems(Tree[] trees) {
		IdGroup idGroup = TreeUtils.getLeafIdGroup(trees[0]);
		CladeSystem[] cladeSystems = new CladeSystem[trees.length];
		for (int i =0; i < cladeSystems.length; i++) {
			cladeSystems[i] = getClades(idGroup, trees[i]);
		}
		return cladeSystems;
	}

	public static void calculateCladeProbabilities(Tree tree, CladeSystem[] cladeSystems) {
	
		CladeSystem cladeSystem = getClades(cladeSystems[0].getIdGroup(), tree);
		
		for (int i =0; i < tree.getInternalNodeCount()-1; i++) {
			Node node = tree.getInternalNode(i);
			boolean[] clade = cladeSystem.getClade(i);
			if (node.isRoot()) throw new RuntimeException("Root node does not have clade probability!");
				
			int cladeCount = 0;
			for (int j = 0; j < cladeSystems.length; j++) {
				if (cladeSystems[j].hasClade(clade)) cladeCount += 1;
			}
			double pr = (double)cladeCount / (double)cladeSystems.length; 
			tree.setAttribute(node, AttributeNode.CLADE_PROBABILITY, new Double(pr));	
		}
	}
	
	/**
	 * creates a clade system from a tree
	 * (using a pre-specified order of sequences)
	 *
	 * @param idGroup  sequence order for the matrix
	 * @param tree
	 */
	public static CladeSystem getClades(IdGroup idGroup, Tree tree) {
		tree.createNodeList();
		
		int size = tree.getInternalNodeCount()-1;
		CladeSystem cladeSystem = new CladeSystem(idGroup, size);
		
		boolean[][] clades = cladeSystem.getCladeArray();
		
		for (int i = 0; i < size; i++) {
			getClade(idGroup, tree.getInternalNode(i), clades[i]);
		}
		
		return cladeSystem;
	}

	/**
	 * creates a clade system from a tree
	 * (using tree-induced order of sequences)
	 *
	 * @param tree
	 */
	public static CladeSystem getClades(Tree tree) {
		IdGroup idGroup = TreeUtils.getLeafIdGroup(tree);
		
		return getClades(idGroup, tree);
	}

	/**
	 * get clade for internal node
	 *
	 * @param idGroup order of labels
	 * @param internalNode Node
	 * @param boolean[] clade
	 */
	public static void getClade(IdGroup idGroup, Node internalNode, boolean[] clade) {
		if (internalNode.isLeaf() || internalNode.isRoot()) {
			throw new IllegalArgumentException("Only internal nodes (and no root) nodes allowed");
		}
		
		// make sure clade is reset
		for (int i = 0; i < clade.length; i++) {
			clade[i] = false;
		}
		
		// mark all leafs downstream of the node
		// AJD removed loop, as doesn't appear to be necessary
		SplitUtils.markNode(idGroup, internalNode, clade);
	}

	/**
	 * checks whether two clades are identical
	 * (assuming they are of the same length
	 * and use the same leaf order)
	 *
	 * @param s1 clade 1
	 * @param s2 clade 2
	 */
	public static boolean isSame(boolean[] s1, boolean[] s2)
	{
		if (s1.length != s2.length) 
			throw new IllegalArgumentException("Clades must be of the same length!");
		
		for (int i = 0; i < s1.length; i++) {
			// clades not identical
			if (s1[i] != s2[i]) return false;
		}	
		return true;
	}
	
	//
	// Private stuff
	//
	
	private IdGroup idGroup;
	private boolean[][] clades;
}
