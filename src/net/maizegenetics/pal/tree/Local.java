package net.maizegenetics.pal.tree;

import net.maizegenetics.pal.io.FormattedOutput;
import net.maizegenetics.pal.math.MersenneTwisterFast;

/**
 * Implements LOCAL (Larget and Simon, 1999) and stochastic NNI moves for unrooted trees.
 * @author Alexei Drummond
 * @version $Id: Local.java,v 1.1 2007/01/12 03:26:17 tcasstevens Exp $
 */
public class Local {
	
	private static MersenneTwisterFast random = new MersenneTwisterFast();
	private static double lambda = 1.0;
	
	public static Tree local(Tree tree) {
		return local(tree, 1.0);
	}
	
	/**
	 * This method does a local interchange on
	 * the given tree: <b>Note</b> the original tree <b>is</b> modified. 
	 * If this behaviour is inappropriate then pass new SimpleTree(tree) 
	 * to this method.
	 * <pre>
	 * Actual		 Logical
	 *     2                 1            4
	 *    /|\                 \          /
	 *   / | \                 \        /
	 *  /  |  \                 \      /       
	 * 1   A   3                 2----3
	 *        / \               /      \
	 *       /   \             /        \
	 *      /     \           /          \
	 *     B       4         A            B
	 *
	 * </pre>
	 * A random internal edge (2,3) is selected and extended in both directions
	 * to create a back bone (1,2,3,4). One of the two internal nodes (2,3) is moved
	 * to a new random position on backbone and the backbone is scale in size.
	 * @param scale determines whether or not the backbone is scaled
	 * @return a perturbation of given tree. 
	 */
	public static Tree local(Tree tree, double scaleFactor) {
		
		if (tree.getRoot().getChildCount() != 3) {
			throw new RuntimeException("Root must have trifurcation!");
		}

		// (node1, node2, node3, node4) is the backbone
		
		//-------------------------------------------------------------
		// select an internal edge (i.e. one not connected to a tip) 
		// uniformly and randomly.
		//-------------------------------------------------------------
		
		// assumes root is last internal node and avoids it
		int pos = random.nextInt(tree.getInternalNodeCount()-1);
		Node node3 = tree.getInternalNode(pos);
		Node node2 = node3.getParent();
		
		//-------------------------------------------------------------
		
		// reroot so that top of edge is root
		TreeUtils.reroot(tree, node2);
		
		int k = random.nextInt(node2.getChildCount());
		
		//System.out.println("getting node1...");
		while (node2.getChild(k) == node3) {
			k = random.nextInt(node2.getChildCount());
		}
		Node node1 = node2.getChild(k);
		Node nodeA = null;
		for (int i =0; i < node2.getChildCount(); i++) {
			if ((node2.getChild(i) != node1) && (node2.getChild(i) != node3)) {
				nodeA = node2.getChild(i);
			}
		}
		
		//System.out.println("getting node4...");
		Node node4, nodeB;
		
		node4 = node3.getChild(0);
		nodeB = node3.getChild(1);
		if (random.nextBoolean()) {
			nodeB = node3.getChild(0);
			node4 = node3.getChild(1);
		}
	
		double backBoneLength = node1.getBranchLength() + node3.getBranchLength() + node4.getBranchLength();
			
		// modify backbone length
		double newLength = backBoneLength * scaleFactor;
		node1.setBranchLength(node1.getBranchLength() * scaleFactor);
		node3.setBranchLength(node3.getBranchLength() * scaleFactor);
		node4.setBranchLength(node4.getBranchLength() * scaleFactor);
		double newpos = random.nextDouble() * newLength;
			
		if (random.nextBoolean()) {
			// detach and reattach A

			double easyLength = node1.getBranchLength() + node3.getBranchLength();
			if (newpos < easyLength) {
				//no topology change
				node1.setBranchLength(newpos);
				node3.setBranchLength(easyLength-newpos);
			} else {
				swapNodes(nodeA, nodeB);
				node1.setBranchLength(easyLength);
				node3.setBranchLength(newpos - easyLength);
				node4.setBranchLength(newLength - newpos);
			}
		} else {
			// detach and reattach B
			double easyLength = node3.getBranchLength() + node4.getBranchLength();
			double hardLength = node1.getBranchLength();
			if (newpos > hardLength) {
				// no topology change
				node3.setBranchLength(newpos - hardLength);
				node4.setBranchLength(newLength - newpos);
			} else {
				swapNodes(node1, node4);
				node1.setBranchLength(newpos);
				node3.setBranchLength(hardLength - newpos);
				node4.setBranchLength(easyLength);
			}
		}
	
		tree.createNodeList();
		NodeUtils.lengths2Heights(tree.getRoot());

		return tree;
	}

	public static Tree stochasticNNI(Tree tree) {
	
		if (tree.getRoot().getChildCount() != 3) {
			throw new RuntimeException("Root must have trifurcation!");
		}
		// (node2, node3) is the backbone
		
		//-------------------------------------------------------------
		// select an internal edge (i.e. one not connected to a tip) 
		// uniformly and randomly.
		//-------------------------------------------------------------
		
		// assumes root is last internal node and avoids it
		int pos = random.nextInt(tree.getInternalNodeCount()-1);
		Node node3 = tree.getInternalNode(pos);
		Node node2 = node3.getParent();
		
		//-------------------------------------------------------------
		
		// reroot so that top of edge is root
		TreeUtils.reroot(tree, node2);
		
		int k = random.nextInt(node2.getChildCount());
		
		while (node2.getChild(k) == node3) {
			k = random.nextInt(node2.getChildCount());
		}
		Node node1 = node2.getChild(k);
		
		Node node4 = node3.getChild(0);
		if (random.nextBoolean()) {
			node4 = node3.getChild(1);
		} 
	
		swapNodes(node1, node4);

		return tree;
	}

	private static void swapNodes(Node n1, Node n2) {
		Node parent1 = n1.getParent();
		Node parent2 = n2.getParent();

		for (int i = 0; i < parent1.getChildCount(); i++) {
			if (parent1.getChild(i) == n1) parent1.removeChild(i);
		}
		for (int i = 0; i < parent2.getChildCount(); i++) {
			if (parent2.getChild(i) == n2) parent2.removeChild(i);
		}
		
		parent1.addChild(n2);
		parent2.addChild(n1);
	}

	public static void print4TaxonTree(Tree tree, java.io.PrintWriter out) {
		
		FormattedOutput fo = FormattedOutput.getInstance();
		
		Node root = tree.getRoot();
		Node taxa1 = null, taxa2 = null, internal1 = null, taxa3 = null, taxa4 = null;
	
		for (int i =0; i < root.getChildCount(); i++) {
			if (root.getChild(i).isLeaf()) {
				if (taxa1 == null) taxa1 = root.getChild(i);
				else taxa2 = root.getChild(i);
			} else internal1 = root.getChild(i);
		}
		taxa3 = internal1.getChild(0);
		taxa4 = internal1.getChild(1);

		displayLabel(out, taxa1.getIdentifier().getName(), 8, true);
		out.print("               ");
		displayLabel(out, taxa3.getIdentifier().getName(), 8, true);
		out.println();
		
		out.println("    \\                /");
		
		out.print("  ");
		fo.displayDecimal(out, taxa1.getBranchLength(), 4);
		out.print("          ");
		fo.displayDecimal(out, taxa3.getBranchLength(), 4);
		out.println();

		out.println("      \\            /");
		displayLabel(out, root.getIdentifier().getName(), 8, false);
		out.print("--");
		fo.displayDecimal(out, internal1.getBranchLength(), 4);
		out.println("--" + internal1.getIdentifier().getName());
		
		out.println("      /            \\");

		out.print("  ");
		fo.displayDecimal(out, taxa2.getBranchLength(), 4);
		out.print("          ");
		fo.displayDecimal(out, taxa4.getBranchLength(), 4);
		out.println();

		out.println("    /                \\");
		
		displayLabel(out, taxa2.getIdentifier().getName(), 8, true);
		out.print("               ");
		displayLabel(out, taxa4.getIdentifier().getName(), 8, true);
		out.println();
	}

	/**
	 * print label with a prespecified length
	 * (label will be shortened or spaces will introduced, if necessary)
	 *
	 * @param out output stream
	 * @param label label to be printed
	 * @param width desired length
	 */
	public static void displayLabel(java.io.PrintWriter out, String label, int width, boolean center)
	{
		int len = label.length();

		if (len == width)
		{
			// Print as is
			out.print(label);
		}
		else if (len < width)
		{
			int first = width-len;
			int second = 0;
		
			if (center) {
				first = first / 2;
				second = first - (width-len);
			}
			// fill rest with spaces
			for (int i = 0; i < first; i++) {
				out.print(' ');
			}
			out.print(label);
			for (int i = 0; i < second; i++) {
				out.print(' ');
			}
		}
		else
		{
			// Print first width characters
			for (int i = 0; i < width; i++)
			{
				out.print(label.charAt(i));
			}			
		}		
	}

	public static final void main(String[] args) {
		
		// create test tree

		Node root = new SimpleNode("I", 0.0);
		Node node3 = new SimpleNode("I", 0.01);
		Node node1 = new SimpleNode("1", 0.01);
		Node nodeA = new SimpleNode("2", 0.01);
		Node node4 = new SimpleNode("3", 0.01);
		Node nodeB = new SimpleNode("4", 0.01);
		
		root.addChild(node1);
		root.addChild(nodeA);
		root.addChild(node3);
		node3.addChild(nodeB);
		node3.addChild(node4);
		
		Tree tree = new SimpleTree(root);
		Tree tree2 = new SimpleTree(tree);
		Tree tree3 = new SimpleTree(tree);
		
		java.io.PrintWriter pw = new java.io.PrintWriter(System.out);

		print4TaxonTree(tree, pw);
		pw.flush();
		System.out.println();
		
		System.out.println("scaled 0.5");
		print4TaxonTree(local(tree, 0.5), pw);
		pw.flush();
		System.out.println();
		System.out.println("scaled 2.0");
		print4TaxonTree(local(tree2, 2.0), pw);
		pw.flush();
		System.out.println();
		System.out.println("NNI");
		print4TaxonTree(stochasticNNI(tree3), pw);
		pw.flush();
		System.out.println();
		
	}	
}
