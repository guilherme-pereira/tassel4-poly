// NodeFactory.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)


package net.maizegenetics.pal.tree;


/**
 * Creates nodes
 * <b>
 * The purpose of this class is to decouple the creation of
 * a class of type "Node" from its actual implementation.  This
 * class should be used instead of calling the constructor
 * of an implementation of "Node"
 * (at the moment "SimpleNode") as it may change in the future.</b><p>
 *
 * Other plans: add features here to recyle old nodes rather than
 * leaving them to the Java garbage collector
 *
 * @author Korbinian Strimmer
 */
import net.maizegenetics.pal.ids.Identifier;
public class NodeFactory
{
	/** create a node */
	public static final Node createNode()
	{
		return new SimpleNode();
	}
	/** create a node, with a specified identifier */
	public static final Node createNode(Identifier id)
	{
		return new SimpleNode(id.getName(),0);
	}
	/** create a node, with a specified identifier */
	public static final Node createNode(Identifier id, double height)
	{
		SimpleNode sn = new SimpleNode(id.getName(),0);
		sn.setNodeHeight(height);
		return sn;
	}
		/** create a node, with a specified identifier */
	public static final Node createNodeBranchLength(double branchLength, Identifier id)
	{
		SimpleNode sn = new SimpleNode(id.getName(),0);
		sn.setBranchLength(branchLength);
		return sn;
	}
	/** constructor used to clone a node and all children */
	public static final Node createNode(Node node)
	{
		return new SimpleNode(node);
	}
	public static final Node createNode(Node[] children) {
		return new SimpleNode(children);
	}
	/**
	 * Create a node with the specified children, and the specified branch height
	 */
	public static final Node createNode(Node[] children, double height) {
		SimpleNode sn = new SimpleNode(children);
		sn.setNodeHeight(height);
		return sn;
	}
	/**
	 * Create a node with the specified children, and the specified branch length
	 */
	public static final Node createNodeBranchLength(double branchLength, Node[] children) {
		SimpleNode sn = new SimpleNode(children);
		sn.setBranchLength(branchLength);
		return sn;
	}
}
