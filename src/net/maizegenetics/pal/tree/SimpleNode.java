// SimpleNode.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)


package net.maizegenetics.pal.tree;

import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.ids.LabelMapping;
import net.maizegenetics.pal.util.BranchLimits;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.Enumeration;
import java.util.Hashtable;


/**
 * data structure for a node (includes branch) in a binary/non-binary
 * rooted/unrooted tree
 *
 * @version $Id: SimpleNode.java,v 1.1 2007/01/12 03:26:17 tcasstevens Exp $
 *
 * @author Korbinian Strimmer
 * @author Alexei Drummond
 */
public class SimpleNode implements AttributeNode {

	/** parent node */
	private Node parent;

	/** number of node as displayed */
	private int number;

	/** sequences associated with node */
	private byte[] sequence;

	/** length of branch to parent node */
	private double length;

	/** standard error of length of branch to parent node */
	private double lengthSE;

	/** height of this node */
	private double height;

	/** identifier of node/associated branch */
	private Identifier identifier;

	/** the attributes associated with this node. */
	private Hashtable attributes = null;

	//
	// Private stuff
	//

	private Node[] child;

	//
	// Serialization Stuff
	//

	static final long serialVersionUID=3472432765038556717L;

	//serialver -classpath ./classes net.maizegenetics.pal.tree.SimpleNode

	/** I like doing things my self! */
	private void writeObject(java.io.ObjectOutputStream out) throws java.io.IOException {
		out.writeByte(3); //Version number
		out.writeObject(parent);
		out.writeInt(number);
		out.writeObject(sequence);
		out.writeDouble(length);
		out.writeDouble(lengthSE);
		out.writeDouble(height);
		out.writeObject(identifier);
		out.writeObject(child);
		out.writeObject(attributes);
	}

	 private void readObject(java.io.ObjectInputStream in) throws IOException, ClassNotFoundException{
			byte version = in.readByte();
			switch(version) {
					case 1: {
						parent = (Node)in.readObject();
						number = in.readInt();
						sequence = (byte[])in.readObject();
						Object partial = (double[][][])in.readObject();
						length = in.readDouble();
						lengthSE = in.readDouble();
						height = in.readDouble();
						in.readDouble();
						identifier = (Identifier)in.readObject();
						child = (Node[])in.readObject();
						break;
					}
					case 2 : {
						//attributes are transient

						parent = (Node)in.readObject();
						number = in.readInt();
						sequence = (byte[])in.readObject();
						length = in.readDouble();
						lengthSE = in.readDouble();
						height = in.readDouble();
						identifier = (Identifier)in.readObject();
						child = (Node[])in.readObject();
						break;
					}
					default : {
						//attributes are not transient!

						parent = (Node)in.readObject();
						number = in.readInt();
						sequence = (byte[])in.readObject();
						length = in.readDouble();
						lengthSE = in.readDouble();
						height = in.readDouble();
						identifier = (Identifier)in.readObject();
						child = (Node[])in.readObject();
						attributes = (Hashtable)in.readObject();
					}
			}
	 }

	// the following constructors should eventually become
	// "friendly" to prevent anyone calling them directly.
	// Instead, the NodeFactory should be used!

	/** constructor default node */
	public SimpleNode()
	{
		parent = null;
		child =null;
		length = 0.0;
		lengthSE = 0.0;
		height = 0.0;
		identifier = Identifier.ANONYMOUS;

		number = 0;
		sequence = null;
	}

	public SimpleNode(String name, double branchLength) {
		this();
		identifier = new Identifier(name);
		length = branchLength;

	}
	/**
	 * Constructor
	 * @param children
	 * @param branchLength
	 * @throws IllegalArgumentException if only one child!
	 */
	protected SimpleNode(Node[] children, double branchLength) {
		this();
		this.child = children;
		if(children.length==1) {
		  throw new IllegalArgumentException("Must have more than one child!");
		}
		for(int i = 0 ; i < child.length ; i++) {
			child[i].setParent(this);
		}
		this.length = branchLength;
	}
	protected SimpleNode(Node[] children) {
		this(children, BranchLimits.MINARC);
	}

	/** constructor used to clone a node and all children */
	public SimpleNode(Node n)
	{
		this(n, true);
	}



	public void reset()
	{
		parent = null;
		child = null;
		length = 0.0;
		lengthSE = 0.0;
		height = 0.0;
		identifier = Identifier.ANONYMOUS;

		number = 0;
		sequence = null;
	}

	public SimpleNode(Node n, boolean keepIds) {
		init(n, keepIds);
		for (int i = 0; i < n.getChildCount(); i++) {
			addChild(new SimpleNode(n.getChild(i), keepIds));
		}
	}

	public SimpleNode(Node n, LabelMapping lm) {
		init(n, true, lm);
		for (int i = 0; i < n.getChildCount(); i++) {
			addChild(new SimpleNode(n.getChild(i), lm));
		}
	}

	protected void init(Node n) {
		init(n, true);
	}
	/**
	 * Initialized node instance variables based on given Node.
	 * children are ignored.
	 */
	protected void init(Node n, boolean keepId) {
		init(n,keepId,null);
	}
	/**
	 * Initialized node instance variables based on given Node.
	 * children are ignored.
	 * @param lm - may be null
	 */
	protected void init(Node n, boolean keepId, LabelMapping lm) {
		parent = null;
		length = n.getBranchLength();
		lengthSE = n.getBranchLengthSE();
		height = n.getNodeHeight();
		if (keepId) {
			if(lm!=null) {
				identifier = lm.getLabelIdentifier(n.getIdentifier());
			} else {
				identifier = n.getIdentifier();
			}
		} else { identifier = Identifier.ANONYMOUS; }

		number = n.getNumber();
		sequence = n.getSequence();

		if (n instanceof AttributeNode) {
			AttributeNode attNode = (AttributeNode)n;
			Enumeration e = attNode.getAttributeNames();
			while ((e != null) && e.hasMoreElements()) {
				String name = (String)e.nextElement();
				setAttribute(name, attNode.getAttribute(name));
			}
		}

		child = null;
	}

	/**
	 * Returns the parent node of this node.
	 */
	public final Node getParent() {
		return parent;
	}

	/** Set the parent node of this node. */
	public void setParent(Node node)
	{
		parent = node;
	}

	/**
	 * removes parent.
	 */
	public final void removeParent() {
		parent = null;
	}

	/**
	 * Returns the sequence at this node, in the form of a String.
	 */
	public String getSequenceString() {
		return new String(sequence);
	}

	/**
	 * Returns the sequence at this node, in the form of an array of bytes.
	 */
	public byte[] getSequence() {
		return sequence;
	}

	/**
	 * Sets the sequence at this node, in the form of an array of bytes.
	 */
	public void setSequence(byte[] s) {
		sequence = s;
	}

	/**
	 * Get the length of the branch attaching this node to its parent.
	 */
	public final double getBranchLength() {
		return length;
	}

	/**
	 * Set the length of the branch attaching this node to its parent.
	 */
	public final void setBranchLength(double value) {
		length = value;
	}

	/**
	 * Get the length SE of the branch attaching this node to its parent.
	 */
	public final double getBranchLengthSE() {
		return lengthSE;
	}

	/**
	 * Set the length SE of the branch attaching this node to its parent.
	 */
	public final void setBranchLengthSE(double value) {
		lengthSE = value;
	}


	/**
	 * Get the height of this node relative to the most recent node.
	 */
	public final double getNodeHeight() {
		return height;
	}

	/**
	 * Set the height of this node relative to the most recent node.
	 * @note corrects children branch lengths
	 */
	public final void setNodeHeight(double value) {
		if(value<0) {
			height = -value;
		} else {
			height = value;
		}
		/*if(child!=null) {
			for(int i = 0 ; i<child.length ; i++) {
				child[i].setBranchLength(height-child[i].getNodeHeight());
			}
		}*/
	}
	/**
	 * Set the height of this node relative to the most recent node.
	 * @param adjustChildBranchLengths if true
	 */
	public final void setNodeHeight(double value,boolean adjustChildBranchLengths) {
		if(value<0) {
			height = -value;
		} else {
			height = value;
		}
		if(adjustChildBranchLengths &&child!=null) {
			for(int i = 0 ; i<child.length ; i++) {
				child[i].setBranchLength(height-child[i].getNodeHeight());
			}
		}
	}

	/**
	 * Returns the identifier for this node.
	 */
	public final Identifier getIdentifier() {
		return identifier;
	}

	/**
	 * Set identifier for this node.
	 */
	public final void setIdentifier(Identifier id) {
		identifier = id;
		//return identifier;
	}

	public void setNumber(int n) {
		number = n;
	}

	public int getNumber() {
		return number;
	}

	/**
	 * get child node
	 *
	 * @param n number of child
	 *
	 * @return child node
	 */
	public Node getChild(int n)
	{
		return child[n];
	}

	/**
	 * set child node
	 *
	 * @param n number
	 * @node node new child node
	 */
	public void setChild(int n, Node node)
	{
		child[n] = node;
		child[n].setParent(this);
	}

	/**
	 * check whether this node is an internal node
	 *
	 * @return result (true or false)
	 */
	public boolean hasChildren()
	{
		return !isLeaf();
	}

	/**
	 * check whether this node is an external node
	 *
	 * @return result (true or false)
	 */
	public boolean isLeaf()	{
		return (getChildCount() == 0);
	}

	/**
	 * check whether this node is a root node
	 *
	 * @return result (true or false)
	 */
	public boolean isRoot()
	{
		if (parent == null)
		{
			return true;
		}
		else
		{
			return false;
		}
	}


	/**
	 * add new child node
	 *
	 * @param n new child node
	 */
	public void addChild(Node n)
	{
		insertChild(n, getChildCount());
	}

	/**
	 * add new child node (insertion at a specific position)
	 *
	 * @param n new child node
	 + @param pos position
	 */
	public void insertChild(Node n, int pos)
	{
		int numChildren = getChildCount();

		Node[] newChild = new Node[numChildren + 1];

		for (int i = 0; i < pos; i++)
		{
			newChild[i] = child[i];
		}
		newChild[pos] = n;
		for (int i = pos; i < numChildren; i++)
		{
			newChild[i+1] = child[i];
		}

		child = newChild;

		n.setParent(this);
	}


	/**
	 * remove child
	 *
	 * @param n number of child to be removed
	 */
	public Node removeChild(int n)
	{
		int numChildren = getChildCount();

		if (n >= numChildren)
		{
			throw new IllegalArgumentException("Nonexistent child");
		}
		Node[] newChild = new Node[numChildren-1];

		for (int i = 0; i < n; i++)
		{
			newChild[i] = child[i];
		}

		for (int i = n; i < numChildren-1; i++)
		{
			newChild[i] = child[i+1];
		}

		Node removed = child[n];

		//remove parent link from removed child!
		removed.setParent(null);

		child = newChild;

		return removed;
	}

	/**
	 * determines the height of this node and its descendants
	 * from branch lengths, assuming contemporaneous tips.
	 */
	public void lengths2HeightsContemp()
	{
		double largestHeight = 0.0;

		if (!isLeaf())
		{
			for (int i = 0; i < getChildCount(); i++)
			{
				NodeUtils.lengths2Heights(getChild(i));

				double newHeight =
					getChild(i).getNodeHeight() + getChild(i).getBranchLength();

				if (newHeight > largestHeight)
				{
					largestHeight = newHeight;
				}
			}
		}

		setNodeHeight(largestHeight);
	}

	/**
	 * Sets a named attribute to the given value.
	 * @param name the name of the attribute
	 * @param value the value to set the attribute
	 */
	public final void setAttribute(String name, Object value) {
		if (attributes == null) attributes = new Hashtable();
		attributes.put(name, value);
	}

	/**
	 * @return the attribute with the given name or null if it doesn't exist.
	 * @param name the name of the attribute.
	 */
	public final Object getAttribute(String name) {
		if (attributes == null) return null;
		return attributes.get(name);
	}

	/**
	 * @return an enumeration of the attributes that this node has or null if the
	 * node has no attributes.
	 */
	public final Enumeration getAttributeNames() {
		if (attributes == null) return null;
		return attributes.keys();
	}

	/**
	 * Returns the number of children this node has.
	 */
	public final int getChildCount() {
		if (child == null) return 0;
		return child.length;
	}

	public String toString() {
		StringWriter sw = new StringWriter();
		NodeUtils.printNH(new PrintWriter(sw), this, true, false, 0, false);
		return sw.toString();
	}
}
