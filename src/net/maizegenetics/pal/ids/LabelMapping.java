// LableMapping.java
//
// (c) 1999-2001 PAL Development Core Team
//
// This package may be distributed under the
// terms of the Lesser GNU General Public License (LGPL)

package net.maizegenetics.pal.ids;

/**
 * Title:        LabelMapping
 * Description:  Allows for the substitution of one label for another
 * @author			 Matthew Goode
 * @version 1.0
 */
import java.util.Enumeration;
import java.util.Hashtable;
public class LabelMapping implements java.io.Serializable {
	Hashtable mappings_ = new Hashtable();

	//
	// Serialization code
	//
	private static final long serialVersionUID=-9217142171228146380L;

	//serialver -classpath ./classes net.maizegenetics.pal.ids.LabelMapping
	private void writeObject(java.io.ObjectOutputStream out) throws java.io.IOException {
		out.writeByte(1); //Version number
		out.writeObject(mappings_);
	}

	private void readObject(java.io.ObjectInputStream in) throws java.io.IOException, ClassNotFoundException{
		byte version = in.readByte();
		switch(version) {
			default : {
				mappings_ = (Hashtable)in.readObject();
				break;
			}
		}
	}
	private LabelMapping(Hashtable mapping) {
		this.mappings_ = mapping;
	}
	private LabelMapping(LabelMapping toCopy) {
		this.mappings_ = (Hashtable)toCopy.mappings_.clone();
	}
	public LabelMapping() { }

	public void addMapping(String id, String label) {
		mappings_.put(id,label);
	}
	public void addMapping(Identifier id, String label) {
		if(id!=null&&id.getName()!=null) {
			mappings_.put(id.getName(),label);
		}
	}
	/**
	 * @param names Names
	 * @param colours associated colours
	 * @note assumes parallel arrays
	 */
	public void addMappings(String[] ids, String[] labels) {
		for(int i = 0 ; i < ids.length ; i++) {
			mappings_.put(ids[i],labels[i]);
		}
	}

	public String getLabel(String id, String defaultLabel) {
		if(id==null||!mappings_.containsKey(id)) {
			return defaultLabel;
		}
		return mappings_.get(id).toString();
	}
	public String getLabel(Identifier id, String defaultLabel) {
		if(id==null) {
			return defaultLabel;
		}
		return getLabel(id.getName(),defaultLabel);
	}
	public String getLabel(Identifier id) {
		return getLabel(id.getName(),id.getName());
	}
	public Identifier getLabelIdentifier(Identifier id) {
		if(id==null) {
			return null;
		}
		return new Identifier(getLabel(id.getName(),id.getName()));
	}
	/**
	 * If a mapping occurs more than once will rename instance to "x 1", "x 2"... and so on where x is the mapping in question
	 */
	public LabelMapping getUniquifiedMappings() {
		Hashtable totals = new Hashtable();
		for(Enumeration e = mappings_.keys() ; e.hasMoreElements() ; ) {
			Object key = e.nextElement();
			Object mapping = mappings_.get(key);
			int count = 1;
			if(totals.containsKey(mapping)) {
				count = ((Integer)totals.get(mapping)).intValue()+1;
			}
			totals.put(mapping,new Integer(count));
		}
		Hashtable counts = new Hashtable();
		Hashtable result = new Hashtable();
		for(Enumeration e = mappings_.keys() ; e.hasMoreElements() ; ) {
			Object key = e.nextElement();
			Object mapping = mappings_.get(key);
			int total = ((Integer)totals.get(mapping)).intValue();
			if(total==1) {
				result.put(key,mapping);
			} else {
				int count = 1;
				if(counts.containsKey(mapping)) {
					count = ((Integer)counts.get(mapping)).intValue()+1;
				}
				counts.put(mapping,new Integer(count));
				result.put(key,mapping+" "+count);
			}
		}
		return new LabelMapping(result);
	}
	public LabelMapping getRelabeled(Relabeller relabeller) {
		Hashtable newMapping = new Hashtable();
		for(Enumeration e = mappings_.keys() ; e.hasMoreElements() ; ) {
			Object key = e.nextElement();
			String old = mappings_.get(key).toString();
			newMapping.put(key,relabeller.getNewLabel(old));
		}
		return new LabelMapping(newMapping);
	}
	public IdGroup getMapped(IdGroup original) {
		String[] oldIDs = IdGroupUtils.getNames(original);
		String[] newIDs = new String[oldIDs.length];
		for(int i = 0 ; i < newIDs.length ; i++) {
			newIDs[i] = getLabel(oldIDs[i],oldIDs[i]);
		}
		return new SimpleIdGroup(newIDs);
	}

	// Static classes

	public static interface Relabeller {
		public String getNewLabel(String oldLabel);
	}
}