package net.maizegenetics.pal.alignment;

import java.io.Serializable;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;


import net.maizegenetics.pal.report.TableReport;

public class GeneticMap implements TableReport, Serializable {
        private static final long serialVersionUID = -5197800047652332969L;
	protected HashMap<String,Integer> markerIndex = new HashMap<String,Integer>();	//a hash map of marker indices for lookup functions
	int[] mapIndex;
	ArrayList<mapFeature> featureList = new ArrayList<mapFeature>();
	boolean hasPhysicalPositions = false;
	ArrayList<mapFeature> physicalFeatureList = new ArrayList<mapFeature>();
	String name;
	
	/**
	 * @param markerID	the marker names
	 * @param locus	the locus, usually a chromosome	
	 * @param geneticPosition the genetic position of the markers, required though values can be NaN
	 * @param physicalPosition the physical positions of the markers, can be null, missing values = -1
	 */
	public GeneticMap(String name, String[] markerID, String[] chromosome, double[] geneticPosition, int[] physicalPosition) {
		this.name = name;
		
		//read all markers with genetic positions in the map array
		int n = markerID.length;
		for (int i = 0; i < n; i++) {
			if (!Double.isNaN(geneticPosition[i])) {
				int pos = -1;
				if (physicalPosition != null) {
					pos = physicalPosition[i];
					hasPhysicalPositions = true;
				}
				mapFeature feature = new mapFeature(markerID[i], chromosome[i], pos, geneticPosition[i]);
				featureList.add(feature);
				if (!Double.isNaN(geneticPosition[i]) && pos > -1) physicalFeatureList.add(feature);
			}
		}
		
		Collections.sort(featureList);
		
		//create the hashmap
		int nFeatures = featureList.size();
		for (int i = 0; i < nFeatures; i++) markerIndex.put(featureList.get(i).id, i);
		
	}
	
	/**
	 * @param markerID	the marker names
	 * @param locus	the locus, usually a chromosome	
	 * @param geneticPosition the genetic position of the markers, required though values can be NaN
	 */
	public GeneticMap(String name, String[] markerID, String[] chromosome, double[] geneticPosition) {
		this(name, markerID, chromosome, geneticPosition, null);
	}
	
	/**
	 * Creates an empty genetic map
	 */
	public GeneticMap(String name) {
		this.name = name;
	}
	
	public void addMarker(String id, String chromosome, double geneticPosition, int physPosition) {
		mapFeature feature = new mapFeature(id, chromosome, physPosition, geneticPosition);
		featureList.add(feature);
		if (physPosition > -1) {
			if (!Double.isNaN(geneticPosition)) physicalFeatureList.add(feature);
			hasPhysicalPositions = true;
		}
		markerIndex.put(id, featureList.size() - 1);
	}

	public void addMarker(String id, String chromosome, double geneticPosition) {
		featureList.add(new mapFeature(id, chromosome, geneticPosition));
	}
	
	public boolean addMarker(String[] mapInfo) {
		String id;
		String chr;
		double genpos;
		int physpos;
		int infoLength = mapInfo.length;
		
		if (infoLength > 2) {
			id = mapInfo[0];
			chr = mapInfo[1];
			try {
				genpos = Double.parseDouble(mapInfo[2]);
			} catch(Exception e) {
				genpos = Double.NaN;
			}
			if (infoLength > 3) {
				try {
					physpos = Integer.parseInt(mapInfo[3]);
				} catch (Exception e) {
					physpos = -1;
				}
			} else {
				physpos = -1;
			}

			addMarker(id, chr, genpos, physpos);
			return true;
		}
		return false;
	}

	/**
	 * @param marker	the name of the marker
	 * @return	the marker index used to retrieve other information about the marker
	 */
	public int getMarkerIndex(String marker) {
		return markerIndex.get(marker).intValue();
	}
	
	/**
	 * @param index	the marker index
	 * @return	the name of the marker at that index position in the map
	 */
	public String getMarkerID(int index) {
		return featureList.get(index).id;
	}
	
	/**
	 * @return	the number of markers in this map
	 */
	public int getNumberOfMarkers() {
		return featureList.size();
	}
	
	/**
	 * @param index	the marker index
	 * @return	the chromosome of the marker, -1 if missing 
	 */
	public String getChromosome(int index){
		return featureList.get(index).chromosome;
	}
	
	/**
	 * @param index	the marker index
	 * @return	the genetic position of the marker, Double.NaN if missing
	 */
	public double getGeneticPosition(int index) {
		return featureList.get(index).geneticPos;
	}
	
	/**
	 * @param index	the marker index
	 * @return	the physical position of the marker, -1 if missing
	 */
	public int getPhysicalPosition(int index) {
		return featureList.get(index).physicalPos;
	}
	
	/**
	 * @param locus
	 * @param position
	 * @return	the genetic position corresponding to the given physical position calculated by interpolating between the map markers that have both genetic and physical positions.
	 */
	public double getGeneticPositionFromPhysical(Locus locus, int position) {
		return getGeneticPositionFromPhysical(locus.getChromosomeName(), position);
	}
	
	/**
	 * @param chromosome
	 * @param position
	 * @return	the genetic position corresponding to the given physical position calculated by interpolating between the map markers that have both genetic and physical positions.
	 */
	public double getGeneticPositionFromPhysical(String chromosome, int position) {
		mapFeature thisFeature = new mapFeature("none", chromosome, position, -1);
		int index = Collections.binarySearch(physicalFeatureList, thisFeature, new Comparator<mapFeature>(){

			@Override
			public int compare(mapFeature pos1, mapFeature pos2) {
				int comp = pos1.chromosome.compareTo(pos2.chromosome);
				if (comp != 0) return comp;
				if (pos1.physicalPos < pos2.physicalPos) return -1;
				if (pos1.physicalPos > pos2.physicalPos) return 1;
				return 0;
			}
			
		});
		
		if (index > -1) return physicalFeatureList.get(index).geneticPos;
		
		int n = physicalFeatureList.size();
		if (index == -1) {
			mapFeature left = physicalFeatureList.get(0);
			mapFeature right = physicalFeatureList.get(1);
			double ratio = (right.geneticPos - left.geneticPos) / (right.physicalPos - left.physicalPos);
			return left.geneticPos - ratio * (left.physicalPos - position); 
		} else if (index == -1 - n) {
			mapFeature left = physicalFeatureList.get(n - 2);
			mapFeature right = physicalFeatureList.get(n - 1);
			double ratio = (right.geneticPos - left.geneticPos) / (right.physicalPos - left.physicalPos);
			return right.geneticPos + ratio * (position - right.physicalPos); 
		} else {
			int insertionPoint = -1 - index;
			mapFeature left = physicalFeatureList.get(insertionPoint - 1);
			mapFeature right = physicalFeatureList.get(insertionPoint);
			double ratio = (right.geneticPos - left.geneticPos) / (right.physicalPos - left.physicalPos);
			return left.geneticPos + ratio * (position - left.physicalPos); 
		}
		
	}
	
	/**
	 * @param locus
	 * @param position
	 * @return	the physical position corresponding to the given genetic position calculated by interpolating between the map markers that have both genetic and physical positions.
	 */
	public int getPhysicalPositionFromGenetic(Locus locus, double position) {
		return getPhysicalPositionFromGenetic(locus.getChromosomeName(), position);
	}
	
	/**
	 * @param chromosome
	 * @param position
	 * @return	the physical position corresponding to the given genetic position calculated by interpolating between the map markers that have both genetic and physical positions.
	 */
	public int getPhysicalPositionFromGenetic(String chromosome, double position) {
		mapFeature thisFeature = new mapFeature("none", chromosome, -1, position);
		int index = Collections.binarySearch(physicalFeatureList, thisFeature, new Comparator<mapFeature>(){

			@Override
			public int compare(mapFeature pos1, mapFeature pos2) {
				int comp = pos1.chromosome.compareTo(pos2.chromosome);
				if (comp != 0) return comp;
				if (pos1.geneticPos < pos2.geneticPos) return -1;
				if (pos1.geneticPos > pos2.geneticPos) return 1;
				return 0;
			}
			
		});
		
		if (index > -1) return physicalFeatureList.get(index).physicalPos;
		
		int n = physicalFeatureList.size();
		if (index == -1) {
			mapFeature left = physicalFeatureList.get(0);
			mapFeature right = physicalFeatureList.get(1);
			double ratio = (right.physicalPos - left.physicalPos) / (right.geneticPos - left.geneticPos);
			return left.physicalPos - (int)(ratio * (left.geneticPos - position)); 
		} else if (index == -1 - n) {
			mapFeature left = physicalFeatureList.get(n - 2);
			mapFeature right = physicalFeatureList.get(n - 1);
			double ratio = (right.physicalPos - left.physicalPos) / (right.geneticPos - left.geneticPos);
			return right.physicalPos + (int)(ratio * (position - right.geneticPos)); 
		} else {
			int insertionPoint = -1 - index;
			mapFeature left = physicalFeatureList.get(insertionPoint - 1);
			mapFeature right = physicalFeatureList.get(insertionPoint);
			double ratio = (right.physicalPos - left.physicalPos) / (right.geneticPos - left.geneticPos);
			return left.physicalPos + (int)(ratio * (position - left.geneticPos)); 
		}
	}

	/**
	 * @return true, if the map has physical positions as well as genetic positions
	 */
	public boolean hasPhysicalPositions() {return hasPhysicalPositions;}
	
	/**
	 * @return	the name of this map
	 */
	public String getName() { return name; }

	@Override
	public int getColumnCount() {
		if (hasPhysicalPositions) return 4;
		return 3;
	}

	public String getColumnName(int col) {
		switch(col) {
		case 0:
			return "Marker";
		case 1:
			return "Chromosome";
		case 2:
			return "Genetic Position";
		case 3:
			return "Phys Position";
		}
		return "";
	}

	@Override
	public int getRowCount() {
		return featureList.size();
	}

	@Override
	public Object getValueAt(int row, int col) {
		switch(col) {
		case 0:
			return featureList.get(row).id;
		case 1:
			return featureList.get(row).chromosome;
		case 2:
			return featureList.get(row).geneticPos;
		case 3:
			return featureList.get(row).physicalPos;
		}
		return "";
	}

	//TableReport functions
	@Override
	public int getElementCount() {
		return getColumnCount() * getRowCount();
	}

	@Override
	public Object[] getRow(int row) {
		int ncol = getColumnCount();
		Object[] thisRow = new Object[ncol];
		for (int i = 0; i < ncol; i++) thisRow[i] = getValueAt(row, i);
		return thisRow;
	}

	@Override
	public Object[] getTableColumnNames() {
		int ncol = getColumnCount();
		Object[] names = new Object[ncol];
		for (int i = 0; i < ncol; i++) names[i] = getColumnName(i);
		return names;
	}

	@Override
	public Object[][] getTableData() {
		int ncols = getColumnCount();
		int nrows = getRowCount();
		Object[][] data = new Object[nrows][ncols];
		for (int r = 0; r < nrows; r++) {
			for (int c = 0; c < ncols; c++) {
				data[r][c] = getValueAt(r, c);
			}
		}
		return null;
	}

	@Override
	public Object[][] getTableData(int start, int end) {
		int ncols = getColumnCount();
		int nrows = end - start + 1;
		Object[][] data = new Object[nrows][ncols];
		for (int r = start; r <= end; r++) {
			for (int c = 0; c < ncols; c++) {
				data[r][c] = getValueAt(r, c);
			}
		}
		return null;
	}

	@Override
	public String getTableTitle() {
		return name;
	}


	class mapFeature implements Comparable<mapFeature>, Serializable {
                private static final long serialVersionUID = -5197800047652332969L;
		String id;
		String chromosome;
		int physicalPos;
		double geneticPos;
		
		public mapFeature(String id, String chromosome, int physical, double genetic) {
			this.id = id;
			this.chromosome = chromosome;
			physicalPos = physical;
			geneticPos = genetic;
		}
		
		public mapFeature(String id, String chromosome, double genetic) {
			this.id = id;
			this.chromosome = chromosome;
			physicalPos = -1;
			geneticPos = genetic;
		}

		@Override
		public boolean equals(Object obj) {
			if (obj instanceof mapFeature) {
				mapFeature feature = (mapFeature) obj;
				if (!chromosome.equals(feature.chromosome)) return false;
				if (!Double.isNaN(physicalPos) && !Double.isNaN(feature.physicalPos) && physicalPos != feature.physicalPos) return false;
				if (geneticPos != feature.geneticPos) return false;
				return true;
			}
			return false;
		}

		@Override
		public int compareTo(mapFeature feature) {
			//if the physical position = Double.Nan it cannot 
			int comp = chromosome.compareTo(feature.chromosome);
			if (comp != 0) return comp;
			if (physicalPos == -1 || feature.physicalPos == -1) {
				if (geneticPos > feature.geneticPos) return 1;
				else if (geneticPos < feature.geneticPos) return -1;
				else if (physicalPos > -1) return -1;
				else if (feature.physicalPos > -1) return 1;
				else return 0;
			} else {
				if (geneticPos > feature.geneticPos && physicalPos > feature.physicalPos) return 1;
				else if (geneticPos < feature.geneticPos && physicalPos < feature.physicalPos) return -1;
				else if (geneticPos == feature.geneticPos && physicalPos > feature.physicalPos) return 1;
				else if (geneticPos == feature.geneticPos && physicalPos < feature.physicalPos) return -1;
				else if (geneticPos == feature.geneticPos && physicalPos == feature.physicalPos) return 0;
				else {
					StringBuilder msg = new StringBuilder("The genetic and physical map positions are inconsistent for marker ");
					msg.append(id).append(" on chr ").append(chromosome).append(" at ");
					msg.append(geneticPos).append(", ").append(physicalPos);
					msg.append(" and marker ");
					msg.append(feature.id).append(" on chr ").append(feature.chromosome).append(" at ");
					msg.append(feature.geneticPos).append(", ").append(feature.physicalPos);
					throw new IllegalArgumentException("The genetic and physical map positions are inconsistent.");
				}
			}
			
		}
		
		
	}

}
