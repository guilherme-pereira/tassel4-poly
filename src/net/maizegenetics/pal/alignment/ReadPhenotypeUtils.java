package net.maizegenetics.pal.alignment;

import java.io.BufferedReader;
import java.io.IOException;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.TreeSet;
import java.util.regex.Pattern;

import net.maizegenetics.pal.ids.SimpleIdGroup;

import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;


import net.maizegenetics.util.Utils;


public class ReadPhenotypeUtils {
	public final static String missing = "\\?|nan|-999";
	private static Pattern missingPattern = Pattern.compile("\\?|nan|-999", Pattern.CASE_INSENSITIVE);
	//to prevent instantiation of this utility class
	private ReadPhenotypeUtils(){}
	
	/**
	 * @param original	the original class or discrete data
	 * @param tonumber	a double array, the same length as the original data, which will hold the class indices
	 * @param name		the trait name
	 * @param type		the trait type
	 * @return			the new Trait
	 */
	public static Trait makeCharacterTrait(String[] original, double[] tonumber, String name, String type) {
		TreeSet<String> dataset = new TreeSet<String>();
		for (String str : original) {
			if (!missingPattern.matcher(str).matches() && !str.contains("?")) dataset.add(str);
		}
		String[] levelNames = new String[dataset.size()];
		dataset.toArray(levelNames);
		int n = original.length;
		for (int i = 0; i < n; i++) {
			tonumber[i] = Arrays.binarySearch(levelNames, original[i]);
			if (tonumber[i] < -0.1) tonumber[i] = Double.NaN;
		}
		Trait trait = new Trait(name, true, type);
		
		trait.setLevelLabels(levelNames);
		return trait;
	}
	
	public static double[] doubleFromCharacterTrait(Trait trait, String[] textdata) {
		TreeSet<String> dataset = new TreeSet<String>();
		int n = textdata.length;
		double[] dbldata = new double[n];
		for (String str : textdata) {
			if (!missingPattern.matcher(str).matches() && !str.contains("?")) dataset.add(str);
		}
		String[] levelNames = new String[dataset.size()];
		dataset.toArray(levelNames);
		for (int i = 0; i < n; i++) {
			dbldata[i] = Arrays.binarySearch(levelNames, textdata[i]);
			if (dbldata[i] < -0.1) dbldata[i] = Double.NaN;
		}
		trait.setProperty(Trait.PROP_LEVEL_LABELS, levelNames);
		return dbldata;
	}

	/**
	 * @param inputFile		the name of the file with the phenotype data
	 * @return				a Phenotype containing the data
	 * @throws IOException
	 */
	public static Phenotype readPhenotypeFile(String inputFile) throws IOException {
		String sep = "\\s+";
		//BufferedReader br = new BufferedReader(new FileReader(inputFile));
        BufferedReader br = Utils.getBufferedReader(inputFile);
		String inputline = br.readLine();
		String[] parsedline = inputline.split(sep);
		br.close();
		
		//determine the input style here
		//does the line start with a number
		//if yes, the input line length 2 or 3?
		//if 2, then it is a phenotype file
		//if 3, then it is a polymorphism file
		
		//if not is it old style, new style annotated, or unannotated
		try {
			Integer.decode(parsedline[0]);
			if (parsedline.length == 3) return readNumericalAlignment(inputFile);
			else if (parsedline.length == 2 ) return readPolymorphismAlignment(inputFile);
			else return readGenericFile(inputFile);
		}
		catch(Exception e){}
		if (parsedline[0].startsWith("<Ann")) return readAnnotatedAlignment(inputFile);
		return readGenericFile(inputFile);
	}
	
	/**
	 * @param inputFile	the input file in TASSEL v2 annotated alignment format
	 * @return			a Phenotype
	 * @throws IOException 
	 */
	public static Phenotype readAnnotatedAlignment(String inputFile) throws IOException {
		String sep = "\\s+";
		//BufferedReader br = new BufferedReader(new FileReader(inputFile));
        BufferedReader br = Utils.getBufferedReader(inputFile);
		String[] parsedline = br.readLine().split(sep);
		String missing = "?";
		
		//possible labels include <TRANSPOSED>, <TAXA_NUMBER>, <LOCUS_NUMBER>,<POLY_TYPE>,<DELIMITED_VALUES>
		int taxaNumber = 0;
		int locusNumber = 0;
		boolean isDelimited = false;
		while (!parsedline[0].equalsIgnoreCase("<taxon_name>")) {
			if (!parsedline[0].equalsIgnoreCase("<Annotated>") ) {
				if((parsedline[0].startsWith("<")==false)||(parsedline[0].endsWith(">")==false)) {
					throw new IllegalArgumentException("Error before or with property: " + parsedline[0]);}
				if((parsedline[1].startsWith("<")==true)||(parsedline[1].endsWith(">")==true)) {
					throw new IllegalArgumentException("Error before or with value: " + parsedline[1]);}
				if (parsedline[0].equalsIgnoreCase("<TRANSPOSED>")) {
					if (parsedline[1].startsWith("N")) 
						throw new IllegalArgumentException("Error in " + inputFile + ": The annotated format only accepts transposed data. Data not imported. ");
				}
				if (parsedline[0].equalsIgnoreCase("<TAXA_NUMBER>")) {
					taxaNumber = Integer.parseInt(parsedline[1]);
				}
				if (parsedline[0].equalsIgnoreCase("<LOCUS_NUMBER>")) {
					locusNumber = Integer.parseInt(parsedline[1]);
				}
				if (parsedline[0].equalsIgnoreCase("<DELIMITED_VALUES>")) {
					if (parsedline[1].startsWith("Y")) isDelimited = true;
				}
			}
            parsedline = br.readLine().split(sep);
		}
		
		//now read taxa list
		String[] taxaNames = new String[taxaNumber];
		if (parsedline[0].equalsIgnoreCase("<taxon_name>")) {
			for (int i = 0; i < taxaNumber; i++) taxaNames[i] = parsedline[i+1].toUpperCase();
		}
		else throw new IllegalArgumentException("No taxon name section in " + inputFile+ " . Data not read.");
		
		//headers
		parsedline = br.readLine().split(sep);
		ArrayList<String> columns = new ArrayList<String>();
		int n = parsedline.length;
		int locusNameIndex = -1;
		int locusPosIndex = -1;
		for (int i = 0; i < n - 1; i++) {
			String colname = parsedline[i];
			if (!colname.startsWith("<") || !colname.endsWith(">")) 
				throw new IllegalArgumentException("Improper format for column names in " + inputFile + ". Data not imported.");
			else if (colname.equalsIgnoreCase("<chromosome_number>")) columns.add(Trait.PROP_CHROMOSOME);
			else if (colname.equalsIgnoreCase("<genetic_position>")) columns.add(Trait.PROP_POSITION);
			else if (colname.equalsIgnoreCase("<locus_name>")) {
				locusNameIndex = i;
				columns.add(Trait.PROP_LOCUS);
			}
			else if (colname.equalsIgnoreCase("<locus_position>")) {
				locusPosIndex = i;
				columns.add(Trait.PROP_LOCUS_POSITION);
			}
		}
		
		//read values
		int colNumber = columns.size();
		int firstValue = colNumber;
		int markerNumber = 1;
		LinkedList<Trait> traitList = new LinkedList<Trait>();
		String[] locusValues = new String[taxaNumber];
		double[] dblValues = new double[taxaNumber];
		DoubleMatrix2D genotypes = DoubleFactory2D.dense.make(taxaNumber, locusNumber);
		for (int loc = 0; loc < locusNumber; loc++) {
			parsedline = br.readLine().split(sep);
			
			if (isDelimited) {
				for (int t = 0; t < taxaNumber; t++) {
					locusValues[t] = parsedline[t + firstValue];
				}
			}
			else {
				for (int t = 0; t < taxaNumber; t++) {
					locusValues[t] = parsedline[firstValue].substring(t, t + 1);
				}
			}
			
			String traitname = "";
			if (locusNameIndex >= 0) traitname = parsedline[locusNameIndex];
			else if (locusPosIndex >= 0) traitname = traitname + "." + parsedline[locusPosIndex];
			else traitname = "m" + markerNumber++;
			Trait trait = makeCharacterTrait(locusValues, dblValues, traitname, Trait.TYPE_MARKER);
			for (int col = 0; col < colNumber; col++) {
				trait.setProperty(columns.get(col), parsedline[col + 1]);
			}
			traitList.add(trait);
			genotypes.viewColumn(loc).assign(dblValues);
		}
		br.close();
		return new SimplePhenotype(new SimpleIdGroup(taxaNames), traitList, genotypes);
	}
	
	/**
	 * @param inputFile	the input file in TASSEL v2 polymorphism format
	 * @return			a Phenotype
	 * @throws IOException 
	 */
	public static Phenotype readPolymorphismAlignment(String inputFile) throws IOException {
		//BufferedReader br = new BufferedReader(new FileReader(inputFile));
        BufferedReader br = Utils.getBufferedReader(inputFile);
		String sep = "[\\s]+";
		String[] parsedline = br.readLine().split("[:\\s]+");
		int taxaNumber = Integer.parseInt(parsedline[0]);
		int locusNumber = Integer.parseInt(parsedline[1]);
		String[] taxaNames = new String[taxaNumber];
		String[] locusNames = new String[locusNumber];
		LinkedList<Trait> traitList = new LinkedList<Trait>();
		DoubleMatrix2D genotypes = DoubleFactory2D.dense.make(taxaNumber, locusNumber);
		
		parsedline = br.readLine().split(sep);
		int firstLocus = 0;
		int n = parsedline.length;
		if (n != locusNumber) {
			if (n == locusNumber + 1) firstLocus = 1;
			else throw new IllegalArgumentException("Number of loci in header not equal to number of locus names in " + inputFile + ". Data not imported.");
		}
		for (int i = 0; i < locusNumber; i++) {
			locusNames[i] = parsedline[i + firstLocus];
		}
		
		String[][] locusText = new String[locusNumber][taxaNumber];
		for (int i = 0; i < taxaNumber; i++) {
			parsedline = br.readLine().split(sep);
			taxaNames[i] = parsedline[0].toUpperCase();
			for (int j = 0; j < locusNumber; j++) {
				locusText[j][i] = parsedline[j + 1];
			}
		}
			
		for (int i = 0; i < locusNumber; i++) {
			double[] dblval = new double[taxaNumber];
			traitList.add(makeCharacterTrait(locusText[i], dblval, locusNames[i], Trait.TYPE_MARKER));
			genotypes.viewColumn(i).assign(dblval);
		}
		
		return new SimplePhenotype(new SimpleIdGroup(taxaNames), traitList, genotypes);
	}

	/**
	 * @param inputFile	the input file in TASSEL v2 numerical format
	 * @return
	 * @throws IOException 
	 */
	public static Phenotype readNumericalAlignment(String inputFile) throws IOException {
		//BufferedReader br = new BufferedReader(new FileReader(inputFile));
        BufferedReader br = Utils.getBufferedReader(inputFile);
		String sep = "\\s+";
		String[] parsedline = br.readLine().trim().split(sep);
		if (parsedline.length != 3) throw new IllegalArgumentException("Improper header in " + inputFile + ". Data not imported.");
		
		int taxaNumber = Integer.parseInt(parsedline[0]);
		int traitNumber = Integer.parseInt(parsedline[1]);
		int headerNumber = Integer.parseInt(parsedline[2]);
		if (headerNumber < 1 || headerNumber > 2) throw new IllegalArgumentException("Numerical format only allows 1 or 2 headers in " + inputFile + ". Data not imported.");
		parsedline = br.readLine().split(sep);
		ArrayList<Trait> traitList = new ArrayList<Trait>();
		if (parsedline.length != traitNumber) throw new IllegalArgumentException("Incorrect number of trait names in " + inputFile + ". Data not imported.");
		for (int i = 0; i < traitNumber; i++) {
			Trait trait = new Trait(parsedline[i], false, Trait.TYPE_DATA);
			traitList.add(trait);
		}
		
		if (headerNumber > 1) {
			parsedline = br.readLine().split(sep);
			if (parsedline.length != traitNumber) throw new IllegalArgumentException("Incorrect number of environment names in " + inputFile + ". Data not imported.");
			for (int i = 0; i < traitNumber; i++) traitList.get(i).addFactor(Trait.FACTOR_ENV, parsedline[i]);
		}
		
		DoubleMatrix2D phenotypes = DoubleFactory2D.dense.make(taxaNumber, traitNumber);
		String[] taxaNames = new String[taxaNumber];
		for (int i = 0; i < taxaNumber; i++) {
			parsedline = br.readLine().split(sep);
//			taxaNames[i] = parsedline[0].toUpperCase();
                        taxaNames[i] = parsedline[0];
			for (int j = 0; j < traitNumber; j++) {
				String strval = parsedline[j+1];
				if (strval.toLowerCase().matches(missing)) {
					phenotypes.setQuick(i, j, Double.NaN);
				}
				else phenotypes.setQuick(i, j, Double.parseDouble(strval));
			}
		}
		
		return new SimplePhenotype(new SimpleIdGroup(taxaNames), traitList, phenotypes);
	}
	
	/**
	 * @param inputFile	the input file with TASSEL v3 annotations or with no input directives
	 * @return
	 * @throws IOException 
	 */
	public static Phenotype readGenericFile(String inputFile) throws IOException {
		
		//BufferedReader br = new BufferedReader(new FileReader(inputFile));
        BufferedReader br = Utils.getBufferedReader(inputFile);
		String inputline = br.readLine();
		Pattern sep = Pattern.compile("\\s+");
		
		String[] traitName = null;
		String[] markerName = null;
		ArrayList<String[]> headerList = new ArrayList<String[]>();
		String[] format = null;
		String[] use = null;
		String[] chromosome = null;
		String[] chrpos = null;
		String[] locus = null;
		String[] locuspos = null;
		boolean isNumeric = false;
		boolean isCharacter = false;
		boolean isData = false;
		boolean isFactor = false;
		boolean isCovariate = false;
		int numberOfColumns = 0;
		
		//process header rows and count the non-blank rows
		int numberOfDataLines = 0;
		while (inputline != null) {
			inputline = inputline.trim();
			String[] parsedline = sep.split(inputline);
			if (parsedline.length > 1 && !inputline.startsWith("<") && !inputline.startsWith("#")) {
				numberOfDataLines++;
			}
			else if (parsedline[0].toUpperCase().equals("<TRAIT>")) {
				traitName = processHeader(numberOfColumns, parsedline, inputFile);
				numberOfColumns = traitName.length;
			}
			else if (parsedline[0].toUpperCase().equals("<MARKER>")) {
				markerName = processHeader(numberOfColumns, parsedline, inputFile);
				numberOfColumns = markerName.length;
			}
			else if (parsedline[0].toUpperCase().equals("<HEADER")) {
				String[] parsedHeader = processHeader(numberOfColumns, parsedline, inputFile);
				headerList.add(parsedHeader);
				numberOfColumns = parsedHeader.length - 1;
			}
			else if (parsedline[0].toUpperCase().equals("<FORMAT>")) {
				format = processHeader(numberOfColumns, parsedline, inputFile);
				numberOfColumns = format.length;
			}
			else if (parsedline[0].toUpperCase().equals("<USE>")) {
				use = processHeader(numberOfColumns, parsedline, inputFile);
				numberOfColumns = use.length;
			}
			else if (parsedline[0].toUpperCase().equals("<CHROMOSOME>")) {
				chromosome = processHeader(numberOfColumns, parsedline, inputFile);
				numberOfColumns = chromosome.length;
			}
			else if (parsedline[0].toUpperCase().equals("<CHROMOSOME_POSITION>")) {
				chrpos = processHeader(numberOfColumns, parsedline, inputFile);
				numberOfColumns = chrpos.length;
			}
			else if (parsedline[0].toUpperCase().equals("<LOCUS>")) {
				locus = processHeader(numberOfColumns, parsedline, inputFile);
				numberOfColumns = locus.length;
			}
			else if (parsedline[0].toUpperCase().equals("<LOCUS_POSITION>")) {
				locuspos = processHeader(numberOfColumns, parsedline, inputFile);
				numberOfColumns = locuspos.length;
			}
			else if (parsedline[0].toUpperCase().equals("<NUMERIC>")) {
				isNumeric = true;
			}
			else if (parsedline[0].toUpperCase().equals("<CHARACTER>")) {
				isCharacter = true;
			}
			else if (parsedline[0].toUpperCase().equals("<DATA>")) {
				isData = true;
			}
			else if (parsedline[0].toUpperCase().equals("<COVARIATE>")) {
				isCovariate = true;
			}
			else if (parsedline[0].toUpperCase().equals("<FACTOR>")) {
				isFactor = true;
			}

			inputline = br.readLine();
		}
		br.close();
		
		//create Traits
		ArrayList<Trait> traitList = new ArrayList<Trait>(numberOfColumns);
		for (int c = 0; c < numberOfColumns; c++) {
			String name;
			boolean discrete;
			String type;
			
			if (traitName == null) {
				if (markerName == null) {
					throw new IllegalArgumentException("Error in " + inputFile + ": neither trait nor marker names are defined.");
				}
				else {
					name = markerName[c];
					type = Trait.TYPE_MARKER;
					if (isCharacter) discrete = true;
					else if (isNumeric) discrete = false;
					else if (format == null) discrete = true;
					else if (format[c].toUpperCase().startsWith("N")) discrete = false;
					else discrete = true;
				}
			}
			else {
				name = traitName[c];
				if (isData) {
					type = Trait.TYPE_DATA;
					discrete = false;
				}
				else if (isFactor) {
					type = Trait.TYPE_FACTOR;
					discrete = true;
				}
				else if (isCovariate) {
					type = Trait.TYPE_COVARIATE;
					discrete = false;
				}
				else if (use == null) {
					type = Trait.TYPE_DATA;
					discrete = false;
				}
				else if (use[c].toUpperCase().startsWith("D")) {
					type = Trait.TYPE_DATA;
					discrete = false;
				}
				else if (use[c].toUpperCase().startsWith("C")) {
					type = Trait.TYPE_COVARIATE;
					discrete = false;
				}
				else if (use[c].toUpperCase().startsWith("F")) {
					type = Trait.TYPE_FACTOR;
					discrete = true;
				}
				else if (use[c].toUpperCase().startsWith("M")) {
					type = Trait.TYPE_MARKER;
					if (isCharacter) discrete = true;
					else if (isNumeric) discrete = false;
					else if (format == null) discrete = true;
					else if (format[c].toUpperCase().startsWith("N")) discrete = false;
					else discrete = true;
				}
				else {
					type = Trait.TYPE_DATA;
					discrete = false;
				}
				
			}
			Trait trait = new Trait(name, discrete, type);
			traitList.add(trait);
			for (String[] header : headerList) {
				trait.addFactor(header[0], header[c + 1]);
			}
			if (chromosome != null) trait.setProperty(Trait.PROP_CHROMOSOME, chromosome[c]);
			if (chrpos != null) trait.setProperty(Trait.PROP_POSITION, Double.parseDouble(chrpos[c]));
			if (locus != null) trait.setProperty(Trait.PROP_LOCUS, locus[c]);
			if (locuspos != null) trait.setProperty(Trait.PROP_LOCUS_POSITION, Integer.parseInt(locuspos[c]));
		}
		
		//process body of data
		String[] taxanames = new String[numberOfDataLines];
		DoubleMatrix2D data	= DoubleFactory2D.dense.make(numberOfDataLines, numberOfColumns);
        br = Utils.getBufferedReader(inputFile);
		inputline = br.readLine();
		int totallines = 0;
		int linecount = 0;
		if (markerName != null && isNumeric) { //if there are a lot of markers saving a String array uses a lot of memory
			while (inputline != null) {
				totallines++;
				inputline = inputline.trim();
				String[] parsedline = sep.split(inputline);
				if (parsedline.length > 1 && !inputline.startsWith("<") && !inputline.startsWith("#")) {
					if (parsedline.length != numberOfColumns + 1) {
						StringBuilder msg = new StringBuilder("Error in ");
						msg.append(inputFile);
						msg.append(" line ").append(totallines);
						msg.append(": Incorrect number of data values.");
						throw new IllegalArgumentException(msg.toString());
					}
					taxanames[linecount] = parsedline[0];
					for (int c = 0; c < numberOfColumns; c++) {
						double val;
						try {
							val = Double.parseDouble(parsedline[c + 1]);
						} catch(NumberFormatException e) {
							val = Double.NaN;
						}
						data.set(linecount, c, val);
					}
					linecount++;
				}
				inputline = br.readLine();

			} 
		} else {
			String[][] textdata = new String[numberOfColumns][numberOfDataLines];
			//br = new BufferedReader(new FileReader(inputFile));
			while (inputline != null) {
				totallines++;
				inputline = inputline.trim();
				String[] parsedline = sep.split(inputline);
				if (parsedline.length > 1 && !inputline.startsWith("<") && !inputline.startsWith("#")) {
					if (parsedline.length != numberOfColumns + 1) {
						StringBuilder msg = new StringBuilder("Error in ");
						msg.append(inputFile);
						msg.append(" line ").append(totallines);
						msg.append(": Incorrect number of data values.");
						throw new IllegalArgumentException(msg.toString());
					}
					taxanames[linecount] = parsedline[0];
					for (int c = 0; c < numberOfColumns; c++) {
						textdata[c][linecount] = parsedline[c + 1];
					}
					linecount++;
				}

				inputline = br.readLine();
			}

			for (int c = 0; c < numberOfColumns; c++) {
				if (traitList.get(c).isDiscrete) {
					data.viewColumn(c).assign(doubleFromCharacterTrait(traitList.get(c), textdata[c]));
				}
				else {
					for (int r = 0; r < numberOfDataLines; r++) {
						String textval = textdata[c][r];
						if (textval.equals("-999")) data.setQuick(r, c, Double.NaN);
						else {
							try {data.setQuick(r, c, Double.parseDouble(textdata[c][r]));}
							catch (Exception e) {data.setQuick(r, c, Double.NaN);}
						}
					}
				}
			}
		}
		
		br.close();
		return new SimplePhenotype(new SimpleIdGroup(taxanames), traitList, data);
	}
	
	private static String[] processHeader(int numberOfColumns, String[] parsedline, String filename) {
		if (parsedline[0].equalsIgnoreCase("<Header")) {
			String headername = parsedline[1].split("[=>\\s]")[1];
			if (!parsedline[1].contains("name=") || headername.length() == 0) {
				StringBuilder msg = new StringBuilder("Error in ");
				msg.append(filename);
				msg.append(": Improperly formatted <Header name=> line.");
				throw new IllegalArgumentException(msg.toString());
			}
			int finalBracketPosition = 0;
			while(!parsedline[finalBracketPosition].contains(">")) finalBracketPosition++;
			if (numberOfColumns == 0) numberOfColumns = parsedline.length - finalBracketPosition - 1;
			else if (numberOfColumns != parsedline.length - finalBracketPosition - 1) {
				StringBuilder msg = new StringBuilder("Error in ");
				msg.append(filename);
				msg.append(": The number of ");
				msg.append(parsedline[0]).append(" ").append(parsedline[1]);
				msg.append(" columns does not match the number of columns in previous header rows");
				throw new IllegalArgumentException(msg.toString());
			}
			
			
			String[] contents = new String[numberOfColumns + 1];
			contents[0] = headername;
			for (int i = 0; i < numberOfColumns; i++) {
				contents[i + 1] = parsedline[i + finalBracketPosition + 1];
			}
			return contents;
		}
		else {
			if (numberOfColumns == 0) numberOfColumns = parsedline.length - 1;
			else if (numberOfColumns != parsedline.length - 1) {
				StringBuilder msg = new StringBuilder("Error in ");
				msg.append(filename);
				msg.append(": The number of ");
				msg.append(parsedline[0]);
				msg.append(" columns does not match the number of columns in previous header rows");
				throw new IllegalArgumentException(msg.toString());
			}
			String[] contents = new String[numberOfColumns];
			for (int i = 0; i < numberOfColumns; i++) {
				contents[i] = parsedline[i + 1];
			}
			return contents;
		}
	}
}
