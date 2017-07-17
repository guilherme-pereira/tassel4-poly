package net.maizegenetics.pal.site;

import com.google.common.base.Preconditions;
import com.google.common.collect.ArrayListMultimap;
import net.maizegenetics.pal.site.Position.Allele;

import java.nio.IntBuffer;
import java.util.*;

/**
 * In memory immutable instance of {@link PositionList}.  Use the {@link PositionArrayList.Builder}
 * to create the list.  This list is sorted by position.
 *
 * @author Ed Buckler
 */
public final class PositionArrayList implements PositionList {
    private final List<Position> mySiteList;
    private final int numPositions;
    private final byte[][] alleles;
    private final Map<Chromosome,ChrOffPos> myChrOffPosTree;
    private final Map<String,Chromosome> myChrNameHash;

    private static class ChrOffPos {
        final int startSiteOff;
        final int endSiteOff;
        final int[] position;
        private ChrOffPos(int startSiteOff, int endSiteOff, int[] position) {
            this.startSiteOff=startSiteOff;
            this.endSiteOff=endSiteOff;
            this.position=position;
        }
    }

    private PositionArrayList(ArrayList<Position> builderList) {
        this.numPositions=builderList.size();
        alleles=new byte[Allele.COUNT][numPositions];
//        refAlleles=new byte[numPositions];
//        majorAlleles=new byte[numPositions];
//        ancAlleles=new byte[numPositions];
//        hiDepAlleles=new byte[numPositions];
        ArrayListMultimap<Chromosome,Integer> pTS=ArrayListMultimap.create();
        mySiteList=new ArrayList<Position>(builderList.size());
        myChrOffPosTree=new TreeMap<Chromosome,ChrOffPos>();
        myChrNameHash=new HashMap();
        int currStart=0;
        Chromosome currChr=builderList.get(0).getChromosome();
        for (int i=0; i<builderList.size(); i++) {
            Position ap=builderList.get(i);
            for (Allele allele : Allele.values()) {
              alleles[allele.index()][i]=ap.getAllele(allele);
            }
            mySiteList.add(ap);
            if((i==(builderList.size()-1))||!ap.getChromosome().equals(currChr)) {
                int end=(i==builderList.size()-1)?i:i-1;
                myChrOffPosTree.put(currChr, new ChrOffPos(currStart, i-1, null));
                currChr=ap.getChromosome();
                currStart=i;
            }
            pTS.put(ap.getChromosome(),ap.getPosition());
        }
        for (Chromosome chr: pTS.keySet()) {
            List<Integer> p=pTS.get(chr);
            int[] intP=new int[p.size()];
            for (int i=0; i<intP.length; i++) {intP[i]=p.get(i);}
            ChrOffPos currOff=myChrOffPosTree.get(chr);
            myChrOffPosTree.put(chr, new ChrOffPos(currOff.startSiteOff, currOff.endSiteOff, intP));
            myChrNameHash.put(chr.getName(),chr);
            }
        pTS=null;
    }

    @Override
    public byte getReferenceAllele(int site) {
        return mySiteList.get(site).getAllele(Allele.REF);
    }
    
    @Override
    public byte[] getReference(int startSite, int endSite) {
        byte[] result = new byte[endSite - startSite];
        System.arraycopy(alleles[Allele.REF.index()],startSite,result,0, result.length);
        return result;
    }

    @Override
    public byte[] getReference() {
        return Arrays.copyOf(alleles[Allele.REF.index()],alleles[Allele.REF.index()].length);
    }

    @Override
    public boolean hasReference() {
        return true;
    }

    @Override
    public String[] getSNPIDs() {
        String[] theIDs=new String[mySiteList.size()];
        for (int i = 0; i < theIDs.length; i++) {
            theIDs[i]=mySiteList.get(i).getSNPID();
        }
        return theIDs;
    }

    @Override
    public String getSNPID(int site) {
        return mySiteList.get(site).getSNPID();
    }

    @Override
    public int getSiteCount() {
        return numPositions;
    }

    @Override
    public int getChromosomeSiteCount(Chromosome chromosome) {
        return myChrOffPosTree.get(chromosome).position.length;
    }

    @Override
    public int[] getStartAndEndOfChromosome(Chromosome chromosome) {
        ChrOffPos cop=myChrOffPosTree.get(chromosome);
        if(cop==null) return null;
        return new int[]{cop.startSiteOff,cop.endSiteOff};
    }

    @Override
    public int getPositionInChromosome(int site) {
        return mySiteList.get(site).getPosition();
    }

    @Override
    public int getSiteOfPhysicalPosition(int physicalPosition, Chromosome chromosome) {
        ChrOffPos cop=myChrOffPosTree.get(chromosome);
        if(cop==null) return Integer.MIN_VALUE;
        int i=Arrays.binarySearch(cop.position, physicalPosition); //AvgPerObj:227.5715ns  for 2million positions
        i+=(i<0)?-cop.startSiteOff:cop.startSiteOff;
        while((i>0)&&(physicalPosition==get(i-1).getPosition())) {i--;} //backup to the first position if there are duplicates
        return i;
    }

    @Override
    public int getSiteOfPhysicalPosition(int physicalPosition, Chromosome chromosome, String snpID) {
        int result=getSiteOfPhysicalPosition(physicalPosition, chromosome);
        if (result < 0) {return result;}
        else {
            if (snpID.equals(getSNPID(result))) {return result;
            } else {
                int index=result;
                while ((index < numPositions) && (getPositionInChromosome(index) == physicalPosition)) {
                    if (snpID.equals(getSNPID(index))) {return index;}
                    result++;
                }
                return -result - 1;
            }
        }
    }

    @Override
    public int[] getPhysicalPositions() {
        int[] result=new int[numPositions];
        IntBuffer ib=IntBuffer.wrap(result);
        for (ChrOffPos cop: myChrOffPosTree.values()) {
            ib.put(cop.position);
        }
        return result;
    }

    @Override
    public String getChromosomeName(int site) {
        return mySiteList.get(site).getChromosome().getName();
    }

    @Override
    public Chromosome getChromosome(int site) {
        return mySiteList.get(site).getChromosome();
    }

    @Override
    public Chromosome getChromosome(String name) {
        return myChrNameHash.get(name);
    }

    @Override
    public Chromosome[] getChromosomes() {
        return myChrOffPosTree.keySet().toArray(new Chromosome[0]);
    }

    @Override
    public int getNumChromosomes() {
        return myChrOffPosTree.size();
    }

    @Override
    public int[] getChromosomesOffsets() {
        int[] result=new int[myChrOffPosTree.size()];
        int index=0;
        for (ChrOffPos cop: myChrOffPosTree.values()) {
            result[index++]=cop.startSiteOff;
        }
        return result;
    }

    @Override
    public int getIndelSize(int site) {
        return mySiteList.get(site).getKnownVariants()[1].length();
    }

    @Override
    public boolean isIndel(int site) {
        return mySiteList.get(site).isIndel();
    }

    @Override
    public String getGenomeAssembly() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean isPositiveStrand(int site) {
        return (1==mySiteList.get(site).getStrand());
    }
    
    // List methods

    @Override
    public int size() {
        return mySiteList.size();
    }

    @Override
    public boolean isEmpty() {
        return mySiteList.isEmpty();
    }

    @Override
    public boolean contains(Object o) {
        return mySiteList.contains(o);
    }

    @Override
    public Iterator<Position> iterator() {
        return mySiteList.iterator();
    }

    @Override
    public Object[] toArray() {
        return mySiteList.toArray();
    }

    @Override
    public <AnnotatedSite> AnnotatedSite[] toArray(AnnotatedSite[] a) {
        return mySiteList.toArray(a);
    }

    /**Not supported immutable class*/
    @Override@Deprecated
    public boolean add(Position e) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    /**Not supported immutable class*/
    @Override@Deprecated
    public boolean remove(Object o) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public boolean containsAll(Collection<?> c) {
        return mySiteList.containsAll(c);
    }

    /**Not supported immutable class*/
    @Override@Deprecated
    public boolean addAll(Collection<? extends Position> c) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    /**Not supported immutable class*/
    @Override@Deprecated
    public boolean addAll(int index, Collection<? extends Position> c) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    /**Not supported immutable class*/
    @Override@Deprecated
    public boolean removeAll(Collection<?> c) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    /**Not supported immutable class*/
    @Override@Deprecated
    public boolean retainAll(Collection<?> c) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    /**Not supported immutable class*/
    @Override@Deprecated
    public void clear() {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public Position get(int index) {
        return mySiteList.get(index);
    }

    /**Not supported immutable class*/
    @Override@Deprecated
    public Position set(int index, Position element) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    /**Not supported immutable class*/
    @Override@Deprecated
    public void add(int index, Position element) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    /**Not supported immutable class*/
    @Override@Deprecated
    public Position remove(int index) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public int indexOf(Object o) {
        return mySiteList.indexOf(o);
    }

    @Override
    public int lastIndexOf(Object o) {
        return mySiteList.lastIndexOf(o);
    }

    @Override
    public ListIterator<Position> listIterator() {
        return listIterator(0);
    }

    @Override
    public ListIterator<Position> listIterator(final int index) {
        return new ListIterator<Position>() {
            private final ListIterator<Position> i= mySiteList.listIterator(index);
            public boolean hasNext()     {return i.hasNext();}
            public Position next()              {return i.next();}
            public boolean hasPrevious() {return i.hasPrevious();}
            public Position previous()          {return i.previous();}
            public int nextIndex()       {return i.nextIndex();}
            public int previousIndex()   {return i.previousIndex();}
            public void remove() {throw new UnsupportedOperationException();}
            public void set(Position e) {throw new UnsupportedOperationException();}
            public void add(Position e) {throw new UnsupportedOperationException();}
        };
    }

    @Override
    public List<Position> subList(int fromIndex, int toIndex) {
        return Collections.unmodifiableList(mySiteList.subList(fromIndex, toIndex));
    }

    /**
     * A builder for creating immutable PositionArrayList
     *
     * <p>Example:
     * <pre>   {@code
     *   PositionArrayList.Builder b=new PositionArrayList.Builder();
     *   for (int i = 0; i <size; i++) {
     *       Position ap=new CoreAnnotatedPosition.Builder(chr[chrIndex[i]],pos[i]).refAllele(refSeq[i]).build();
     *       b.add(ap);
     *       }
     *   instance=b.build();}
     * <p></p>
     * If being built separately from the genotypes, then use validate ordering to make sure sites are added in the
     * indended order.  This list WILL be sorted.
     * <p>Builder instances can be reused - it is safe to call {@link #build}
     * multiple times to build multiple lists in series. Each new list
     * contains the one created before it.
     */
    public static class Builder {
        private final ArrayList<Position> contents = new ArrayList<Position>();

        /**
         * Creates a new builder. The returned builder is equivalent to the builder
         * generated by {@link }.
         */
        public Builder() {}

        /**
         * Adds {@code element} to the {@code PositionList}.
         *
         * @param element the element to add
         * @return this {@code Builder} object
         * @throws NullPointerException if {@code element} is null
         */
        public Builder add(Position element) {
            Preconditions.checkNotNull(element, "element cannot be null");
            contents.add(element);
            return this;
        }

        /**
         * Adds each element of {@code elements} to the {@code PositionList}.
         *
         * @param elements the {@code Iterable} to add to the {@code PositionList}
         * @return this {@code Builder} object
         * @throws NullPointerException if {@code elements} is or contains null
         */
        public Builder addAll(Iterable<? extends Position> elements) {
            if (elements instanceof Collection) {
                @SuppressWarnings("unchecked")
                Collection<? extends Position> collection = (Collection<? extends Position>) elements;
                contents.ensureCapacity(contents.size() + collection.size());
            }
            for (Position elem : elements) {
                Preconditions.checkNotNull(elem, "elements contains a null");
                contents.add(elem);
            }
            return this;
        }

        /*
        Returns whether List is already ordered.  Important to check this if genotype and sites are separately built, as the
         PositionArrayList must be sorted, and will be with build.
         */
        public boolean validateOrdering() {
            boolean result=true;
            Position startAP=contents.get(0);
            for (Position ap:contents) {
              if(ap.compareTo(startAP)<0) return false;
            }
            return result;
        }

        /**
         * Returns a newly-created {@code ImmutableList} based on the contents of
         * the {@code Builder}.
         */
        public PositionList build() {
            if(!validateOrdering()) {
                System.out.println("Beginning Sort of Position List");
                Collections.sort(contents);
                System.out.println("Finished Sort of Position List");
            }
            return new PositionArrayList(contents);
        }
    }

}
