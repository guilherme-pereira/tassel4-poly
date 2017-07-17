package net.maizegenetics.pal.taxa;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;

import java.util.*;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.prefs.TasselPrefs;
import org.apache.log4j.Logger;

/**
 * In memory immutable instance of {@link TaxaList}. Basic list of taxa
 * (samples) that are used in Alignments and other purposes.
 *
 * Use {@link TaxaListBuilder} to instantiate.
 *
 * @author Ed Buckler
 *
 */
public class TaxaArrayList implements TaxaList {

    private static final Logger myLogger = Logger.getLogger(TaxaArrayList.class);
    private final List<AnnotatedTaxon> myTaxaList;
    private final int myNumTaxa;
    private final Multimap<String, Integer> myNameToIndex;

    TaxaArrayList(TaxaListBuilder builder) {

        List<AnnotatedTaxon> srcList = builder.getImmutableList();
        myTaxaList = new ArrayList<AnnotatedTaxon>(srcList.size());
        myNumTaxa = srcList.size();
        myNameToIndex = HashMultimap.create(srcList.size() * 2, 1);
        int index = 0;
        for (AnnotatedTaxon annotatedTaxon : srcList) {
            myTaxaList.add(annotatedTaxon);
            if (myNameToIndex.containsKey(annotatedTaxon.getFullName())) {
                myLogger.warn("init: Taxa name is duplicated :" + annotatedTaxon.getFullName());
            }
            myNameToIndex.put(annotatedTaxon.getFullName(), index);

            // Ed, we need to talk about this. -Terry
            //if (!annotatedTaxon.getFullName().equals(annotatedTaxon.getName())) {
            //    myNameToIndex.put(annotatedTaxon.getName(), index);
            //}

            index++;
        }
    }

    @Override
    public int getTaxaCount() {
        return myNumTaxa;
    }

    @Override
    public String getTaxaName(int index) {
        return myTaxaList.get(index).getName();
    }

    @Override
    public String getFullTaxaName(int index) {
        return myTaxaList.get(index).getFullName();
    }

    @Override
    public int size() {
        return myNumTaxa;
    }

    @Override
    public List<Integer> getIndicesMatchingTaxon(String name) {

        TasselPrefs.TASSEL_IDENTIFIER_JOIN_TYPES type = TasselPrefs.getIDJoinStrict();

        if (type == TasselPrefs.TASSEL_IDENTIFIER_JOIN_TYPES.Strict) {
            return new ArrayList<Integer>(myNameToIndex.get(name));
        }

        List<Integer> result = new ArrayList<Integer>(1);
        for (int i = 0, n = getTaxaCount(); i < n; i++) {
            if (get(i).equals(name)) {
                result.add(i);
            }
        }

        return result;

    }

    @Override
    public List<Integer> getIndicesMatchingTaxon(Identifier taxon) {

        TasselPrefs.TASSEL_IDENTIFIER_JOIN_TYPES type = TasselPrefs.getIDJoinStrict();

        if (type == TasselPrefs.TASSEL_IDENTIFIER_JOIN_TYPES.Strict) {
            return new ArrayList<Integer>(myNameToIndex.get(taxon.getFullName()));
        }

        List<Integer> result = new ArrayList<Integer>(1);
        for (int i = 0, n = getTaxaCount(); i < n; i++) {
            if (get(i).equals(taxon)) {
                result.add(i);
            }
        }

        return result;

    }

    @Override
    public boolean isEmpty() {
        return myTaxaList.isEmpty();
    }

    @Override
    public boolean contains(Object o) {
        return myTaxaList.contains(o);
    }

    @Override
    public Iterator<AnnotatedTaxon> iterator() {
        return myTaxaList.iterator();
    }

    @Override
    public Object[] toArray() {
        return myTaxaList.toArray();
    }

    @Override
    public <T> T[] toArray(T[] a) {
        return myTaxaList.toArray(a);
    }

    @Override
    public boolean add(AnnotatedTaxon annotatedTaxon) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public boolean remove(Object o) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public boolean containsAll(Collection<?> c) {
        return myTaxaList.containsAll(c);
    }

    @Override
    public boolean addAll(Collection<? extends AnnotatedTaxon> c) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public boolean addAll(int index, Collection<? extends AnnotatedTaxon> c) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public boolean removeAll(Collection<?> c) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public boolean retainAll(Collection<?> c) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public void clear() {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public AnnotatedTaxon get(int index) {
        return myTaxaList.get(index);
    }

    @Override
    public AnnotatedTaxon set(int index, AnnotatedTaxon element) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public void add(int index, AnnotatedTaxon element) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public AnnotatedTaxon remove(int index) {
        throw new UnsupportedOperationException("This Class is Immutable.");
    }

    @Override
    public int indexOf(Object o) {
        return myTaxaList.indexOf(o);
    }

    @Override
    public int lastIndexOf(Object o) {
        return myTaxaList.lastIndexOf(o);
    }

    @Override
    public ListIterator<AnnotatedTaxon> listIterator() {
        return listIterator(0);
    }

    @Override
    public ListIterator<AnnotatedTaxon> listIterator(final int index) {
        return new ListIterator<AnnotatedTaxon>() {
            private final ListIterator<AnnotatedTaxon> i = myTaxaList.listIterator(index);

            public boolean hasNext() {
                return i.hasNext();
            }

            public AnnotatedTaxon next() {
                return i.next();
            }

            public boolean hasPrevious() {
                return i.hasPrevious();
            }

            public AnnotatedTaxon previous() {
                return i.previous();
            }

            public int nextIndex() {
                return i.nextIndex();
            }

            public int previousIndex() {
                return i.previousIndex();
            }

            public void remove() {
                throw new UnsupportedOperationException();
            }

            public void set(AnnotatedTaxon e) {
                throw new UnsupportedOperationException();
            }

            public void add(AnnotatedTaxon e) {
                throw new UnsupportedOperationException();
            }
        };
    }

    @Override
    public List<AnnotatedTaxon> subList(int fromIndex, int toIndex) {
        return Collections.unmodifiableList(myTaxaList.subList(fromIndex, toIndex));
    }
}
