package net.maizegenetics.pal.taxa;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.HDF5LinkInformation;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.HapMapHDF5Constants;
import net.maizegenetics.pal.ids.IdGroup;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * A builder for creating immutable {@link TaxaList} instances.
 * Example:<pre>   {@code
 *   TaxaListBuilder tlb=new TaxaListBuilder();
 *   for (int i = 0; i < 10; i++) {
 *      AnnotatedTaxon at= new AnnotatedTaxon.Builder("Z"+i+":Line:mays:Zea")
 *           .inbreedF(0.99f)
 *           .parents("B73","B97")
 *           .pedigree("(B73xB97)S6I1")
 *           .build();
 *       tlb.add(at);
 *       }
 *   TaxaList tl=tlb.build();}</pre>
 *   <p></p>
 *   If building from HDF5:<pre>
 *   {@code
 *   TaxaList tl=new TaxaListBuilder().buildFromHDF5(testMutFile);
 *   }</pre>
 *
 * @author Ed Buckler
 */
public class TaxaListBuilder {
    private final List<AnnotatedTaxon> myTaxaList;

    public TaxaListBuilder() {
        myTaxaList=new ArrayList<AnnotatedTaxon>();
    }

    public TaxaListBuilder add(AnnotatedTaxon taxon) {
        myTaxaList.add(taxon);
        return this;
    }

    public TaxaListBuilder addAll(List<AnnotatedTaxon> taxa) {
        for (AnnotatedTaxon annotatedTaxon : taxa) {
            myTaxaList.add(annotatedTaxon);
        }
        return this;
    }

    public TaxaListBuilder addAll(Alignment a) {
        //TODO change this over to get TaxaList from the Alignment, and then add the annotated alignments
        IdGroup ids=a.getIdGroup();
        for (int i = 0; i <ids.getIdCount() ; i++) {
            myTaxaList.add(new AnnotatedTaxon.Builder(ids.getIdentifier(i)).build());
        }
        return this;
    }

    /*Sort the taxa by their natural order (alphabetically by name)*/
    public TaxaListBuilder sort() {
        Collections.sort(myTaxaList);
        return this;
    }

    public TaxaList build() {
        return new TaxaArrayList(this);
    }

    public TaxaList buildFromHDF5(String hdf5FileName) {
        IHDF5Reader reader= HDF5Factory.openForReading(hdf5FileName);
        myTaxaList.clear();
        List<HDF5LinkInformation> fields=reader.getAllGroupMemberInformation(HapMapHDF5Constants.GENOTYPES, true);
        for (HDF5LinkInformation is : fields) {
            if(is.isDataSet()==false) continue;
            myTaxaList.add(new AnnotatedTaxon.Builder(is.getName()).build());
        }
        sort();
        return build();
    }

    //Default package private method to hand the list to the instance
    List<AnnotatedTaxon> getImmutableList() {
         return Collections.unmodifiableList(myTaxaList);
    }



}
