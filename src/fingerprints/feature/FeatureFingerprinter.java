/* Copyright (C) 2013       Syed Asad Rahman <S9asad@gmail.com>
 *               
 *
 * Contact: cdk-devel@lists.sourceforge.net
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * All we ask is that proper credit is given for our work, which includes
 * - but is not limited to - adding the above copyright notice to the beginning
 * of your source code files, and to any copyright notice that you may distribute
 * with programs based on this work.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package fingerprints.feature;

import static fingerprints.helper.RandomNumber.generateMersenneTwisterRandomNumber;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.lang3.builder.HashCodeBuilder;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.RingSet;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.exception.InvalidSmilesException;
import org.openscience.cdk.fingerprint.BitSetFingerprint;
import org.openscience.cdk.fingerprint.Fingerprinter;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.fingerprint.ICountFingerprint;
import org.openscience.cdk.fingerprint.IFingerprinter;
import org.openscience.cdk.graph.PathTools;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.ringsearch.AllRingsFinder;
import org.openscience.cdk.ringsearch.RingSearch;
import org.openscience.cdk.ringsearch.SSSRFinder;
import org.openscience.cdk.similarity.Tanimoto;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.RingSetManipulator;
import org.openscience.cdk.tools.periodictable.PeriodicTable;

/**
 * Chemical Fingerprint based on the atom neighbourhood information Generates a
 * fingerprint for a given {@link IAtomContainer}. Fingerprints are
 * one-dimensional bit arrays, where bits are set according to a the occurrence
 * of a particular structural feature (See for example the Daylight inc. theory
 * manual for more information). Fingerprints allow for a fast screening step to
 * exclude candidates for a substructure search in a database. They are also a
 * means for determining the similarity of chemical structures.
 *
 * <p>
 * A fingerprint is generated for an AtomContainer with this code:
 * <pre>
 * Molecule molecule = new Molecule();
 * IFingerprinter fingerprinter =
 * new FeatureFingerprinter();
 * BitSet fingerprint = fingerprinter.getFingerprint(molecule);
 * fingerprint.size(); // returns 512 by default
 * fingerprint.length(); // returns the highest set bit
 * </pre></p>
 *
 * <p>
 * The FingerPrinter assumes that hydrogens are explicitly given! Furthermore,
 * if pseudo atoms or atoms with malformed symbols are present, their atomic
 * number is taken as one more than the last element currently supported in
 * {@link PeriodicTable}.
 *
 * <p>
 * Unlike the {@link Fingerprinter}, this fingerprinter does not take into
 * account aromaticity but it differentiates between ring and non ring system.
 *
 * <p>
 * The integer hashing is done using the CRC32 algorithm, using the Java CRC32
 * class, which is the same formula/parameters as used by PNG files, and
 * described in:</p>
 *
 * <a href="http://www.w3.org/TR/PNG/#D-CRCAppendix">http://www.w3.org/TR/PNG/#D-CRCAppendix</a>
 *
 *
 * @cdk.keyword fingerprint
 * @cdk.keyword similarity
 * @cdk.module standard
 * @cdk.githash
 */
public class FeatureFingerprinter implements IFingerprinter {

    private static final long serialVersionUID = 96896986897971L;

    public static void main(String[] args) {
        String s1 = "N[C@@H](CC1=CC=C(O)C=C1)C(O)=O";
////        String s1 = "N[C@@H](CC1=CC=C(O)C=C1)C(*)=O";
        String s2 = "OC(=O)C(=O)CC1=CC=C(O)C=C1";

//        /*
////         Ring-1
////         */
//        String s1 = "C1CCCCC1";
////        /*
////         Ring-2
////         */
//        String s2 = "C1CCC2CCCCC2C1";
//////        /*
//////         Linear
//////         */
//        String s2 = "CCCCCC";
        final SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());

        IBitFingerprint bitFingerprint1 = null;
        try {
            IAtomContainer parseSmilesQ = sp.parseSmiles(s1);
            FeatureFingerprinter featureFingerprint = new FeatureFingerprinter();
            bitFingerprint1 = featureFingerprint.getBitFingerprint(parseSmilesQ);
            System.out.println("features " + bitFingerprint1.asBitSet());
        } catch (InvalidSmilesException ex) {
            Logger.getLogger(FeatureFingerprinter.class.getName()).log(Level.SEVERE, null, ex);
        } catch (CDKException ex) {
            Logger.getLogger(FeatureFingerprinter.class.getName()).log(Level.SEVERE, null, ex);
        }

        IBitFingerprint bitFingerprint2 = null;
        try {
            IAtomContainer parseSmilesQ = sp.parseSmiles(s2);
            FeatureFingerprinter featureFingerprint = new FeatureFingerprinter();
            bitFingerprint2 = featureFingerprint.getBitFingerprint(parseSmilesQ);
            System.out.println("features " + bitFingerprint2.asBitSet());
        } catch (InvalidSmilesException ex) {
            Logger.getLogger(FeatureFingerprinter.class.getName()).log(Level.SEVERE, null, ex);
        } catch (CDKException ex) {
            Logger.getLogger(FeatureFingerprinter.class.getName()).log(Level.SEVERE, null, ex);
        }

        try {
            if (bitFingerprint1 != null && bitFingerprint2 != null) {
                double calculate = Tanimoto.calculate(bitFingerprint1, bitFingerprint2);
                System.out.println("Tanimoto: " + calculate);
            }
        } catch (Exception ex) {
            Logger.getLogger(FeatureFingerprinter.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    /**
     * The default length of created fingerprints.
     */
    public final static int DEFAULT_SIZE = 512;

    private int fingerprintSize;
    private final boolean DEBUG = false;
    public final AllRingsFinder arf;

    /**
     * Constructs a fingerprint generator that creates fingerprints of the given
     * size, using a generation algorithm with the given search depth.
     *
     * @param size The desired size of the fingerprint
     */
    public FeatureFingerprinter(int size) {
        this.fingerprintSize = size;
        this.arf = new AllRingsFinder();
    }

    /**
     * Creates a fingerprint generator of length <code>DEFAULT_SIZE</code>.
     */
    public FeatureFingerprinter() {
        this(DEFAULT_SIZE);
    }

    @Override
    public IBitFingerprint getBitFingerprint(IAtomContainer container) throws CDKException {
        BitSet bitSet = new BitSet(fingerprintSize);
        for (int i = 0; i < bitSet.length(); i++) {
            bitSet.set(i, false);
        }

        try {
            IAtomContainer clonedContainer = container.clone();

            try {
                AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(clonedContainer);
            } catch (CDKException ex) {
                Logger.getLogger(FeatureGenerator.class.getName()).log(Level.SEVERE, null, ex);
            }

            /*
             Generate hashing information for atoms using connectivity information etc.
             */
            Map<IAtom, String> atomInvariantsMap = FeatureGenerator.getAtomInvariants(clonedContainer);
            /*
             Store all the atoms
             */
            for (IAtom a : atomInvariantsMap.keySet()) {
                int hashCode = a.getSymbol().hashCode();
                long b = hashCode >= 0 ? hashCode : ((hashCode & 0x7FFFFFFF) | (1L << 31));
                bitSet.set((int) (b % fingerprintSize));
            }
            /*
             Store all the atom invariants
             */
            for (String invariant : atomInvariantsMap.values()) {
                int hashCode = invariant.hashCode();
                if (DEBUG) {
                    System.out.println("invariant " + invariant);
                }
                long b = hashCode >= 0 ? hashCode : ((hashCode & 0x7FFFFFFF) | (1L << 31));
                bitSet.set((int) (b % fingerprintSize));
            }

            IRingSet rings = new RingSet();
            // sets SSSR information
            SSSRFinder finder = new SSSRFinder(container);
            IRingSet sssr = finder.findEssentialRings();
            rings.add(sssr);
            RingSetManipulator.sort(rings);
            setRingBits(bitSet, rings);

        } catch (CloneNotSupportedException exception) {
            throw new CDKException(
                    "Exception while cloning the input: " + exception.getMessage(),
                    exception);
        }
        return new BitSetFingerprint(bitSet);
    }

    private void setRingBits(BitSet bitSet, IRingSet rings) {
        Map<String, Integer> ringMap = new HashMap<>();
        for (IAtomContainer ring : rings.atomContainers()) {
            List<String> list = new ArrayList<>();
            for (IAtom a : ring.atoms()) {
                list.add(a.getSymbol());
            }
            Collections.sort(list, new NaturalOrderComparator());
            StringBuilder s = new StringBuilder();
            for (String f : list) {
                s.append(f);
            }
            s.trimToSize();
            String pattern = s.toString();
            if (ringMap.containsKey(pattern)) {
                Integer counter = ringMap.get(pattern) + 1;
                ringMap.put(pattern, counter);
            } else {
                ringMap.put(pattern, 1);
            }
        }

        for (String p : ringMap.keySet()) {
            int counter = ringMap.get(p);
            while (counter-- > 0) {
                int hashCode = p.concat(String.valueOf(counter)).hashCode();
                long b = hashCode >= 0 ? hashCode : ((hashCode & 0x7FFFFFFF) | (1L << 31));
                bitSet.set((int) (b % fingerprintSize));
            }
        }
    }

    @Override
    public ICountFingerprint getCountFingerprint(IAtomContainer iac) throws CDKException {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public Map<String, Integer> getRawFingerprint(IAtomContainer iac) throws CDKException {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public int getSize() {
        return fingerprintSize;
    }

}
