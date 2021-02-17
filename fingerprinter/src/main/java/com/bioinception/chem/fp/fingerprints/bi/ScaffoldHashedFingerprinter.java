/* 
 * Copyright (C) 2002-2007  Christoph Steinbeck <steinbeck@users.sf.net>
 *               2020-2021  Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>

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
package com.bioinception.chem.fp.fingerprints.bi;

import static com.bioinception.chem.fp.fingerprints.bi.PathEncoder.encodePaths;
import static com.bioinception.chem.fp.fingerprints.bi.PathEncoder.isPseudoAtom;
import static com.bioinception.chem.fp.fingerprints.helper.RandomNumber.generateMersenneTwisterRandomNumber;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.ringsearch.AllRingsFinder;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.util.AbstractMap.SimpleImmutableEntry;
import java.util.Arrays;
import java.util.BitSet;
import java.util.List;
import java.util.Map;
import org.apache.commons.lang3.builder.HashCodeBuilder;
import org.openscience.cdk.fingerprint.AbstractFingerprinter;
import org.openscience.cdk.fingerprint.BitSetFingerprint;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.fingerprint.ICountFingerprint;
import org.openscience.cdk.exception.Intractable;
import org.openscience.cdk.fingerprint.IFingerprinter;
import org.openscience.cdk.graph.CycleFinder;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.tools.manipulator.RingSetManipulator;

/**
 * Generates a fingerprint for a given AtomContainer. Fingerprints are
 * one-dimensional bit arrays, where bits are set according to a the occurrence
 * of a particular structural feature (See for example the Daylight inc. theory
 * manual for more information). Fingerprints allow for a fast screening step to
 * exclude candidates for a substructure search in a database. They are also a
 * means for determining the similarity of chemical structures.
 *
 * This sets
 *
 * - 16 rings mcb/sssr
 *
 * - 128 rings path
 *
 * - 128 path 0-3
 *
 * - 256 path 0-5
 *
 * - 624 path 0-7
 *
 * <p>
 *
 * A fingerprint is generated for an AtomContainer with this code:
 * <pre>
 * Molecule molecule = new Molecule();
 * IFingerprinter fingerprinter = new ScaffoldHashedFingerprinter();
 * IBitFingerprint fingerprint = fingerprinter.getBitFingerprint(molecule);
 * fingerprint.size(); // returns 1024 by default
 * fingerprint.length(); // returns the highest set bit
 * </pre>
 * <p>
 *
 * The FingerPrinter assumes that hydrogen's are explicitly given! Furthermore,
 * if pseudo atoms or atoms with malformed symbols are present, their atomic
 * number is taken as one more than the last element currently supported in
 * {@link org.openscience.cdk.tools.periodictable.PeriodicTable}.
 *
 * <font color="#FF0000">Warning: The aromaticity detection for this
 * FingerPrinter relies on AllRingsFinder, which is known to take very long for
 * some molecules with many cycles or special cyclic topologies. Thus, the
 * AllRingsFinder has a built-in timeout of 5 seconds after which it aborts and
 * throws an Exception. If you want your SMILES generated at any expense, you
 * need to create your own AllRingsFinder, set the timeout to a higher value,
 * and assign it to this FingerPrinter. In the vast majority of cases, however,
 * the defaults will be fine. </font>
 * <p>
 *
 * <font color="#FF0000">Another Warning : The daylight manual says:
 * "Fingerprints are not so definite: if a fingerprint indicates a pattern is
 * missing then it certainly is, but it can only indicate a pattern's presence
 * with some probability." In the case of very small molecules, the probability
 * that you get the same fingerprint for different molecules is high. </font>
 * </p>
 *
 */
public class ScaffoldHashedFingerprinter extends AbstractFingerprinter implements IFingerprinter {

    /**
     * Throw an exception if too many paths (per atom) are generated.
     */
    private final static int DEFAULT_PATH_LIMIT = 42000;

    /**
     * The default length of created fingerprints.
     */
    public final static int DEFAULT_SIZE = 1024;
    /**
     * The default search depth used to create the fingerprints.
     */
    public final static int DEFAULT_SEARCH_DEPTH = 7;

    private int size;
    private int searchDepth;
    private int pathLimit = DEFAULT_PATH_LIMIT;

    private boolean hashPseudoAtoms = false;

    static int debugCounter = 0;

    private static ILoggingTool logger = LoggingToolFactory
            .createLoggingTool(ScaffoldHashedFingerprinter.class);

    /**
     * Creates a fingerprint generator of length <code>DEFAULT_SIZE</code> and
     * with a search depth of <code>DEFAULT_SEARCH_DEPTH</code>.
     */
    public ScaffoldHashedFingerprinter() {
        this(DEFAULT_SIZE, DEFAULT_SEARCH_DEPTH);
    }

    public ScaffoldHashedFingerprinter(int size) {
        this(size, DEFAULT_SEARCH_DEPTH);
    }

    /**
     * Constructs a fingerprint generator that creates fingerprints of the given
     * size, using a generation algorithm with the given search depth.
     *
     * @param size The desired size of the fingerprint
     * @param searchDepth The desired depth of search (number of bonds)
     */
    public ScaffoldHashedFingerprinter(int size, int searchDepth) {
        this.size = size;
        this.searchDepth = searchDepth;

    }

    @Override
    protected List<Map.Entry<String, String>> getParameters() {
        return Arrays.<Map.Entry<String, String>>asList(
                new SimpleImmutableEntry<>("searchDepth", Integer.toString(searchDepth)),
                new SimpleImmutableEntry<>("pathLimit", Integer.toString(pathLimit)),
                new SimpleImmutableEntry<>("hashPseudoAtoms", Boolean.toString(hashPseudoAtoms))
        );
    }

    /**
     * Generates a fingerprint of the default size for the given AtomContainer.
     *
     * @param container The AtomContainer for which a Fingerprint is generated
     * @param ringFinder An instance of
     * {@link org.openscience.cdk.ringsearch.AllRingsFinder}
     * @exception CDKException if there is a timeout in ring or aromaticity
     * perception
     * @return A {@link BitSet} representing the fingerprint
     */
    public IBitFingerprint getBitFingerprint(IAtomContainer container, AllRingsFinder ringFinder) throws CDKException {
        logger.debug("Entering Fingerprinter");
        logger.debug("Starting Aromaticity Detection");
        long before = System.currentTimeMillis();

        if (!hasPseudoAtom(container.atoms())) {
            AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(container);
            Aromaticity.cdkLegacy().apply(container);
        }
        long after = System.currentTimeMillis();
        logger.debug("time for aromaticity calculation: " + (after - before) + " milliseconds");
        logger.debug("Finished Aromaticity Detection");

        // all cycles or relevant or essential
//        CycleFinder cf = Cycles.or(Cycles.all(),
//                Cycles.or(Cycles.relevant(),
//                        Cycles.essential()));
        //minimal cycle basis
        CycleFinder cf = Cycles.mcb();

        IRingSet rings = null;
        try {
            Cycles cycles = cf.find(container);
            rings = cycles.toRingSet();
            RingSetManipulator.markAromaticRings(rings);
            RingSetManipulator.sort(rings);
        } catch (Intractable e) {
            // ignore error - edge short cycles do not check tractability
        }
        /*
         * Encode Rings
         */
        int size0 = 16;
        BitSet bitSet0 = new BitSet(size0);
        if (rings != null) {
            setRingBits(bitSet0, rings, size0);
        }
//        System.out.println("BitSet - 0 " + bitSet0);

        /*
         * Encode Rings Path
         */
        int size1 = 128;
        BitSet bitSet1 = new BitSet(size1);
        if (rings != null) {
            for (int i = 0; i < rings.getAtomContainerCount(); i++) {
                IAtomContainer ring = rings.getAtomContainer(i);
                encodePaths(ring, 0, searchDepth, bitSet1, size1, pathLimit, hashPseudoAtoms);
            }
        }

//         System.out.println("BitSet - 1 " + bitSet1);

        /*
         * Encode Atoms
         */
        int size2 = 128;
        BitSet bitSet2 = new BitSet(size2);
        encodePaths(container, 0, 3, bitSet2, size2, pathLimit, hashPseudoAtoms);
//        System.out.println("BitSet - 2 " + bitSet2);
        /*
         * Encode Paths 1 to 5 length
         */
        int size3 = 256;
        BitSet bitSet3 = new BitSet(size3);
        encodePaths(container, 0, 5, bitSet3, size3, pathLimit, hashPseudoAtoms);
//        System.out.println("BitSet - 3 " + bitSet3);

        /*
         * Encode Paths 1 to search depth length
         */
        int size4 = size - (size3 + size2 + size1 + size0);
        BitSet bitSet4 = new BitSet(size4);
        encodePaths(container, 0, searchDepth, bitSet4, size4, pathLimit, hashPseudoAtoms);
//        System.out.println("BitSet - 4 " + bitSet4);

        /*
         * Set all bits
         */
        BitSet concatenate_vectors = concatenate_vectors(bitSet1, bitSet0);
        concatenate_vectors = concatenate_vectors(bitSet2, concatenate_vectors);
        concatenate_vectors = concatenate_vectors(bitSet3, concatenate_vectors);
        concatenate_vectors = concatenate_vectors(bitSet4, concatenate_vectors);
//        System.out.println("Concat BitSet " + concatenate_vectors);

        BitSet bitSet = new BitSet(size);
        bitSet.or(concatenate_vectors);
//        System.out.println("BitSet: " + bitSet);
        return new BitSetFingerprint(bitSet);
    }

    private void setRingBits(BitSet bitset, IRingSet rings, int maxRingSize) {
//        System.out.println("Rings " + rings.getAtomContainerCount());
        int ringSize = 0;
        for (IAtomContainer ring : rings.atomContainers()) {
            int atomCount = ring.getAtomCount();
//            System.out.println("Ring size " + atomCount);
            if (atomCount < maxRingSize) {
                if (ringSize < atomCount) {
//                    System.out.println(ringSize + ", Ring size " + atomCount);
                    int toHashCode = new HashCodeBuilder(17, 37).append(atomCount).toHashCode();
                    int ringPosition = (int) generateMersenneTwisterRandomNumber(maxRingSize, toHashCode);
                    bitset.set(ringPosition);
                    ringSize++;
                }
            }
        }
    }

    /**
     * Generates a fingerprint of the default size for the given AtomContainer.
     *
     * @param container The AtomContainer for which a Fingerprint is generated
     * @return
     * @throws org.openscience.cdk.exception.CDKException
     */
    @Override
    public IBitFingerprint getBitFingerprint(IAtomContainer container) throws CDKException {
        return getBitFingerprint(container, null);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public Map<String, Integer> getRawFingerprint(IAtomContainer iAtomContainer) throws CDKException {
        throw new UnsupportedOperationException();
    }

    /**
     *
     * @param limit
     */
    public void setPathLimit(int limit) {
        this.pathLimit = limit;
    }

    /**
     *
     * @param value
     */
    public void setHashPseudoAtoms(boolean value) {
        this.hashPseudoAtoms = value;
    }

    /**
     *
     * @return
     */
    public int getSearchDepth() {
        return searchDepth;
    }

    @Override
    public int getSize() {
        return size;
    }

    @Override
    public ICountFingerprint getCountFingerprint(IAtomContainer container) throws CDKException {
        throw new UnsupportedOperationException();
    }

    private static boolean hasPseudoAtom(Iterable<IAtom> path) {
        for (IAtom atom : path) {
            if (isPseudoAtom(atom)) {
                return true;
            }
        }
        return false;
    }

    /**
     * Concatenates two bitSet objects
     *
     * @param vector_1_in
     * @param vector_2_in
     * @return
     */
    public static BitSet concatenate_vectors(BitSet vector_1_in, BitSet vector_2_in) {
        BitSet vector_1_in_clone = (BitSet) vector_1_in.clone();
        BitSet vector_2_in_clone = (BitSet) vector_2_in.clone();

        BitSet vectorConcat = new BitSet(vector_1_in.length() + vector_2_in.length());

        int index = -1;
        while (index < (vector_2_in_clone.length() - 1)) {
            index = vector_2_in_clone.nextSetBit((index + 1));
            vectorConcat.set(index);
        }

        index = -1;
        while (index < (vector_1_in_clone.length() - 1)) {
            index = vector_1_in_clone.nextSetBit((index + 1));
            vectorConcat.set(vector_2_in_clone.length() + index);
        }

        return vectorConcat;
    }

}
