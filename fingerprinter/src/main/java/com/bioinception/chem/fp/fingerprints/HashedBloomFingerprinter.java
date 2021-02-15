/* $Revision$ $Author$ $Date$
 *
 * Copyright (C) 2002-2007  Christoph Steinbeck <steinbeck@users.sf.net>
 *               2011-2012  Syed Asad Rahman <asad@ebi.ac.uk>
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
package com.bioinception.chem.fp.fingerprints;

import com.bioinception.chem.fp.fingerprints.helper.BloomFilter;
import com.bioinception.chem.fp.fingerprints.helper.MoleculeWalker;
import com.bioinception.chem.fp.fingerprints.helper.RandomNumber;
import static com.bioinception.chem.fp.fingerprints.helper.RandomNumber.generateMersenneTwisterRandomNumber;
import com.bioinception.chem.fp.fingerprints.interfaces.IFingerprinter;
import com.bioinception.chem.fp.fingerprints.interfaces.IWalker;
import java.util.BitSet;
import java.util.Map;
import org.apache.commons.lang3.builder.HashCodeBuilder;
import org.openscience.cdk.RingSet;
import org.openscience.cdk.aromaticity.CDKHueckelAromaticityDetector;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.ringsearch.AllRingsFinder;
import org.openscience.cdk.ringsearch.SSSRFinder;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.RingSetManipulator;
import org.openscience.cdk.fingerprint.BitSetFingerprint;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.fingerprint.ICountFingerprint;

/*
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk> 2007-2011
 *
 */
public class HashedBloomFingerprinter extends RandomNumber implements IFingerprinter {

    /**
     * The default length of created fingerprints.
     */
    public final static int DEFAULT_SIZE = 1024;
    /**
     * The default search depth used to create the fingerprints.
     */
    public final static int DEFAULT_SEARCH_DEPTH = 8;
    private int fingerPrintSize;
    private BloomFilter<String> bloomFilter;
    private int ringBitCount;
    private boolean respectRingMatches;
    private int searchDepth;
    static int debugCounter = 0;
    // do all ring perception
    private AllRingsFinder arf;
    private static ILoggingTool logger
            = LoggingToolFactory.createLoggingTool(HashedBloomFingerprinter.class);

    /**
     * Creates a fingerprint generator of length <code>DEFAULT_SIZE</code> and
     * with a search depth of <code>DEFAULT_SEARCH_DEPTH</code>.
     */
    public HashedBloomFingerprinter() {
        this(DEFAULT_SIZE, DEFAULT_SEARCH_DEPTH);
    }

    public HashedBloomFingerprinter(int size) {
        this(size, DEFAULT_SEARCH_DEPTH);
    }

    /**
     * Constructs a fingerprint generator that creates fingerprints of the given
     * fingerPrintSize, using a generation algorithm with the given search
     * depth.
     *
     * @param fingerPrintSize The desired fingerPrintSize of the fingerprint
     * @param searchDepth The desired depth of search
     */
    public HashedBloomFingerprinter(int fingerPrintSize, int searchDepth) {
        this.searchDepth = searchDepth;
        this.respectRingMatches = false;
        this.ringBitCount = 10;
        this.arf = new AllRingsFinder();
        setFingerprintLength(fingerPrintSize);
    }

    public int getRingBitCount() {
        return ringBitCount;
    }

    public void reserveRingBits(int reserveBits) {
        if (reserveBits > bloomFilter.getBitArraySize()) {
            throw new IllegalStateException("Attempting to set more ring bits than total bits.");
        }
        this.ringBitCount = reserveBits;
        this.bloomFilter = new BloomFilter<>(fingerPrintSize - getRingBitCount());
    }

    private void setFingerprintLength(int length) {
        if (length < ringBitCount) {
            throw new IllegalStateException("Attempting to set fingerprint length below reserved ring bit count " + ringBitCount);
        }
        this.fingerPrintSize = length;
        this.bloomFilter = new BloomFilter<>(fingerPrintSize - getRingBitCount());
    }

    public int getFingerprintLength() {
        return bloomFilter.getBitArraySize() + ringBitCount;
    }

    /**
     * Generates a fingerprint of the default fingerPrintSize for the given
     * AtomContainer.
     *
     * @TestMethod("testGetFingerprint_IAtomContainer")
     * @param container The AtomContainer for which a Fingerprint is generated
     * @param ringFinder An instance of
     * {@link org.openscience.cdk.ringsearch.AllRingsFinder}
     * @exception CDKException if there is a timeout in ring or aromaticity
     * perception
     * @return A {@link BitSet} representing the fingerprint
     */
    public IBitFingerprint getFingerprint(IAtomContainer container,
            AllRingsFinder ringFinder)
            throws CDKException {
        if (ringFinder != null) {
            this.arf = ringFinder;
        }
        logger.debug("Entering Fingerprinter");
        logger.debug("Starting Aromaticity Detection");
        long before = System.currentTimeMillis();
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(container);
        CDKHueckelAromaticityDetector.detectAromaticity(container);
        long after = System.currentTimeMillis();
        logger.debug("time for aromaticity calculation: "
                + (after - before) + " milliseconds");
        logger.debug("Finished Aromaticity Detection");
        findPaths(container, searchDepth);
        return new BitSetFingerprint(generateFingerprint(container));
    }

    /**
     * Generates a fingerprint of the default fingerPrintSize for the given
     * AtomContainer.
     *
     * @TestMethod("testGetFingerprint_IAtomContainer")
     * @param container The AtomContainer for which a Fingerprint is generated
     * @return
     * @throws CDKException
     */
    @Override
    public IBitFingerprint getBitFingerprint(IAtomContainer container)
            throws CDKException {
        return getFingerprint(container, null);
    }

    /**
     * {@inheritDoc}
     *
     * @param atomContainer
     * @return
     * @throws CDKException
     */
    @Override
    public Map<String, Integer> getRawFingerprint(IAtomContainer atomContainer) throws CDKException {
        throw new UnsupportedOperationException();
    }

    /**
     * Get all paths of lengths 0 to the specified length.
     *
     * This method will find all paths upto length N starting from each atom in
     * the molecule and return the unique set of such paths.
     *
     * @param container The molecule to search
     * @param searchDepth The maximum path length desired
     */
    protected void findPaths(IAtomContainer container, int searchDepth) {
        IWalker walker = new MoleculeWalker(searchDepth, container);
        // convert paths to BitSet
        bloomFilter.addAll(walker.getPaths());
    }

    /**
     * @TestMethod("testGetSize")
     * @TestMethod("testGetSearchDepth")
     * @return
     */
    @Override
    public int getSearchDepth() {
        return searchDepth;
    }

    @Override
    public int getSize() {
        return fingerPrintSize;
    }

    /**
     * Should match rings to rings and non-rings to non-rings
     *
     * @return the respect ring matches
     */
    @Override
    public boolean isRespectRingMatches() {
        return respectRingMatches;
    }

    /**
     * Ring matches are allowed and non-ring to ring matches are discarded
     *
     * @param respectRingMatches respect the ring-to-ring matches and discard
     * non-ring to ring matches
     */
    @Override
    public void setRespectRingMatches(boolean respectRingMatches) {
        this.respectRingMatches = respectRingMatches;
    }

    private BitSet generateFingerprint(IAtomContainer container) {
        BitSet walkBits = bloomFilter.toBitSet();
        BitSet result = new BitSet(getFingerprintLength());
        result.or(walkBits);
        if (isRespectRingMatches()) {
            IRingSet rings = new RingSet();
            IRingSet allRings;
            try {
                allRings = arf.findAllRings(container);
                rings.add(allRings);
            } catch (CDKException e) {
                logger.debug(e.toString());
            }

            // sets SSSR information
            SSSRFinder finder = new SSSRFinder(container);
            IRingSet sssr = finder.findEssentialRings();
            rings.add(sssr);
            RingSetManipulator.markAromaticRings(rings);
            RingSetManipulator.sort(rings);
            setRingBits(result, rings);
        }
        bloomFilter.clear();
        return result;
    }

    private void setRingBits(BitSet bitset, IRingSet rings) {
        int ringSize = 0;
        for (IAtomContainer ring : rings.atomContainers()) {
            int atomCount = ring.getAtomCount();
            if (atomCount < ringBitCount) {
                if (ringSize < atomCount) {
                    int toHashCode = new HashCodeBuilder(17, 37).append(atomCount).toHashCode();
                    int ringPosition = (int) generateMersenneTwisterRandomNumber(ringBitCount, toHashCode);
                    int index = bloomFilter.getBitArraySize() + (ringPosition - 2);
                    if (index < getFingerprintLength()) {
                        bitset.set(index);
                    }
                    ringSize++;
                }
            }
        }
    }

    @Override
    public IBitFingerprint getBitFingerprint(IAtomContainer container, AllRingsFinder ringFinder) throws CDKException {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public ICountFingerprint getCountFingerprint(IAtomContainer iac) throws CDKException {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public boolean isRespectFormalCharges() {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public void setRespectFormalCharges(boolean respectFormalCharges) {
        throw new UnsupportedOperationException("Not supported yet.");
    }

    @Override
    public String getVersionDescription() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public BitSet getFingerprint(IAtomContainer iac) throws CDKException {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
}
