/* $Revision$ $Author$ $Date$
 *
 * Copyright (C) 2002-2007  Christoph Steinbeck <steinbeck@users.sf.net>
 *               2011       Syed Asad Rahman <asad@ebi.ac.uk>
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
package fingerprints;

import fingerprints.interfaces.IFingerprinter;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import org.openscience.cdk.annotations.TestClass;
import org.openscience.cdk.annotations.TestMethod;
import org.openscience.cdk.aromaticity.CDKHueckelAromaticityDetector;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IRingSet;
import org.openscience.cdk.ringsearch.AllRingsFinder;
import org.openscience.cdk.ringsearch.SSSRFinder;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.RingSetManipulator;
import fingerprints.helper.MoleculeWalker;
import fingerprints.helper.RandomNumber;
import fingerprints.interfaces.IWalker;
import java.util.*;
import org.openscience.cdk.fingerprint.BitSetFingerprint;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.fingerprint.ICountFingerprint;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.interfaces.IAtomContainerSet;

/**
 * Generates a fingerprint for a given AtomContainer. Fingerprints are one-dimensional bit arrays, where bits are set
 * according to a the occurrence of a particular structural feature (See for example the Daylight inc. theory manual for
 * more information). Fingerprints allow for a fast screening step to exclude candidates for a substructure search in a
 * database. They are also a means for determining the similarity of chemical structures. <p>
 *
 * A fingerprint is generated for an AtomContainer with this code:
 * <pre>
 *   AtomContainer molecule = new AtomContainer();
 *   IFingerprinter fingerprinter = new HashedFingerprinter();
 *   BitSet fingerprint = fingerprinter.getFingerprint(molecule);
 *   This will match ring system with rings.
 *   fingerprinter.setRespectRingMatches(true);
 *   fingerprint.fingerprintLength(); // returns 1024 by default
 *   fingerprint.length(); // returns the highest set bit
 * </pre> <p>
 *
 * The FingerPrinter assumes that hydrogens are explicitly given!
 *
 * <font color="#FF0000">Warning: The aromaticity detection for this FingerPrinter relies on AllRingsFinder, which is
 * known to take very long for some molecules with many cycles or special cyclic topologies. Thus, the AllRingsFinder
 * has a built-in timeout of 5 seconds after which it aborts and throws an Exception. If you want your SMILES generated
 * at any expense, you need to create your own AllRingsFinder, set the timeout to a higher value, and assign it to this
 * FingerPrinter. In the vast majority of cases, however, the defaults will be fine. </font> <p>
 *
 * <font color="#FF0000">Another Warning : The daylight manual says: "Fingerprints are not so definite: if a fingerprint
 * indicates a pattern is missing then it certainly is, but it can only indicate a pattern's presence with some
 * probability." In the case of very small molecules, the probability that you get the same fingerprint for different
 * molecules is high. </font> </p>
 *
 * @author Syed Asad Rahman (2011), Christoph Steinbeck (2002-2007) @cdk.created 07-11-2011 @cdk.keyword fingerprint
 * @cdk.keyword similarity @cdk.module standard @cdk.githash
 */
@TestClass("org.openscience.cdk.fingerprint.FingerprinterTest")
public class HashedFingerprinter extends RandomNumber implements IFingerprinter {

    /**
     * The default length of created fingerprints.
     */
    public final static int DEFAULT_SIZE = 1024;
    /**
     * The default search depth used to create the fingerprints.
     */
    public final static int DEFAULT_SEARCH_DEPTH = 8;
    private int fingerprintLength;
    private boolean respectRingMatches;
    private boolean respectFormalCharges;
    private boolean respectStereoAssignments;
    private int searchDepth;
    private static ILoggingTool logger =
            LoggingToolFactory.createLoggingTool(HashedFingerprinter.class);
    private AllRingsFinder arf;

    /**
     * Creates a fingerprint generator of length
     * <code>DEFAULT_SIZE</code> and with a search depth of
     * <code>DEFAULT_SEARCH_DEPTH</code>.
     */
    public HashedFingerprinter() {
        this(DEFAULT_SIZE, DEFAULT_SEARCH_DEPTH);
    }

    public HashedFingerprinter(int size) {
        this(size, DEFAULT_SEARCH_DEPTH);
    }

    /**
     * Constructs a fingerprint generator that creates fingerprints of the given fingerprintLength, using a generation
     * algorithm with the given search depth.
     *
     * @param fingerprintLength The desired fingerprintLength of the fingerprint
     * @param searchDepth The desired depth of search
     */
    public HashedFingerprinter(int fingerprintLength, int searchDepth) {
        this.fingerprintLength = fingerprintLength;
        this.searchDepth = searchDepth;
        this.respectRingMatches = false;
        this.respectFormalCharges = false;
        this.respectStereoAssignments = false;
        this.arf = new AllRingsFinder();
    }

    /**
     * Generates a fingerprint of the default fingerprintLength for the given AtomContainer.
     *
     * @param container The AtomContainer for which a Fingerprint is generated
     * @param ringFinder An instance of
     *                   {@link org.openscience.cdk.ringsearch.AllRingsFinder}
     * @exception CDKException if there is a timeout in ring or aromaticity perception
     * @return A {@link BitSet} representing the fingerprint
     */
    @TestMethod("testGetFingerprint_IAtomContainer")
    public IBitFingerprint getBitFingerprint(
            IAtomContainer container,
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
        BitSet bitSet = new BitSet(fingerprintLength);

        if (!ConnectivityChecker.isConnected(container)) {
            IAtomContainerSet partitionedMolecules = ConnectivityChecker.partitionIntoMolecules(container);
            for (IAtomContainer atomContainer : partitionedMolecules.atomContainers()) {
                addUniquePath(atomContainer, bitSet);
            }
        } else {
            addUniquePath(container, bitSet);
        }
        return new BitSetFingerprint(bitSet);
    }

    private void addUniquePath(IAtomContainer container, BitSet bitSet) {
        Integer[] hashes = findPaths(container, searchDepth);
        for (Integer hash : hashes) {
            int position = (int) generateMersenneTwisterRandomNumber(fingerprintLength, hash.intValue());
            bitSet.set(position);
        }
    }

    /**
     * Generates a fingerprint of the default fingerprintLength for the given AtomContainer.
     *
     * @param container The AtomContainer for which a Fingerprint is generated
     * @return
     * @throws CDKException
     */
    @TestMethod("testGetFingerprint_IAtomContainer")
    @Override
    public IBitFingerprint getBitFingerprint(IAtomContainer container)
            throws CDKException {
        return getBitFingerprint(container, null);
    }

    @Override
    public Map<String, Integer> getRawFingerprint(IAtomContainer atomContainer) throws CDKException {
        Map<String, Integer> uniquePaths = new TreeMap<String, Integer>();
        if (!ConnectivityChecker.isConnected(atomContainer)) {
            IAtomContainerSet partitionedMolecules = ConnectivityChecker.partitionIntoMolecules(atomContainer);
            for (IAtomContainer container : partitionedMolecules.atomContainers()) {
                addUniquePaths(container, uniquePaths);
            }
        } else {
            addUniquePaths(atomContainer, uniquePaths);
        }
        return uniquePaths;
    }

    private void addUniquePaths(IAtomContainer atomContainer, Map<String, Integer> uniquePaths) {
        Integer[] hashes = findPaths(atomContainer, searchDepth);
        for (Integer hash : hashes) {
            int position = (int) generateMersenneTwisterRandomNumber(fingerprintLength, hash.intValue());
            uniquePaths.put(new Integer(position).toString(), hash);
        }
    }

    /**
     * Get all paths of lengths 0 to the specified length.
     *
     * This method will find all paths upto length N starting from each atom in the molecule and return the unique set
     * of such paths.
     *
     * @param container The molecule to search
     * @param searchDepth The maximum path length desired
     * @return A map of path strings, keyed on themselves
     */
    protected Integer[] findPaths(IAtomContainer container, int searchDepth) {

        IWalker walker = new MoleculeWalker(searchDepth, container);
        // convert paths to hashes
        List<Integer> paths = new ArrayList<Integer>();
        int patternIndex = 0;

        for (String s : walker.getPaths()) {
            int toHashCode = s.hashCode();
            paths.add(patternIndex, toHashCode);
            patternIndex++;
        }

        if (isRespectRingMatches()) {
            SSSRFinder finder = new SSSRFinder(container);
            IRingSet sssr = finder.findEssentialRings();
            RingSetManipulator.sort(sssr);
            int ringCounter = sssr.getAtomContainerCount();
            for (Iterator<IAtomContainer> it = sssr.atomContainers().iterator(); it.hasNext();) {
                IAtomContainer ring = it.next();
                int toHashCode = String.valueOf(ringCounter * ring.getAtomCount()).hashCode();
                paths.add(patternIndex, toHashCode);
                patternIndex++;
                ringCounter--;
            }
        }

        if (isRespectFormalCharges()) {
            List<String> l = new ArrayList<String>();
            for (Iterator<IAtom> it = container.atoms().iterator(); it.hasNext();) {
                IAtom atom = it.next();
                int charge = atom.getFormalCharge() == null ? 0 : atom.getFormalCharge().intValue();
                if (charge != 0) {
                    l.add(atom.getSymbol().concat(String.valueOf(charge)));
                }
            }
            Collections.sort(l);
            int toHashCode = l.hashCode();
            paths.add(patternIndex, toHashCode);
            patternIndex++;
        }

        if (isRespectStereoAssignments()) {
            List<String> l = new ArrayList<String>();
            /*
             * atom stereo parity
             */
            for (Iterator<IAtom> it = container.atoms().iterator(); it.hasNext();) {
                IAtom atom = it.next();
                int st = atom.getStereoParity() == null ? 0 : atom.getStereoParity().intValue();
                if (st != 0) {
                    l.add(atom.getSymbol().concat(String.valueOf(st)));
                }
            }
            Collections.sort(l);
            int toHashCode = l.hashCode();
            paths.add(patternIndex, toHashCode);
            patternIndex++;
        }

        if (container.getSingleElectronCount() > 0) {
            StringBuilder radicalInformation = new StringBuilder();
            radicalInformation.append("RAD: ".concat(String.valueOf(container.getSingleElectronCount())));
            paths.add(patternIndex, radicalInformation.toString().hashCode());
            patternIndex++;
        }
        if (container.getLonePairCount() > 0) {
            StringBuilder lpInformation = new StringBuilder();
            lpInformation.append("LP: ".concat(String.valueOf(container.getLonePairCount())));
            paths.add(patternIndex, lpInformation.toString().hashCode());
            patternIndex++;
        }
        return paths.toArray(new Integer[paths.size()]);
    }

    /**
     *
     * @return
     */
    @TestMethod("testGetSearchDepth")
    @Override
    public int getSearchDepth() {
        return searchDepth;
    }

    @TestMethod("testGetSize")
    @Override
    public int getSize() {
        return fingerprintLength;
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
     * @param respectRingMatches respect the ring-to-ring matches and discard non-ring to ring matches
     */
    @Override
    public void setRespectRingMatches(boolean respectRingMatches) {
        this.respectRingMatches = respectRingMatches;
    }

    /**
     * @return the respectFormalCharges
     */
    @Override
    public boolean isRespectFormalCharges() {
        return respectFormalCharges;
    }

    /**
     * @param respectFormalCharges the flag to set if formal charge is checked
     */
    @Override
    public void setRespectFormalCharges(boolean respectFormalCharges) {
        this.respectFormalCharges = respectFormalCharges;
    }

    /**
     * @return the respectStereoAssignments
     */
    public boolean isRespectStereoAssignments() {
        return respectStereoAssignments;
    }

    /**
     * @param respectStereoAssignments the respectStereoAssignments to set
     */
    public void setRespectStereoAssignments(boolean respectStereoAssignments) {
        this.respectStereoAssignments = respectStereoAssignments;
    }

    @Override
    public ICountFingerprint getCountFingerprint(IAtomContainer iac) throws CDKException {
        throw new UnsupportedOperationException("Not supported yet.");
    }
}
