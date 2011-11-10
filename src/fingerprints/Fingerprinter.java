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

import helper.RandomNumber;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import org.apache.commons.lang3.builder.HashCodeBuilder;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.annotations.TestClass;
import org.openscience.cdk.annotations.TestMethod;
import org.openscience.cdk.aromaticity.CDKHueckelAromaticityDetector;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.PathTools;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.ringsearch.AllRingsFinder;
import org.openscience.cdk.tools.ILoggingTool;
import org.openscience.cdk.tools.LoggingToolFactory;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.periodictable.PeriodicTable;

/**
 *  Generates a fingerprint for a given AtomContainer. Fingerprints are
 *  one-dimensional bit arrays, where bits are set according to a the 
 *  occurrence of a particular structural feature (See for example the 
 *  Daylight inc. theory manual for more information). Fingerprints allow for 
 *  a fast screening step to exclude candidates for a substructure search in a 
 *  database. They are also a means for determining the similarity of chemical 
 *  structures. <p>
 *
 *  A fingerprint is generated for an AtomContainer with this code: <pre>
 *   Molecule molecule = new Molecule();
 *   IFingerprinter fingerprinter = new Fingerprinter();
 *   BitSet fingerprint = fingerprinter.getFingerprint(molecule);
 *   This will match ring system with rings.
 *   fingerprinter.setRespectRingMatches(true);
 *   fingerprint.size(); // returns 1024 by default
 *   fingerprint.length(); // returns the highest set bit
 * </pre> <p>
 *
 *  The FingerPrinter assumes that hydrogens are explicitly given!
 *
 *  <font color="#FF0000">Warning: The aromaticity detection for this
 *  FingerPrinter relies on AllRingsFinder, which is known to take very long
 *  for some molecules with many cycles or special cyclic topologies. Thus, 
 *  the AllRingsFinder has a built-in timeout of 5 seconds after which it 
 *  aborts and throws an Exception. If you want your SMILES generated at any 
 *  expense, you need to create your own AllRingsFinder, set the timeout to a 
 *  higher value, and assign it to this FingerPrinter. In the vast majority of 
 *  cases, however, the defaults will be fine. </font> <p>
 *
 *  <font color="#FF0000">Another Warning : The daylight manual says:
 *  "Fingerprints are not so definite: if a fingerprint indicates a pattern is
 *  missing then it certainly is, but it can only indicate a pattern's presence
 *  with some probability." In the case of very small molecules, the 
 *  probability that you get the same fingerprint for different molecules is 
 *  high. </font>
 *  </p>
 *
 * @author         Syed Asad Rahman (2011), Christoph Steinbeck (2002-2007)
 * @cdk.created    07-11-2011
 * @cdk.keyword    fingerprint
 * @cdk.keyword    similarity
 * @cdk.module     standard
 * @cdk.githash
 */
@TestClass("org.openscience.cdk.fingerprint.FingerprinterTest")
public class Fingerprinter implements IFingerprinter {

    private int size;
    private boolean respectRingMatches;
    private int searchDepth;
    static int debugCounter = 0;
    private static int[] hashes;
    private static ILoggingTool logger =
            LoggingToolFactory.createLoggingTool(Fingerprinter.class);
    private static final Map<String, String> patterns = new HashMap<String, String>() {

        private static final long serialVersionUID = 3348458944893841L;

        {
            put("R", "*");
            put("X", "**");
        }
    };

    /**
     * Creates a fingerprint generator of length <code>DEFAULT_SIZE</code>
     * and with a search depth of <code>DEFAULT_SEARCH_DEPTH</code>.
     */
    public Fingerprinter() {
        this(DEFAULT_SIZE, DEFAULT_SEARCH_DEPTH);
    }

    public Fingerprinter(int size) {
        this(size, DEFAULT_SEARCH_DEPTH);
    }

    /**
     * Constructs a fingerprint generator that creates fingerprints of
     * the given size, using a generation algorithm with the given search
     * depth.
     *
     * @param  size        The desired size of the fingerprint
     * @param  searchDepth The desired depth of search
     */
    public Fingerprinter(int size, int searchDepth) {
        this.size = size;
        this.searchDepth = searchDepth;
        this.respectRingMatches = false;
    }

    /**
     * Generates a fingerprint of the default size for the given AtomContainer.
     *
     * @param container The AtomContainer for which a Fingerprint is generated
     * @param ringFinder An instance of 
     *                   {@link org.openscience.cdk.ringsearch.AllRingsFinder}
     * @exception CDKException if there is a timeout in ring or aromaticity 
     *                         perception
     * @return A {@link BitSet} representing the fingerprint
     */
    @TestMethod("testGetFingerprint_IAtomContainer")
    @Override
    public BitSet getFingerprint(IAtomContainer container,
            AllRingsFinder ringFinder)
            throws CDKException {

        int position = -1;
        logger.debug("Entering Fingerprinter");
        logger.debug("Starting Aromaticity Detection");
        long before = System.currentTimeMillis();
        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(container);
        CDKHueckelAromaticityDetector.detectAromaticity(container);
        long after = System.currentTimeMillis();
        logger.debug("time for aromaticity calculation: "
                + (after - before) + " milliseconds");
        logger.debug("Finished Aromaticity Detection");
        BitSet bitSet = new BitSet(size);

        hashes = findPaths(container, searchDepth);
        for (int hash : hashes) {
            position = (int) RandomNumber.generateMersenneTwisterRandomNumber(size, hash);
            bitSet.set(position);
        }

        return bitSet;
    }

    /**
     * Generates a fingerprint of the default size for the given AtomContainer.
     *
     *@param container The AtomContainer for which a Fingerprint is generated
     * @return
     * @throws CDKException  
     */
    @TestMethod("testGetFingerprint_IAtomContainer")
    @Override
    public BitSet getFingerprint(IAtomContainer container)
            throws CDKException {
        return getFingerprint(container, null);
    }

    /** {@inheritDoc}
     * @param atomContainer
     * @return 
     * @throws CDKException  
     */
    @Override
    public Map<String, Integer> getRawFingerprint(IAtomContainer atomContainer) throws CDKException {
        Map<String, Integer> uniquePaths = new TreeMap<String, Integer>();
        for (int hash : hashes) {
            int position = (int) RandomNumber.generateMersenneTwisterRandomNumber(size, hash);
            uniquePaths.put(new Integer(position).toString(), hash);
        }
        return uniquePaths;
    }

    /**
     * Get all paths of lengths 0 to the specified length.
     *
     * This method will find all paths upto length N starting from each
     * atom in the molecule and return the unique set of such paths.
     *
     * @param container The molecule to search
     * @param searchDepth The maximum path length desired
     * @return A map of path strings, keyed on themselves
     */
    protected int[] findPaths(IAtomContainer container, int searchDepth) {

        List<String> pseudoAtoms = new ArrayList<String>();
        int pseduoAtomCounter = 0;

        List<StringBuffer> allPaths = new ArrayList<StringBuffer>();

        Map<IAtom, Map<IAtom, IBond>> cache = new HashMap<IAtom, Map<IAtom, IBond>>();
        for (IAtom sourceAtom : container.atoms()) {
            List<List<IAtom>> pathsOfLengthUpto = PathTools.getPathsOfLengthUpto(container, sourceAtom, DEFAULT_SEARCH_DEPTH);
            for (List<IAtom> path : pathsOfLengthUpto) {
                StringBuffer sb = new StringBuffer();
                IAtom x = path.get(0);

                if (x instanceof IPseudoAtom) {
                    if (!pseudoAtoms.contains(x.getSymbol())) {
                        pseudoAtoms.add(pseduoAtomCounter++, x.getSymbol());
                    }
                    sb.append((char) (PeriodicTable.getElementCount()
                            + pseudoAtoms.indexOf(x.getSymbol()) + 1));
                } else {
                    Integer atnum = PeriodicTable.getAtomicNumber(x.getSymbol());
                    if (atnum != null) {
                        sb.append(toAtomPattern(x));
                    } else {
                        sb.append((char) PeriodicTable.getElementCount() + 1);
                    }
                }

                for (int i = 1; i < path.size(); i++) {
                    final IAtom[] y = {path.get(i)};
                    Map<IAtom, IBond> m = cache.get(x);
                    final IBond[] b = {m != null ? m.get(y[0]) : null};
                    if (b[0] == null) {
                        b[0] = container.getBond(x, y[0]);
                        cache.put(x,
                                new HashMap<IAtom, IBond>() {

                                    {
                                        put(y[0], b[0]);
                                    }
                                    private static final long serialVersionUID = 0xb3a7a32449fL;
                                });
                    }
                    sb.append(getBondSymbol(b[0]));
                    sb.append(toAtomPattern(y[0]));
                    x = y[0];
                }

                // we store the lexicographically lower one of the
                // string and its reverse
                StringBuffer revForm = new StringBuffer(sb);
                revForm.reverse();
                if (sb.toString().compareTo(revForm.toString()) <= 0) {
                    allPaths.add(sb);
                } else {
                    allPaths.add(revForm);
                }
            }
        }

        pseudoAtoms.clear();

        // now lets clean stuff up
        Set<String> cleanPath = new HashSet<String>();
        for (StringBuffer s : allPaths) {
            String s1 = s.toString().trim();
            if (s1.equals("")) {
                continue;
            }
            if (cleanPath.contains(s1)) {
                continue;
            }
            String s2 = s.reverse().toString().trim();
            if (cleanPath.contains(s2)) {
                continue;
            }
            cleanPath.add(s2);
        }

        // convert paths to hashes
        int[] hashes = new int[cleanPath.size()];
        int i = 0;
        for (String s : cleanPath) {
            hashes[i++] = new HashCodeBuilder(17, 37).append(s).toHashCode();
        }
        return hashes;
    }

    private String toAtomPattern(IAtom atom) {
        Double charge = atom.getCharge() == null ? 0. : atom.getCharge();
        Double stereoParity = atom.getStereoParity() == null ? 0. : atom.getStereoParity();
        Integer atomNum = atom.getAtomicNumber() == null ? 0 : atom.getAtomicNumber();
        int isRingAtom = 0;
        if (isRespectRingMatches()) {
            isRingAtom = atom.getFlag(CDKConstants.ISINRING) ? 1 : 0;
        }

        String atomConfiguration = atom.getSymbol()
                + ":" + charge.toString()
                + ":" + stereoParity.toString()
                + ":" + atomNum
                + ":" + isRingAtom;

        if (!patterns.containsKey(atomConfiguration)) {
            String generatedPattern = generateNewPattern();
            patterns.put(atomConfiguration, generatedPattern);
        }
        return patterns.get(atomConfiguration);
    }

    /**
     *  Gets the bondSymbol attribute of the Fingerprinter class
     *
     *@param  bond  Description of the Parameter
     *@return       The bondSymbol value
     */
    protected String getBondSymbol(IBond bond) {
        String bondSymbol = "";
        if (bond.getFlag(CDKConstants.ISAROMATIC)) {
            bondSymbol += "@";
        } else if (bond.getOrder() == IBond.Order.SINGLE) {
            bondSymbol += "-";
        } else if (bond.getOrder() == IBond.Order.DOUBLE) {
            bondSymbol += "=";
        } else if (bond.getOrder() == IBond.Order.TRIPLE) {
            bondSymbol += "#";
        }

        return bondSymbol;
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
        return size;
    }

    private String generateNewPattern() {
        int patternSize = patterns.size() + 1;
        StringBuilder st = new StringBuilder(patternSize);
        for (int i = 0; i < patternSize; i++) {
            st.append('*');
        }
        return st.toString();
    }

    /**
     * Should match rings to rings and non-rings to non-rings
     * @return the respect ring matches
     */
    @Override
    public boolean isRespectRingMatches() {
        return respectRingMatches;
    }

    /**
     * Ring matches are allowed and non-ring to ring matches are
     * discarded
     * @param respectRingMatches respect the ring-to-ring matches 
     * and discard non-ring to ring matches
     */
    @Override
    public void setRespectRingMatches(boolean respectRingMatches) {
        this.respectRingMatches = respectRingMatches;
    }
}
