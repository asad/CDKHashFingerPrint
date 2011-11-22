/* $Revision$ $Author$ $Date$
 *
 * Copyright (C) 2011       Syed Asad Rahman <asad@ebi.ac.uk>
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
package fingerprints.helper;

import fingerprints.interfaces.IWalker;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.graph.PathTools;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.tools.periodictable.PeriodicTable;

/**
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk> 2007-2011
 */
public class MoleculeWalker implements IWalker, Serializable {

    private static final long serialVersionUID = 0x3b728f46;
    private int maximumDepth;
    private final IAtomContainer atomContainer;
    private Set<String> cleanPath;
    private static final Map<String, String> patterns = new HashMap<String, String>() {

        private static final long serialVersionUID = 3348458944893841L;

        {
            put("R", "*");
            put("X", "**");
        }
    };

    /**
     * 
     * @param maximumDepth
     * @param atomContainer
     */
    public MoleculeWalker(int maximumDepth, IAtomContainer atomContainer) {
        this.cleanPath = new HashSet<String>();
        this.atomContainer = atomContainer;
        this.maximumDepth = maximumDepth;
        findPaths();
    }

    /**
     * 
     * @param atomContainer
     */
    public MoleculeWalker(IAtomContainer atomContainer) {
        this(atomContainer.getAtomCount(), atomContainer);
    }

    /**
     * @return the maximumDepth
     */
    @Override
    public int getMaximumDepth() {
        return maximumDepth;
    }

    /**
     * @param maximumDepth 
     */
    @Override
    public void setMaximumDepth(int maximumDepth) {
        this.maximumDepth = maximumDepth;
    }

    /**
     * @return the cleanPath
     */
    @Override
    public Set<String> getPaths() {
        return cleanPath;
    }

    /**
     * @return the cleanPath
     */
    @Override
    public int getPathCount() {
        return cleanPath.size();
    }

    private void findPaths() {
        List<String> pseudoAtoms = new ArrayList<String>();
        int pseduoAtomCounter = 0;

        List<StringBuffer> allPaths = new ArrayList<StringBuffer>();

        Map<IAtom, Map<IAtom, IBond>> cache = new HashMap<IAtom, Map<IAtom, IBond>>();
        for (IAtom sourceAtom : atomContainer.atoms()) {
            List<List<IAtom>> pathsOfLengthUpto = PathTools.getPathsOfLengthUpto(atomContainer, sourceAtom, getMaximumDepth());
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
                        b[0] = atomContainer.getBond(x, y[0]);
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
    }

    private String toAtomPattern(IAtom atom) {
        Double stereoParity = atom.getStereoParity() == null ? 0. : atom.getStereoParity();
        Integer atomicNumber = atom.getAtomicNumber() == null ? 0 : atom.getAtomicNumber();
        
        String atomConfiguration = atom.getSymbol()
                + ":" + stereoParity.toString()
                + ":" + atomicNumber;

        if (!patterns.containsKey(atomConfiguration)) {
            String generatedPattern = generateNewPattern();
            patterns.put(atomConfiguration, generatedPattern);
        }
        return patterns.get(atomConfiguration);
    }

    /**
     *  Gets the bondSymbol attribute of the HashedFingerprinter class
     *
     *@param  bond  Description of the Parameter
     *@return       The bondSymbol value
     */
    private String getBondSymbol(IBond bond) {
        String bondSymbol = "";
        if (bond.getFlag(CDKConstants.ISAROMATIC)) {
            bondSymbol += "@";
        } else if (bond.getOrder() == IBond.Order.SINGLE) {
            bondSymbol += "-";
        } else if (bond.getOrder() == IBond.Order.DOUBLE) {
            bondSymbol += "=";
        } else if (bond.getOrder() == IBond.Order.TRIPLE) {
            bondSymbol += "#";
        } else if (bond.getOrder() == IBond.Order.QUADRUPLE) {
            return "&";
        }

        return bondSymbol;
    }

    private String generateNewPattern() {
        int patternSize = patterns.size() + 1;
        StringBuilder st = new StringBuilder(patternSize);
        for (int i = 0; i < patternSize; i++) {
            st.append('*');
        }
        return st.toString();
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (String path : cleanPath) {
            sb.append(path).append("->");
        }
        return sb.toString();
    }
}
