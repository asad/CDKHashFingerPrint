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

import java.io.Serializable;
import java.util.*;
import org.openscience.cdk.graph.PathTools;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IPseudoAtom;
import org.openscience.cdk.tools.periodictable.PeriodicTable;
import fingerprints.interfaces.IWalker;
import org.openscience.cdk.interfaces.*;

/**
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk> 2007-2011
 */
public class MoleculeWalker implements IWalker, Serializable {

    private static final long serialVersionUID = 0x3b728f46;
    private final IAtomContainer atomContainer;
    private final Set<String> cleanPath;
    private final List<String> pseudoAtoms;
    private int pseduoAtomCounter;
    private final List<StringBuffer> allPaths;
    private final Map<IAtom, Map<IAtom, IBond>> cache;
    private final int maximumDepth;

    /**
     *
     * @param maximumDepth
     * @param atomContainer
     */
    public MoleculeWalker(int maximumDepth, IAtomContainer atomContainer) {
        this.cleanPath = new HashSet<String>();
        this.atomContainer = atomContainer;
        this.maximumDepth = maximumDepth;
        this.pseudoAtoms = new ArrayList<String>();
        this.pseduoAtomCounter = 0;
        this.allPaths = new ArrayList<StringBuffer>();
        this.cache = new HashMap<IAtom, Map<IAtom, IBond>>();
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
     * @return the cleanPath
     */
    @Override
    public Set<String> getPaths() {
        return Collections.unmodifiableSet(cleanPath);
    }

    /**
     * @return the cleanPath
     */
    @Override
    public int getPathCount() {
        return cleanPath.size();
    }

    private void findPaths() {
        pseudoAtoms.clear();
        getBondPaths();
        for (StringBuffer s : allPaths) {
            String s1 = s.toString().trim();
            if (s1.equals("")) {
                continue;
            }
            cleanPath.add(s1);
        }
    }

    private void getBondPaths() {
        for (IAtom sourceAtom : atomContainer.atoms()) {
            List<List<IAtom>> pathsOfLengthUpto = PathTools.getPathsOfLengthUpto(atomContainer, sourceAtom, getMaximumDepth());
            for (List<IAtom> path : pathsOfLengthUpto) {
                StringBuffer sb = new StringBuffer();
                IAtom x = path.get(0);

                if (x instanceof IPseudoAtom) {
                    if (!pseudoAtoms.contains(x.getSymbol())) {
                        pseudoAtoms.add(pseduoAtomCounter, x.getSymbol());
                        pseduoAtomCounter += 1;
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
                    if (!allPaths.contains(sb)) {
                        allPaths.add(sb);
                    }
                } else {
                    if (!allPaths.contains(revForm)) {
                        allPaths.add(revForm);
                    }
                }
            }
        }
    }

    private String toAtomPattern(IAtom atom) {
        return atom.getSymbol();
    }

    /**
     * Gets the bondSymbol attribute of the HashedFingerprinter class
     *
     * @param bond Description of the Parameter
     * @return The bondSymbol value
     */
    private char getBondSymbol(IBond bond) {
        if (isSP2Bond(bond)) {
            return '@';
        } else {
            switch (bond.getOrder()) {
                case SINGLE:
                    return '1';
                case DOUBLE:
                    return '2';
                case TRIPLE:
                    return '3';
                case QUADRUPLE:
                    return '4';
                default:
                    return '5';
            }
        }
    }

    /**
     * Returns true if the bond binds two atoms, and both atoms are SP2.
     */
    private boolean isSP2Bond(IBond bond) {
        if (bond.getAtomCount() == 2
                && bond.getAtom(0).getHybridization() == IAtomType.Hybridization.SP2
                && bond.getAtom(1).getHybridization() == IAtomType.Hybridization.SP2) {
            return true;
        }
        return false;
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
