/* $Revision$ $Author$ $Date$
 *
 * Copyright (C) 2002-2007  Christoph Steinbeck <steinbeck@users.sf.net>
 *               2021  Syed Asad Rahman <asad@ebi.ac.uk>
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
package com.bioinception.chem.fp.fingerprints.bi;

import java.util.BitSet;
import java.util.List;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

/**
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
public class PathEncoder {

    static void encodePaths(IAtomContainer mol, int minDepth, int maxDepth, BitSet fp, int size, int pathLimit, boolean hashPseudoAtoms) throws CDKException {
        State state = new State(mol, fp, size, minDepth + 1, maxDepth + 1);
        for (IAtom atom : mol.atoms()) {
            state.setNumPaths(0);
            state.visit(atom);
            traversePaths(state, atom, null, pathLimit, hashPseudoAtoms);
            state.unvisit(atom);
        }
    }

    /**
     *
     * @param state
     * @param beg
     * @param prev
     * @param pathLimit
     * @param hashPseudoAtoms
     * @throws CDKException
     */
    static void traversePaths(State state, IAtom beg, IBond prev, int pathLimit, boolean hashPseudoAtoms) throws CDKException {
        if (!hashPseudoAtoms && isPseudoAtom(beg)) {
            return;
        }
        state.push(beg, prev);
        state.addHash(encodeUniquePath(state.getAtomPath(), state.getBondPath(), state.buffer));
        if (state.getNumPaths() > pathLimit) {
            throw new CDKException("Too many paths! Structure is likely a cage, reduce path length or increase path limit");
        }
//        System.out.println("state.getMinDepth() "+state.getMinDepth());
//        System.out.println("state.getMaxDepth() "+state.getMaxDepth());
        if (state.getAtomPath().size() >= state.getMinDepth()
                && state.getAtomPath().size() < state.getMaxDepth()) {
//            System.out.println("state.getAtomPath().size() "+state.getAtomPath().size());
            for (IBond bond : state.getBonds(beg)) {
                if (bond.equals(prev)) {
                    continue;
                }
                final IAtom nbr = bond.getOther(beg);
                if (state.visit(nbr)) {
                    traversePaths(state, nbr, bond, pathLimit, hashPseudoAtoms);
                    state.unvisit(nbr); // traverse all paths
                }
            }
        }
        state.pop();
    }

    /**
     *
     * @param apath
     * @param bpath
     * @param buffer
     * @return
     */
    static int encodeUniquePath(List<IAtom> apath, List<IBond> bpath, StringBuilder buffer) {
        if (bpath.isEmpty()) {
            return getAtomSymbol(apath.get(0)).hashCode();
        }
        final int x;
        if (compare(apath, bpath) >= 0) {
            x = hashPath(apath, bpath);
        } else {
            x = hashRevPath(apath, bpath);
        }
        return x;
    }

    static int hashPath(List<IAtom> apath, List<IBond> bpath) {
        int hash = 0;
        hash = appendHash(hash, getAtomSymbol(apath.get(0)));
        for (int i = 1; i < apath.size(); i++) {
            final IAtom next = apath.get(i);
            final IBond bond = bpath.get(i - 1);
            hash = appendHash(hash, getBondSymbol(bond));
            hash = appendHash(hash, getAtomSymbol(next));
        }
        return hash;
    }

    static int hashRevPath(List<IAtom> apath, List<IBond> bpath) {
        int hash = 0;
        int last = apath.size() - 1;
        hash = appendHash(hash, getAtomSymbol(apath.get(last)));
        for (int i = last - 1; i >= 0; i--) {
            final IAtom next = apath.get(i);
            final IBond bond = bpath.get(i);
            hash = appendHash(hash, getBondSymbol(bond));
            hash = appendHash(hash, getAtomSymbol(next));
        }
        return hash;
    }

    static int appendHash(int hash, String str) {
        int len = str.length();
        for (int i = 0; i < len; i++) {
            hash = 31 * hash + str.charAt(0);
        }
        return hash;
    }

    /**
     * Gets the bondSymbol attribute of the Fingerprinter class
     *
     * @param bond Description of the Parameter
     * @return The bondSymbol value
     */
    static String getBondSymbol(IBond bond) {
        if (bond.isAromatic()) {
            return ":";
        }
        switch (bond.getOrder()) {
            case SINGLE:
                return "-";
            case DOUBLE:
                return "=";
            case TRIPLE:
                return "#";
            default:
                return "";
        }
    }

    static String getAtomSymbol(IAtom atom) {
        // XXX: backwards compatibility
        // This is completely random, I believe the intention is because
        // paths were reversed with string manipulation to de-duplicate
        // (only the lowest lexicographically is stored) however this
        // doesn't work with multiple atom symbols:
        // e.g. Fe-C => C-eF vs C-Fe => eF-C
        // A dirty hack is to replace "common" symbols with single letter
        // equivalents so the reversing is less wrong
        switch (getElem(atom)) {
            case 0:  // *
                return "*";
            case 6:  // C
                return "C";
            case 7:  // N
                return "N";
            case 8:  // O
                return "O";
            case 17: // Cl
                return "X";
            case 35: // Br
                return "Z";
            case 14: // Si
                return "Y";
            case 33: // As
                return "D";
            case 3: // Li
                return "L";
            case 34: // Se
                return "E";
            case 11:  // Na
                return "G";
            case 20:  // Ca
                return "J";
            case 13:  // Al
                return "A";
        }
        return atom.getSymbol();
    }

    static int getElem(IAtom atom) {
        Integer elem = atom.getAtomicNumber();
        if (elem == null) {
            elem = 0;
        }
        return elem;
    }

    /**
     * Compares a path of atoms with it's self to give the lexicographically
     * lowest traversal (forwards or backwards).
     *
     * @param apath path of atoms
     * @param bpath path of bonds
     * @return &lt;0 forward is lower &gt;0 reverse is lower
     */
    static int compare(List<IAtom> apath, List<IBond> bpath) {
        int i = 0;
        int len = apath.size();
        int j = len - 1;
        int cmp = compare(apath.get(i), apath.get(j));
        if (cmp != 0) {
            return cmp;
        }
        i++;
        j--;
        while (j != 0) {
            cmp = compare(bpath.get(i - 1), bpath.get(j));
            if (cmp != 0) {
                return cmp;
            }
            cmp = compare(apath.get(i), apath.get(j));
            if (cmp != 0) {
                return cmp;
            }
            i++;
            j--;
        }
        return 0;
    }

    /**
     * Compares bonds symbols lexicographical
     *
     * @param a bond a
     * @param b bond b
     * @return comparison &lt;0 a is less than b, &gt;0 a is more than b
     */
    static int compare(IBond a, IBond b) {
        return getBondSymbol(a).compareTo(getBondSymbol(b));
    }

    /**
     * Compares atom symbols lexicographical
     *
     * @param a atom a
     * @param b atom b
     * @return comparison &lt;0 a is less than b, &gt;0 a is more than b
     */
    static int compare(IAtom a, IAtom b) {
        final int elemA = getElem(a);
        final int elemB = getElem(b);
        if (elemA == elemB) {
            return 0;
        }
        return getAtomSymbol(a).compareTo(getAtomSymbol(b));
    }

    /**
     *
     * @param a
     * @return
     */
    static boolean isPseudoAtom(IAtom a) {
        return getElem(a) == 0;
    }
}
