/* $Revision$ $Author$ $Date$
 *
 * Copyright (C) 2002-2007  Christoph Steinbeck <steinbeck@users.sf.net>
 *               2020-2021  Syed Asad Rahman <asad@ebi.ac.uk>
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

import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashSet;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;

/**
 *
 * @author Syed Asad Rahman <asad.rahman@bioinceptionlabs.com>
 */
class State {

    /**
     * @return the numPaths
     */
    public int getNumPaths() {
        return numPaths;
    }

    /**
     * @param numPaths the numPaths to set
     */
    public void setNumPaths(int numPaths) {
        this.numPaths = numPaths;
    }

    /**
     * @return the maxDepth
     */
    public int getMaxDepth() {
        return maxDepth;
    }

    /**
     * @return the minDepth
     */
    public int getMinDepth() {
        return minDepth;
    }

    /**
     * @return the bpath
     */
    public List<IBond> getBondPath() {
        return bpath;
    }

    /**
     * @return the apath
     */
    public List<IAtom> getAtomPath() {
        return apath;
    }

    private int numPaths = 0;
    private final Random rand = new Random();
    private final BitSet fp;
    private final IAtomContainer mol;
    private final Set<IAtom> visited = new HashSet<>();
    private final List<IAtom> apath = new ArrayList<>();
    private final List<IBond> bpath = new ArrayList<>();
    private final int maxDepth;
    private final int minDepth;
    private final int fpsize;
    private final Map<IAtom, List<IBond>> cache = new IdentityHashMap<>();
    public StringBuilder buffer = new StringBuilder();

    public State(IAtomContainer mol, BitSet fp, int fpsize, int minDepth, int maxDepth) {
        this.mol = mol;
        this.fp = fp;
        this.fpsize = fpsize;
        this.minDepth = minDepth;
        this.maxDepth = maxDepth;
    }

    List<IBond> getBonds(IAtom atom) {
        List<IBond> bonds = cache.get(atom);
        if (bonds == null) {
            bonds = mol.getConnectedBondsList(atom);
            cache.put(atom, bonds);
        }
        return bonds;
    }

    boolean visit(IAtom a) {
        return visited.add(a);
    }

    boolean unvisit(IAtom a) {
        return visited.remove(a);
    }

    void push(IAtom atom, IBond bond) {
        getAtomPath().add(atom);
        if (bond != null) {
            getBondPath().add(bond);
        }
    }

    void pop() {
        if (!apath.isEmpty()) {
            getAtomPath().remove(getAtomPath().size() - 1);
        }
        if (!bpath.isEmpty()) {
            getBondPath().remove(getBondPath().size() - 1);
        }
    }

    void addHash(int x) {
        setNumPaths(getNumPaths() + 1);
        rand.setSeed(x);
        // XXX: fp.set(x % size); would work just as well but would encode a
        //      different bit
        fp.set(rand.nextInt(fpsize));
    }
}
