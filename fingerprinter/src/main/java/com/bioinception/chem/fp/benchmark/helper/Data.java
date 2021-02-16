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
package com.bioinception.chem.fp.benchmark.helper;

import java.util.BitSet;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public class Data {

    /**
     * @return the fpName
     */
    public String getFPName() {
        return fpName;
    }

    /**
     * @param fpName the fpName to set
     */
    public void setFPName(String fpName) {
        this.fpName = fpName;
    }

    private String fpName;
    private BitSet fingerprint;
    private final IAtomContainer atomContainer;

    /**
     * Store the fingerprint and its structure
     *
     * @param fingerprint
     * @param atomContainer
     */
    public Data(BitSet fingerprint, IAtomContainer atomContainer) {
        this.fingerprint = fingerprint;
        this.atomContainer = atomContainer;
    }

    /**
     * Store the fingerprint and its structure
     *
     * @param atomContainer
     */
    public Data(IAtomContainer atomContainer) {
        this.atomContainer = atomContainer;
    }

    /**
     * @return the fingerprint
     */
    public BitSet getFingerprint() {
        return fingerprint;
    }

    /**
     * @return the atomContainer
     */
    public IAtomContainer getAtomContainer() {
        return atomContainer;
    }
}
