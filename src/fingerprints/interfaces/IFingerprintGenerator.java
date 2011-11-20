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
package fingerprints.interfaces;

import java.util.BitSet;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk> 2007-2011
 */
public interface IFingerprintGenerator {

    BitSet getEStateFingerprinter(IAtomContainer mol) throws CDKException;

    int getEStateFingerprinterSize();

    BitSet getExtPubChemFingerprint(IAtomContainer mol) throws CDKException;

    int getExtPubChemFingerprintSize();

    BitSet getExtendedFingerprint(IAtomContainer mol) throws CDKException;

    int getExtendedFingerprintSize();

    BitSet getGraphOnlyFingerprinter(IAtomContainer mol) throws CDKException;

    int getGraphOnlyFingerprinterSize();

    BitSet getHashedFingerprint(IAtomContainer mol) throws CDKException;

    int getHashedFingerprintDepth();

    int getHashedFingerprintSize();

    BitSet getHybridizationFingerprinter(IAtomContainer mol) throws Exception;

    int getHybridizationFingerprinterSize() throws Exception;

    BitSet getMACCSFingerprinter(IAtomContainer mol) throws CDKException;

    int getMACCSFingerprinterSize();

    BitSet getPubChemFingerprint(IAtomContainer mol) throws CDKException;

    int getPubChemFingerprintSize();

    BitSet getStandardFingerprint(IAtomContainer mol) throws CDKException;

    int getStandardFingerprintDepth();

    int getStandardFingerprintSize();

    BitSet getSubstructureFingerprinter(IAtomContainer mol) throws CDKException;

    int getSubstructureFingerprinterSize();
}
