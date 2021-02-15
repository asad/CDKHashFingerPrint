/*
 * Copyright (C) 2003-2021 Syed Asad Rahman <asad @ ebi.ac.uk>.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301  USA
 */
package com.bioinception.chem.fp.fingerprints.interfaces;

import java.util.BitSet;
import java.util.Map;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.fingerprint.IBitFingerprint;
import org.openscience.cdk.fingerprint.ICountFingerprint;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.ringsearch.AllRingsFinder;

/**
 *
 * @author Asad
 */
public interface IFingerprinter extends org.openscience.cdk.fingerprint.IFingerprinter {

    /**
     * Generates a fingerprint of the default fingerprintLength for the given
     * AtomContainer.
     *
     * @TestMethod(value = "testGetFingerprint_IAtomContainer")
     * @param container The AtomContainer for which a Fingerprint is generated
     * @param ringFinder An instance of
     * {@link org.openscience.cdk.ringsearch.AllRingsFinder}
     * @exception CDKException if there is a timeout in ring or aromaticity
     * perception
     * @return A {@link BitSet} representing the fingerprint
     */
    IBitFingerprint getBitFingerprint(IAtomContainer container, AllRingsFinder ringFinder) throws CDKException;

    /**
     * Generates a fingerprint of the default fingerprintLength for the given
     * AtomContainer.
     *
     * @TestMethod(value = "testGetFingerprint_IAtomContainer")
     * @param container The AtomContainer for which a Fingerprint is generated
     * @return
     * @throws CDKException
     */
    @Override
    IBitFingerprint getBitFingerprint(IAtomContainer container) throws CDKException;

    @Override
    ICountFingerprint getCountFingerprint(IAtomContainer iac) throws CDKException;

    @Override
    Map<String, Integer> getRawFingerprint(IAtomContainer atomContainer) throws CDKException;

    /**
     * @TestMethod(value = "testGetSearchDepth")
     * @return
     */
    int getSearchDepth();

    /*
    @TestMethod(value = "testGetSize")
     */
    @Override
    int getSize();

    /**
     * @return the respectFormalCharges
     */
    boolean isRespectFormalCharges();

    /**
     * Should match rings to rings and non-rings to non-rings
     *
     * @return the respect ring matches
     */
    boolean isRespectRingMatches();

    /**
     * @param respectFormalCharges the flag to set if formal charge is checked
     */
    void setRespectFormalCharges(boolean respectFormalCharges);

    /**
     * Ring matches are allowed and non-ring to ring matches are discarded
     *
     * @param respectRingMatches respect the ring-to-ring matches and discard
     * non-ring to ring matches
     */
    void setRespectRingMatches(boolean respectRingMatches);
}
