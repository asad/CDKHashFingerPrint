/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package fingerprints;

import java.util.BitSet;
import java.util.Map;
import org.openscience.cdk.annotations.TestMethod;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.ringsearch.AllRingsFinder;

/**
 *
 * @author Asad
 */
public interface IFingerprinter extends org.openscience.cdk.fingerprint.IFingerprinter {
    /**
     * The default search depth used to create the fingerprints.
     */
    int DEFAULT_SEARCH_DEPTH = 8;
    /**
     * The default length of created fingerprints.
     */
    int DEFAULT_SIZE = 2048;

    /**
     * Generates a fingerprint of the default size for the given AtomContainer.
     *
     * @param container The AtomContainer for which a Fingerprint is generated
     * @param ringFinder An instance of
     * {@link org.openscience.cdk.ringsearch.AllRingsFinder}
     * @exception CDKException if there is a timeout in ring or aromaticity
     * perception
     * @return A {@link BitSet} representing the fingerprint
     */
    @TestMethod(value = "testGetFingerprint_IAtomContainer")
    BitSet getFingerprint(IAtomContainer container, AllRingsFinder ringFinder) throws CDKException;

    /**
     * Generates a fingerprint of the default size for the given AtomContainer.
     *
     * @param container The AtomContainer for which a Fingerprint is generated
     * @return
     * @throws CDKException
     */
    @TestMethod(value = "testGetFingerprint_IAtomContainer")
    @Override
    BitSet getFingerprint(IAtomContainer container) throws CDKException;

    /**
     * {@inheritDoc}
     * @param atomContainer
     * @return
     * @throws CDKException
     */
    @Override
    Map<String, Integer> getRawFingerprint(IAtomContainer atomContainer) throws CDKException;

    @TestMethod(value = "testGetSearchDepth")
    int getSearchDepth();

    /**
     * 
     * @return
     */
    @TestMethod(value = "testGetSize")
    @Override
    int getSize();

    /**
     * Should match rings to rings and non-rings to non-rings
     * @return the respect ring matches
     */
    boolean isRespectRingMatches();

    /**
     * Ring matches are allowed and non-ring to ring matches are
     * discarded
     * @param respectRingMatches respect the ring-to-ring matches
     * and discard non-ring to ring matches
     */
    void setRespectRingMatches(boolean respectRingMatches);
    
}
