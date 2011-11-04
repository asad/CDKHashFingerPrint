/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package helper;

import java.util.BitSet;
import org.openscience.cdk.interfaces.IAtomContainer;

/**
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public class Data {

    private BitSet fingerprint;
    private IAtomContainer atomContainer;

    /**
     * Store the fingerprint and its structure
     * @param fingerprint
     * @param atomContainer
     */
    public Data(BitSet fingerprint, IAtomContainer atomContainer) {
        this.fingerprint = fingerprint;
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
