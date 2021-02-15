/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package com.bioinception.chem.fp.fingerprints.interfaces;

import java.util.Set;

/**
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk> 2007-2011
 */
public interface IWalker {

    /**
     * @return the maximumDepth
     */
    int getMaximumDepth();

    /**
     * @return the cleanPath
     */
    int getPathCount();

    /**
     * @return the cleanPath
     */
    Set<String> getPaths();
}
