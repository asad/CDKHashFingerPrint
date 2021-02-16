/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.bioinception.chem.fp.fingerprints.helper;

import java.util.BitSet;

/**
 *
 * @author asad
 */
public class FingerprinterTool {

    /**
     *
     * @param bs1
     * @param bs2
     * @return
     */
    public static boolean isSubset(BitSet bs1, BitSet bs2) {
//        System.out.println("bs1 " + bs1.cardinality());
//        System.out.println("bs2 " + bs2.cardinality());
        BitSet clone = (BitSet) bs1.clone();
        clone.and(bs2);
        return bs2.cardinality() > 0 && clone.equals(bs2);
    }
}
