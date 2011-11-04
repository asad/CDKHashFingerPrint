/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package helper;

import org.apache.commons.math.random.MersenneTwister;
import org.apache.commons.math.random.RandomAdaptor;
import org.apache.commons.math.random.RandomGenerator;

/**
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public class RandomNumber {

    /**
     * Mersenne Twister Random Number
     * @param maximum
     * @param hashCode
     * @return
     */
    public static long generateMersenneTwisterRandomNumber(int maximum, long hashCode) {
        RandomGenerator rg = new RandomAdaptor(new MersenneTwister(hashCode));
        return rg.nextInt(maximum);
    }
    
        /**
     * Mersenne Twister Random Number
     * @param maximum
     * @return
     */
    public static long generateMersenneTwisterRandomNumber(int maximum) {
        RandomGenerator rg = new RandomAdaptor(new MersenneTwister());
        return rg.nextInt(maximum);
    }
}
