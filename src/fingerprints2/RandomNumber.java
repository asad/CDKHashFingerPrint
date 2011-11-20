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
package fingerprints2;

import org.apache.commons.math.random.MersenneTwister;
import org.apache.commons.math.random.RandomAdaptor;
import org.apache.commons.math.random.RandomGenerator;

/**
 * @author Syed Asad Rahman <asad@ebi.ac.uk> 2007-2011
 */
public class RandomNumber {

    /**
     * Mersenne Twister Random Number 
     * for a hashcode within a range between 0 to maximum 
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
